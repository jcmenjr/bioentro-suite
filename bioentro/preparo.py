#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preparo v0.1.0 — Genome Preparation Utility for bioentro-suite

Prepares genomic files for downstream analysis with Bakta, Panaroo,
and bioentro-suite. Handles naming inconsistencies from NCBI downloads,
reference genomes, and locally sequenced assemblies.

Problem this solves
───────────────────
NCBI genomes arrive with names like:
    GCF_000006765.1_ASM676v1_genomic.fna

Bakta uses the filename stem as the locus tag prefix, producing:
    GCF_000006765.1_ASM676v1_genomic_00001  ← unusable as an ID

Panaroo uses the GFF filename as the sample name in its output matrix,
so inconsistent names propagate into all downstream analyses.

This script ensures clean, consistent, biologically meaningful names
are set BEFORE Bakta runs, so all downstream IDs are correct.

Subcommands
───────────
  detect   : Scan a directory of genome files, auto-detect naming patterns
             (NCBI, assembler output, custom), and generate a sample_map.tsv
             draft for manual review before renaming.

  rename   : Apply a sample_map.tsv to rename genome files. Optionally
             rewrites FASTA headers to match the new name. Non-destructive:
             always writes to a new output directory.

  check    : Verify Bakta output completeness and ID consistency across
             all samples. Reports missing files and ID mismatches.

  panaroo  : Generate a ready-to-run Panaroo shell script from a directory
             of Bakta GFF outputs.

Supported source types
──────────────────────
  ncbi      : GCF_*/GCA_* accession-based names (RefSeq / GenBank)
  assembler : SPAdes / Flye / Unicycler output (NODE_*, contig_*, scaffold_*)
  custom    : Any other pattern — falls back to sequential numbering

Usage:
    preparo detect  -i genomes_dir/ -o sample_map.tsv
    preparo rename  -i genomes_dir/ -m sample_map.tsv -o renamed_dir/
    preparo check   -i bakta_outputs/
    preparo panaroo -i bakta_outputs/ -o run_panaroo.sh

Requirements:
    biopython >= 1.79

Author:  Julio Méndez
License: MIT
"""

from __future__ import annotations

__version__ = "0.1.0"

# ---------------------------------------------------------------------------
# Standard library
# ---------------------------------------------------------------------------
import argparse
import csv
import logging
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

# ---------------------------------------------------------------------------
# Third-party
# ---------------------------------------------------------------------------
from Bio import SeqIO

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# FASTA extensions we recognize as genome files
GENOME_EXTENSIONS: Tuple[str, ...] = (".fna", ".fa", ".fasta")

# GFF extensions from Bakta
GFF_EXTENSIONS: Tuple[str, ...] = (".gff", ".gff3")

# All file types Bakta produces per sample — used for completeness check
BAKTA_EXPECTED_SUFFIXES: Tuple[str, ...] = (
    ".gff3",    # annotation + sequence (input to Panaroo)
    ".faa",     # protein sequences (input to bioentro protein mode)
    ".ffn",     # nucleotide gene sequences
    ".fna",     # genome sequence
    ".tsv",     # annotation table
    ".txt",     # summary stats
    ".log",     # run log
)

# Regex patterns for source-type detection
_NCBI_PATTERN     = re.compile(r"^(GCF|GCA)_\d{9}\.\d+", re.IGNORECASE)
_ASSEMBLER_PATTERN = re.compile(
    r"^(NODE|contig|scaffold|sequence|ctg|scf)[\._\-]", re.IGNORECASE
)

# Characters not allowed in locus tag prefixes (INSDC standard)
_INVALID_LOCUS_CHARS = re.compile(r"[^A-Za-z0-9_]")

# Maximum recommended locus tag prefix length (INSDC guideline: ≤ 12 chars)
MAX_LOCUS_TAG_LEN: int = 12

# ---------------------------------------------------------------------------
# Logger
# ---------------------------------------------------------------------------
logger = logging.getLogger(__name__)


def _configure_logging(verbose: bool = False) -> None:
    """Configure root logger. Called only from main()."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s — %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


# ---------------------------------------------------------------------------
# Custom exceptions
# ---------------------------------------------------------------------------

class PreparoError(Exception):
    """Base class for preparo errors."""

class InputDirectoryError(PreparoError):
    """Missing or empty input directory."""

class SampleMapError(PreparoError):
    """Malformed or missing sample map file."""

class NamingConflictError(PreparoError):
    """Duplicate or conflicting names in sample map."""


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SampleEntry:
    """One row in sample_map.tsv."""
    original_file: str      # original filename (without directory)
    new_name: str           # clean sample name (used as locus tag prefix)
    source_type: str        # ncbi | assembler | custom
    notes: str = ""         # auto-generated hints for the user


@dataclass
class CheckResult:
    """Result of checking one Bakta sample directory."""
    sample_name: str
    gff_found: bool
    faa_found: bool
    missing_files: List[str] = field(default_factory=list)
    id_prefix_ok: bool = True
    id_prefix_found: str = ""
    warnings: List[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.gff_found and self.faa_found and not self.missing_files and self.id_prefix_ok


# ---------------------------------------------------------------------------
# Source-type detection helpers
# ---------------------------------------------------------------------------

def detect_source_type(stem: str) -> str:
    """Classify a filename stem as 'ncbi', 'assembler', or 'custom'.

    Args:
        stem: Filename without extension (e.g. 'GCF_000006765.1_ASM676v1_genomic').

    Returns:
        One of: 'ncbi', 'assembler', 'custom'.
    """
    if _NCBI_PATTERN.match(stem):
        return "ncbi"
    if _ASSEMBLER_PATTERN.match(stem):
        return "assembler"
    return "custom"


def _suggest_name_from_ncbi(stem: str, fna_path: Path, index: int) -> Tuple[str, str]:
    """Suggest a clean name for an NCBI genome file.

    Strategy:
      1. Try to extract organism info from the FASTA header [organism=...].
      2. Fall back to the GCF/GCA accession (short form).

    Args:
        stem:     Filename stem.
        fna_path: Path to the genome file (read first record header).
        index:    1-based index for disambiguation.

    Returns:
        (suggested_name, note) tuple.
    """
    # Try to read the first FASTA header
    try:
        record = next(SeqIO.parse(fna_path, "fasta"))
        desc = record.description

        # [organism=Pseudomonas aeruginosa PAO1] → Pseudomonas_aeruginosa_PAO1
        match = re.search(r"\[organism=([^\]]+)\]", desc)
        if match:
            org = match.group(1).strip()
            # Abbreviate: "Pseudomonas aeruginosa PAO1" → "Pae_PAO1"
            parts = org.split()
            if len(parts) >= 2:
                abbreviated = f"{parts[0][0]}{parts[1][:2]}_{('_'.join(parts[2:]) or str(index))}"
            else:
                abbreviated = org.replace(" ", "_")
            clean = _sanitize_name(abbreviated)
            return clean, f"organism='{org}' from FASTA header"

        # [strain=PAO1] fallback
        match = re.search(r"\[strain=([^\]]+)\]", desc)
        if match:
            strain = _sanitize_name(match.group(1).strip())
            return strain, f"strain='{match.group(1)}' from FASTA header"

    except (StopIteration, Exception) as exc:
        logger.debug("Could not read FASTA header for '%s': %s", fna_path.name, exc)

    # Last resort: short accession
    acc_match = re.match(r"(GCF|GCA)_(\d+)\.\d+", stem, re.IGNORECASE)
    if acc_match:
        short = f"{acc_match.group(1)}{acc_match.group(2)[:6]}_{index:02d}"
        return _sanitize_name(short), "accession-based fallback"

    return f"genome_{index:03d}", "no organism info found — please edit"


def _suggest_name_from_assembler(stem: str, index: int) -> Tuple[str, str]:
    """Suggest a name for an assembler-output genome.

    These files typically have no useful organism info in the name.
    We generate a sequential placeholder the user should edit.
    """
    return f"sample_{index:03d}", "assembler output — edit with your sample ID"


def _sanitize_name(name: str) -> str:
    """Convert a string to a valid locus tag prefix.

    Rules (INSDC standard):
      - Only alphanumeric + underscore
      - No leading digits
      - Maximum MAX_LOCUS_TAG_LEN characters
      - Uppercase
    """
    clean = _INVALID_LOCUS_CHARS.sub("_", name).strip("_")
    # Remove consecutive underscores
    clean = re.sub(r"_+", "_", clean)
    # Cannot start with a digit
    if clean and clean[0].isdigit():
        clean = "S" + clean
    return clean[:MAX_LOCUS_TAG_LEN].upper()


# ---------------------------------------------------------------------------
# detect subcommand
# ---------------------------------------------------------------------------

def cmd_detect(
    input_dir: Path,
    output_map: Path,
    extensions: Tuple[str, ...] = GENOME_EXTENSIONS,
) -> List[SampleEntry]:
    """Scan a directory and generate a sample_map.tsv draft.

    The generated file is a starting point for manual editing —
    especially for assembler outputs where no organism info is available.

    Args:
        input_dir:  Directory containing genome FASTA files.
        output_map: Path where sample_map.tsv will be written.
        extensions: File extensions to consider as genome files.

    Returns:
        List of SampleEntry objects (also written to output_map).

    Raises:
        InputDirectoryError: If input_dir does not exist or is empty.
    """
    if not input_dir.is_dir():
        raise InputDirectoryError(f"Not a directory: '{input_dir}'")

    files = sorted(
        f for f in input_dir.iterdir()
        if f.is_file() and f.suffix.lower() in extensions
    )
    if not files:
        raise InputDirectoryError(
            f"No genome files found in '{input_dir}' "
            f"(looking for {extensions})."
        )

    logger.info("Found %d genome file(s) in '%s'.", len(files), input_dir)

    entries: List[SampleEntry] = []
    seen_names: Dict[str, int] = {}   # name → count, for deduplication

    for i, fna_path in enumerate(files, start=1):
        stem = fna_path.stem
        source = detect_source_type(stem)

        if source == "ncbi":
            name, note = _suggest_name_from_ncbi(stem, fna_path, i)
        elif source == "assembler":
            name, note = _suggest_name_from_assembler(stem, i)
        else:
            # Custom: sanitize the stem directly
            name = _sanitize_name(stem)
            note = "custom filename — verify the suggested name"

        # Deduplicate: if name already seen, append index
        if name in seen_names:
            seen_names[name] += 1
            name = f"{name}_{seen_names[name]:02d}"
            note += " (deduplicated)"
        else:
            seen_names[name] = 1

        entries.append(SampleEntry(
            original_file=fna_path.name,
            new_name=name,
            source_type=source,
            notes=note,
        ))

    # Write TSV
    output_map.parent.mkdir(parents=True, exist_ok=True)
    with output_map.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["original_file", "new_name", "source_type", "notes"])
        for e in entries:
            writer.writerow([e.original_file, e.new_name, e.source_type, e.notes])

    logger.info(
        "Sample map written to '%s'. REVIEW AND EDIT before running 'rename'.",
        output_map,
    )

    # Print summary to stdout
    ncbi_n = sum(1 for e in entries if e.source_type == "ncbi")
    asm_n  = sum(1 for e in entries if e.source_type == "assembler")
    cus_n  = sum(1 for e in entries if e.source_type == "custom")
    print(f"\nDetected {len(entries)} genome file(s):")
    print(f"  NCBI accession : {ncbi_n}")
    print(f"  Assembler output: {asm_n}")
    print(f"  Custom          : {cus_n}")
    print(f"\nSample map draft → '{output_map}'")
    print("  ⚠  Review and edit 'new_name' column before running 'rename'.")
    print("  ⚠  Names must be unique, ≤12 chars, alphanumeric + underscore.")
    print("  ⚠  These names will become Bakta locus tag prefixes.")

    return entries


# ---------------------------------------------------------------------------
# rename subcommand
# ---------------------------------------------------------------------------

def _load_sample_map(map_path: Path) -> List[SampleEntry]:
    """Load and validate a sample_map.tsv file.

    Args:
        map_path: Path to the TSV generated by detect (and edited by user).

    Returns:
        List of validated SampleEntry objects.

    Raises:
        SampleMapError:     If file is missing or malformed.
        NamingConflictError: If new_name values are not unique.
    """
    if not map_path.exists():
        raise SampleMapError(f"Sample map not found: '{map_path}'")

    entries: List[SampleEntry] = []
    try:
        with map_path.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"original_file", "new_name"}
            if not required.issubset(set(reader.fieldnames or [])):
                raise SampleMapError(
                    f"Sample map must have columns: {required}. "
                    f"Found: {reader.fieldnames}"
                )
            for lineno, row in enumerate(reader, start=2):
                orig = row["original_file"].strip()
                name = row["new_name"].strip()
                if not orig or not name:
                    raise SampleMapError(
                        f"Line {lineno}: 'original_file' and 'new_name' cannot be empty."
                    )
                sanitized = _sanitize_name(name)
                if sanitized != name.upper():
                    logger.warning(
                        "Line %d: name '%s' will be sanitized to '%s'.",
                        lineno, name, sanitized,
                    )
                entries.append(SampleEntry(
                    original_file=orig,
                    new_name=sanitized,
                    source_type=row.get("source_type", "custom").strip(),
                    notes=row.get("notes", "").strip(),
                ))
    except (KeyError, csv.Error) as exc:
        raise SampleMapError(f"Cannot parse sample map: {exc}") from exc

    # Check uniqueness of new_name
    seen: Dict[str, str] = {}
    for e in entries:
        if e.new_name in seen:
            raise NamingConflictError(
                f"Duplicate new_name '{e.new_name}' for files "
                f"'{seen[e.new_name]}' and '{e.original_file}'. "
                "All new names must be unique."
            )
        seen[e.new_name] = e.original_file

    return entries


def _rewrite_fasta_headers(src: Path, dst: Path, new_name: str) -> int:
    """Copy a FASTA file rewriting headers to include the new sample name.

    New header format: >{new_name}_ctg{N}  (original description appended)
    This ensures contig IDs are traceable to their sample even when files
    are merged downstream.

    Args:
        src:      Source FASTA path.
        dst:      Destination path.
        new_name: Clean sample name.

    Returns:
        Number of records written.
    """
    records = list(SeqIO.parse(src, "fasta"))
    for i, rec in enumerate(records, start=1):
        original_desc = rec.description
        rec.id = f"{new_name}_ctg{i:04d}"
        rec.description = f"{rec.id} | {original_desc}"
        rec.name = rec.id
    SeqIO.write(records, dst, "fasta")
    return len(records)


def cmd_rename(
    input_dir: Path,
    map_path: Path,
    output_dir: Path,
    rewrite_headers: bool = True,
) -> None:
    """Rename genome files according to sample_map.tsv.

    Non-destructive: always writes to a new output directory.
    Optionally rewrites FASTA headers to embed the sample name.

    Args:
        input_dir:       Directory containing original genome files.
        map_path:        Validated sample_map.tsv.
        output_dir:      Destination directory for renamed files.
        rewrite_headers: If True, rewrite FASTA headers with new_name prefix.

    Raises:
        InputDirectoryError: If input_dir is missing.
        SampleMapError:      If map is malformed.
        NamingConflictError: If duplicate names exist.
    """
    if not input_dir.is_dir():
        raise InputDirectoryError(f"Not a directory: '{input_dir}'")

    entries = _load_sample_map(map_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    renamed = 0
    skipped = 0
    warnings: List[str] = []

    for entry in entries:
        src = input_dir / entry.original_file
        if not src.exists():
            warnings.append(
                f"'{entry.original_file}' not found in '{input_dir}' — skipped."
            )
            skipped += 1
            continue

        # Preserve original extension
        dst = output_dir / f"{entry.new_name}{src.suffix.lower()}"

        if dst.exists():
            warnings.append(
                f"Destination '{dst.name}' already exists — skipped "
                f"(delete output_dir to force re-run)."
            )
            skipped += 1
            continue

        if rewrite_headers and src.suffix.lower() in GENOME_EXTENSIONS:
            n_records = _rewrite_fasta_headers(src, dst, entry.new_name)
            logger.debug(
                "Renamed + rewrote headers: '%s' → '%s' (%d contigs).",
                src.name, dst.name, n_records,
            )
        else:
            shutil.copy2(src, dst)
            logger.debug("Copied: '%s' → '%s'.", src.name, dst.name)

        renamed += 1

    if warnings:
        for w in warnings:
            logger.warning(w)

    print(f"\nRename complete:")
    print(f"  Renamed : {renamed}")
    print(f"  Skipped : {skipped}")
    print(f"  Output  : '{output_dir}'")
    if renamed > 0:
        print(
            f"\nNext step: run Bakta on each file using its filename stem as --locus-tag.\n"
            f"Example:\n"
            f"  for f in {output_dir}/*.fna; do\n"
            f"    name=$(basename \"$f\" .fna)\n"
            f"    bakta --db /path/to/bakta_db --prefix \"$name\" "
            f"--locus-tag \"$name\" --output bakta_out/\"$name\" \"$f\"\n"
            f"  done"
        )


# ---------------------------------------------------------------------------
# check subcommand
# ---------------------------------------------------------------------------

def cmd_check(bakta_dir: Path) -> List[CheckResult]:
    """Verify completeness and ID consistency of Bakta outputs.

    Scans a directory where each subdirectory is one sample's Bakta output.
    Also accepts a flat directory where all Bakta files share a common prefix.

    Args:
        bakta_dir: Root directory of Bakta outputs.

    Returns:
        List of CheckResult, one per sample.

    Raises:
        InputDirectoryError: If bakta_dir does not exist.
    """
    if not bakta_dir.is_dir():
        raise InputDirectoryError(f"Not a directory: '{bakta_dir}'")

    # Detect layout: per-sample subdirs or flat
    subdirs = [d for d in bakta_dir.iterdir() if d.is_dir()]
    gff_files_flat = list(bakta_dir.glob("*.gff3")) + list(bakta_dir.glob("*.gff"))

    results: List[CheckResult] = []

    if subdirs:
        logger.info(
            "Detected per-sample subdirectory layout (%d subdirs).", len(subdirs)
        )
        for sample_dir in sorted(subdirs):
            result = _check_sample_dir(sample_dir, sample_dir.name)
            results.append(result)
    elif gff_files_flat:
        logger.info(
            "Detected flat layout (%d GFF files in root).", len(gff_files_flat)
        )
        for gff in sorted(gff_files_flat):
            sample_name = gff.stem
            result = _check_sample_dir(bakta_dir, sample_name, prefix=sample_name)
            results.append(result)
    else:
        raise InputDirectoryError(
            f"No Bakta output found in '{bakta_dir}'. "
            "Expected either per-sample subdirectories or .gff3 files."
        )

    # Print report
    ok_count = sum(1 for r in results if r.ok)
    print(f"\nBakta output check: {ok_count}/{len(results)} samples OK\n")
    print(f"{'Sample':<25} {'GFF':>5} {'FAA':>5} {'IDs':>5} {'Status'}")
    print("─" * 55)
    for r in results:
        status = "✓ OK" if r.ok else "✗ ISSUES"
        id_ok  = "✓" if r.id_prefix_ok else "✗"
        print(
            f"{r.sample_name:<25} "
            f"{'✓' if r.gff_found else '✗':>5} "
            f"{'✓' if r.faa_found else '✗':>5} "
            f"{id_ok:>5}  {status}"
        )
        for w in r.warnings:
            print(f"  ⚠  {w}")

    if ok_count < len(results):
        print(
            f"\n{len(results) - ok_count} sample(s) need attention. "
            "Fix issues before running Panaroo."
        )
    else:
        print("\nAll samples look good. You can proceed with Panaroo.")

    return results


def _check_sample_dir(
    directory: Path,
    sample_name: str,
    prefix: Optional[str] = None,
) -> CheckResult:
    """Check one sample's Bakta output for completeness and ID consistency."""
    pfx = prefix or sample_name
    result = CheckResult(sample_name=sample_name, gff_found=False, faa_found=False)

    # Check for required files
    for suffix in BAKTA_EXPECTED_SUFFIXES:
        candidate = directory / f"{pfx}{suffix}"
        if not candidate.exists():
            if suffix in (".gff3",):
                result.gff_found = False
                result.missing_files.append(f"{pfx}.gff3")
            elif suffix in (".faa",):
                result.faa_found = False
                result.missing_files.append(f"{pfx}.faa")
            else:
                result.warnings.append(f"Optional file missing: {pfx}{suffix}")
        else:
            if suffix == ".gff3":
                result.gff_found = True
            elif suffix == ".faa":
                result.faa_found = True

    # Check ID prefix consistency in the FAA file
    faa_path = directory / f"{pfx}.faa"
    if faa_path.exists():
        try:
            first_rec = next(SeqIO.parse(faa_path, "fasta"))
            result.id_prefix_found = first_rec.id.rsplit("_", 1)[0]
            if result.id_prefix_found.upper() != pfx.upper():
                result.id_prefix_ok = False
                result.warnings.append(
                    f"ID prefix mismatch: expected '{pfx}', "
                    f"found '{result.id_prefix_found}' in {pfx}.faa. "
                    "Did Bakta run with --locus-tag set to the sample name?"
                )
        except StopIteration:
            result.warnings.append(f"{pfx}.faa is empty.")
        except Exception as exc:
            result.warnings.append(f"Could not read {pfx}.faa: {exc}")

    return result


# ---------------------------------------------------------------------------
# panaroo subcommand
# ---------------------------------------------------------------------------

def cmd_panaroo(
    bakta_dir: Path,
    output_script: Path,
    threads: int = 8,
    mode: str = "strict",
    panaroo_outdir: Optional[Path] = None,
) -> None:
    """Generate a ready-to-run Panaroo shell script.

    Collects all .gff3 files from bakta_dir (flat or subdirectory layout)
    and writes a shell script with the complete Panaroo command.

    Args:
        bakta_dir:      Root directory of Bakta outputs.
        output_script:  Path for the generated .sh file.
        threads:        Number of threads for Panaroo (default: 8).
        mode:           Panaroo QC mode: 'strict', 'moderate', 'sensitive'.
        panaroo_outdir: Panaroo output directory. Defaults to ./panaroo_output.

    Raises:
        InputDirectoryError: If no GFF3 files are found.
    """
    if not bakta_dir.is_dir():
        raise InputDirectoryError(f"Not a directory: '{bakta_dir}'")

    # Collect GFF files from subdirs or flat
    gff_files: List[Path] = []
    for subdir in sorted(bakta_dir.iterdir()):
        if subdir.is_dir():
            gff_files.extend(sorted(subdir.glob("*.gff3")))
    if not gff_files:
        gff_files = sorted(bakta_dir.glob("*.gff3"))
    if not gff_files:
        raise InputDirectoryError(
            f"No .gff3 files found in '{bakta_dir}' or its subdirectories."
        )

    outdir = panaroo_outdir or Path("panaroo_output")
    sample_names = [f.stem for f in gff_files]

    script_lines = [
        "#!/usr/bin/env bash",
        f"# Generated by preparo v{__version__}",
        f"# {len(gff_files)} sample(s) detected",
        "#",
        "# Review parameters before running:",
        f"#   --clean-mode {mode}   (strict | moderate | sensitive)",
        f"#   -t {threads}              (threads)",
        "#",
        "# Activate your conda environment first:",
        "#   conda activate panaroo",
        "",
        "set -euo pipefail",
        "",
        f"mkdir -p {outdir}",
        "",
        "panaroo \\",
        f"    --clean-mode {mode} \\",
        f"    -t {threads} \\",
        f"    -o {outdir} \\",
        "    --remove-invalid-genes \\",
    ]

    # Add each GFF file as input
    for i, gff in enumerate(gff_files):
        sep = " \\" if i < len(gff_files) - 1 else ""
        script_lines.append(f"    {gff}{sep}")

    script_lines += [
        "",
        "echo 'Panaroo finished. Output in: " + str(outdir) + "'",
        "",
        "# Expected outputs for bioentro-suite:",
        f"#   {outdir}/gene_presence_absence.csv  ← pangenome matrix",
        f"#   {outdir}/pan_genome_reference.fa     ← representative sequences",
        f"#   {outdir}/core_gene_alignment.aln      ← core genome alignment",
    ]

    output_script.parent.mkdir(parents=True, exist_ok=True)
    output_script.write_text("\n".join(script_lines) + "\n", encoding="utf-8")
    output_script.chmod(0o755)   # make executable

    print(f"\nParanoo command script written to '{output_script}'")
    print(f"  Samples   : {len(gff_files)}")
    print(f"  Mode      : {mode}")
    print(f"  Threads   : {threads}")
    print(f"  Output dir: {outdir}")
    print(f"\nSample names that will appear in Panaroo matrix:")
    for name in sample_names:
        print(f"  {name}")
    print(f"\nRun with:  bash {output_script}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="preparo",
        description=(
            f"preparo v{__version__} — Genome Preparation Utility\n\n"
            "Prepares genome files for Bakta → Panaroo → bioentro-suite pipeline.\n"
            "Handles NCBI downloads, assembler outputs, and custom genomes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Typical workflow:\n"
            "  1. preparo detect  -i genomes_dir/ -o sample_map.tsv\n"
            "  2. Edit sample_map.tsv  (set meaningful new_name values)\n"
            "  3. preparo rename  -i genomes_dir/ -m sample_map.tsv -o renamed/\n"
            "  4. Run Bakta on each file in renamed/ (see rename output for example)\n"
            "  5. preparo check   -i bakta_outputs/\n"
            "  6. preparo panaroo -i bakta_outputs/ -o run_panaroo.sh\n"
            "  7. bash run_panaroo.sh\n\n"
            "Run 'preparo <subcommand> --help' for subcommand-specific options."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )

    sub = parser.add_subparsers(dest="subcommand", metavar="SUBCOMMAND")
    sub.required = True

    def _add_shared(p: argparse.ArgumentParser) -> None:
        p.add_argument("-v", "--verbose", action="store_true",
                       help="Enable debug-level logging")

    # ── detect ────────────────────────────────────────────────────────────────
    det = sub.add_parser(
        "detect",
        help="Scan genomes directory and generate sample_map.tsv draft",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  preparo detect -i genomes/ -o sample_map.tsv\n"
            "  preparo detect -i ncbi_downloads/ -o map.tsv\n"
        ),
    )
    _add_shared(det)
    det.add_argument("-i", "--input", required=True,
                     help="Directory containing genome FASTA files")
    det.add_argument("-o", "--output", default="sample_map.tsv",
                     help="Output sample map TSV (default: sample_map.tsv)")

    # ── rename ────────────────────────────────────────────────────────────────
    ren = sub.add_parser(
        "rename",
        help="Rename genome files using sample_map.tsv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  preparo rename -i genomes/ -m sample_map.tsv -o renamed/\n"
            "  preparo rename -i genomes/ -m sample_map.tsv -o renamed/ "
            "--no-rewrite-headers\n"
        ),
    )
    _add_shared(ren)
    ren.add_argument("-i", "--input", required=True,
                     help="Directory containing original genome files")
    ren.add_argument("-m", "--map", required=True,
                     help="sample_map.tsv (from 'detect', edited by user)")
    ren.add_argument("-o", "--output", required=True,
                     help="Output directory for renamed files")
    ren.add_argument("--no-rewrite-headers", action="store_true",
                     help="Copy files without rewriting FASTA headers")

    # ── check ─────────────────────────────────────────────────────────────────
    chk = sub.add_parser(
        "check",
        help="Verify Bakta output completeness and ID consistency",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  preparo check -i bakta_outputs/\n"
            "  preparo check -i bakta_outputs/ -v  # show debug detail\n"
        ),
    )
    _add_shared(chk)
    chk.add_argument("-i", "--input", required=True,
                     help="Directory with Bakta outputs (subdirs or flat)")

    # ── panaroo ───────────────────────────────────────────────────────────────
    pan = sub.add_parser(
        "panaroo",
        help="Generate a ready-to-run Panaroo shell script",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  preparo panaroo -i bakta_outputs/ -o run_panaroo.sh\n"
            "  preparo panaroo -i bakta_outputs/ -o run.sh -t 16 --mode moderate\n"
        ),
    )
    _add_shared(pan)
    pan.add_argument("-i", "--input", required=True,
                     help="Directory with Bakta GFF3 outputs")
    pan.add_argument("-o", "--output", default="run_panaroo.sh",
                     help="Output shell script path (default: run_panaroo.sh)")
    pan.add_argument("-t", "--threads", type=int, default=8,
                     help="Number of Panaroo threads (default: 8)")
    pan.add_argument("--mode", default="strict",
                     choices=["strict", "moderate", "sensitive"],
                     help="Panaroo QC mode (default: strict)")
    pan.add_argument("--outdir", default=None,
                     help="Panaroo output directory (default: panaroo_output/)")

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point.

    Returns:
        0 on success, 1 on handled error, 130 on keyboard interrupt.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _configure_logging(verbose=args.verbose)

    try:
        if args.subcommand == "detect":
            cmd_detect(
                input_dir=Path(args.input),
                output_map=Path(args.output),
            )

        elif args.subcommand == "rename":
            cmd_rename(
                input_dir=Path(args.input),
                map_path=Path(args.map),
                output_dir=Path(args.output),
                rewrite_headers=not args.no_rewrite_headers,
            )

        elif args.subcommand == "check":
            results = cmd_check(bakta_dir=Path(args.input))
            # Exit code 1 if any sample failed
            if any(not r.ok for r in results):
                return 1

        elif args.subcommand == "panaroo":
            cmd_panaroo(
                bakta_dir=Path(args.input),
                output_script=Path(args.output),
                threads=args.threads,
                mode=args.mode,
                panaroo_outdir=Path(args.outdir) if args.outdir else None,
            )

        return 0

    except PreparoError as exc:
        logger.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logger.warning("Interrupted by user.")
        return 130


if __name__ == "__main__":
    sys.exit(main())