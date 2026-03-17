#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import logging
import math
import os
import zlib
import re
from collections import Counter
from Bio import SeqIO

# 1. LOGGING CONFIGURATION
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# 2. CORE FUNCTIONS
def calculate_shannon_kl(sequence, k, is_protein=True):
    """Calculates Shannon Entropy and KL Divergence."""
    alphabet_size = 20 if is_protein else 4
    h_max = math.log2(alphabet_size**k)
    p_uniform = 1 / (alphabet_size**k)
    
    if len(sequence) < k:
        return 0.0, 0.0, h_max
    
    invalid_char = 'X' if is_protein else 'N'
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1) 
             if invalid_char not in sequence[i:i+k]]
    
    if not kmers:
        return 0.0, 0.0, h_max
        
    counts = Counter(kmers)
    total = sum(counts.values())
    
    h_real = 0.0
    kl_divergence = 0.0
    for c in counts.values():
        p_i = c / total
        h_real -= p_i * math.log2(p_i)
        kl_divergence += p_i * math.log2(p_i / p_uniform)
        
    return h_real, kl_divergence, h_max

def calculate_kolmogorov(sequence):
    """Approximates Kolmogorov Complexity using zlib."""
    if not sequence:
        return 0.0
    try:
        bytes_seq = sequence.encode('utf-8')
        compressed = zlib.compress(bytes_seq, level=9)
        return round(len(compressed) / len(bytes_seq), 4)
    except Exception:
        return 0.0

def process_dna_mode(fna_path, gff_path, kmer_size):
    """Processes genome and extracts CLEAN Organism Name."""
    try:
        records = list(SeqIO.parse(fna_path, "fasta"))
        if not records:
            return None
            
        # Extracts only the clean name [organism=...]
        full_desc = records[0].description
        match = re.search(r'\[organism=([^\]]+)\]', full_desc)
        org_name = match.group(1) if match else full_desc.split(' ')[0]
        
        raw_seq = "".join(str(r.seq).upper() for r in records)
        genome_seq = "".join(b for b in raw_seq if b in "ACGNT")
        
        size_mb = len(genome_seq) / 1e6
        h_real, kl_div, h_max = calculate_shannon_kl(genome_seq, kmer_size, is_protein=False)
        kolmo = calculate_kolmogorov(genome_seq)
        
        efficiency = (h_real / h_max) * 100 if h_max > 0 else 0
        gc_pct = ((genome_seq.count('G') + genome_seq.count('C')) / len(genome_seq)) * 100 if genome_seq else 0

        total_cds = 0
        with open(gff_path, 'r') as gff:
            for line in gff:
                if not line.startswith("#") and "\tCDS\t" in line:
                    total_cds += 1
        
        return {
            "organism": org_name, "size_mb": size_mb, "gc": gc_pct, "h_max": h_max, 
            "h_real": h_real, "kl_div": kl_div, "efficiency": efficiency,
            "kolmogorov": kolmo, "total_cds": total_cds, "density": total_cds / size_mb if size_mb > 0 else 0
        }
    except Exception as e:
        logging.error(f"Error en DNA mode: {e}"); sys.exit(1)

def process_protein_mode(faa_path, kmer_size):
    """Processes proteome and cleans up Protein Labels."""
    results = []
    threshold = (20**kmer_size) / 2 
    try:
        for record in SeqIO.parse(faa_path, "fasta"):
            seq_str = str(record.seq).upper()
            
            # Clean product name by removing ID and any bracketed info
            product = record.description.replace(record.id, "").strip()
            product = re.sub(r'\[[^\]]+\]', '', product).strip() # Deep cleaning
            if not product: product = "hypothetical protein"
            
            h_real, kl_div, h_max = calculate_shannon_kl(seq_str, kmer_size, is_protein=True)
            kolmo = calculate_kolmogorov(seq_str)
            
            results.append({
                "ID": record.id, "Product": product, "Len": len(seq_str), 
                "H_max": h_max, "H_real": h_real, "KL": kl_div, 
                "Eff": (h_real/h_max)*100 if h_max > 0 else 0,
                "Kolmo": kolmo, "Status": "Reliable" if len(seq_str) >= threshold else "Low_Sample"
            })
        return results
    except Exception as e:
        logging.error(f"Error en Protein Mode: {e}"); sys.exit(1)

# 3. CLI ARGUMENTS
def parse_args():
    parser = argparse.ArgumentParser(description="bioentro v1.1.0: Informational Complexity Suite")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fasta", help="DNA sequences (.fna)")
    group.add_argument("-a", "--aminoacids", help="Protein sequences (.faa)")
    parser.add_argument("-g", "--gff", help="GFF annotation")
    parser.add_argument("-k", "--kmer", type=int, default=None, help="K-mer size")
    parser.add_argument("-o", "--output", required=True, help="Output .tsv file")
    return parser.parse_args()

# 4. MAIN
def main():
    args = parse_args()
    k_eff = args.kmer if args.kmer is not None else (6 if args.fasta else 1)

    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.aminoacids:
        logging.info(f"Processing Proteins (k={k_eff})...")
        data = process_protein_mode(args.aminoacids, k_eff)
        with open(args.output, 'w') as out:
            out.write("Prot_ID\tProduct\tLen_AA\tH_max\tH_real\tKL_Div\tEfficiency\tKolmogorov\tStatus\n")
            for p in data:
                out.write(f"{p['ID']}\t{p['Product']}\t{p['Len']}\t{p['H_max']:.4f}\t{p['H_real']:.4f}\t"
                          f"{p['KL']:.4f}\t{p['Eff']:.2f}\t{p['Kolmo']:.4f}\t{p['Status']}\n")
    else:
        if not args.gff:
            logging.error("DNA mode requires a GFF file (-g)."); sys.exit(1)
        logging.info(f"Processing DNA (k={k_eff})...")
        res = process_dna_mode(args.fasta, args.gff, k_eff)
        with open(args.output, 'w') as out:
            out.write("Organism\tSize_Mb\tGC_Pct\tH_max\tH_real\tKL_Div\tEfficiency\tKolmogorov\tCDS\tDensity\n")
            out.write(f"{res['organism']}\t{res['size_mb']:.4f}\t{res['gc']:.2f}\t"
                      f"{res['h_max']:.4f}\t{res['h_real']:.4f}\t{res['kl_div']:.4f}\t"
                      f"{res['efficiency']:.2f}\t{res['kolmogorov']:.4f}\t{res['total_cds']}\t{res['density']:.2f}\n")
    
    logging.info(f"Success! Output: {args.output}")

if __name__ == "__main__":
    main()