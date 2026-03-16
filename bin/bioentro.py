#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import logging
import math
import os
from collections import Counter
from Bio import SeqIO

# 1. LOGGING CONFIGURATION
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# 2. CALC FUNCTIONS
def get_metrics(sequence, k, is_protein=False):
    """ Calculate metrics adapting the alphabet (4 for DNA, 20 for Protein)."""
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

def process_dna_mode(fna_path, gff_path, kmer_size):
    """Process the complete genome (Nucleotides)."""
    try:
        records = list(SeqIO.parse(fna_path, "fasta"))
        raw_seq = "".join(str(r.seq).upper() for r in records)
        genome_seq = "".join(b for b in raw_seq if b in "ACGNT")
        
        size_mb = len(genome_seq) / 1e6
        h_real, kl_div, h_max = get_metrics(genome_seq, kmer_size, is_protein=False)
        efficiency = (h_real / h_max) * 100 if h_max > 0 else 0
        
        gc_count = genome_seq.count('G') + genome_seq.count('C')
        gc_pct = (gc_count / len(genome_seq)) * 100 if genome_seq else 0

        total_cds = 0
        with open(gff_path, 'r') as gff:
            for line in gff:
                if not line.startswith("#") and "\tCDS\t" in line:
                    total_cds += 1
        
        return {
            "size_mb": size_mb, "gc": gc_pct, "h_max": h_max, 
            "h_real": h_real, "kl_div": kl_div, "efficiency": efficiency,
            "total_cds": total_cds, "density": total_cds / size_mb if size_mb > 0 else 0
        }
    except Exception as e:
        logging.error(f"Error in DNA mode: {e}")
        sys.exit(1)

def process_protein_mode(faa_path, kmer_size):
    """Process proteins with a statistical reliability semaphore."""
    results = []
    # Threshold for reliability: at least half of the possible k-mers should be present in the sequence
    threshold = (20**kmer_size) / 2 
    
    try:
        for record in SeqIO.parse(faa_path, "fasta"):
            seq_len = len(record.seq)
            h_real, kl_div, h_max = get_metrics(str(record.seq).upper(), kmer_size, is_protein=True)
            efficiency = (h_real / h_max) * 100 if h_max > 0 else 0
            status = "Reliable" if seq_len >= threshold else "Low_Sample"
            
            results.append({
                "ID": record.id,
                "Len_AA": seq_len,
                "H_max": h_max,
                "H_real": h_real,
                "KL_Div": kl_div,
                "Efficiency": efficiency,
                "Status": status
            })
        return results
    except Exception as e:
        logging.error(f"Error in Protein Mode: {e}")
        sys.exit(1)

# 3. CLI ARGUMENTS
def parse_args():
    parser = argparse.ArgumentParser(description="bioentro v1.0.0: DNA and Proteins")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--fasta", help="fasta file with DNA sequences (.fna)")
    group.add_argument("-a", "--aminoacids", help="fasta file with protein sequences (.faa)")
    parser.add_argument("-g", "--gff", help="GFF annotation (required for DNA mode)")
    parser.add_argument("-k", "--kmer", type=int, default=2, help="K-mer size (usually DNA=6, AA=2)")
    parser.add_argument("-o", "--output", required=True, help="Output .tsv file")
    return parser.parse_args()

# 4. ACCESS POINT
def main():
    args = parse_args()
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.aminoacids:
        cepa_id = os.path.basename(args.aminoacids).split('.')[0]
        logging.info(f"Protein mode: Processing {cepa_id}")
        data = process_protein_mode(args.aminoacids, args.kmer)
        
        # WRITING FOR PROTEIN MODE
        with open(args.output, 'w') as out:
            # Header with all columns
            out.write("Prot_ID\tLen_AA\tH_max\tH_real\tKL_Div\tEfficiency\tStatus\n")
            for p in data:
                out.write(f"{p['ID']}\t{p['Len_AA']}\t{p['H_max']:.4f}\t{p['H_real']:.4f}\t"
                          f"{p['KL_Div']:.4f}\t{p['Efficiency']:.2f}\t{p['Status']}\n")
    else:
        if not args.gff:
            logging.error("The DNA mode requires a GFF file (-g).")
            sys.exit(1)
        cepa_id = os.path.basename(args.fasta).split('.')[0]
        logging.info(f"DNA mode: Processing {cepa_id}")
        res = process_dna_mode(args.fasta, args.gff, args.kmer)
        
        # WRITING FOR DNA MODE
        with open(args.output, 'w') as out:
            out.write("ID\tSize_Mb\tGC_Pct\tH_max\tH_real\tKL_Div\tEfficiency\tCDS\tDensity\n")
            out.write(f"{cepa_id}\t{res['size_mb']:.4f}\t{res['gc']:.2f}\t{res['h_max']:.4f}\t"
                      f"{res['h_real']:.4f}\t{res['kl_div']:.4f}\t{res['efficiency']:.2f}\t"
                      f"{res['total_cds']}\t{res['density']:.2f}\n")
    
    logging.info(f"Success. File generated at: {args.output}")

if __name__ == "__main__":
    main()