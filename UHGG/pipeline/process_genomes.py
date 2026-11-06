
# =======
# IMPORTS
# =======

import os
import gzip
import time
import argparse
from itertools import product

import numpy as np
import pandas as pd


# ==========================
# AMINO ACID & CODON MAPPING
# ==========================

# Alphabetic order of amino acids (for the output array)
amino_acids = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]

# Codon mapping dictionary
codon_table = {
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "AAT": "Asn", "AAC": "Asn",
    "GAT": "Asp", "GAC": "Asp",
    "TGT": "Cys", "TGC": "Cys",
    "CAA": "Gln", "CAG": "Gln",
    "GAA": "Glu", "GAG": "Glu",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "CAT": "His", "CAC": "His",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "AAA": "Lys", "AAG": "Lys",
    "ATG": "Met",
    "TTT": "Phe", "TTC": "Phe",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "TGG": "Trp",
    "TAT": "Tyr", "TAC": "Tyr",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TAA": "STOP", "TAG": "STOP", "TGA": "STOP"
}


# ================
# HELPER FUNCTIONS
# ================

# Creates an ordered amino acid / peptide list and index (based on k)
def create_index(k):
    peptides = list(product(amino_acids, repeat=k))
    index = {p: i for i, p in enumerate(peptides)} 
    columns = ["".join(p) for p in peptides]
    return index, columns

# Translation table for reversed coding sequences
complement = str.maketrans("ATCG", "TAGC")
def reverse_complement(seq):
    return seq.translate(complement)[::-1]

# Extracts the CDS from a contig
def get_sequence(seq, start, end, strand, cds_id):
    assert seq is not None, f"{cds_id}: no matching contig found"
    assert strand in ("+", "-"), f"{cds_id}: not a valid strand value '{strand}'"
    cds_seq = seq[start-1:end]
    assert len(cds_seq) == end - start + 1, f"{cds_id}: length does not match expected"
    if strand == "-":
        cds_seq = reverse_complement(cds_seq)
    return cds_seq

# Creates a list of codons for a given sequence
def get_codons(seq, phase, cds_id):
    seq = seq[phase:]
    assert len(seq) % 3 == 0, f"{cds_id}: length not divisible by 3"
    codons = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if "N" not in codon:
            codons.append(codon)
    return codons

# Maps the list of codons to their corresponding amino acids
def translate(codons):
    aa_seq = []
    for codon in codons:
        aa = codon_table.get(codon)
        if aa == "STOP":
            break
        aa_seq.append(aa)
    return aa_seq

# Produces counts of amino acids / peptides for each sequence
def count_peptides(amino_acids, k, aa_index, size):
    counts = np.zeros(size, dtype=np.int32)
    for i in range(len(amino_acids) - k + 1):
        peptide = tuple(amino_acids[i:i+k])
        if peptide in aa_index:
            counts[aa_index[peptide]] += 1
    return counts

# Parses the .gff file
def parse_gff(filepath):
    features = []
    sequences = {}
    current_seq = []
    in_fasta = False
    current_contig = None
    
    with gzip.open(filepath, "rt") as f:
        for line in f:
            line = line.strip()
            if in_fasta:
                if line.startswith(">"):
                    if current_contig:
                        sequences[current_contig] = "".join(current_seq)
                    current_contig = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
            elif line == "##FASTA":
                in_fasta = True
            elif line.startswith("#"):
                continue
            else:
                fields = line.split("\t")
                if fields[2] == "CDS":
                    contig_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    phase = int(fields[7]) if fields[7] != "." else 0
                    cds_id = f"{contig_id}_{start}_{end}"
                    features.append((contig_id, cds_id, start, end, strand, phase))
        if current_contig:
            sequences[current_contig] = "".join(current_seq)
    
    return features, sequences


# =================
# GENOME PROCESSING
# =================

# Returns dataframe with one row per CDS per genome
def process_sequences(genome_dir, genome_id, k=1):
    filepath = os.path.join(genome_dir, f"{genome_id}.gff.gz")
    features, sequences = parse_gff(filepath)
    index, columns = create_index(k)
    size = len(index)
    rows = []
    
    for contig_id, cds_id, start, end, strand, phase in features:
        contig_seq = sequences.get(contig_id)
        cds_seq = get_sequence(contig_seq, start, end, strand, cds_id)
        codons = get_codons(cds_seq, phase, cds_id)
        aa_seq = translate(codons)
        counts = count_peptides(aa_seq, k, index, size)
        cds_length = end - start + 1
        rows.append([genome_id, cds_id, cds_length] + counts.tolist())
    
    return pd.DataFrame(rows, columns=["genome","cds_id","cds_length"] + columns)

# Returns dataframe agregated by genome
def process_genome(genome_dir, genome_id, k=1):
    filepath = os.path.join(genome_dir, f"{genome_id}.gff.gz")
    features, sequences = parse_gff(filepath)
    index, columns = create_index(k)
    size = len(index)
    total_counts = np.zeros(size, dtype=np.int32)
    total_length = 0
    
    for contig_id, cds_id, start, end, strand, phase in features:
        contig_seq = sequences.get(contig_id)
        cds_seq = get_sequence(contig_seq, start, end, strand, cds_id)
        codons = get_codons(cds_seq, phase, cds_id)
        aa_seq = translate(codons)
        counts = count_peptides(aa_seq, k, index, size)
        total_counts += counts
        total_length += end - start + 1
    
    row = [genome_id, total_length] + total_counts.tolist()
    return pd.DataFrame([row], columns=["genome","cum_cds_length"] + columns)


# ================
# BATCH PROCESSING
# ================

genome_dir = "/insomnia001/depts/pmg/users/ic2465/Projects/MANU_combine/ComBinE/data/UHGG/genomes"

# Creates a list of genome names from a batch file
def read_batch(batch_file):
    with open(batch_file) as f:
        genome_ids = [line.strip() for line in f if line.strip() != ""]
    return genome_ids

# Parses parameters, runs read_batch, runs process_genome, and outputs as a .parquet file
def process_batch():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch_file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--k", type=int, default=1)
    args = parser.parse_args()

    genome_ids = read_batch(args.batch_file)
    dfs = []
    start_time = time.time()

    for genome_id in genome_ids:
        df = process_genome(genome_dir, genome_id, k=args.k)
        dfs.append(df)

    elapsed = time.time() - start_time
    result = pd.concat(dfs, ignore_index=True)
    result.to_parquet(args.output)

    print(f"Batch completed in {elapsed:.2f} seconds.")


# ==================
# RUNNING THE SCRIPT
# ==================

if __name__ == "__main__":
    process_batch()
