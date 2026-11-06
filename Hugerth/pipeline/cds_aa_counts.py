
import os
import argparse
import numpy as np
import pandas as pd

amino_acids = ["Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
    "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"]

codon_table = {
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg","AGA":"Arg","AGG":"Arg",
    "AAT":"Asn","AAC":"Asn",
    "GAT":"Asp","GAC":"Asp",
    "TGT":"Cys","TGC":"Cys",
    "CAA":"Gln","CAG":"Gln",
    "GAA":"Glu","GAG":"Glu",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
    "CAT":"His","CAC":"His",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile",
    "TTA":"Leu","TTG":"Leu","CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "AAA":"Lys","AAG":"Lys",
    "ATG":"Met",
    "TTT":"Phe","TTC":"Phe",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser","AGT":"Ser","AGC":"Ser",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "TGG":"Trp",
    "TAT":"Tyr","TAC":"Tyr",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TAA":"STOP","TAG":"STOP","TGA":"STOP"
}

aa_index = {aa: i for i, aa in enumerate(amino_acids)}
complement = str.maketrans("ATCG", "TAGC")

def reverse_complement(seq):
    return seq.translate(complement)[::-1]

def parse_attributes(field):
    attrs = {}
    for item in field.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = val
    return attrs

def parse_gff(filepath):
    features = []
    sequences = {}
    current_seq = []
    in_fasta = False
    current_contig = None

    with open(filepath) as f:
        for line in f:
            line = line.rstrip("\n")
            if in_fasta:
                if line.startswith(">"):
                    if current_contig:
                        sequences[current_contig] = "".join(current_seq)
                    current_contig = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            elif line == "##FASTA":
                in_fasta = True
            elif line.startswith("#") or line.strip() == "":
                continue
            else:
                fields = line.split("\t")
                if len(fields) < 9 or fields[2] != "CDS":
                    continue
                attrs = parse_attributes(fields[8])
                contig_id = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                phase = int(fields[7]) if fields[7] != "." else 0
                cds_id = attrs.get("locus_tag") or attrs.get("ID") or f"{contig_id}_{start}_{end}"
                features.append((contig_id, cds_id, start, end, strand, phase))
        if current_contig:
            sequences[current_contig] = "".join(current_seq)

    return features, sequences

def count_cds(seq, phase):
    seq = seq[phase:]
    counts = np.zeros(len(amino_acids), dtype=np.int32)
    skipped = 0
    n_full = len(seq) // 3
    if len(seq) % 3 != 0:
        skipped += 1 
    for i in range(n_full):
        codon = seq[3*i:3*i+3]
        aa = codon_table.get(codon)
        if aa is None:
            skipped += 1
            continue
        if aa == "STOP":
            break
        counts[aa_index[aa]] += 1
    return counts, skipped

def process_sample(gff_path, sample_id):
    features, sequences = parse_gff(gff_path)
    rows = []
    for contig_id, cds_id, start, end, strand, phase in features:
        contig_seq = sequences.get(contig_id)
        if contig_seq is None:
            continue
        cds_seq = contig_seq[start-1:end]
        if strand == "-":
            cds_seq = reverse_complement(cds_seq)
        counts, skipped = count_cds(cds_seq, phase)
        rows.append([sample_id, contig_id, cds_id, end-start+1] + counts.tolist() + [skipped])
    columns = ["sample_id","contig_id","cds_id","cds_length"] + amino_acids + ["skipped_codons"]
    return pd.DataFrame(rows, columns=columns)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--sample_id", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = process_sample(args.gff, args.sample_id)
    df.to_parquet(args.output, index=False)
    print(f"{args.sample_id}: {len(df)} CDS rows")

if __name__ == "__main__":
    main()