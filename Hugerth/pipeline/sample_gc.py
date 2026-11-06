import os
import gzip
import argparse
import pandas as pd

SRC = "/insomnia001/depts/pmg/KoremLab/Projects/vagDB/public_samples/assemblies"

def gc_of(path):
    gc = 0
    at = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                continue
            s = line.strip().upper()
            gc += s.count("G") + s.count("C")
            at += s.count("A") + s.count("T")
    total = gc + at
    return 100*gc/total if total > 0 else float("nan")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.samples) as f:
        ids = [line.strip() for line in f if line.strip()]

    rows = []
    for sid in ids:
        path = os.path.join(SRC, f"{sid}.contigs.fasta.gz")
        rows.append([sid, gc_of(path)])
        print(f"{sid} done")

    pd.DataFrame(rows, columns=["sample_id","GC_content"]).to_parquet(args.output, index=False)

if __name__ == "__main__":
    main()