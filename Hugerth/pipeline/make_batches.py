
import os
import argparse
import pandas as pd

SRC = "/insomnia001/depts/pmg/KoremLab/Projects/vagDB/public_samples/assemblies"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True)
    parser.add_argument("--column", default="file_name")
    parser.add_argument("--batch_dir", required=True)
    parser.add_argument("--n_batches", type=int, default=58)
    parser.add_argument("--missing_log", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.csv, encoding="utf-8-sig")
    prefixes = [str(p).strip() for p in df[args.column] if str(p).strip() != ""]

    found = []
    missing = []
    for p in prefixes:
        path = os.path.join(SRC, f"{p}.contigs.fasta.gz")
        if os.path.isfile(path):
            found.append((p, os.path.getsize(path)))
        else:
            missing.append(p)

    with open(args.missing_log, "w") as f:
        for p in missing:
            f.write(p + "\n")

    found.sort(key=lambda x: x[1], reverse=True)
    batches = [[] for _ in range(args.n_batches)]
    for i, (p, size) in enumerate(found):
        batches[i % args.n_batches].append(p)

    os.makedirs(args.batch_dir, exist_ok=True)
    for i, batch in enumerate(batches):
        with open(os.path.join(args.batch_dir, f"batch_{i:03d}.txt"), "w") as f:
            for p in batch:
                f.write(p + "\n")

    print(f"found {len(found)}, missing {len(missing)}, wrote {args.n_batches} batches")

if __name__ == "__main__":
    main()