
import os
import glob
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet_dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.parquet_dir, "*.parquet")))
    dfs = [pd.read_parquet(f) for f in files]
    result = pd.concat(dfs, ignore_index=True)
    result.to_parquet(args.output, index=False)
    print(f"merged {len(files)} samples, {len(result)} rows")

if __name__ == "__main__":
    main()