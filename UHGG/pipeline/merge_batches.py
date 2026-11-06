
import glob
import pandas as pd

files = sorted(glob.glob("output/batch_*.parquet"))
dfs = [pd.read_parquet(f) for f in files]
result = pd.concat(dfs, ignore_index=True)
result.to_parquet("/insomnia001/depts/pmg/users/gw2508/Projects/UHGG/final_aa_profiles.parquet")

print(result.shape)