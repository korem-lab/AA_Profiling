import pandas as pd

parquet_path = "/insomnia001/depts/pmg/users/gw2508/Projects/UHGG/final_aa_profiles.parquet"
metadata_path = "/insomnia001/depts/pmg/users/ic2465/Projects/MANU_combine/ComBinE/data/UHGG/UHGGv2.0.2_metadata.csv"
output_path = "/insomnia001/depts/pmg/users/gw2508/Projects/UHGG/rep_species_profiles.parquet"

aa_columns = [
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val"
]

df = pd.read_parquet(parquet_path)
meta = pd.read_csv(metadata_path)
meta = meta.rename(columns={"Genome": "genome"})

meta_reps = meta[meta["genome"] == meta["Species_rep"]].copy()
meta_reps = meta_reps[["genome", "GC_content", "Genome_type", "Lineage", "Completeness", "Contamination"]]
df = df.merge(meta_reps, on="genome", how="inner")

df["total_aa"] = df[aa_columns].sum(axis=1)
for aa in aa_columns:
    df[aa + "_pct"] = df[aa] / df["total_aa"]

df.to_parquet(output_path, index=False)
print(f"Saved {len(df)} species reps to {output_path}")