
import pandas as pd

parquet_path = "/insomnia001/depts/pmg/users/gw2508/Projects/UHGG/final_aa_profiles.parquet"
metadata_path = "/insomnia001/depts/pmg/users/ic2465/Projects/MANU_combine/ComBinE/data/UHGG/UHGGv2.0.2_metadata.csv"
output_path = "/insomnia001/depts/pmg/users/gw2508/Projects/UHGG/top_species_profiles.parquet"

aa_columns = [
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val"
]

df = pd.read_parquet(parquet_path)
meta = pd.read_csv(metadata_path)
meta = meta.rename(columns={"Genome": "genome"})

meta = meta[["genome", "Species_rep", "GC_content", "Genome_type", "Lineage", "Completeness", "Contamination"]]
df = df.merge(meta, on="genome", how="inner")

top_species = df["Species_rep"].value_counts().head(3).index
df = df[df["Species_rep"].isin(top_species)].copy()

df["total_aa"] = df[aa_columns].sum(axis=1)
for aa in aa_columns:
    df[aa + "_pct"] = df[aa] / df["total_aa"]

df.to_parquet(output_path, index=False)
print(df["Species_rep"].value_counts().loc[top_species])
print(f"Saved {len(df)} genomes from top 3 species to {output_path}")