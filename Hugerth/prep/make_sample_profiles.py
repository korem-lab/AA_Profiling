import pandas as pd

aa_columns = ["Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile",
    "Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val"]
pct_columns = [aa + "_pct" for aa in aa_columns]

df = pd.read_parquet("hugerth_cds_aa_counts.parquet")

samples = df.groupby("sample_id")[aa_columns].sum()
samples["n_cds"] = df.groupby("sample_id").size()
samples["total_aa"] = samples[aa_columns].sum(axis=1)
for aa in aa_columns:
    samples[aa + "_pct"] = 100*samples[aa]/samples["total_aa"]
samples = samples.reset_index()

gc = pd.read_parquet("sample_gc.parquet")
samples = samples.merge(gc, on="sample_id", how="left")

samples.to_parquet("hugerth_sample_profiles.parquet", index=False)
print(samples.shape)
print(samples["GC_content"].isna().sum())
print(samples[["sample_id","n_cds","total_aa","GC_content","Ala_pct","Lys_pct"]].head())