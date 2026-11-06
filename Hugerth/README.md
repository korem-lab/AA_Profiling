# Hugerth Vaginal Metagenome Profiling

Annotates the Hugerth et al. (2023) vaginal metagenome assemblies with Prokka, computes the amino acid composition of each sample from its coding sequences (CDS), and tests whether composition differs across CSTs.

## Structure

```
hugerth/
├── README.md
├── pipeline/
│   ├── make_batches.py
│   ├── run_prokka.sh
│   ├── cds_aa_counts.py
│   ├── merge_parquets.py
│   └── sample_gc.py
├── prep/
│   └── make_sample_profiles.py
├── analysis/
│   ├── hugerth_analysis.qmd
│   └── single_day_analysis.qmd
└── reports/
    ├── hugerth_analysis.pdf
    └── single_day_analysis.pdf
```

## Pipeline

`make_batches.py` checks which samples have an assembly, logs any that are missing, and splits the rest into 58 batches balanced by file size.

This runs as a Slurm array job (58 tasks). Each task calls `run_prokka.sh` on one batch of samples. For each sample:
- annotates the assembly with `prokka --metagenome`
- extracts each CDS from the resulting GFF
- translates up to the first stop codon, skipping ambiguous codons
- writes one row per CDS with the 20 amino acid counts

Next, the per-sample outputs are merged with `merge_parquets.py`, and `sample_gc.py` computes assembly-wide GC content for each sample.
- `make_sample_profiles.py` sums CDS counts to one row per sample, then adds CDS count, total amino acids, per-sample percentages, and GC content.

## Analysis

- `hugerth_analysis.qmd`: amino acid distributions, GC content vs. amino acid composition, PCA, and t-SNE across all samples. PERMANOVA and PERMDISP are then run.
- `single_day_analysis.qmd`: the same core analysis restricted to the most-sampled visit, so no subject appears more than once.

## Inputs (not included but available on Insomnia)

- Hugerth assemblies
- Sample metadata
- Ravel CST assignments

Outputs are written to `/insomnia001/depts/pmg/users/gw2508/Projects/Hugerth/`.
