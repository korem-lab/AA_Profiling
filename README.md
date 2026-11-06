# AA_Profiling

Amino acid composition profiling of prokaryotic genomes and vaginal metagenomes, done with the Korem Lab (PMG, Columbia University). Thank you to Izaak Coleman for his guidance on this project.

## Projects

| Folder | Description |
|---|---|
| [uhgg/](uhgg/) | Amino acid profiles by genome across the UHGG catalog, with dimensionality reduction by phylum and within the most-sampled species. |
| [hugerth/](hugerth/) | Amino acid profiles by sample from the Hugerth et al. (2023) vaginal metagenomes, annotated with Prokka and compared across community state types (CSTs). |

Each project folder has its own README covering the pipeline, run order, and required inputs.

## Layout

Both projects follow the same structure:

- pipeline/: Slurm jobs and scripts run on the Insomnia HPC cluster
- prep/: scripts that build analysis-ready tables from the pipeline output
- analysis/: Quarto notebooks
- reports/: Notebooks rendered to PDF

## Data

No data is included in this repository. Input genomes, assemblies, metadata, and output files (.parquet) are on Insomnia.
