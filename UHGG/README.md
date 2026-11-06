# UHGG Amino Acid Profiling

Computes the amino acid composition of every genome in the UHGG v2.0.2 catalog from its coding sequences (CDS), then examines how composition varies with GC content and taxonomy.

## Structure

```
uhgg/
├── README.md
├── pipeline/
│   ├── execute.sh
│   ├── process_genomes.py
│   └── merge_batches.py
├── prep/
│   ├── rep_species_prep.py
│   └── top_species_prep.py
├── analysis/
│   ├── rep_species_clustering.qmd
│   └── top_species_clustering.qmd
└── reports/
    ├── rep_species_clustering.pdf
    └── top_species_clustering.pdf
```

## Pipeline

This runs as a Slurm array job (97 tasks). Each task calls `process_genomes.py` on one batch of genome IDs. For each genome:
- parses the `.gff.gz` annotation
- extracts each CDS
- translates up to the first stop codon, skipping codons that contain `N`
- sums amino acid counts across all CDS

The result is one row per genome with the 20 amino acid counts and the cumulative CDS length.

Next, UHGG metadata is merged (GC content, genome type, lineage, completeness, contamination) and per-genome percentages are calculated.
- `rep_species_prep.py` keeps one representative genome per species.
- `top_species_prep.py` keeps every genome from the three species with the most genomes.

## Analysis

- `rep_species_clustering.qmd`: GC content vs. amino acid composition, PCA, and t-SNE across species representatives, colored by phylum.
- `top_species_clustering.qmd`: the same plots within each of the top three species, colored by genome type (Isolate vs. MAG).

## Inputs (not included but available on Insomnia)

- UHGG genomes
- UHGG metadata
- Batch lists

Outputs are written to `/insomnia001/depts/pmg/users/gw2508/Projects/get_profiles/`.
