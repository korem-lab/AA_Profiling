#!/bin/bash
#SBATCH --job-name=aa_profiling
#SBATCH --account=pmg
#SBATCH --array=0-96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --time=00:45:00
#SBATCH --output=logs/batch_%a.out
#SBATCH --error=logs/batch_%a.err

module load anaconda/2023.09
source /insomnia001/shared/apps/anaconda/2023.09/etc/profile.d/conda.sh
conda activate /insomnia001/depts/pmg/users/gw2508/envs/bio_stats

BATCH_ID=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

python process_genomes.py --batch_file batches/batch_${BATCH_ID}.txt --output output/batch_${BATCH_ID}.parquet