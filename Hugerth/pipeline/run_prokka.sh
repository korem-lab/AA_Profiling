#!/bin/bash
#SBATCH --job-name=hugerth_prokka
#SBATCH --account=pmg
#SBATCH --array=0-57
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/batch_%a.out
#SBATCH --error=logs/batch_%a.err

set -u

source /insomnia001/shared/apps/anaconda/2023.09/etc/profile.d/conda.sh

PROJ=/insomnia001/depts/pmg/users/gw2508/Projects/Hugerth
SRC=/insomnia001/depts/pmg/KoremLab/Projects/vagDB/public_samples/assemblies
TASK=$(printf "%03d" ${SLURM_ARRAY_TASK_ID})
BATCH=${PROJ}/batches/batch_${TASK}.txt
WORK=${PROJ}/tmp/task_${TASK}

mkdir -p ${WORK}

while read SAMPLE; do
    [ -z "${SAMPLE}" ] && continue
    OUT=${PROJ}/parquets/${SAMPLE}.parquet
    if [ -f "${OUT}" ]; then
        echo "skipping ${SAMPLE}, parquet exists"
        continue
    fi

    zcat ${SRC}/${SAMPLE}.contigs.fasta.gz > ${WORK}/${SAMPLE}.fa
    if [ ! -s "${WORK}/${SAMPLE}.fa" ]; then
        echo "empty assembly ${SAMPLE}, skipping"
        rm -f ${WORK}/${SAMPLE}.fa
        continue
    fi

    conda activate prokka
    prokka --metagenome --cpus 4 --force \
        --outdir ${WORK}/prokka_${SAMPLE} \
        --prefix ${SAMPLE} \
        --locustag ${SAMPLE} \
        ${WORK}/${SAMPLE}.fa
    conda deactivate

    if [ -f "${WORK}/prokka_${SAMPLE}/${SAMPLE}.gff" ]; then
        conda activate bio_stats
        python ${PROJ}/cds_aa_counts.py \
            --gff ${WORK}/prokka_${SAMPLE}/${SAMPLE}.gff \
            --sample_id ${SAMPLE} \
            --output ${OUT}
        conda deactivate
    else
        echo "prokka produced no gff for ${SAMPLE}"
    fi

    rm -rf ${WORK}/prokka_${SAMPLE} ${WORK}/${SAMPLE}.fa
done < ${BATCH}

rmdir ${WORK} 2>/dev/null
echo "task ${TASK} finished"