#!/bin/bash
#SBATCH --job-name=gsamixer_step1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --time=5:00:00
#SBATCH --output=logs/gsamixer_step1.log

cd /vast/scratch/users/gerring.z/icc4/scripts/genesets/gsamixer

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

# Get the current file name based on array index
${MIXER_PY} split_sumstats \
  --trait1-file /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANNABIS_EVR_250919_gsamixer.sumstats.gz \
  --out /vast/scratch/users/gerring.z/icc4/data/gwas/gsamixer/MA_CANNABIS_EVR_250919.chr@.sumstats.gz

${MIXER_PY} split_sumstats \
  --trait1-file /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANN_FACTOR_250919_gsamixer.sumstats.gz \
  --out /vast/scratch/users/gerring.z/icc4/data/gwas/gsamixer/MA_CANN_FACTOR_250919.chr@.sumstats.gz

${MIXER_PY} split_sumstats \
  --trait1-file /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919_gsamixer.sumstats.gz \
  --out /vast/scratch/users/gerring.z/icc4/data/gwas/gsamixer/MA_FREQ_ALL_250919.chr@.sumstats.gz


