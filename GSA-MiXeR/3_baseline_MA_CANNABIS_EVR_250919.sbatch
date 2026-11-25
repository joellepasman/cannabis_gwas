#!/bin/bash
#SBATCH --job-name=gsamixer_MA_CANNABIS_EVR_250919_step3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=25gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/gsamixer_MA_CANNABIS_EVR_250919_step3.log

cd /vast/scratch/users/gerring.z/icc4/scripts/genesets/gsamixer 

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} plsa --gsa-base \
  --trait1-file /vast/scratch/users/gerring.z/icc4/data/gwas/gsamixer/MA_CANNABIS_EVR_250919.chr@.sumstats.gz \
  --use-complete-tag-indices \
  --bim-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --loadlib-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
  --annot-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz \
  --go-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/gsa-mixer-baseline-annot_10mar2023.txt \
  --seed 1000 --exclude-ranges MHC --hardprune-r2 0.6 --threads 8 \
  --out /vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANNABIS_EVR_250919_base


