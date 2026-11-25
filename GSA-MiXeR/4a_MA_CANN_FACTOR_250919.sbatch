#!/bin/bash
#SBATCH --job-name=gsamixer_MA_CANN_FACTOR_250919_step4a
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/gsamixer_MA_CANN_FACTOR_250919_step4a_%a.log

cd /vast/scratch/users/gerring.z/MA_CANN_FACTOR_250919/scripts/genesets/gsamixer

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

${MIXER_PY} plsa --gsa-full \
  --trait1-file /vast/scratch/users/gerring.z/icc4/data/gwas/gsamixer/MA_CANN_FACTOR_250919.chr@.sumstats.gz \
  --out /vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANN_FACTOR_250919_full \
  --bim-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
  --use-complete-tag-indices \
  --loadlib-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
  --go-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/gsa-mixer-gene-annot_10mar2023.csv \
  --go-file-test /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/gsa-mixer-hybridLOO-annot_10mar2023.csv \
  --annot-file /stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/baseline_v2.2_1000G.EUR.QC.@.annot.gz \
  --load-params-file /vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANN_FACTOR_250919_base.json \
  --seed 1000 --exclude-ranges MHC --hardprune-r2 0.6 --threads 8

