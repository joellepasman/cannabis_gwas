#!/bin/bash
#SBATCH --job-name=gsamixer_step2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=24:00:00
#SBATCH --output=gsamixer_step2.log

cd /vast/scratch/users/gerring.z/gsamixer

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

for chri in {1..22}; do ${MIXER_PY} plsa \
      --bim-file ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --use-complete-tag-indices \
      --chr2use ${chri} \
      --savelib-file ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bin \
      --out ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@; done

