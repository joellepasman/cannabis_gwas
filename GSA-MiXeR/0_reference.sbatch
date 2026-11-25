#!/bin/bash
#SBATCH --job-name=gsamixer_reference
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=20gb
#SBATCH --time=24:00:00
#SBATCH --output=gsamixer_reference.log

cd /vast/scratch/users/gerring.z/gsamixer

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"

for chri in {1..22}; do apptainer exec -B /vast gsa-mixer_latest.sif python /tools/mixer/precimed/mixer.py ld --bfile 1000G.EUR.QC.$chri --r2min 0.05 --ldscore-r2min 0.01 --out 1000G.EUR.QC.$chri.ld --ld-window-kb 10000; done

${MIXER_PY} snps --bim-file ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --ld-file ./reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.ld --chr2use 21-22 --r2 0.6 --maf 0.05 --subset 20000 --out 1000G.EUR.QC.@.snps --seed 123



