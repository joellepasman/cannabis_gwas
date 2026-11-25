#!/bin/bash
#SBATCH --job-name=magma_MA_FREQ_ALL_250919
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --time=1:00:00
#SBATCH --output=magma_MA_FREQ_ALL_250919.log

cd /vast/scratch/users/gerring.z/icc4/scripts/gbt/magma_v110

module load magma-gpu/2.7.0-11.3.1

./magma --bfile /stornext/Bioinf/data/lab_bahlo/users/gerring.z/magma/g1000_eur --pval /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919.ma ncol=N --gene-annot /vast/scratch/users/gerring.z/icc4/data/magma/g1000_ensembl.genes.annot --out /vast/scratch/users/gerring.z/icc4/results/gbt/magma/magma_MA_FREQ_ALL_250919

