#!/bin/bash
#SBATCH --job-name=mbat_MA_CANN_FACTOR_250919
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH --output=mbat_MA_CANN_FACTOR_250919.log

cd /vast/scratch/users/gerring.z/icc4/scripts/gbt/gcta-1.94.1-linux-kernel-3-x86_64

./gcta64 --bfile /stornext/Bioinf/data/lab_bahlo/users/gerring.z/smr/ld/1000G.EUR.merged --mBAT-combo /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANN_FACTOR_250919.ma --mBAT-gene-list /vast/scratch/users/gerring.z/icc4/scripts/gbt/gcta-1.94.1-linux-kernel-3-x86_64/mBAT_cpp_manual_sample/glist_ensgid_hg19_v40.txt --out /vast/scratch/users/gerring.z/icc4/results/gbt/mbat/mbat_MA_CANN_FACTOR_250919 --thread-num 10
