#!/bin/bash
#SBATCH --job-name=ldsc_munge_sumstats
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=24:00:00
#SBATCH --output=munge_MA_FREQ_ALL_250919_%j.log

cd /vast/scratch/users/gerring.z/icc4/scripts/partitioned_ldsc 

module load miniconda3/latest
conda activate ldsc

python ./ldsc/munge_sumstats.py --sumstats /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919.ma --out /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919_pldsc --snp SNP --N-col n --a1 A1 --a2 A2 --frq freq --merge-alleles w_hm3.snplist --a1-inc --chunksize 500000 --signed-sumstats b,0
 
