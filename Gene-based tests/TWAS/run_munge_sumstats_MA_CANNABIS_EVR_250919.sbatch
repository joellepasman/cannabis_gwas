#!/bin/bash
#SBATCH --job-name=ldsc_munge_sumstats
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=munge_%j.log

cd /home/users/allstaff/gerring.z/projects/icc4/scripts/fusion_twas-master/ldsc

pwd; hostname; date

module load miniconda3/latest
conda activate ldsc

python munge_sumstats.py --sumstats /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANNABIS_EVR_250919.ma --out /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANNABIS_EVR_250919 --snp SNP --N-col n --a1 A1 --a2 A2 --frq freq --signed-sumstats b,0
 
