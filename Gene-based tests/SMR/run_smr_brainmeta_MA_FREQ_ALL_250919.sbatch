#!/bin/bash
#SBATCH --job-name=MA_FREQ_ALL_250919
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=02:00:00
#SBATCH --array=1-22
#SBATCH --output=MA_FREQ_ALL_250919-%A_%a.out

cd /home/users/allstaff/gerring.z/projects/icc4/scripts/smr/smr-1.3.1-linux-x86_64

./smr --bfile /stornext/Bioinf/data/lab_bahlo/users/gerring.z/twas/ld/LDREF/1000G.EUR.$SLURM_ARRAY_TASK_ID --gwas-summary /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919.ma --beqtl-summary /stornext/Bioinf/data/lab_bahlo/users/gerring.z/smr/BrainMeta_cis_eqtl_summary/BrainMeta_cis_eQTL_chr$SLURM_ARRAY_TASK_ID --out /vast/scratch/users/gerring.z/icc4/results/gbt/smr/smr_brainmeta_MA_FREQ_ALL_250919_chr$SLURM_ARRAY_TASK_ID --thread-num 10 
