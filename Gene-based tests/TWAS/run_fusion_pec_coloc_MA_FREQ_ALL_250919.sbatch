#!/bin/bash
#SBATCH --job-name=MA_FREQ_ALL_250919
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --array=1-22
#SBATCH --output=MA_FREQ_ALL_250919-%A_%a.log

cd /home/users/allstaff/gerring.z/projects/icc4/scripts/fusion_twas-master

module load R/4.4.1

Rscript FUSION.assoc_test.R --sumstats /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919.sumstats.gz --weights /stornext/Bioinf/data/lab_bahlo/users/gerring.z/twas/weights/PEC.BRAIN.RNASEQ.SUBSET_N.pos --weights_dir /stornext/Bioinf/data/lab_bahlo/users/gerring.z/twas/weights/ --ref_ld_chr /stornext/Bioinf/data/lab_bahlo/users/gerring.z/twas/ld/1000G_EUR_Phase3_plink/1000G.EUR.QC. --chr $SLURM_ARRAY_TASK_ID --out /vast/scratch/users/gerring.z/icc4/results/gbt/twas/twas_pec_coloc_MA_FREQ_ALL_250919_chr$SLURM_ARRAY_TASK_ID.dat  --coloc_P 0.999999 --GWASN 462303 --PANELN weight_n_panel.txt 
