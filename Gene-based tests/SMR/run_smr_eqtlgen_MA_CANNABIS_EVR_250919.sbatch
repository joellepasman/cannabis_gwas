#!/bin/bash
#SBATCH --job-name=MA_CANNABIS_EVR_250919
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=8gb
#SBATCH --time=05:00:00
#SBATCH --output=MA_CANNABIS_EVR_250919.log

cd /home/users/allstaff/gerring.z/projects/icc4/scripts/smr/smr-1.3.1-linux-x86_64

./smr --bfile /stornext/Bioinf/data/lab_bahlo/users/gerring.z/smr/ld/1000G.EUR.merged --gwas-summary /vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANNABIS_EVR_250919.ma --beqtl-summary /stornext/Bioinf/data/lab_bahlo/users/gerring.z/smr/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --out /vast/scratch/users/gerring.z/icc4/results/gbt/smr/smr_eqtlgen_MA_CANNABIS_EVR_250919 --thread-num 10 
