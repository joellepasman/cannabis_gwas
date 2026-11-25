#!/bin/bash
#SBATCH --job-name=MA_CANNABIS_EVR_250919_step4b
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=36gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/MA_CANNABIS_EVR_250919_step4b_%j.log

cd /vast/scratch/users/gerring.z/icc4/scripts/genesets/gsamixer

module load apptainer

# Setup container and tools
export MIXER_SIF=gsa-mixer_latest.sif
export MIXER_PY="singularity exec -B /vast -B /stornext ${MIXER_SIF} python /tools/mixer/precimed/mixer.py"
export MAGMA="singularity exec -B /vast -B /stornext ${MIXER_SIF} magma"

# Shared paths
BFILE=/stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/1000G_EUR_Phase3_plink/merged/1000G.EUR.QC
SNP_LOC=${BFILE}.bim
GENE_LOC=/stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/magma-gene-annot_10mar2023.csv
GENESET_ANNOT=/stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/magma-geneset-annot_10mar2023.csv
GENE_ANNOT=/stornext/Bioinf/data/lab_bahlo/users/gerring.z/mixer/reference/magma_annotation.step1.genes.annot

# ensure log and result directories exist
mkdir -p logs
mkdir -p /vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer

# Step 1: SNP-to-gene annotation (only once)
#$MAGMA --snp-loc "$SNP_LOC" \
#       --gene-loc "$GENE_LOC" \
#       --out "$ANNOT_OUT" \
#       --annotate window=0

# step 2: Gene-level analysis
$MAGMA --pval "/vast/scratch/users/gerring.z/icc4/data/gwas/MA_CANNABIS_EVR_250919.ma" snp-id=SNP pval=P ncol=N \
       --bfile "$BFILE" \
       --gene-annot "$GENE_ANNOT" \
       --out "/vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANNABIS_EVR_250919_magma.step2"

# step 3: geneset-level analysis
$MAGMA --gene-results "/vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANNABIS_EVR_250919_magma.step2.genes.raw" \
       --set-annot "$GENESET_ANNOT" \
       --out "/vast/scratch/users/gerring.z/icc4/results/genesets/gsamixer/MA_CANNABIS_EVR_250919_magma"

