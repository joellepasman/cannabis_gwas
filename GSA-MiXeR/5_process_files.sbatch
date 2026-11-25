#!/bin/bash
#SBATCH --job-name=gsamixer_step5
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --time=01:00:00
#SBATCH --output=gsamixer_step5.log

cd /vast/scratch/users/gerring.z/icc4/scripts/genesets/gsamixer

module load apptainer

export MIXER_SIF=gsa-mixer_latest.sif
export PYTHON="singularity exec -B /vast -B /stornext $BIND ${MIXER_SIF} python"

$PYTHON process_gsa_mixer_output_FDR.py 

