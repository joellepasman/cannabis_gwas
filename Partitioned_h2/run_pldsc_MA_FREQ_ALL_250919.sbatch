#!/bin/bash
#SBATCH --job-name=sldsc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gerring.z@wehi.edu.au
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --time=12:00:00
#SBATCH --output=pldsc_MA_FREQ_ALL_250919_%j.log

cd /vast/scratch/users/gerring.z/icc4/scripts/partitioned_ldsc

module load miniconda3/latest
conda activate ldsc

declare -a celltypes=("Amygdala_excitatory" "Astrocyte" "Bergmann_glia" "Cerebellar_inhibitory" "CGE_interneuron" "Choroid_plexus" "Committed_oligodendrocyte_precursor" "Deep_layer_corticothalamic_and_6b" "Deep_layer_intratelencephalic" "Deep_layer_near_projecting" "Eccentric_medium_spiny_neuron" "Ependymal" "Fibroblast" "Hippocampal_CA1_3" "Hippocampal_CA4" "Hippocampal_dentate_gyrus" "LAMP5_LHX6_and_Chandelier" "Lower_rhombic_lip" "Mammillary_body" "Medium_spiny_neuron" "MGE_interneuron" "Microglia" "Midbrain_derived_inhibitory" "Miscellaneous" "Oligodendrocyte" "Oligodendrocyte_precursor" "Splatter" "Thalamic_excitatory" "Upper_layer_intratelencephalic" "Upper_rhombic_lip" "Vascular")

arraylength=${#celltypes[@]}

for (( i=0; i<${arraylength}; i++ ));
do
python ./ldsc/ldsc.py --h2 /vast/scratch/users/gerring.z/icc4/data/gwas/MA_FREQ_ALL_250919_pldsc.sumstats.gz --ref-ld-chr /stornext/Bioinf/data/lab_bahlo/users/gerring.z/ldsc/baseline_v1.2/baseline.,/stornext/Bioinf/data/lab_bahlo/users/gerring.z/ldsc/superclusters/${celltypes[$i]}/baseline. --w-ld-chr /stornext/Bioinf/data/lab_bahlo/users/gerring.z/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --overlap-annot --frqfile-chr /stornext/Bioinf/data/lab_bahlo/users/gerring.z/ldsc/1000G_Phase3_frq/1000G.EUR.QC. --print-coefficients --out /vast/scratch/users/gerring.z/icc4/results/partitioned_heritability/MA_FREQ_ALL_250919_${celltypes[$i]}
done
