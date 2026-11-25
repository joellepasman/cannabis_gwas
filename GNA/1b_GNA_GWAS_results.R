##################################################
# Cannabis GWAS                                  #
# - Get results                                  #
##################################################

require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)

working_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GNA"


#=============================#
# CHECK AND COMBINE RESULTS   #
#=============================#

# combine all batches
setwd(file.path(working_dir,"SNPeffects"))

system(paste0("
head -n1 GNA_gwas1_1_* > networkGWAS1_combined.csv
head -n1 GNA_gwas2_1_* > networkGWAS2_combined.csv
for file in GNA_gwas1_*.csv; do tail -n +2 $file >> networkGWAS1_combined.csv; done
for file in GNA_gwas2_*.csv; do tail -n +2 $file >> networkGWAS2_combined.csv; done
"))

# read in results
files <- c("networkGWAS1_combined.csv","networkGWAS2_combined.csv")
snp_results <- lapply(files,fread)


### estimate effective N of factors NOTE very approximate estimate. Can estimate more accurately if needed
# apply MAF filter (as low and high MAF can bias the result)
#neff <- lapply(snp_results, function(x){
#  x %>%
#    filter(between(MAF, 0.1, 0.4)) %>%
#    mutate(neff = 1/((2*MAF*(1-MAF))*SE^2)) %>%
#    pull(neff) %>%
#    mean() %>%
#    floor()
#})

### cleaned sumstats
freq <- fread("/working/lab_esked/jacksonT/reference_files/GenomicSEM/reference.1000G.maf.0.005.eur_A1frq.txt")

snp_results <- lapply(1:length(snp_results), function(x){
  snp_results[[x]] %>%
    rename(BETA = est, Z = Z_Estimate, P = Pval_Estimate) %>%
    inner_join(freq) %>%
    #mutate(NEFF = neff[[x]]) %>%
    select(SNP,CHR,BP,A1,A2,A1FRQ,BETA,SE,Z,P,NEFF)
})



# quick check of manhattan
library(qqman)

tmp <- snp_results[[2]]

manhattan(filter(tmp, p_CUD < 0.01), p = "p_CUD")
manhattan(filter(tmp, p_EVER < 0.01), p = "p_EVER")
manhattan(filter(tmp, p_FREQ_woPast < 0.01), p = "p_FREQ_woPast")





### write final file
### just taking common1 - which is fix_measurement=TRUE [the default in GSEM]
setwd(working_dir)

fwrite(snp_results[[1]], file = "Cann_factor_assoc.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


## clump results 
dir.create(file.path(working_dir,"clumped"))
setwd(file.path(working_dir,"clumped"))

pbs_plink <- paste0("#!/bin/bash
#PBS -N plink_clump
#PBS -l mem=15gb,walltime=02:00:00,ncpus=5
#PBS -o ",working_dir,"/clumped
#PBS -e ",working_dir,"/clumped
module load plink/1.90b6.8

cd ",working_dir,"

plink --bfile /reference/data/UKBB_500k/versions/lab_esked/LDref/50kRandomWhiteBrit/chr1-22_ukb_50k_LD_ref_INFO_0.8_MAF_0.01 \\
--clump Cann_factor_assoc.tsv \\
--clump-p1 5e-8 \\
--clump-p2 0.05 \\
--clump-r2 0.1 \\
--clump-kb 250 \\
--out ",file.path(working_dir,"clumped/Internalising_clumped_r0.1")," \\
--threads 5

")

write(pbs_plink,file = "clump_plink.pbs")
system("qsub clump_plink.pbs")

