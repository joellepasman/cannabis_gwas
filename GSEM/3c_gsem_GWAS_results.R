##################################################
# Cannabis GWAS                                  #
# - Get results                                  #
##################################################

require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)
require(GenomicSEM)

working_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM"


#=============================#
# CHECK AND COMBINE RESULTS   #
#=============================#

# combine all batches
setwd(file.path(working_dir,"SNPeffects"))

system(paste0("
head -n1 GSEM_common1_1_* > Common1_GSEM_combined.csv
head -n1 GSEM_common2_1_* > Common2_GSEM_combined.csv
for file in GSEM_common1_*.csv; do tail -n +2 $file >> Common1_GSEM_combined.csv; done
for file in GSEM_common2_*.csv; do tail -n +2 $file >> Common2_GSEM_combined.csv; done
"))

# read in results
files <- c("Common1_GSEM_combined.csv","Common2_GSEM_combined.csv")
snp_results <- lapply(files,fread)

# check for SNPs with errors to inspect (zero)
fail <- lapply(snp_results, filter, error != 0 | warning != 0)


### estimate effective N of factors NOTE very approximate estimate.
# apply MAF filter (as low and high MAF can bias the result)
#neff <- lapply(snp_results, function(x){
#  x %>%
#    filter(between(MAF, 0.1, 0.4)) %>%
#    mutate(neff = 1/((2*MAF*(1-MAF))*SE^2)) %>%
#    pull(neff) %>%
#    mean() %>%
#    floor()
#})

### Effective N - used METAL in order to get the effective N accounting for sample overlap
files <- c("/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM/Metal/META_CANN_effN1.txt", "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM/Metal/META_CANN_effN_noPast1.txt")

nlist <- lapply(files, function(x){
  fread(x) %>%
    rename(SNP = MarkerName) %>%
    select(SNP, N)
})


### cleaned sumstats
freq <- fread("/working/lab_esked/jacksonT/reference_files/GenomicSEM/reference.1000G.maf.0.005.eur_A1frq.txt")

snp_results <- lapply(1:length(snp_results), function(x){
  snp_results[[x]] %>%
    rename(BETA = est, Z = Z_Estimate, P = Pval_Estimate) %>%
    inner_join(freq) %>%
    inner_join(nlist[[x]]) %>%
    rename(EAF = A1FRQ, NEFF = N) %>%
    select(SNP,CHR,BP,A1,A2,EAF,BETA,SE,Z,P,NEFF)
})



# quick check of manhattan
library(qqman)
lapply(snp_results, function(x) x %>% filter(P < 0.001) %>% manhattan())




### write final file
### just taking common1 - which is fix_measurement=TRUE [the default in GSEM]
setwd(working_dir)

fwrite(snp_results[[1]], file = "ICC4_CANN_FACTOR_250212.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
fwrite(snp_results[[2]], file = "ICC4_CANN_FACTOR_noPast_250212.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


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
for file in ICC4_CANN_*.txt
do
plink --bfile /reference/data/UKBB_500k/versions/lab_esked/LDref/50kRandomWhiteBrit/chr1-22_ukb_50k_LD_ref_INFO_0.8_MAF_0.01 \\
--clump $file \\
--clump-p1 5e-8 \\
--clump-p2 0.05 \\
--clump-r2 0.1 \\
--clump-kb 250 \\
--out ",file.path(working_dir,"clumped/${file}_clumped_r0.1")," \\
--threads 5
done
")

write(pbs_plink,file = "clump_plink.pbs")
system("qsub clump_plink.pbs")



#### ldsc - h2snp etc.

setwd(working_dir)

# munge
files = c("ICC4_CANN_FACTOR_250212.txt", "ICC4_CANN_FACTOR_noPast_250212.txt")
trait.names = c("CANN_wPast","CANN_noPast")

munge(files=files, hm3="/working/lab_esked/jacksonT/reference_files/ldsc/w_hm3.noMHC.snplist", trait.names=trait.names)

# ldsc
ld="/working/lab_esked/jacksonT/reference_files/ldsc/eur_w_ld_chr/"
wld="/working/lab_esked/jacksonT/reference_files/ldsc/eur_w_ld_chr/"
traits = paste0(trait.names,".sumstats.gz")
sample.prev <- c(NA,NA)
population.prev <- c(NA,NA)

LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, trait.names=trait.names, ld=ld, wld=wld, stand = TRUE)
saveRDS(LDSCoutput, file="LDSCoutput_CANN_FACTOR_PastvsNoPast.RDS")
