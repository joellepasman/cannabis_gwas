#####################################
# Cannabis Use GSEM - multi LDSC    #
# ICC                               #
#####################################

require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)
require(GenomicSEM)

working_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM"
munged_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM/munged_sumstats"



#=====================#
# Munge summary stats #
#=====================#
setwd(munged_dir)

files = c("/working/lab_esked/jacksonT/cannabisUseGWAS/ICC4_meta/sumstats/MA_CANNABIS_EVR_250212_withNEFF.txt",
          "/working/lab_esked/jacksonT/cannabisUseGWAS/ICC4_meta/meta_freq_metal/ICC4_FREQ_META_250212.txt",
          "/working/lab_esked/jacksonT/cannabisUseGWAS/ICC4_meta/meta_freq_metal_noPast/ICC4_FREQ_META_noPast_250212.txt",
          "/working/lab_esked/jacksonT/cannabisUseGWAS/ICC4_meta/sumstats/CUD_clean.txt")
trait.names = c("EVER","FREQ_wPast","FREQ_woPast","CUD")

munge(files=files, hm3="/working/lab_esked/jacksonT/reference_files/ldsc/w_hm3.noMHC.snplist", trait.names=trait.names)


#====================#
# Multivariable LDSC #
#====================#
setwd(working_dir)

ld="/working/lab_esked/jacksonT/reference_files/ldsc/eur_w_ld_chr/"
wld="/working/lab_esked/jacksonT/reference_files/ldsc/eur_w_ld_chr/"

## with Past
trait.names = c("EVER","FREQ_wPast","CUD")
traits = file.path(munged_dir, paste0(trait.names,".sumstats.gz"))
sample.prev <- c(0.5,NA,0.5)
population.prev <- c(0.43,NA,0.02)
LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, trait.names=trait.names, ld=ld, wld=wld, stand = TRUE)
saveRDS(LDSCoutput, file="LDSCoutput_cannabis_ALL.RDS")

## without Past
trait.names = c("EVER","FREQ_woPast","CUD")
traits = file.path(munged_dir, paste0(trait.names,".sumstats.gz"))
sample.prev <- c(0.5,NA,0.5)
population.prev <- c(0.43,NA,0.02)
LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, trait.names=trait.names, ld=ld, wld=wld, stand = TRUE)
saveRDS(LDSCoutput, file="LDSCoutput_cannabis_ALL_noPast.RDS")


#===========================#
# Prepare sumstats for GWAS #
#===========================#

setwd(working_dir)
trait.names = c("EVER","FREQ_wPast","FREQ_woPast","CUD")
se.logit <- c(T,F,F,T)
OLS <- c(F,T,T,F)
linprob <- c(F,F,F,F)
ref <- "/working/lab_esked/jacksonT/reference_files/GenomicSEM/reference.1000G.maf.0.005.eur.txt"
sumstats_gsem <- sumstats(files=files, ref=ref, trait.names=trait.names, se.logit=se.logit, OLS=OLS, linprob=linprob, info.filter=.6, maf.filter=0.01, parallel=FALSE)

# write sumstats - with past
sumstats_gsem %>%
  select(!c(beta.FREQ_woPast, se.FREQ_woPast)) %>%
  fwrite(., file = "sumstats_cannabis_ALL.tsv", sep = "\t")

# write sumstats - without past
sumstats_gsem %>%
  select(!c(beta.FREQ_wPast, se.FREQ_wPast)) %>%
  fwrite(., file = "sumstats_cannabis_ALL_noPast.tsv", sep = "\t")
