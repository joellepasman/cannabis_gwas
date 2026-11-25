require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)
require(GenomicSEM)

setwd("/working/lab_esked/jacksonT/cannabisUseGWAS/ICC4_meta/sumstats/")


#================#
# Clean sumstats #
#================#

# lifetime use
cohorts <- fread("cohorts_ICC_ever.csv")
cohorts <- cohorts %>% mutate(Neff = 4*(1/(1/cases + 1/controls)))
neff <- sum(cohorts$Neff)
prev <- sum(cohorts$cases)/(sum(cohorts$cases)+sum(cohorts$controls))

ever <- fread("MA_CANNABIS_EVR_250212.txt")
ever <- ever %>% mutate(N_EFF = neff) %>%
  select(!c("N"))
fwrite(ever, file = "MA_CANNABIS_EVR_250212_withNEFF.txt", sep = " ", quote = FALSE)

# CUD (Levey 2023)
cud <- fread("dbGAP_EUR_CanMVP_META.gz")
cud <- cud %>%
  select(!MarkerName) %>%
  rename(P = P.value) %>%
  mutate(N_EFF = 161053,  #reported in paper
         A1 = case_match(A1, "g"~"G", "c"~"C", "t"~"T", "a"~"A", "d"~"D", "i"~"I"),
         A2 = case_match(A2, "g"~"G", "c"~"C", "t"~"T", "a"~"A", "d"~"D", "i"~"I")) %>%
  na.omit()
# add EAF column from ever phenotype
cud <- cud %>%
  inner_join(select(ever,SNP,A1,A2,EAF), by = "SNP") %>%
  filter(((A1.x==A1.y) & (A2.x==A2.y)) | ((A1.x==A2.y) & (A2.x==A1.y))) %>%
  mutate(FREQ = if_else(A1.x == A1.y, EAF, 1-EAF)) %>%
  select(SNP,CHR,BP,A1.x,A2.x,FREQ,BETA,SE,P,N_EFF) %>%
  rename(A1=A1.x, A2=A2.x, EAF=FREQ)
fwrite(cud, file = "CUD_clean.txt.gz", sep = " ", quote = FALSE)

# FREQ already cleaned - just copied from original to cleaned