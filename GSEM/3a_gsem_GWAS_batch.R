#####################################
# Cannabis Use common factor GWAS   #
# ICC                               #
#####################################

require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)
require(GenomicSEM)

working_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM"

#================================#
# Run model WITHOUT SNP effects  #
#================================#

setwd(working_dir)

LDSCoutput1 <- readRDS("LDSCoutput_cannabis_ALL.RDS")
LDSCoutput2 <- readRDS("LDSCoutput_cannabis_ALL_noPast.RDS")

model1 <- "F =~ EVER + FREQ_wPast + CUD
FREQ_wPast~~a*FREQ_wPast
a>0.001"
model2 <- "F =~ EVER + FREQ_woPast + CUD
FREQ_woPast~~a*FREQ_woPast
a>0.001"

output1 <- usermodel(LDSCoutput1, estimation = "DWLS", model = model1)
output1
output2 <- usermodel(LDSCoutput2, estimation = "DWLS", model = model2)
output2 

fwrite(output1[[2]], file = "Common_factor_model_results.csv", quote = F)
fwrite(output2[[2]], file = "Common_factor_model_results_noPast.csv", quote = F)

#==========================================#
# SNP effects - create pbs batch job files #
#==========================================#

##### specify below arguments and paths
ldsc_out1 = "LDSCoutput_cannabis_ALL.RDS" #path to file of multivariate LDSC R object
ldsc_out2 = "LDSCoutput_cannabis_ALL_noPast.RDS"
sumstats1 = "sumstats_cannabis_ALL.tsv"      #path to sumstats (output from Genomic SEM sumstats() function)
sumstats2 = "sumstats_cannabis_ALL_noPast.tsv"

num_SNPs = 7254575 # total number of SNPs in sumstats
batch_size = 3000 # number of SNPs to run per batch
memory = "3500MB" # specify memory per batch
walltime = "02:00:00" # specify time per patch

model1 <- "F =~ EVER + FREQ_wPast + CUD
FREQ_wPast~~a*FREQ_wPast
a>0.001
F~SNP
SNP~~SNP"

model2 <- "F =~ EVER + FREQ_woPast + CUD
FREQ_woPast~~a*FREQ_woPast
a>0.001
F~SNP
SNP~~SNP"

##### RUN BELOW TO GENERATE SCRIPTS AND FILES

setwd(working_dir)
dir.create("SNPeffects")
dir.create("pbs")


### 1) array values file 

num_batch <- ceiling(num_SNPs/batch_size)
values <- data.frame(matrix(ncol = 2, nrow = num_batch))
values[1] <- seq(from=1, to=num_SNPs, by=batch_size)
values[2] <- seq(from=batch_size, by = batch_size, length.out = num_batch)
values[num_batch,2] <- num_SNPs
write.table(values, file = "run_gsemGWAS_array.values", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



### 2) pbs script

pbs <- paste0("#!/bin/bash
#PBS -N GenomicSEM
#PBS -l mem=",memory,",walltime=",walltime,",ncpus=1
#PBS -l avx2=True
#PBS -J 1-",num_batch,"
#PBS -o ",working_dir,"/pbs
#PBS -e ",working_dir,"/pbs
module load R/4.3.1
cd $PBS_O_WORKDIR

VARIABLES=($(awk -v var=$PBS_ARRAY_INDEX 'NR==var' run_gsemGWAS_array.values))
currStart=${VARIABLES[0]}
currEnd=${VARIABLES[1]}

Rscript run_gsemGWAS_array.R $currStart $currEnd")

write(pbs,file="run_gsemGWAS_array.pbs")



### 3) Genomic SEM script

GSEM <- paste0("args = commandArgs(trailingOnly=TRUE)
n_start <- args[1]
n_stop <- args[2]


require(data.table)
require(GenomicSEM)
LDSCoutput1 <- readRDS('",ldsc_out1,"')
LDSCoutput2 <- readRDS('",ldsc_out2,"')
sumstats1 <- fread('",sumstats1,"', data.table = FALSE)
sumstats2 <- fread('",sumstats2,"', data.table = FALSE)
subsetSNPs1 <- sumstats1[paste0(n_start):paste0(n_stop),]
subsetSNPs2 <- sumstats2[paste0(n_start):paste0(n_stop),]
model1 <- '",model1,"'
model2 <- '",model2,"'

print('Starting the common factor model GWAS!')
results1 <- userGWAS(covstruc=LDSCoutput1,SNPs=subsetSNPs1,estimation='DWLS',model=model1,sub='F~SNP',parallel=FALSE,GC='none')
results2 <- userGWAS(covstruc=LDSCoutput2,SNPs=subsetSNPs2,estimation='DWLS',model=model2,sub='F~SNP',parallel=FALSE,GC='none')

print('Finished the common factor model GWAS! Writing results to file.')
write.csv(results1, file = paste0('",working_dir,"/SNPeffects/GSEM_common1_',n_start,'_',n_stop,'.csv'))
write.csv(results2, file = paste0('",working_dir,"/SNPeffects/GSEM_common2_',n_start,'_',n_stop,'.csv'))

")
write(GSEM,file="run_gsemGWAS_array.R")



#### submit job
system("qsub run_gsemGWAS_array.pbs")

