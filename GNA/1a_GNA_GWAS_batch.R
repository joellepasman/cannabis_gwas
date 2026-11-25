#####################################
# Cannabis Use network GWAS         #
# ICC                               #
#####################################

require(data.table)
options(datatable.fread.datatable=FALSE)
require(tidyverse)
require(GNA)

working_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GNA"
data_dir = "/working/lab_esked/jacksonT/cannabisUseGWAS/GSEM"

#================================#
# Run model WITHOUT SNP effects  #
#================================#

setwd(working_dir)

LDSCoutput1 <- readRDS(file.path(data_dir,"LDSCoutput_cannabis_ALL.RDS"))
LDSCoutput2 <- readRDS(file.path(data_dir,"LDSCoutput_cannabis_ALL_noPast.RDS"))

output1 <- traitNET(LDSCoutput1)
output2 <- traitNET(LDSCoutput2)

saveRDS(output1$model_results$sparse$omega, file = "sparse_omega1.RDS")
saveRDS(output2$model_results$sparse$omega, file = "sparse_omega2.RDS")


#==========================================#
# SNP effects - create pbs batch job files #
#==========================================#

##### specify below arguments and paths
ldsc_out1 = file.path(data_dir,"LDSCoutput_cannabis_ALL.RDS") #path to file of multivariate LDSC R object
ldsc_out2 = file.path(data_dir,"LDSCoutput_cannabis_ALL_noPast.RDS")
sumstats1 = file.path(data_dir,"sumstats_cannabis_ALL.tsv")      #path to sumstats (output from Genomic SEM sumstats() function)
sumstats2 = file.path(data_dir,"sumstats_cannabis_ALL_noPast.tsv")
omega1 = file.path(working_dir,"sparse_omega1.RDS")
omega2 = file.path(working_dir,"sparse_omega2.RDS")

num_SNPs = 7254575 # total number of SNPs in sumstats
batch_size = 3000 # number of SNPs to run per batch
memory = "3500MB" # specify memory per batch
walltime = "02:00:00" # specify time per patch


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
write.table(values, file = "run_networkGWAS_array.values", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



### 2) pbs script

pbs <- paste0("#!/bin/bash
#PBS -N GNA_networkGWAS
#PBS -l mem=",memory,",walltime=",walltime,",ncpus=1
#PBS -l avx2=True
#PBS -J 1-",num_batch,"
#PBS -o ",working_dir,"/pbs
#PBS -e ",working_dir,"/pbs
module load R/4.3.1
cd $PBS_O_WORKDIR

VARIABLES=($(awk -v var=$PBS_ARRAY_INDEX 'NR==var' run_networkGWAS_array.values))
currStart=${VARIABLES[0]}
currEnd=${VARIABLES[1]}

Rscript run_networkGWAS_array.R $currStart $currEnd")

write(pbs,file="run_networkGWAS_array.pbs")



### 3) GNA script

GNA <- paste0("args = commandArgs(trailingOnly=TRUE)
n_start <- args[1]
n_stop <- args[2]


require(data.table)
require(GNA)
LDSCoutput1 <- readRDS('",ldsc_out1,"')
LDSCoutput2 <- readRDS('",ldsc_out2,"')
sumstats1 <- fread('",sumstats1,"', data.table = FALSE)
sumstats2 <- fread('",sumstats2,"', data.table = FALSE)
subsetSNPs1 <- sumstats1[paste0(n_start):paste0(n_stop),]
subsetSNPs2 <- sumstats2[paste0(n_start):paste0(n_stop),]
omega1 <- readRDS('",omega1,"')
omega2 <- readRDS('",omega2,"')

print('Starting the network GWAS!')
results1 <- gwasNET(covstruc=LDSCoutput1,SNPs=subsetSNPs1,fix_omega=omega1,parallel=FALSE)
results2 <- gwasNET(covstruc=LDSCoutput2,SNPs=subsetSNPs2,fix_omega=omega2,parallel=FALSE)

print('Finished the network GWAS! Writing results to file.')
write.csv(results1, file = paste0('",working_dir,"/SNPeffects/GNA_gwas1_',n_start,'_',n_stop,'.csv'))
write.csv(results2, file = paste0('",working_dir,"/SNPeffects/GNA_gwas2_',n_start,'_',n_stop,'.csv'))

")
write(GNA,file="run_networkGWAS_array.R")



#### submit job
system("qsub run_networkGWAS_array.pbs")

