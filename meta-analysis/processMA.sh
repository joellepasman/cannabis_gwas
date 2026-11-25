#!/sw/arch/RHEL8/EB_production/2022/software/R/4.2.1-foss-2022a/bin/Rscript
#SBATCH --job-name=processMA
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=rome
#SBATCH --output=slurmout/processMA_%j.out


################################################################################
## Format

library(tidyverse)
library(data.table)
library(tidyGWAS)
library(R.utils)

# Format and tidy

formatsumstats <- function(trait) {

	df <- fread(paste0("MA_", trait, "1.txt"))

	# Add Z, rename columns 

	if (!("Zscore" %in% colnames(df))) {
	df <- df %>% 
		mutate(Zscore = Effect/ StdErr) %>%
		rename(B = Effect, SE = StdErr)
	}

	sumstats <- df %>%
		distinct(MarkerName, .keep_all = TRUE) %>%
		mutate(CHR = str_extract(MarkerName, "^[^_]*"),
		POS = str_extract(MarkerName, "(?<=_)[^_]*")) %>%
		select(!MarkerName) %>%
		rename(EffectAllele = Allele1,
			OtherAllele = Allele2,
			EAF = Freq1,
			P = `P-value`,
			Z = Zscore) 

	# Using tidygwas to add RS number and beta/SE

	cleaned <- tidyGWAS(
	 tbl = sumstats, 
	 output_format = "csv",
	 logfile = TRUE,
	  dbsnp_path = "/gpfs/work5/0/prjs0778/dbSNP155_v09/",
	  repair_cols = T,	
	 output_dir = paste0("tidy_", trait, "250919")) 

	# Filter MAF and remove duplicate RS IDs 
	cleaned <- cleaned %>%
		filter(EAF > .01 & EAF < .99) %>%
		distinct(RSID, .keep_all = TRUE)

	# inner_join with HRC to retain high-quality SNPs
	hrc <- fread("/gpfs/home5/joepas/sumstats/reference/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.rsid_map.gz") 

	cleaned <- cleaned %>%
		inner_join(hrc, by=c("RSID" = "rsid")) %>%
		select(!c(chr, pos))

	cleaned <- cleaned %>%
		rename(SNP = RSID, BP = POS_37, A1 = EffectAllele, A2 = OtherAllele, BETA =B) %>% 
		select(SNP, CHR, BP, A1, A2, EAF, BETA, SE, P, N) %>%
		drop_na()

	fwrite(cleaned, paste0("MA_", trait, "_250919.txt"), sep = " ", quote = F)

}

traits <- c("CANNABIS_EVR", "FREQ_LIFETIME", "FREQ_MAX", "FREQ_PAST")
	 
map(traits, formatsumstats)

system("gzip MA_*_250919.txt")

