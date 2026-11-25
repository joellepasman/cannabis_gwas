rm(list=ls())
library(data.table)
library(gprofiler2)
library(tidyr)
library(dplyr)
library(biomaRt)
library(tools)
library(ggplot2)
library(forcats)
library(stringr)
library(readxl)
library(purrr)
library(patchwork)
library(scales)
library(grid)
library(cowplot)
library(patchwork)
library(ggvenn)
library(openxlsx)
library(ggplotify)
library(eulerr)


##########################################
# IMPORT DATASETS                        #
##########################################

# define working directories
smr_dir <- "/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/gbt/smr"
twas_dir <- "/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/gbt/twas"
magma_dir <- "/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/gbt/magma"
mbat_dir <- "/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/gbt/mbat"

## smr

# brainmeta
initiation_brainmeta_files <- list.files(smr_dir,"smr_brainmeta_MA_CANNABIS_EVR_250919",full.names = TRUE)
frequency_brainmeta_files <- list.files(smr_dir,"smr_brainmeta_MA_FREQ_ALL_250919",full.names = TRUE)
factor_brainmeta_files <- list.files(smr_dir,"smr_brainmeta_MA_CANN_FACTOR_250919",full.names = TRUE)

initiation_eSMR_brainmeta <- lapply(initiation_brainmeta_files,read.delim)
initiation_eSMR_brainmeta <- do.call(rbind.data.frame,initiation_eSMR_brainmeta)

frequency_eSMR_brainmeta <- lapply(frequency_brainmeta_files,read.delim)
frequency_eSMR_brainmeta <- do.call(rbind.data.frame,frequency_eSMR_brainmeta)

factor_eSMR_brainmeta <- lapply(factor_brainmeta_files,read.delim)
factor_eSMR_brainmeta <- do.call(rbind.data.frame,factor_eSMR_brainmeta)

# eqtlgen
initiation_eSMR_eqtlgen <- read.delim(paste0(smr_dir,"/smr_eqtlgen_MA_CANNABIS_EVR_250919.smr"))
frequency_eSMR_eqtlgen <- read.delim(paste0(smr_dir,"/smr_eqtlgen_MA_FREQ_ALL_250919.smr"))
factor_eSMR_eqtlgen <- read.delim(paste0(smr_dir,"/smr_eqtlgen_MA_CANN_FACTOR_250919.smr"))

# build lists for brainmeta and eqtlgen results 
smr_brainmeta_res <- list(initiation_eSMR_brainmeta,frequency_eSMR_brainmeta,factor_eSMR_brainmeta)
names(smr_brainmeta_res) <- c("Initiation","Frequency","Factor")

smr_eqtlgen_res = list(initiation_eSMR_eqtlgen,frequency_eSMR_eqtlgen,factor_eSMR_eqtlgen)
names(smr_eqtlgen_res) <- c("Initiation","Frequency","Factor")

## twas

# initiation
initiation_twas_files <- list.files(twas_dir,"CANNABIS_EVR_250919",full.names = T)
initiation_twas_list <- lapply(initiation_twas_files,read.delim)
initiation_twas <- do.call(rbind.data.frame,initiation_twas_list)

# frequency
frequency_twas_files <- list.files(twas_dir,"FREQ_ALL_250919",full.names = T)
frequency_twas_list <- lapply(frequency_twas_files,read.delim)
frequency_twas <- do.call(rbind.data.frame,frequency_twas_list)

# factor
factor_twas_files <- list.files(twas_dir,"CANN_FACTOR_250919",full.names = T)
factor_twas_list <- lapply(factor_twas_files,read.delim)
factor_twas <- do.call(rbind.data.frame,factor_twas_list)

# combine data
twas_res <- list(initiation_twas,frequency_twas,factor_twas)
names(twas_res) <- c("Initiation","Frequency","Factor")

# define function to convert columns to numeric
convert_to_numeric <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- as.numeric(data[[col]])
  }
  return(data)
}

# convert columns to numeric
columns_to_convert <- colnames(initiation_twas[,c(19:25)])
for(i in 1:length(twas_res)){
  twas_res[[i]] <- convert_to_numeric(twas_res[[i]], columns_to_convert)
}

# remove missing z scores
for(i in 1:length(twas_res)){
  twas_res[[i]] <- twas_res[[i]][!is.na(twas_res[[i]]$TWAS.Z),]
}

# split pec and gtex
twas_pec_res <- lapply(twas_res,function(x) subset(x,x$PANEL=="PEC.BRAIN.RNASEQ.RSID"))
twas_gtex_res <- lapply(twas_res,function(x) subset(x,x$PANEL!="PEC.BRAIN.RNASEQ.RSID"))
names(twas_pec_res) <- c("Initiation","Frequency","Factor")
names(twas_gtex_res) <- c("Initiation","Frequency","Factor")

# create ensg for twas pec and gtex
twas_pec_res <- lapply(twas_pec_res, function(df) {
  df$ensgid <- sub("\\..*$", "", sapply(strsplit(df$FILE, "/"), "[[", 12))
  df
})

# create ensg for twas pec and gtex
twas_gtex_res <- lapply(twas_gtex_res, function(df) {
  df$ensgid <- sub("\\..*$", "", sapply(strsplit(df$FILE, "/"), "[[", 11))
  df
})


## magma ##

# import magma results
magma_initiation <- as.data.frame(fread(paste0(magma_dir,"/magma_MA_CANNABIS_EVR_250919.genes.out")))
magma_frequency <- as.data.frame(fread(paste0(magma_dir,"/magma_MA_FREQ_ALL_250919.genes.out")))
magma_factor <- as.data.frame(fread(paste0(magma_dir,"/magma_MA_CANN_FACTOR_250919.genes.out")))

# create list
magma_res <- list(magma_initiation,magma_frequency,magma_factor)
names(magma_res) <- c("Initiation","Frequency","Factor")

## mbat

# import results
mbat_initiation <- read.delim(paste0(mbat_dir,"/mbat_MA_CANNABIS_EVR_250919.gene.assoc.mbat"))
mbat_frequency <- read.delim(paste0(mbat_dir,"/mbat_MA_FREQ_ALL_250919.gene.assoc.mbat"))
mbat_factor <- read.delim(paste0(mbat_dir,"/mbat_MA_CANN_FACTOR_250919.gene.assoc.mbat"))

# create list
mbat_res <- list(mbat_initiation,mbat_frequency,mbat_factor)
names(mbat_res) <- c("Initiation","Frequency","Factor")

## calculate fdr and bonferroni thresholds

# magma
magma_res <- lapply(magma_res, function(file) {
  file <- file[order(file$P),]
  file$fdr <- p.adjust(file$P, "fdr")
  file$bonferroni <- p.adjust(file$P, "bonferroni")
  return(file)
})

# mbat-combo
mbat_res <- lapply(mbat_res, function(file) {
  file <- file[order(file$P_mBATcombo),]
  file$fdr <- p.adjust(file$P_mBATcombo, "fdr")
  file$bonferroni <- p.adjust(file$P_mBATcombo, "bonferroni")
  return(file)
})

# smr eqtlgen
smr_eqtlgen_res <- lapply(smr_eqtlgen_res, function(file) {
  file <- file[order(file$p_SMR),]
  file$fdr <- p.adjust(file$p_SMR, "fdr")
  file$bonferroni <- p.adjust(file$p_SMR, "bonferroni")
  return(file)
})

# smr eqtlgen
smr_brainmeta_res <- lapply(smr_brainmeta_res, function(file) {
  file <- file[order(file$p_SMR),]
  file$fdr <- p.adjust(file$p_SMR, "fdr")
  file$bonferroni <- p.adjust(file$p_SMR, "bonferroni")
  return(file)
})

# twas pec
twas_pec_res <- lapply(twas_pec_res, function(file) {
  file <- file[order(file$TWAS.P),]
  file$fdr <- p.adjust(file$TWAS.P, "fdr")
  file$bonferroni <- p.adjust(file$TWAS.P, "bonferroni")
  return(file)
})

# twas gtex
twas_gtex_res <- lapply(twas_gtex_res, function(file) {
  file <- file[order(file$TWAS.P),]
  file$fdr <- p.adjust(file$TWAS.P, "fdr")
  file$bonferroni <- p.adjust(file$TWAS.P, "bonferroni")
  return(file)
})

# select tissue with smallest p value for twas gtex
twas_gtex_res <- lapply(twas_gtex_res, function(df) {
  df$TISSUE <- sapply(strsplit(df$PANEL, "\\."), "[[", 3)
  df
})

twas_gtex_res_minP <- lapply(twas_gtex_res, function(df) {
  as.data.frame(
    df %>%
      group_by(ID) %>%
      slice_min(TWAS.P, n = 1, with_ties = FALSE) %>%
      ungroup()
  )
})

twas_gtex_res_minP <- lapply(twas_gtex_res_minP, function(df) {
  df[order(df$TWAS.P), ]
})

for(i in 1:length(smr_brainmeta_res)){
  smr_brainmeta_res[[i]]$probeID=sapply(strsplit(smr_brainmeta_res[[i]]$probeID,"\\."),"[[",1)
}

# create a 'parent list' of all gene-based results
parent_list <- list(magma_res,mbat_res,smr_brainmeta_res,smr_eqtlgen_res,twas_gtex_res_minP,twas_pec_res)


# harmonize gene names
renamed_list <- lapply(parent_list, function(inner_list) {
  lapply(inner_list, function(df) {
    if ("gene_id" %in% names(df)) {
      return(df)
    }
    old_names <- c("GENE","probeID", "ensgid", "Gene")
    matched <- which(names(df) %in% old_names)
    if (length(matched) >= 1) {
      names(df)[matched[1]] <- "gene_id"
    }
    return(df)
  })
})

# harmonize p values
renamed_list <- lapply(renamed_list, function(inner_list) {
  lapply(inner_list, function(df) {
    old_names <- c("P", "P_mBATcombo", "p_SMR", "TWAS.P")
    matched <- which(names(df) %in% old_names)
    if (length(matched) == 1) {
      names(df)[matched] <- "pvalue"
    }
    return(df)
  })
})

names(renamed_list) <- c("MAGMA","mBAT","SMR_BrainMeta", "SMR_eQTLgen","TWAS_GTEx", "TWAS_PEC")

# select FDR < 0.05
filtered_nested_list <- lapply(renamed_list, function(inner_list) {
  lapply(inner_list, function(file) {
    file[which(file$fdr < 0.05), "gene_id"]
  })
})

# identify unique genes with FDR < 0.05
fdr_genes <- unique(filtered_nested_list)
flat_list <- unlist(filtered_nested_list, recursive = FALSE)
unique_gene_ids <- unique(unlist(flat_list))

##########################################
# CREATE GENE SPINE FOR MASTER GENE LIST #
##########################################

# collect gene coordinates for genes with FDR < 0.05
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")

gene_info <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id","chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = unique_gene_ids,
  mart = ensembl
)

# create a gene spine for merging of gene-based association statistics
gene_spine = gene_info[which(gene_info$chromosome_name %in% seq(1:22)),]
names(gene_spine) = c("hgnc_symbol","gene_id","chromosome_name","start_position","end_position")

# Ensure gene_id column is character for safe joins
gene_spine$gene_id <- as.character(gene_spine$gene_id)

# Define the columns to keep from each nested data frame
cols_to_keep <- c("gene_id", "pvalue", "fdr", "bonferroni")

# Perform the merge across the nested list
merged_list <- lapply(renamed_list, function(inner_list) {
  lapply(inner_list, function(df) {
    df <- df[, intersect(cols_to_keep, names(df)), drop = FALSE]
    if ("gene_id" %in% names(df)) {
      merge(gene_spine, df, by = "gene_id", all.x = TRUE)
    } else {
      df
    }
  })
})

names(merged_list) <- names(renamed_list)

# Define the columns to rename
cols_to_rename <- c("pvalue", "fdr", "bonferroni")

# Apply renaming to each method × trait data frame
master_list <- lapply(merged_list, function(trait_list) {
  lapply(names(trait_list), function(trait_name) {
    df <- trait_list[[trait_name]]
    for (col in cols_to_rename) {
      if (col %in% names(df)) {
        names(df)[names(df) == col] <- paste0(trait_name, "_", col)
      }
    }
    return(df)
  }) |> setNames(names(trait_list)) 
})

# Merge all 12 trait data frames per method
merge_keys <- c("gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")

merged_by_measure <- lapply(master_list, function(trait_list) {
  Reduce(function(df1, df2) {
    full_join(df1, df2, by = merge_keys)
  }, trait_list)
})

# Columns to merge on
merge_keys <- c("gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")

# rename all non-key columns to include method name
prefixed_dfs <- mapply(function(df, method_name) {
  non_key_cols <- setdiff(names(df), merge_keys)
  names(df)[names(df) %in% non_key_cols] <- paste0(method_name, "_", non_key_cols)
  return(df)
}, merged_by_measure, names(merged_by_measure), SIMPLIFY = FALSE)

# merge all 7 method-level data frames by merge_keys
final_merged <- Reduce(function(x, y) {
  full_join(x, y, by = merge_keys)
}, prefixed_dfs)

# remove rows with NA across all columns
final_merged_clean <- final_merged[rowSums(!is.na(final_merged[, -c(1:5)])) > 0, ]

# count the number of significant results
final_merged_clean %>%
  summarise(
    across(matches("(_bonferroni|_fdr)$"), ~ sum(.x < 0.05, na.rm = TRUE))
  ) %>%
  pivot_longer(
    everything(),
    names_to = "trait",
    values_to = "signif_count"
  ) %>%
  mutate(
    type = case_when(
      str_ends(trait, "_bonferroni") ~ "bonf_signif_count",
      str_ends(trait, "_fdr") ~ "fdr_signif_count"
    ),
    trait = str_remove(trait, "_bonferroni$|_fdr$")
  ) %>%
  pivot_wider(
    names_from = type,
    values_from = signif_count
  )

##########################################
# CREATE TIERED LIST OF GENES            #  
##########################################

# Function to rename and merge
rename_and_merge <- function(base_df, new_data, name_prefix, by_x = "gene_id", by_y = "ensgid", cols = c("COLOC.PP0", "COLOC.PP4")) {
  renamed <- new_data[, c(by_y, cols[cols %in% colnames(new_data)])] %>%
    rename_with(~ paste0(name_prefix, "_", .x), .cols = -all_of(by_y))
  merge(base_df, renamed, by.x = by_x, by.y = by_y, all.x = TRUE)
}

gene_rank <- final_merged_clean

# TWAS_PEC
for (i in seq_along(twas_pec_res)) {
  prefix <- paste0("TWAS_PEC_", names(twas_pec_res)[i])
  gene_rank <- rename_and_merge(gene_rank, twas_pec_res[[i]], prefix, cols = "COLOC.PP4")
}

# TWAS_GTEx (only COLOC.PP4)
for (i in seq_along(twas_gtex_res_minP)) {
  prefix <- paste0("TWAS_GTEx_", names(twas_gtex_res_minP)[i])
  gene_rank <- rename_and_merge(gene_rank, twas_gtex_res_minP[[i]], prefix, cols = "COLOC.PP4")
}

# SMR BrainMeta (uses gene_id and p_HEIDI)
for (i in seq_along(smr_brainmeta_res)) {
  prefix <- paste0("SMR_BrainMeta_", names(smr_brainmeta_res)[i])
  renamed <- smr_brainmeta_res[[i]][, c("probeID", "p_HEIDI")] %>%
    rename_with(~ paste0(prefix, "_", .x), .cols = -probeID)
  gene_rank <- merge(gene_rank, renamed, by.x = "gene_id", by.y = "probeID", all.x = TRUE)
}

# SMR eQTLgen (uses gene_id and p_HEIDI)
for (i in seq_along(smr_eqtlgen_res)) {
  prefix <- paste0("SMR_eQTLgen_", names(smr_eqtlgen_res)[i])
  renamed <- smr_eqtlgen_res[[i]][, c("probeID", "p_HEIDI")] %>%
    rename_with(~ paste0(prefix, "_", .x), .cols = -probeID)
  gene_rank <- merge(gene_rank, renamed, by.x = "gene_id", by.y = "probeID", all.x = TRUE)
}

###########################################################################
# COUNT NUMBER OF TIMES GENE IS SIGNIFICANT ACROSS ALL PHENOTYPES/METHODS #
###########################################################################

# identify all columns that match "_bonferroni"
bonf_cols <- grep("_bonferroni$", names(gene_rank), value = TRUE)
fdr_cols <- grep("_fdr$", names(gene_rank), value = TRUE)

# count number of times bonferroni < 0.05 per row
gene_rank$bonf_all_n <- rowSums(gene_rank[bonf_cols] < 0.05, na.rm = TRUE)
gene_rank$fdr_all_n <- rowSums(gene_rank[fdr_cols] < 0.05, na.rm = TRUE)

######################
# SPLIT BY PHENOTYPE #
######################

# Define phenotype names
phenotypes <- c("Initiation", "Frequency", "Factor")

# Columns to always keep
base_cols <- c("gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position")
sum_cols <- c("bonf_all_n", "fdr_all_n")

# Create a named list of data frames split by phenotype, removing suffixes
gene_rank_by_pheno <- lapply(phenotypes, function(pheno) {
  pheno_cols <- grep(paste0("_", pheno, "_"), names(gene_rank), value = TRUE)
  selected_cols <- c(base_cols, pheno_cols, sum_cols)
  df <- gene_rank[, selected_cols, drop = FALSE]
  
  # Remove _Initiation / _Frequency / _Factor from column names
  names(df) <- gsub(paste0("_", pheno, "_"), "_", names(df))
  
  return(df)
})

# Name the list by phenotype
names(gene_rank_by_pheno) <- phenotypes

# where hgnc_symbol is missing, replace with gene_id
gene_rank_by_pheno <- lapply(gene_rank_by_pheno, function(df) {
  df %>%
    mutate(hgnc_symbol = if_else(is.na(hgnc_symbol) | hgnc_symbol == "", gene_id, hgnc_symbol))
})



###############################################################
# COUNT NUMBER OF TIMES GENE IS SIGNIFICANT WITHIS PHENOTYPES #            
###############################################################

gene_rank_by_pheno <- lapply(gene_rank_by_pheno, function(df) {
  bonf_cols <- grep("_bonferroni$", names(df), value = TRUE)
  fdr_cols <- grep("_fdr$", names(df), value = TRUE)
  df$bonf_n <- rowSums(df[bonf_cols] < 0.05, na.rm = TRUE)
  df$fdr_n <- rowSums(df[fdr_cols] < 0.05, na.rm = TRUE)
  return(df)
})

# arrange by bonferroni n
gene_rank_by_pheno <- lapply(gene_rank_by_pheno, function(df) {
  df %>%
    arrange(desc(bonf_n))
})

########################
# TAG TWAS COLOC > 0.8 #
########################

# TWAS methods and phenotypes
twas_methods <- c("TWAS_GTEx", "TWAS_PEC")
phenotypes <- names(gene_rank_by_pheno) 

# Loop over each phenotype-specific data frame
for (pheno in phenotypes) {
  df <- gene_rank_by_pheno[[pheno]]
  
  for (method in twas_methods) {
    
    # Define column names
    bonf_col <- paste0(method, "_bonferroni")
    coloc_col <- paste0(method, "_COLOC.PP4")
    
    # Check if both columns are present
    if (bonf_col %in% names(df) && coloc_col %in% names(df)) {
      
      # Create indicator: Bonferroni < 0.05 and COLOC.PP4 > 0.8
      indicator_col <- paste0("bonf_and_coloc_", method)
      df[[indicator_col]] <- df[[bonf_col]] < 0.05 & df[[coloc_col]] > 0.8
      
    } else {
      warning(paste("Missing required columns for", method, pheno))
    }
  }
  
  # Save back to the list
  gene_rank_by_pheno[[pheno]] <- df
}

########################
# TAG SMR HEIDI > 0.05 #
########################

smr_methods <- c("SMR_BrainMeta", "SMR_eQTLgen")
phenotypes <- names(gene_rank_by_pheno) 

for (pheno in phenotypes) {
  df <- gene_rank_by_pheno[[pheno]]
  
  for (method in smr_methods) {
    
    # Define column names
    bonf_col <- paste0(method, "_bonferroni")
    heidi_col <- paste0(method, "_p_HEIDI")
    
    # Check if both columns exist
    if (bonf_col %in% names(df) && heidi_col %in% names(df)) {
      indicator_col <- paste0("bonf_and_heidi_", method)
      df[[indicator_col]] <- df[[bonf_col]] < 0.05 & df[[heidi_col]] > 0.05
      
    } else {
      warning(paste("Missing required columns for", method))
    }
  }
  
  # Save updated df back to the list
  gene_rank_by_pheno[[pheno]] <- df
}

##########################################################################
# CALCULATE MEAN -LOG10(P) AND ORDER BY NUMBER OF ASSOCIATIONS, LOG MEAN #
##########################################################################

for (pheno in phenotypes) {
  df <- gene_rank_by_pheno[[pheno]]
  pval_cols <- grep("_pvalue$", names(df), value = TRUE)
  
  if (length(pval_cols) > 0) {
    # Extract p-value matrix
    pval_matrix <- df[, pval_cols]
    log10_matrix <- -log10(pval_matrix)
    
    # Count non-missing values per row
    valid_counts <- rowSums(!is.na(log10_matrix))
    
    # Compute rowmean only if at least 3 non-missing values
    log10p_rowmean <- ifelse(
      valid_counts >= 3,
      rowMeans(log10_matrix, na.rm = TRUE),
      NA_real_
    )
    
    df$log10p_rowmean <- log10p_rowmean
  } else {
    warning(paste("No _pvalue columns found for", pheno))
  }
  
  gene_rank_by_pheno[[pheno]] <- df
}

# export supplementary information

# export each table as supplementary information
wb <- createWorkbook()
for (phenotype in names(gene_rank_by_pheno)) {
    df <- gene_rank_by_pheno[[phenotype]]
    sheet_name <- phenotype
    sheet_name <- substr(sheet_name, 1, 31)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = df)
}

# Save the workbook
saveWorkbook(wb, file = "/vast/scratch/users/gerring.z/icc4/results/suppl_tab-gene_res.xlsx", overwrite = TRUE)


##########################################
# DEFINE HIGH CONFIDENCE GENES           #  
##########################################

prioritised_gene_rank <- lapply(gene_rank_by_pheno, function(df) {
  df %>%
    mutate(
      ## --- Raw Bonferroni significance per method ---
      MAGMA_bonf_sig        = MAGMA_bonferroni        < 0.05,
      mBAT_bonf_sig         = mBAT_bonferroni         < 0.05,
      TWAS_PEC_bonf_sig     = TWAS_PEC_bonferroni     < 0.05,
      TWAS_GTEx_bonf_sig    = TWAS_GTEx_bonferroni    < 0.05,
      SMR_BrainMeta_bonf_sig= SMR_BrainMeta_bonferroni< 0.05,
      SMR_eQTLgen_bonf_sig  = SMR_eQTLgen_bonferroni  < 0.05,
      
      ## Any Bonferroni significance in any method
      any_bonf_sig = as.integer(
        MAGMA_bonf_sig | mBAT_bonf_sig |
          TWAS_PEC_bonf_sig | TWAS_GTEx_bonf_sig |
          SMR_BrainMeta_bonf_sig | SMR_eQTLgen_bonf_sig
      ),
      
      ## --- Evidence flags (Bonferroni + QC) ---
      
      ## Gene-based GWAS evidence
      Evidence_GeneGWAS = as.integer(MAGMA_bonf_sig | mBAT_bonf_sig),
      
      ## TWAS brain (PEC)
      Evidence_TWAS_brain = as.integer(
        TWAS_PEC_bonf_sig & TWAS_PEC_COLOC.PP4 >= 0.8
      ),
      
      ## TWAS blood (GTEx)
      Evidence_TWAS_blood = as.integer(
        TWAS_GTEx_bonf_sig & TWAS_GTEx_COLOC.PP4 >= 0.8
      ),
      
      ## SMR brain (BrainMeta)
      Evidence_SMR_brain = as.integer(
        SMR_BrainMeta_bonf_sig & SMR_BrainMeta_p_HEIDI > 0.01
      ),
      
      ## SMR blood (eQTLgen)
      Evidence_SMR_blood = as.integer(
        SMR_eQTLgen_bonf_sig & SMR_eQTLgen_p_HEIDI > 0.01
      ),
      
      ## Collapse expression evidence
      n_expr_blood  = Evidence_TWAS_blood + Evidence_SMR_blood,
      n_expr_brain  = Evidence_TWAS_brain + Evidence_SMR_brain,
      n_expr_total  = n_expr_blood + n_expr_brain,
      n_evidence_total = Evidence_GeneGWAS + n_expr_total,
      
      Tier = case_when(
        Evidence_GeneGWAS == 1 & n_expr_total >= 1 ~ "Tier1_high_conf",
        
        Evidence_GeneGWAS == 1 & n_expr_total == 0 ~ "Tier2_gene_only",
        
        Evidence_GeneGWAS == 0 & n_expr_total >= 2 ~ "Tier2_expr_only",
        
        any_bonf_sig == 1 ~ "Tier3_suggestive",
        
        ## Tier 0: nothing convincing
        TRUE ~ "Tier0_none"
      )
    )
})

tier_order <- c(
  "Tier1_high_conf",
  "Tier2_expr_only",
  "Tier2_gene_only",
  "Tier3_suggestive",
  "Tier0_none"
)

# Apply ordering to every data frame in your nested list

sorted_list <- lapply(prioritised_gene_rank, function(df) {
  df %>%
    mutate(
      Tier = factor(Tier, levels = tier_order)
    ) %>%
    arrange(Tier)
})

wb <- createWorkbook()

for (ph in names(sorted_list)) {
  
  # Select only columns 1:27 and 43:53
  df_out <- sorted_list[[ph]][ , c(1:27, 43:53)]
  
  addWorksheet(wb, sheetName = ph)
  writeData(wb, sheet = ph, x = df_out)
}

# Save workbook
saveWorkbook(wb, "/vast/scratch/users/gerring.z/icc4/results/gene_prioritisation_results.xlsx", overwrite = TRUE)

## Unique Bonferroni genes across all methods
long_df <- bind_rows(prioritised_gene_rank, .id = "phenotype")

gene_based_sig_genes <- long_df %>%
  filter(
    MAGMA_bonferroni        < 0.05 |
      mBAT_bonferroni         < 0.05 |
      TWAS_PEC_bonferroni     < 0.05 |
      TWAS_GTEx_bonferroni    < 0.05 |
      SMR_BrainMeta_bonferroni < 0.05 |
      SMR_eQTLgen_bonferroni   < 0.05
  ) %>%
  distinct(hgnc_symbol) %>%
  pull(hgnc_symbol)

n_gene_based_sig <- length(gene_based_sig_genes)
n_gene_based_sig


# function to count unique gene-based significant genes within one phenotype
count_gene_sig_by_pheno <- function(ph) {
  long_df %>%
    filter(
      phenotype == ph,
      MAGMA_bonferroni        < 0.05 |
        mBAT_bonferroni         < 0.05 |
        TWAS_PEC_bonferroni     < 0.05 |
        TWAS_GTEx_bonferroni    < 0.05 |
        SMR_BrainMeta_bonferroni < 0.05 |
        SMR_eQTLgen_bonferroni   < 0.05
    ) %>%
    distinct(hgnc_symbol) %>%
    nrow()
}

# Apply to each phenotype
phenos <- unique(long_df$phenotype)

bonf_sig_by_pheno <- setNames(
  sapply(phenos, count_gene_sig_by_pheno),
  phenos
)

bonf_sig_by_pheno


# A gene is Bonferroni-significant in any gene-based test if ANY of these are true:
gene_sig_filter <- expression(
  MAGMA_bonferroni        < 0.05 |
    mBAT_bonferroni         < 0.05 |
    TWAS_PEC_bonferroni     < 0.05 |
    TWAS_GTEx_bonferroni    < 0.05 |
    SMR_BrainMeta_bonferroni < 0.05 |
    SMR_eQTLgen_bonferroni   < 0.05
)

# Count unique genes per phenotype × tier
bonf_by_pheno_tier <- long_df %>%
  filter(eval(gene_sig_filter)) %>%
  group_by(phenotype, Tier) %>%
  summarise(
    n_unique_genes = n_distinct(hgnc_symbol),
    genes = list(unique(hgnc_symbol)),
    .groups = "drop"
  )

bonf_by_pheno_tier


tier1_by_pheno <- long_df %>%
  filter(Tier == "Tier1_high_conf") %>%
  group_by(phenotype) %>%
  summarise(genes = list(unique(hgnc_symbol)),
            .groups = "drop")

tier1_by_pheno

g_init  <- tier1_by_pheno$genes[[which(tier1_by_pheno$phenotype=="Initiation")]]
g_fact  <- tier1_by_pheno$genes[[which(tier1_by_pheno$phenotype=="Factor")]]
g_freq  <- tier1_by_pheno$genes[[which(tier1_by_pheno$phenotype=="Frequency")]]

overlap_init_factor <- intersect(g_init, g_fact)
overlap_init_freq   <- intersect(g_init, g_freq)
overlap_factor_freq <- intersect(g_fact, g_freq)

overlap_all3 <- Reduce(intersect, list(g_init, g_fact, g_freq))
length(overlap_all3)

tier1_union <- unique(c(g_init, g_fact, g_freq))
length(tier1_union)

tier1_overlap_summary <- tibble(
  comparison = c("Init ∩ Factor", "Init ∩ Frequency", "Factor ∩ Frequency", "All 3"),
  n_overlap  = c(
    length(overlap_init_factor),
    length(overlap_init_freq),
    length(overlap_factor_freq),
    length(overlap_all3)
  ),
  genes = list(
    overlap_init_factor,
    overlap_init_freq,
    overlap_factor_freq,
    overlap_all3
  )
)

tier1_overlap_summary


## rank -- use for plotting

ranked_top20 <- long_df %>%
  filter(Tier %in% c("Tier1_high_conf", "Tier2_gene_only", "Tier2_expr_only")) %>%
  arrange(desc(log10p_rowmean)) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  dplyr::select(hgnc_symbol, Tier, log10p_rowmean)

ranked_top20

gene_rank_by_pheno <- lapply(gene_rank_by_pheno, function(df) {
  df %>%
    filter(!grepl("^ENSG", hgnc_symbol)) %>%   # remove ENSG-prefixed gene symbols
    arrange(desc(log10p_rowmean))
})

gene_rank_by_pheno_top20 <- lapply(gene_rank_by_pheno, function(df) {
  df[which(df$hgnc_symbol %in% ranked_top20$hgnc_symbol),]
})

# export genes for annotation with chatGPT agent

write.table(gene_rank_by_pheno_top20$Initiation[,c(2:5)],file="/vast/scratch/users/gerring.z/icc4/results/genes_to_annotate.txt",row.names = F,col.names = T,sep="\t", quote=F)

##########################################
# PLOT DATA                              #  
##########################################

# Expecting names like "Initiation", "Frequency", "Factor" (case-insensitive OK)
pheno_map <- c(
  Factor     = grep("^Fact", names(gene_rank_by_pheno), ignore.case = TRUE, value = TRUE)[1],
  Initiation = grep("^Init", names(gene_rank_by_pheno), ignore.case = TRUE, value = TRUE)[1],
  Frequency  = grep("^Freq", names(gene_rank_by_pheno), ignore.case = TRUE, value = TRUE)[1]
)
stopifnot(all(!is.na(pheno_map)))

get_sig_genes <- function(df) {
  stopifnot(all(c("gene_id","bonf_n") %in% names(df)))
  unique(df$gene_id[df$bonf_n >= 1])
}

sig_sets <- list(
  "All Cannabis" = get_sig_genes(prioritised_gene_rank[[pheno_map["Factor"]]]),
  "Ever Use"     = get_sig_genes(prioritised_gene_rank[[pheno_map["Initiation"]]]),
  "All Frequency"= get_sig_genes(prioritised_gene_rank[[pheno_map["Frequency"]]])
)

### venn diagram scaled

library(eulerr)

fit <- euler(sig_sets)

labs <- fit$original   # or fitted(fit) if you prefer the fitted values

# Replace 0-count regions with blank labels
labs[labs == 0] <- ""

# smaller Venn directly from eulerr
venn_small <- as.ggplot(
  plot(
    fit,
    venn = list(
      lwd     = 2,
      lty     = 2,
      col     = "black",
      scaling = 0.5 
    ),
    fills = c("#d95f02", "#66a61e", "#7570b3"),
    alpha = 0.6,
    edges = TRUE,
    labels = list(
      font = 2,
      cex  = 1.2,
      pos  = c(2, 1, 2)
    ),
    quantities = list(
      labels = labs,
      cex    = 1.2
    )
  )
) +
  theme_void() +
  theme(
    # Enough headroom for the top label, but not ridiculous
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# Initialize output list
plot_data_by_pheno <- list()

for (pheno in names(gene_rank_by_pheno_top20)) {
  df <- gene_rank_by_pheno_top20[[pheno]]
  
  # Identify relevant columns for pivoting (no phenotype in name)
  reshape_cols <- grep(
    "^(MAGMA_|mBAT_|SMR_BrainMeta_|SMR_eQTLgen_|TWAS_GTEx_|TWAS_PEC_).*(pvalue|fdr|bonferroni|p_HEIDI|COLOC\\.PP4)$",
    names(df), value = TRUE
  )
  
  if (length(reshape_cols) == 0) {
    warning(paste("No matching columns found for", pheno, "- skipping."))
    next
  }
  
  # Reshape long (only method and statistic)
  df_long <- df %>%
    pivot_longer(
      cols = all_of(reshape_cols),
      names_to = c("method", "statistic"),
      names_pattern = "^(.*)_(pvalue|fdr|bonferroni|p_HEIDI|COLOC\\.PP4)$",
      values_to = "value"
    )
  
  # Pivot wider
  df_wide <- df_long %>%
    pivot_wider(
      names_from = statistic,
      values_from = value
    )
  
  # Add annotations
  df_wide <- df_wide %>%
    mutate(
      sig_label = case_when(
        bonferroni < 0.05 ~ "**",
        fdr < 0.05 ~ "*",
        TRUE ~ ""
      ),
      log10p = -log10(pvalue),
      chromosome_name = factor(chromosome_name, levels = c(as.character(1:22), "X", "Y")),
      hgnc_symbol = fct_reorder2(hgnc_symbol, chromosome_name, start_position),
      method = recode(method,
                      "MAGMA" = "MAGMA",
                      "mBAT" = "mBAT-combo",
                      "SMR_BrainMeta" = "SMR Brain",
                      "SMR_eQTLgen" = "SMR Blood",
                      "TWAS_GTEx" = "TWAS GTEx",
                      "TWAS_PEC" = "TWAS PEC")
    ) %>%
    arrange(as.numeric(chromosome_name), start_position) %>%
    mutate(hgnc_symbol = factor(hgnc_symbol, levels = unique(hgnc_symbol)))
  
  # Save to list
  plot_data_by_pheno[[pheno]] <- df_wide
}

# where hgnc_symbol is missing, replace with gene_id
plot_data_by_pheno <- lapply(plot_data_by_pheno, function(df) {
  df %>%
    mutate(hgnc_symbol = if_else(is.na(hgnc_symbol) | hgnc_symbol == "", gene_id, hgnc_symbol))
})

################################
# PLOT WITH THEMES #
################################

# Define chromosome factor levels for correct numeric ordering
plot_data_by_pheno$Factor <- plot_data_by_pheno$Factor %>%
  mutate(
    chromosome_name = factor(chromosome_name, levels = c(as.character(1:22), "X", "Y")),
    hgnc_symbol = fct_reorder2(hgnc_symbol, chromosome_name, start_position)
  )

# Custom theme (as before)
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks = element_blank(),
    panel.spacing.x = unit(0.2, "lines"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()
  )

# Subset data to Factor phenotype
plot_data_by_pheno <- lapply(plot_data_by_pheno, function(df) {
  subset(df, !is.na(pvalue))
})

# Identify rows meeting COLOC or HEIDI thresholds
highlight_df <- list()
for(i in 1:length(plot_data_by_pheno)) {
highlight_df[[i]] <- plot_data_by_pheno[[i]] %>%
  filter(
    (grepl("TWAS", method) & COLOC.PP4 > 0.8 & bonferroni < 0.05) |
      (grepl("SMR", method)  & p_HEIDI > 0.05 & bonferroni < 0.05)
  ) %>%
  mutate(
    border_color = "black",
    border_linetype = case_when(
      grepl("TWAS", method) ~ "dashed",
      grepl("SMR", method)  ~ "dashed"
    )
  )
}
names(highlight_df) <- phenotypes

##########################################
# PLOT DATA WITH SINGLE GREYSCALE LEGEND #
##########################################

all_genes_df <- bind_rows(plot_data_by_pheno, .id = "phenotype") %>%
  distinct(hgnc_symbol, chromosome_name, start_position)

# Order by chromosome and position
gene_order <- all_genes_df %>%
  arrange(
    as.numeric(chromosome_name),
    start_position
  ) %>%
  pull(hgnc_symbol)


# Step 1: Your three colored plots (no legend shown)
p1 <- ggplot(plot_data_by_pheno$Factor, aes(x = method, y = hgnc_symbol, fill = log10p)) +
  geom_tile(color = NA) +
  geom_tile(data = highlight_df$Factor, aes(color = border_color, linetype = border_linetype), fill = NA, size = 0.6) +
  geom_text(aes(label = sig_label), color = "black", size = 4.5, vjust = 0.8) +
  scale_y_discrete(limits = rev(gene_order)) +
  scale_fill_gradientn(
    colours = c("white", "#d95f02", "black"),
    values = rescale(c(0, 15, 30)),
    limits = c(0, 30),
    na.value = "white",
    oob = squish,
    guide = "none"
  ) +
  scale_color_identity() +
  scale_linetype_identity() +
  labs(title = "All Cannabis", x = NULL, y = "") +
  custom_theme

p2 <- ggplot(plot_data_by_pheno$Initiation, aes(x = method, y = hgnc_symbol, fill = log10p)) +
  geom_tile(color = NA) +
  geom_tile(data = highlight_df$Initiation, aes(color = border_color, linetype = border_linetype), fill = NA, size = 0.6) +
  geom_text(aes(label = sig_label), color = "black", size = 4.5, vjust = 0.8) +
  scale_y_discrete(limits = rev(gene_order)) +
  scale_fill_gradientn(
    colours = c("white", "#66a61e", "black"),
    values = rescale(c(0, 15, 30)),
    limits = c(0, 30),
    na.value = "white",
    oob = squish,
    guide = "none"
  ) +
  scale_color_identity() +
  scale_linetype_identity() +
  labs(title = "Ever Use", x = NULL, y = "") +
  custom_theme + 
  theme(axis.text.y = element_blank(),
                      axis.ticks.y = element_blank())

p3 <- ggplot(plot_data_by_pheno$Frequency, aes(x = method, y = hgnc_symbol, fill = log10p)) +
  geom_tile(color = NA) +
  geom_tile(data = highlight_df$Frequency, aes(color = border_color, linetype = border_linetype), fill = NA, size = 0.6) +
  geom_text(aes(label = sig_label), color = "black", size = 4.5, vjust = 0.8) +
  scale_y_discrete(limits = rev(gene_order)) +
  scale_fill_gradientn(
    colours = c("white", "#7570b3", "black"),
    values = rescale(c(0, 15, 30)),
    limits = c(0, 30),
    na.value = "white",
    oob = squish,
    guide = "none"
  ) +
  scale_color_identity() +
  scale_linetype_identity() +
  labs(title = "All Frequency", x = NULL, y = "") +
  custom_theme + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())


# create shared legend plot
legend_plot <- ggplot(data.frame(x = 1, y = 1, log10p = 15), aes(x, y, fill = log10p)) +
  geom_tile(alpha = 0, size = 0) +  # invisible tile
  scale_fill_gradientn(
    colours = c("white", "grey50", "black"),
    values = rescale(c(0, 15, 30)),
    limits = c(0, 30),
    name = expression(-log[10](P)),
    oob = squish,
    na.value = "white"
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.key.width = unit(2.8, "cm"),
    legend.key.height = unit(0.75, "cm"),
    legend.title = element_text(size = 14, hjust = 0, vjust = 0.85),
    legend.text = element_text(size = 12),
    legend.box.margin = margin(0, 0, 0, 10),
    plot.margin = margin(0, 0, 0, 0)
  )

# combine all plots
final_plot <- (p1 + p2 + p3 + wrap_elements(legend_plot)) +
  plot_layout(ncol = 3, guides = "collect", heights = c(10, 1)) &
  theme(
    legend.position = "bottom",
    legend.justification = "right",
    plot.title = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 5)
  )

# format shared legend
legend_shared <- cowplot::get_legend(
  p1 +
    scale_fill_gradientn(
      colours = c("white", "grey60", "black"),
      limits = c(0, 30),
      name = expression(-log[10](P)),
      oob = scales::squish,
      na.value = "white"
    ) +
    guides(
      fill = guide_colourbar(
        title.theme = element_text(size = 14, face = "bold"),
        label.theme = element_text(size = 12),
        barwidth = unit(0.8, "cm"),
        barheight = unit(5, "cm"),
        label.position = "left"
      ),
      color = "none",
      linetype = "none"
    ) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.key.height = unit(3, "cm"),
      legend.key.width = unit(1, "cm"),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(-30, 0, 0, 0)
    )
)

# remove legends from all plots and center titles
p1 <- p1 + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

p2 <- p2 + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

p3 <- p3 + 
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# align top row plots (horizontal alignment, left axis)
aligned <- align_plots(p1, p2, p3, align = "h", axis = "l")
main_plot <- plot_grid(plotlist = aligned, ncol = 3, rel_widths = c(1, 1, 1))

combined_plot <- plot_grid(
  main_plot,
  legend_shared,
  rel_widths = c(1, 0.1),
  ncol = 2
)

# print
combined_plot

# remove margins around each plot
p1 <- p1 + theme(plot.margin = margin(t=0, r=0, b=0, l=0))
p2 <- p2 + theme(plot.margin = margin(t=0, r=0, b=0, l=0))
p3 <- p3 + theme(plot.margin = margin(t=0, r=0, b=0, l=0))

# align top row plots (horizontal alignment, left axis)
aligned_top <- align_plots(p1, p2, p3, align = "h", axis = "l")

##############################################x
# ADD GSA MiXER plots                        #
##############################################

setwd("/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/genesets/gsamixer")

# Get all supplementary tables in the current directory
excel_files <- list.files(pattern = "*_SupplementaryTables\\.xlsx$")

# Read each Excel file and all its sheets into a nested list
nested_tables <- map(excel_files, function(file) {
  sheets <- excel_sheets(file)
  sheet_data <- map(sheets, ~ read_excel(file, sheet = .x))
  names(sheet_data) <- sheets
  sheet_data
})

# Name the top-level elements after the Excel filenames
names(nested_tables) <- c("FACTOR","EVR","FREQ_MA")

# collapse data frames
main_res <- map2_dfr(
  nested_tables,
  names(nested_tables),
  ~ .x[["ST5"]] %>%
    mutate(across(everything(), as.character)) %>%
    mutate(PHENOTYPE = .y)
) %>%
  arrange(as.numeric(MAGMA_P))


exploratory_res <- map2_dfr(
  nested_tables,
  names(nested_tables),
  ~ .x[["ST7"]] %>%
    mutate(across(everything(), as.character)) %>%
    mutate(PHENOTYPE = .y)
) %>%
  arrange(as.numeric(MAGMA_FDR))


#####################
# plot main results #
#####################

# Prepare and clean data
collapsed_st5 <- main_res %>%
  mutate(
    enrich = as.numeric(enrich),
    se_enrich = as.numeric(se_enrich),
    MAGMA_P = as.numeric(MAGMA_P),
    MAGMA_FDR = as.numeric(MAGMA_FDR)
  )

# keep highest enrichment per GO
collapsed_st5_filtered <- collapsed_st5 %>%
  group_by(GO) %>%
  slice_max(order_by = enrich, n = 1, with_ties = FALSE) %>%
  ungroup()

# define common GO term ordering (descending by enrichment)
go_levels <- collapsed_st5_filtered %>%
  arrange(enrich) %>%
  pull(GO)

# reverse y-axis order
collapsed_st5_filtered <- collapsed_st5_filtered %>%
  mutate(
    GO = factor(GO, levels = go_levels),
    log10_MAGMA_P = -log10(MAGMA_P)
  )

# Define FDR threshold line for MAGMA_P
fdr_thresh <- max(collapsed_st5_filtered$MAGMA_P[collapsed_st5_filtered$MAGMA_FDR < 0.05], na.rm = TRUE)
fdr_cutoff <- ifelse(is.finite(fdr_thresh), -log10(fdr_thresh), NA)

# Define common theme
clean_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(b=0, l=0, t=10, r=0),
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14)
  )

# Plot 1: GSA-MiXeR enrichment with error bars
collapsed_st5_filtered <- collapsed_st5_filtered %>%
  mutate(
    PHENOTYPE_CLEAN = case_when(
      grepl("FACTOR", PHENOTYPE) ~ "All Cannabis",
      grepl("EVR", PHENOTYPE) ~ "Ever Use",
      grepl("FREQ_MA", PHENOTYPE) ~ "All Frequency",
      TRUE ~ PHENOTYPE
    )
  )

# order names
collapsed_st5_filtered$PHENOTYPE_CLEAN <- factor(
  collapsed_st5_filtered$PHENOTYPE_CLEAN,
  levels = c("All Frequency", "Ever Use", "All Cannabis")
)


phenotype_colors <- c(
  "All Cannabis" = "#d95f02",
  "Ever Use" = "#66a61e",
  "All Frequency" = "#7570b3"
)

# clean pathway title string
collapsed_st5_filtered <- collapsed_st5_filtered %>%
  mutate(
    GO = GO %>%
      str_replace_all(
        fixed("GOBP_NEGATIVE_REGULATION_OF_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY"),
        "GOBP_ADENYLATE_CYCLASE_ACTIVATING_GPCR_PATHWAY"
      ) %>%
      str_replace_all(
        regex("GOBP_REGULATION_OF_ADENYLATE_CYCLASE_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY", ignore_case = TRUE),
        ""
      ) %>%
      str_remove(regex("^(?:GO:(?:BP|MF|CC)|GOBP|GOMF|GOCC|SYNGO):?\\s*[_\\-]*\\s*", ignore_case = TRUE)) %>%
      str_remove("^.*?:\\s*") %>%
      str_replace_all("_", " ") %>%
      str_replace_all(
        regex("ADENYLATE CYCLASE ACTIVATING G\\s*PROTEIN COUPLED RECEPTOR", ignore_case = TRUE),
        "adenylate cyclase-activating GPCR"
      ) %>%
      str_to_lower() %>%
      str_replace("^([a-z])", ~ toupper(.x)) %>%
      str_squish() %>%
      str_replace_all(regex("\\bgpcr\\b", ignore_case = TRUE), "GPCR") %>%
      sapply(function(x) if (nchar(x) > 60) str_wrap(x, width = 60) else x) %>%
      as.character()
  ) %>%
  filter(GO != "")


p4 <- ggplot(
  collapsed_st5_filtered %>%
    mutate(GO = fct_reorder(GO, enrich, .desc = FALSE)),  # reorder by enrichment
  aes(x = enrich, y = GO)
) +
  geom_point(aes(color = PHENOTYPE_CLEAN), size = 2, show.legend = TRUE) +
  geom_errorbarh(aes(xmin = enrich - se_enrich, xmax = enrich + se_enrich), 
                 height = 0.2, show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0, 4), breaks = seq(0, 4, 1), expand = c(0, 0)) +
  scale_color_manual(
    name = "Cannabis use\nphenotype",
    values = phenotype_colors,
    breaks = c("All Cannabis","Ever Use", "All Frequency")
  ) +
  labs(x = "Fold Enrichment", y = "", title = "") +
  clean_theme +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(colour = "grey90", size = 0.4),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2)
  )



p5 <- ggplot(
  collapsed_st5_filtered %>%
    mutate(GO = fct_reorder(GO, enrich, .desc = FALSE)),  # reorder by enrichment
  aes(x = log10_MAGMA_P, y = GO, color = PHENOTYPE_CLEAN)
) + 
  geom_point(size = 2) +
  geom_vline(xintercept = fdr_cutoff, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_color_manual(
    name = "Cannabis use\nphenotype",
    values = phenotype_colors,
    breaks = c("All Cannabis","Ever Use", "All Frequency")
  ) +
  labs(x = expression(-log[10](MAGMA~italic(P))), y = NULL, title = "MAGMA Gene Set Association") +
  clean_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

(p4 + plot_spacer() + p5) +
  plot_layout(ncol = 3, widths = c(1, 0.1, 1)) +
  plot_annotation() &
  theme(legend.position = "right")


right_col <- (p4 + plot_spacer() + p5) +
  plot_layout(ncol = 3, widths = c(1, 0.1, 1)) +
  plot_annotation(tag_levels = NULL) &
  theme(
    legend.position = "right",
    plot.margin = margin(t = 10, r = 5, b = 5, l = 5)
  )

# Assemble top and bottom rows with specified relative widths
aplot <- venn_small
bplot <- plot_grid(aligned_top[[1]], aligned_top[[2]], aligned_top[[3]], ncol = 3, rel_widths = c(1.2, 0.95, 0.95))

main_plot <- plot_grid(
  aplot,
  bplot,
  legend_shared,
  ncol = 3,
  rel_widths = c(0.7, 1.2, 0.2),
  labels = c("A", "B"),
  label_size = 20,
  label_x = 0,
  label_y = 1,
  rel_heights = c(0.7, 1),  # or adjust as needed
  align = "v",
  axis = "lr"
)

## add partitioned heritability results ##

# definite results path
res_path = "/stornext/Bioinf/data/lab_bahlo/users/gerring.z/projects/icc4/results/partitioned_heritability/"
out_path = "/vast/scratch/users/gerring.z/icc4/results/partitioned_heritability/cleaned/"

## factor

# import results
icc4_factor = list.files(res_path,"^MA_CANN_FACTOR_250919_.*\\.results$")
icc4_factor_ls = lapply(paste0(res_path,icc4_factor),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("MA_CANN_FACTOR_250919_|\\.results", "", x)
}
names(icc4_factor_ls) = lapply(icc4_factor, clean_string)

# collapse results to data frame
icc4_factor_df = do.call(rbind.data.frame,icc4_factor_ls)
icc4_factor_df$cell_type = sapply(strsplit(rownames(icc4_factor_df),"\\."),"[[",1)
rownames(icc4_factor_df) = NULL
icc4_factor_df = icc4_factor_df[order(icc4_factor_df$Enrichment_p),]

## export cell type results for each dataset
icc4_factor_df$Category = gsub("L2_0","",icc4_factor_df$Category)
icc4_factor_celltype = subset(icc4_factor_df,Category=="L2_1")

# calcualte fdr and bonferroni
icc4_factor_celltype$fdr = p.adjust(icc4_factor_celltype$Enrichment_p,"fdr")
icc4_factor_celltype$bonferroni = p.adjust(icc4_factor_celltype$Enrichment_p,"bonferroni")

## initiation

# import results
icc4_initiation = list.files(res_path,"^MA_CANNABIS_EVR_250919_.*\\.results$")
icc4_initiation_ls = lapply(paste0(res_path,icc4_initiation),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("MA_CANNABIS_EVR_250919_|\\.results", "", x)
}
names(icc4_initiation_ls) = lapply(icc4_initiation, clean_string)

# collapse results to data frame
icc4_initiation_df = do.call(rbind.data.frame,icc4_initiation_ls)
icc4_initiation_df$cell_type = sapply(strsplit(rownames(icc4_initiation_df),"\\."),"[[",1)
rownames(icc4_initiation_df) = NULL
icc4_initiation_df = icc4_initiation_df[order(icc4_initiation_df$Enrichment_p),]

## export cell type results for each dataset
icc4_initiation_df$Category = gsub("L2_0","",icc4_initiation_df$Category)
icc4_initiation_celltype = subset(icc4_initiation_df,Category=="L2_1")

# calcualte fdr and bonferroni
icc4_initiation_celltype$fdr = p.adjust(icc4_initiation_celltype$Enrichment_p,"fdr")
icc4_initiation_celltype$bonferroni = p.adjust(icc4_initiation_celltype$Enrichment_p,"bonferroni")

## frequency

# import results
icc4_frequency = list.files(res_path,"^MA_FREQ_ALL_250919_.*\\.results$")
icc4_frequency_ls = lapply(paste0(res_path,icc4_frequency),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("MA_FREQ_ALL_250919_|\\.results", "", x)
}
names(icc4_frequency_ls) = lapply(icc4_frequency, clean_string)

# collapse results to data frame
icc4_frequency_df = do.call(rbind.data.frame,icc4_frequency_ls)
icc4_frequency_df$cell_type = sapply(strsplit(rownames(icc4_frequency_df),"\\."),"[[",1)
rownames(icc4_frequency_df) = NULL
icc4_frequency_df = icc4_frequency_df[order(icc4_frequency_df$Enrichment_p),]

## export cell type results for each dataset
icc4_frequency_df$Category = gsub("L2_0","",icc4_frequency_df$Category)
icc4_frequency_celltype = subset(icc4_frequency_df,Category=="L2_1")

# calcualte fdr and bonferroni
icc4_frequency_celltype$fdr = p.adjust(icc4_frequency_celltype$Enrichment_p,"fdr")
icc4_frequency_celltype$bonferroni = p.adjust(icc4_frequency_celltype$Enrichment_p,"bonferroni")

# create list of phenotypes
pldsc_res = list(icc4_initiation_celltype,icc4_frequency_celltype,icc4_factor_celltype)
names(pldsc_res) = c("Initiation","Frequency","Factor")

# merge results (dataset, tissue, enrichment, enrichment p value) and reshape for heatmap type analysis
pldsc_res_df = do.call(rbind.data.frame,pldsc_res)
pldsc_res_df$dataset = sapply(strsplit(rownames(pldsc_res_df),"\\."),"[[",1)
rownames(pldsc_res_df) = NULL

# plot data
plot_dat = pldsc_res_df[,c("dataset","cell_type","Enrichment","Enrichment_p","fdr","bonferroni")]

# define cell type names and dataset
plot_dat$cell_type = gsub("_"," ",plot_dat$cell_type)

# control dataset order
plot_dat$dataset <- factor(plot_dat$dataset,
                           levels = c("Initiation", "Frequency", "Factor"))

# create a categorical variable for significance 
plot_dat$significance = ifelse(plot_dat$bonferroni<0.05 & plot_dat$Enrichment>1, "Enriched, Bonferroni", 
                               ifelse(plot_dat$fdr<0.05 & plot_dat$Enrichment>1, "Enriched, FDR", 
                                      ifelse(plot_dat$bonferroni<0.05 & plot_dat$Enrichment<1, "Depleted, Bonferroni", 
                                             ifelse(plot_dat$fdr<0.05 & plot_dat$bonferroni>0.05 & plot_dat$Enrichment<1, "Depleted, FDR", "None"))))

#plot_dat$significance = ifelse(plot_dat$fdr<0.05 & plot_dat$Enrichment>1, "Enriched", 
#                               ifelse(plot_dat$fdr<0.05 & plot_dat$Enrichment<1, "Depleted", "None"))

# recode categorical variable
plot_dat$significance_ordered <- factor(plot_dat$significance, levels=c("Enriched, Bonferroni", "Enriched, FDR", "Depleted, Bonferroni","Depleted, FDR", "None"), ordered=TRUE)

# plot results on continuous scale
p6 <- ggplot(plot_dat, aes(x = cell_type, y = dataset)) +
  geom_tile(aes(fill = Enrichment), color = "black", size = 0.3, width = 0.9, height = 0.9) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1) +
  coord_flip() +
  scale_y_discrete(position = "right")  +
  scale_x_discrete(limits=rev) +
  labs(fill = expression("Partitioned h"^2*SNP)) +
  geom_text(data = plot_dat[plot_dat$fdr < 0.05, ], 
            aes(label = "*"), 
            color = "black", 
            size = 8,
            vjust= 0.75) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(size = 12, color="black", angle=45, hjust=0.02),
    axis.text.y = element_text(size = 12, color="black", hjust = 1),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12))

# master combine

bottom_row <- plot_grid(
  p4 + theme(plot.margin = margin(0, -10, 0, 0)),
  p6 + theme(plot.margin = margin(0, 0, -10, 0)),
  ncol = 2,
  rel_widths = c(1.2, 0.8),
  labels = c("C", "D"),
  label_size = 20,
  label_x = 0,
  label_y = 1,
  align = "hv",   # tightly align both horizontally and vertically
  axis = "tblr"   # align all axes (top, bottom, left, right)
)

# combine with gsa mixer results
stacked_plot <- plot_grid(
  main_plot,
  bottom_row,
  #ncol = 1,
  ncol = 1,
  rel_heights = c(1, 1.2),
  rel_widths = c(1, 1),
  labels = c(" ", " "),
  label_size = 24,
  label_fontface = "bold",
  label_x = 0,
  label_y = 1
)








