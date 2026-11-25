rm(list=ls())
library(qvalue)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# definite results path
res_path = "/home/users/allstaff/gerring.z/projects/icc4/results/pldsc/"
out_path = "/home/users/allstaff/gerring.z/projects/icc4/results/export/"

## initiation

# import results
icc4_initiation = list.files(res_path,"^ICC4_CANNABIS_EVR_240812_.*\\.results$")
icc4_initiation_ls = lapply(paste0(res_path,icc4_initiation),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("ICC4_CANNABIS_EVR_240812_|\\.results", "", x)
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

# export
write.csv(icc4_initiation_celltype,file=paste0(out_path,"initation_pldsc.csv"),row.names=F)

## frequency

# import results
icc4_frequency = list.files(res_path,"^ICC4_FREQ_META_241127_.*\\.results$")
icc4_frequency_ls = lapply(paste0(res_path,icc4_frequency),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("ICC4_FREQ_META_241127_|\\.results", "", x)
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

# export
write.csv(icc4_frequency_celltype,file=paste0(out_path,"frequency_pldsc.csv"),row.names=F)

## factor

# import results
icc4_factor = list.files(res_path,"^icc4_cann_factor_.*\\.results$")
icc4_factor_ls = lapply(paste0(res_path,icc4_factor),read.delim)

# define function to clean string and name tissues/cell types
clean_string = function(x) {
  gsub("icc4_cann_factor_|\\.results", "", x)
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

# export
write.csv(icc4_factor_celltype,file=paste0(out_path,"factor_pldsc.csv"),row.names=F)

####

###

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
ggplot(plot_dat, aes(x = cell_type, y = dataset)) +
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
    axis.ticks.y = element_blank())


# Create the plot - discrete
#ggplot(plot_dat, aes(x = cell_type, y = dataset)) +
#  geom_tile(aes(fill = significance_ordered), color = "black", size = 0.3, width = 0.9, height = 0.9) +
#  scale_fill_manual(values = c("#0072B2","#56B4E9","#D55E00","white")) +
#  coord_flip() +
    #  scale_y_discrete(position = "right")  +
    #  scale_x_discrete(limits=rev) +
    #  labs(fill = expression("Partitioned h"^2*SNP)) +
    #  theme(
    #    axis.title.x = element_blank(),
    #    axis.title.y = element_blank(),
    #    panel.background = element_rect(fill = "white", color = NA),  
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(), 
    #    axis.text.x = element_text(size = 12, color="black", angle=45, hjust=0.02),
    #    axis.text.y = element_text(size = 12, color="black", hjust = 1),
    #    axis.ticks.x = element_blank(),
    #    axis.line = element_blank(),
    #    axis.ticks.y = element_blank())


