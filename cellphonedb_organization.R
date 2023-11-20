rm(list = ls())
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggtext))
suppressMessages(library(gridtext))
suppressMessages(library(viridis))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(circlize))
suppressMessages(library(ComplexHeatmap))

cpdb_out_dir <- '/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/6_cci_results/'
mean_df <- read.table(paste(cpdb_out_dir, 'means.txt', sep = ''), sep = '\t', check.names = FALSE, header = TRUE)
sig_mean_df <- read.table(paste(cpdb_out_dir, 'significant_means.txt', sep = ''), sep = '\t', check.names = FALSE, header = TRUE)
pval_df <- read.table(paste(cpdb_out_dir, 'pvalues.txt', sep = ''), sep = '\t', check.names = FALSE, header = TRUE)
print(dim(mean_df))
print(dim(pval_df))
print(dim(sig_mean_df))
mean_df$uid <- str_c(mean_df[,1], 'xxx', mean_df[,2])
pval_df$uid <- str_c(pval_df[,1], 'xxx', pval_df[,2])
print(table(pval_df$receptor_a, pval_df$receptor_b))
sig_mean_df$uid <- str_c(sig_mean_df[,1], 'xxx', sig_mean_df[,2])

pair_col <- colnames(mean_df)[str_detect(colnames(mean_df), '\\|')]

gath_mean_df <- gather(mean_df, 'cell_pairs', 'mean', pair_col)
rownames(gath_mean_df) <- str_c(gath_mean_df$uid, '+++', gath_mean_df$cell_pairs)
gath_mean_df$cell_1 <- str_split_fixed(gath_mean_df$cell_pairs, '\\|', n = 2)[,1]
gath_mean_df$cell_2 <- str_split_fixed(gath_mean_df$cell_pairs, '\\|', n = 2)[,2]

gath_pval_df <- gather(pval_df, 'cell_pairs', 'pval', pair_col)
rownames(gath_pval_df) <- str_c(gath_pval_df$uid, '+++', gath_pval_df$cell_pairs)

gath_sigavg_df <- gather(sig_mean_df, 'cell_pairs', 'sig_means', pair_col)
rownames(gath_sigavg_df) <- str_c(gath_sigavg_df$uid, '+++', gath_sigavg_df$cell_pairs)

clean_sigavg_df <- gath_sigavg_df[!is.na(gath_sigavg_df$sig_means),]
clean_sigavg_df <- clean_sigavg_df[str_detect(clean_sigavg_df$cell_pairs, 'TAM'),]
self_inter_pairs <- c("PDL1+TAM|PDL1+TAM", "PDL1-TAM|PDL1-TAM", "PDL1+TAM|PDL1-TAM", "PDL1-TAM|PDL1+TAM")
clean_sigavg_df<- clean_sigavg_df[!str_detect(clean_sigavg_df$cell_pairs, 'GC'),]
clean_sigavg_df <- clean_sigavg_df[!str_detect(clean_sigavg_df$cell_pairs, 'NonTAM'),]
clean_sigavg_df <- clean_sigavg_df[str_detect(clean_sigavg_df$cell_pairs, 'TAM'),]
clean_sigavg_df <- clean_sigavg_df[!(clean_sigavg_df$cell_pairs %in% self_inter_pairs),]
clean_sigavg_df[,c('cell_1', 'cell_2')] <- str_split_fixed(clean_sigavg_df$cell_pairs, '\\|', n=2)
clean_sigavg_df$receptor_cell <- ifelse(clean_sigavg_df$receptor_a == 'True', clean_sigavg_df$cell_1, clean_sigavg_df$cell_2)
clean_sigavg_df$ligand_cell <- ifelse(clean_sigavg_df$receptor_a == 'False', clean_sigavg_df$cell_1, clean_sigavg_df$cell_2)
clean_sigavg_df$TAM <- ifelse(str_detect(clean_sigavg_df$cell_pairs, 'PDL1\\+'), 'PDL1+', 'PDL1-')
clean_sigavg_df$PAIR <- ifelse(str_detect(clean_sigavg_df$cell_pairs, '\\|PDL1'), clean_sigavg_df$cell_1, clean_sigavg_df$cell_2)
clean_sigavg_df$TAM_EXPR <- ifelse(str_detect(clean_sigavg_df$receptor_cell, 'TAM'), 'Receptor', 'Ligand')
clean_sigavg_df$pval <- gath_pval_df[rownames(clean_sigavg_df), 'pval']
clean_sigavg_df$pflag <- ifelse(clean_sigavg_df$pval <= 0.01, 'p<=0.01', ifelse(clean_sigavg_df$pval <= 0.05, '0.01<p<=0.05', 'ns'))
clean_sigavg_df[,c('ip_a', 'ip_b')] <- str_split_fixed(clean_sigavg_df$interacting_pair, '_', n = 2)
clean_sigavg_df$col_gene_a <- ifelse(clean_sigavg_df$receptor_a == 'True', 
				     str_c("<span style='color:dodgerblue'>", clean_sigavg_df$ip_a,"</span>"),
				     str_c("<span style='color:firebrick'>", clean_sigavg_df$ip_a, "</span>"))
clean_sigavg_df$col_gene_b <- ifelse(clean_sigavg_df$receptor_b == 'True', 
				     str_c("<span style='color:dodgerblue'>", clean_sigavg_df$ip_b, "</span>"),
				     str_c("<span style='color:firebrick'>", clean_sigavg_df$ip_b, "</span>"))
clean_sigavg_df$col_inter_pairs <- str_c(clean_sigavg_df$col_gene_a, "-", clean_sigavg_df$col_gene_b)
print(head(clean_sigavg_df))

clean_spr_sigavg_df <- spread(clean_sigavg_df[,c('TAM', 'sig_means', 'PAIR', 'interacting_pair', 'TAM_EXPR', 'secreted')], 
			      'TAM', 'sig_means')
write.csv(clean_spr_sigavg_df, paste(cpdb_out_dir, 'sig_means_tam_df.csv', sep = ''))

