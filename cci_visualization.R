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
suppressMessages(library(readxl))
suppressMessages(library(scales))
suppressMessages(library(ComplexHeatmap))

input_dir <- '/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/6_cci_results/'
cci_xlsx <- 'selected_cci_v1.xlsx'

cci_df <- read_excel(paste(input_dir, cci_xlsx, sep = ""))
cci_df$col_ligand <- str_c("<span style='color:firebrick'>", cci_df$ligand, "</span>")
cci_df$col_receptor <- str_c("<span style='color:dodgerblue'>", cci_df$receptor, "</span>")
cci_df$col_lr <- str_c(cci_df$col_ligand, "-", cci_df$col_receptor)

cci_df$col_ligandc <- str_c("<span style='color:firebrick4'>", cci_df$ligand_cell, "</span>")
cci_df$col_receptorc <- str_c("<span style='color:dodgerblue4'>", cci_df$receptor_cell, "</span>")
cci_df$col_lrc <- str_c(cci_df$col_ligandc, "-", cci_df$col_receptorc)

cci_df$order_col <- ifelse(str_detect(cci_df$ligand_cell, "TAM"), cci_df$receptor_cell, cci_df$ligand_cell)


tmp_cci_df <- cci_df[cci_df$pipeline=='CellPhoneDB',]
tmp_cci_df$raw_value <- as.numeric(tmp_cci_df$value)
tmp_cci_df$value <- rescale(tmp_cci_df$raw_value, to = c(0,100))
facet_mean_dot <- ggplot(tmp_cci_df, aes(x = reorder(col_lrc, order_col), y = col_lr, color = value)) +
  geom_point(size=6) +
  facet_wrap(~order_col, scales = 'free', nrow = 1) +
  scale_color_viridis_c(option = 'plasma') +
  labs(color = 'Significant Means\nCellPhoneDB(v2)', 
       y = "<span style='color:firebrick'>Ligand</span>-<span style='color:dodgerblue'>Receptor</span> Pairs",
       x = "<span style='color:firebrick4'>Cell Expressing Ligand</span>-<span style='color:dodgerblue4'>Cell Expressing Receptor</span>") +
  #	scale_color_distiller(palette = 'Spectral') +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10),#element_text(angle = 45, hjust = 1),
        axis.text.y = element_markdown(size = 10),
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12))
ggsave(paste(input_dir, 'cpdb_compare_significant_interaction_mean_facet.png', sep = ''), dpi = 300, width = 12, height = 4)

tmp_cci_df <- cci_df[cci_df$pipeline=='ICELLNET',]
tmp_cci_df$value <- as.numeric(tmp_cci_df$value)
facet_mean_dot <- ggplot(tmp_cci_df, aes(x = reorder(col_lrc, order_col), y = col_lr, color = value)) +
  geom_point(size=6) +
  facet_wrap(~order_col, scales = 'free', nrow = 1) +
  scale_color_viridis_c(option = 'plasma') +
  labs(color = 'Comm. Scores\nICELLNET', 
       y = "<span style='color:firebrick'>Ligand</span>-<span style='color:dodgerblue'>Receptor</span> Pairs",
       x = "<span style='color:firebrick4'>Cell Expressing Ligand</span>-<span style='color:dodgerblue4'>Cell Expressing Receptor</span>") +
  #	scale_color_distiller(palette = 'Spectral') +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1, size = 10),#element_text(angle = 45, hjust = 1),
        axis.text.y = element_markdown(size = 10),
        axis.title.x = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12))
ggsave(paste(input_dir, 'icn_compare_significant_interaction_mean_facet.png', sep = ''), dpi = 300, width = 16, height = 9)
