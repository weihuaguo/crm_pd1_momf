rm(list = ls())
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(ggVennDiagram))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))


input_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub/3_signature"
sig_file <- "organized_tam_signatures_231025.xlsx"

sig_df <- read_excel(paste(input_dir, sig_file, sep = "/"))
sig_df$uid <- str_c(sig_df$tam_type, " from ", sig_df$resource)
print(head(sig_df))

cts_df <- sig_df %>%
	group_by(uid, gene) %>%
	summarize(n = n())
print(as.data.frame(cts_df))

cat("Compare validation cohorts\n")
pdl1_df <- sig_df[sig_df$resource != "Collected" & sig_df$resource != "in-house signature",]
val_df <- pdl1_df[pdl1_df$tam_type == "PD-L1+",]
y <- list()
for (irs in unique(val_df$resource)) {
	y[[irs]] <- val_df$gene[val_df$resource == irs]
}

vdgg <-ggVennDiagram(y, label_alpha = 0) +
	scale_fill_distiller(palette = "YlOrRd", direction = 1) +
	scale_color_brewer(palette = "Dark2") +
	labs(title = "PD-L1+ validation", fill = "Gene number")
ggsave(paste(input_dir, "/pdl1_pos_validation_venn_diagram.png", sep = ""), vdgg, dpi = 300, width = 9, height = 6)

val_df <- pdl1_df[pdl1_df$tam_type == "PD-L1-",]
y <- list()
for (irs in unique(val_df$resource)) {
	y[[irs]] <- val_df$gene[val_df$resource == irs]
}

vdgg <-ggVennDiagram(y, label_alpha = 0) +
	scale_fill_distiller(palette = "YlOrRd", direction = 1) +
	scale_color_brewer(palette = "Dark2") +
	labs(title = "PD-L1- validation", fill = "Gene number")
ggsave(paste(input_dir, "/pdl1_neg_validation_venn_diagram.png", sep = ""), vdgg, dpi = 300, width = 9, height = 6)

cat("Compare PD-L1+/- vs M1/2\n")
all_rs <- unique(sig_df$resource)
all_rs <- all_rs[all_rs != "Collected"]
for (irs in all_rs) {
	cat(irs, "\n")
	tmp_sig_df <- sig_df[sig_df$resource == irs | sig_df$resource == "Collected",]
#	print(head(tmp_sig_df))
	x <- list()
	x[['M1']] <- tmp_sig_df$gene[tmp_sig_df$tam_type == "M1"]
	x[['M2']] <- tmp_sig_df$gene[tmp_sig_df$tam_type == "M2"]
	x[['PD-L1+']] <- tmp_sig_df$gene[tmp_sig_df$tam_type == "PD-L1+"]
	x[['PD-L1-']] <- tmp_sig_df$gene[tmp_sig_df$tam_type == "PD-L1-"]


	vdgg <-ggVennDiagram(x, label_alpha = 0) +
		scale_fill_distiller(palette = "YlOrRd", direction = 1) +
		scale_color_manual(values = c("M1" = "gold1", "M2" = "mediumpurple1", "PD-L1-" = "mediumpurple4", "PD-L1+" = "gold4")) +
		scale_color_brewer(palette = "Paired") +
		labs(title = irs, fill = "Gene number")
	ggsave(paste(input_dir, "/", irs, "_to_m12_venn_diagram.png", sep = ""), vdgg, dpi = 300, width = 9, height = 6)
#	print(x)

	spr_df <- spread(tmp_sig_df[,c("tam_type", "gene", "weight")], "tam_type", "weight")
	spr_df[is.na(spr_df)] <- 0.0
	print(colnames(spr_df))

	hm_df <- as.data.frame(spr_df)
	rownames(hm_df) <- hm_df$gene
	hm_df$gene <- NULL

	hm <- Heatmap(t(hm_df),
		      name = "Weight",
		      column_title = irs,
		      cluster_rows = TRUE,
		      cluster_columns = TRUE,
		      column_names_gp = gpar(fontsize = 7), 
		      row_names_gp = gpar(fontsize = 11),
		      show_column_names = TRUE,
		      heatmap_legend_param = list(direction = "horizontal", 
						  labels_gp = gpar(fontsize = 12),
						  title_gp = gpar(fontsize = 16),
						  legend_width = unit(3.6, "cm")))
	png(paste(input_dir, "/", irs, "_signature_weight_clustered_heatmap.png", sep = ""), res = 300, units = 'in', width = 2+0.1*length(unique(tmp_sig_df$gene)), height = 3)
	draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	hm <- Heatmap(t(hm_df),
		      name = "Weight",
		      column_title = irs,
		      cluster_rows = FALSE,
		      cluster_columns = TRUE,
		      column_names_gp = gpar(fontsize = 7), 
		      row_names_gp = gpar(fontsize = 11),
		      show_column_names = TRUE,
		      heatmap_legend_param = list(direction = "horizontal", 
						  labels_gp = gpar(fontsize = 12),
						  title_gp = gpar(fontsize = 16),
						  legend_width = unit(3.6, "cm")))
	png(paste(input_dir, "/", irs,  "_signature_weight_ordered_heatmap.png", sep = ""), res = 300, units = 'in', width = 2+0.1*length(unique(tmp_sig_df$gene)), height = 3)
	draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()
}

cat("Heatmap\n")
spr_df <- spread(sig_df[,c("uid", "gene", "weight")], "uid", "weight")
spr_df[is.na(spr_df)] <- 0.0
print(colnames(spr_df))

hm_df <- as.data.frame(spr_df)
rownames(hm_df) <- hm_df$gene
hm_df$gene <- NULL

col_ann <- sig_df[!duplicated(sig_df$uid),]
rownames(col_ann) <- col_ann$uid
col_ann <- col_ann[colnames(hm_df),]
print(col_ann)

col_ann$tam_type <- as.factor(col_ann$tam_type)
tam_col <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(col_ann$tam_type)))
names(tam_col) <- levels(col_ann$tam_type)

col_ann$resource <- as.factor(col_ann$resource)
rs_col <- colorRampPalette(brewer.pal(12,"Set1"))(length(levels(col_ann$resource)))
names(rs_col) <- levels(col_ann$resource)

smp_hm_ann <- rowAnnotation(`TAM type` = col_ann$tam_type,
				Resource = col_ann$resource,
				col = list(`TAM type` = tam_col,
					   Resource = rs_col
					   ),
				annotation_legend_param = list(`TAM type` = list(nrow = 2, labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16)),
							       Resource = list(nrow = 2, labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 16))
				)
)

hm <- Heatmap(t(hm_df),
	      name = "Weight",
	      cluster_rows = TRUE,
	      cluster_columns = TRUE,
	      column_names_gp = gpar(fontsize = 7), 
	      row_names_gp = gpar(fontsize = 11),
	      show_column_names = TRUE,
	      row_split = col_ann$tam_type,
	      left_annotation = smp_hm_ann,
	      heatmap_legend_param = list(direction = "horizontal", 
					  labels_gp = gpar(fontsize = 12),
					  title_gp = gpar(fontsize = 16),
					  legend_width = unit(3.6, "cm")))
png(paste(input_dir, "all_signature_weight_clustered_heatmap.png", sep = "/"), res = 300, units = 'in', width = 36, height = 6)
draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
gar <- dev.off()

hm <- Heatmap(t(hm_df),
	      name = "Weight",
	      cluster_rows = FALSE,
	      cluster_columns = TRUE,
	      column_names_gp = gpar(fontsize = 7), 
	      row_names_gp = gpar(fontsize = 11),
	      show_column_names = TRUE,
	      row_order = order(col_ann$resource),
	      left_annotation = smp_hm_ann,
	      heatmap_legend_param = list(direction = "horizontal", 
					  labels_gp = gpar(fontsize = 12),
					  title_gp = gpar(fontsize = 16),
					  legend_width = unit(3.6, "cm")))
png(paste(input_dir, "all_signature_weight_ordered_heatmap.png", sep = "/"), res = 300, units = 'in', width = 36, height = 6)
draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
gar <- dev.off()
