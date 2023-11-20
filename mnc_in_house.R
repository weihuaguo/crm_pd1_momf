rm(list = ls())

subMonocle2 <- function(cts, metadata, plotPf, alignBatch, group, tifres = 300, rdsSave = FALSE, deNovoFlag = TRUE, pseudoTime = FALSE) {
	cat("Start to analyze with Monocle 2...\n")
	mncSt <- Sys.time()
	pd <- new("AnnotatedDataFrame", data = metadata)
	mncCDS <- newCellDataSet(cts, phenoData = pd, expressionFamily=negbinomial.size())
	mncCDS <- estimateSizeFactors(mncCDS)
	cat("\tFind the genes used for ordering cells...\n")
	tmpT <- Sys.time()
	seurat_gene <- differentialGeneTest(mncCDS, fullModelFormulaStr = "~PDL1") # NOTE: try use PDL1/seurat_clusters
	cat("\t")
	print(Sys.time()-tmpT)
	tmpT <- Sys.time()
	cat("\tStart to calculate pseudotime...\n")
	ordering_genes <- row.names(subset(seurat_gene, qval < 0.01))[order(seurat_gene$qval)][1:2000] # NOTE: choose different top n genes 
	# NOTE: top 1000 + seurat_clusters/PDL1 works
	mncCDS <- setOrderingFilter(mncCDS, ordering_genes)
	mncCDS <- reduceDimension(mncCDS, max_components = 2, method = "DDRTree")
	mncCDS <- orderCells(mncCDS)
	print(Sys.time()-tmpT)

	## TODO: refine the dot size and transparency for better visualization
	clusterGG <- plot_cell_trajectory(mncCDS, color_by = "seurat_clusters")
	ggsave(plot = clusterGG, filename = paste(plotPf, "traj_ddrtree_with_clusters.png", sep = ""), 
	       dpi = tifres, width = 9, height = 6)
	clusterGG <- plot_cell_trajectory(mncCDS, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 2)
	ggsave(plot = clusterGG, filename = paste(plotPf, "traj_ddrtree_with_clusters.png", sep = ""), 
	       dpi = tifres, width = 12, height = 6)

	groupGG <- plot_cell_trajectory(mncCDS, color_by = group)
	ggsave(plot = groupGG, filename = paste(plotPf, "traj_ddrtree_with_", group, ".png", sep = ""), 
	       dpi = tifres, width = 9, height = 6)

	groupGG <- plot_cell_trajectory(mncCDS, color_by = group) + facet_wrap(as.formula(paste("~", group)), nrow=1)
	ggsave(plot = groupGG, filename = paste(plotPf, "traj_ddrtree_with_facet_", group,".png", sep = ""), 
	       dpi = tifres, width = 9, height = 6)

	# TODO: add the heatmap of gene expression ordered by pseudotime
	pseudoGG <- plot_cell_trajectory(mncCDS, color_by = "Pseudotime")
	ggsave(plot = pseudoGG, filename = paste(plotPf, "traj_ddrtree_with_pseudotime.png", sep = ""),
	       dpi = tifres, width = 9, height = 6)
	print(head(seurat_gene))
	if (rdsSave) {saveRDS(mncCDS, paste(plotPf, "monocle2_object_with_trajectory.RDS", sep = ""))}
	cat("Time cost\t")
	print(Sys.time()-mncSt)
	return(mncCDS)
}

cat("Loading packages...\n")
#suppressMessages(library(Signac))
suppressMessages(library(Seurat))
#suppressMessages(library(SeuratWrappers))
suppressMessages(library(readxl))
suppressMessages(library(Matrix))
suppressMessages(library(ggplot2))
suppressMessages(library(viridis))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(patchwork))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(monocle))

cat("Setting up environment...\n")
workDir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub"
dataDir <- paste(workDir, "1_inhouse_scrna_seq", sep = "/") # Please change this to repeat the work
resDir <- paste(workDir, "2_inhouse_scrna_seq_result", sep = "/")

expID <- "inhouse_3prime_tumor_final"
expDir <- paste(resDir, expID, sep = "/")

m2Flag <- FALSE
m2BinFlag <- TRUE


use_res <- 0.5
subFolder <- "myeloid_sub_cluster"
subDir <- paste(expDir, subFolder, sep = "/")
resFolder <- paste("tam_subcluster_res", use_res, sep = "_")
tamPf <- paste(expDir, "/", subFolder, "/", resFolder, "/", resFolder, "_", sep = "")

if (m2Flag) {
	tamObj <- readRDS(paste(tamPf, "Seurat_Objects_Clustered.RDS", sep = ""))

	print(tamObj)
	cellAnnDf <- as.data.frame(read_excel(paste(tamPf, "ALL_POS_MARKERS.xlsx", sep = "")))
	rownames(cellAnnDf) <- cellAnnDf[,1]
	cellAnnDf[,1] <- NULL
	cellMarkerDf <- cellAnnDf
	cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$myeloid_type),]
	print(head(cellAnnDf))

	cols <- brewer.pal(length(unique(cellAnnDf$myeloid_type)), "Set2")
	names(cols) <- unique(cellAnnDf$myeloid_type)
	cols['PD-L1- TAM'] <- "dodgerblue"
	cols['PD-L1+ TAM'] <- "firebrick"

	pdl1_cols <- c("PD-L1- TAM" = "dodgerblue", "PD-L1+ TAM" = "firebrick")

	tamObj@meta.data$PDL1 <- "NA"
	for (ic in cellAnnDf$cluster) {
		tamObj@meta.data$PDL1[tamObj@meta.data$seurat_clusters == ic] <- cellAnnDf$myeloid_type[cellAnnDf$cluster == ic]
	}
	tamObj@meta.data$tam_yn <- ifelse(str_detect(tamObj@meta.data$PDL1, "TAM"), "TAM", "Non-TAM")
	proObj <- subset(tamObj, subset = tam_yn == "TAM")
	print(unique(proObj@meta.data$PDL1))

	tmpCtsMat <- as.matrix(proObj[["RNA"]]@counts)
	tmpMeta <- proObj@meta.data
	tmpMeta$seurat_clusters <- factor(tmpMeta$seurat_clusters)
	print(proObj)
	print(head(tmpMeta))

	tmpMNC2Obj <- subMonocle2(cts = tmpCtsMat, metadata = tmpMeta, plotPf = tamPf, alignBatch = "sampleID", group = 'PDL1',
				  rdsSave = TRUE, deNovoFlag = FALSE, pseudoTime = TRUE)
}

if (m2BinFlag) {
	cds <- readRDS(paste(tamPf, "monocle2_object_with_trajectory.RDS", sep = ""))
	print(cds)

	pdl1GG <- plot_cell_trajectory(cds, color_by = "PDL1") + 
		facet_wrap(~PDL1, nrow=1) + 
		scale_color_manual(values = c("dodgerblue", "firebrick"))
	ggsave(plot = pdl1GG, filename = paste(tamPf, "traj_ddrtree_with_facet_pdl1_refine.png", sep = ""), 
	       dpi = 600, width = 9, height = 5)


	mnc_pdf <- pData(cds)
	mnc_pdf$seurat_clusters <- factor(mnc_pdf$seurat_clusters)
	mnc_expr <- exprs(cds)
	print(mnc_expr[1:9,1:6])
	print(dim(mnc_expr))
	mnc_pdf$CD14 <- t(mnc_expr)[,"CD14"]
	print(head(mnc_pdf))
	sct_gg <- ggplot(mnc_pdf, aes(y = CD14, x = Pseudotime, color = PDL1)) +
		geom_point() + 
		scale_color_manual(values = c("dodgerblue", "firebrick")) +
		theme_classic()
	ggsave(paste(tamPf, "Response_scatter_cd14_pseudotime.png", sep = ""), sct_gg, dpi = 300, width = 6, height = 6)

	cat("Find the genes associated with pseudotime...\n")
	pt_de_res <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 8)
	write.csv(pt_de_res, paste(tamPf, "mnc2_pseudotime_differential_genes.csv", sep = ""))
	sig_gene_names <- row.names(subset(pt_de_res, qval < 0.1))
	top_pt_names <- row.names(subset(pt_de_res, qval < 0.1))[order(pt_de_res$qval)][1:20]
	print(top_pt_names)

	png(paste(tamPf, "mnc2_pseudotime_gene_heatmap.png", sep = ""), res = 300, width = 6, height = 9, units = 'in')
	print(plot_pseudotime_heatmap(cds[top_pt_names,], num_clusters = 3, cores = 8, show_rownames = T))
	gar <- dev.off()

	cat("Analyzing branches in single-cell trajectories...\n")
	BEAM_res <- BEAM(cds, branch_point = 1, cores = 16)
	BEAM_res <- BEAM_res[order(BEAM_res$qval),]
	write.csv(BEAM_res, paste(tamPf, "mnc2_branch_differential_genes.csv", sep = ""))
	top_beam_genes <- row.names(subset(BEAM_res, qval < 0.1))[order(BEAM_res$qval)][1:20]
	print(top_beam_genes)

	png(paste(tamPf, "mnc2_BEAM_gene_heatmap.png", sep = ""), res = 300, width = 6, height = 12, units = 'in')
	print(plot_genes_branched_heatmap(cds[top_beam_genes,], branch_point = 1, num_clusters = 4, 
					  cores = 8, use_gene_short_name = T, show_rownames = T))
	gar <- dev.off()

	top_beam_genes <- row.names(subset(BEAM_res, qval < 0.1))[order(BEAM_res$qval)][1:10]
	beam_gg <- plot_genes_branched_pseudotime(cds[top_beam_genes,], branch_point = 1, color_by = "PDL1", ncol = 2) +
		scale_color_manual(values = c("dodgerblue", "firebrick"))
	png(paste(tamPf, "beam_branched_pseudotime_gene_expression.png", sep = ""), res = 300, width = 6, height = 10, units = 'in')
	print(beam_gg)
	gar <- dev.off()



	mnc_pdf <- pData(cds)
	mnc_pdf$seurat_clusters <- factor(mnc_pdf$seurat_clusters)
	print(head(mnc_pdf))

	n <- 10
	## Manual histogram
	mnc_pdf$pt_bin <- 0
	step <- (max(mnc_pdf$Pseudotime) - min(mnc_pdf$Pseudotime))/(n-1)
	pbk <- min(mnc_pdf$Pseudotime)
	for (i in 1:n) {
		bk <- min(mnc_pdf$Pseudotime) + i*step
		tmp_mask <- mnc_pdf$Pseudotime < bk & mnc_pdf$Pseudotime >= pbk
		mnc_pdf$pt_bin[tmp_mask]<- i
		pbk <- bk
	}
	write.csv(mnc_pdf, paste(tamPf, "mnc2_pseudotime_binned_phenodata.csv", sep = ""))
	bin_cts_df <- as.data.frame(table(mnc_pdf$pt_bin, mnc_pdf$PDL1)) # NOTE: change the "PDL1"
	colnames(bin_cts_df) <- c("ptb", "Group", "cts")
	print(head(bin_cts_df))

	bin_cts_df <- bin_cts_df %>%
		group_by(Group, ptb) %>%
		mutate(Freq_per_group = cts/sum(cts))
	write.csv(bin_cts_df, paste(tamPf, "mnc2_pseudotime_binned_cts.csv", sep = ""))

	bar_gg <- ggplot(bin_cts_df, aes(x = ptb, y = cts, fill = Group, color = Group)) +
		geom_bar(stat = "identity", position = "fill", alpha = 0.8) +
		labs(x = "Pseudotime bins", y = "PD-L1+/-TAM in pseudotime (%)") +
		scale_color_manual(values = c("dodgerblue", "firebrick")) +
		scale_fill_manual(values = c("dodgerblue", "firebrick")) +
		theme_bw() +
		theme(axis.text.x = element_blank())
	ggsave(paste(tamPf, "mnc2_pt_bin_bar_group.png", sep = ""), dpi = 600, height = 4, width = 4.5)
}

