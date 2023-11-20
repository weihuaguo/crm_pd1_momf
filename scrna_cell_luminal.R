# Customized scripts for PD-L1 TAM
# Breast cancer scRNA-seq data
# Weihua Guo, Ph.D.
# 11/15/21

rm(list = ls())
options(future.globals.maxSize=171798691840)

srQCPrePro <- function(srsc, plotPf, res = 300) {
	cat("Start to QC and pre-processing Seurat object...\n")
	mito.genes <- grep(pattern = "^MT-", x = rownames(srsc@assays[["RNA"]]), value = TRUE)
	srsc[["percent.mt"]] <- PercentageFeatureSet(srsc, pattern = "^MT-")

	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6)

	qcsc1 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "percent.mt")
	qcsc2 <- FeatureScatter(srsc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	cqcgg <- CombinePlots(plots = list(qcsc1, qcsc2))
	ggsave(plot = cqcgg, filename = paste(plotPf, "QCPlot2.png", sep = ""), dpi = res, width = 9, height = 6)
	srsc <- subset(srsc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

	###################
	srsc <- NormalizeData(srsc, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

	################
	srsc <- FindVariableFeatures(srsc, selection.method = "vst", nfeatures = 2000)
	top10Features <- head(VariableFeatures(srsc), 18)
	varGG <- VariableFeaturePlot(srsc)
	labVargg <- LabelPoints(plot = varGG, points = top10Features, repel = TRUE)
	ggsave(plot = labVargg, filename = paste(plotPf, "variable_features.png", sep = ""), dpi = res, width = 9, height = 6)
	return(srsc)
}

inteSrPro <- function(srsc, plotPf, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE) {
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) 

	srsc <- ScaleData(srsc)

	#############
	srsc <- RunPCA(srsc)
	if (FALSE) { #FPCA
	pcaGG1 = VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
	ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), dpi = res, width = 9, height = 9)
	
	pcaGG2 = DimPlot(srsc, reduction = "pca")
	ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), dpi = res, width = 9, height = 6)

	png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
	DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
	gar = dev.off()
	} #FPCA

	if (jsFlag) {
		jst = Sys.time()
		cat("Running Jack Straw method to determine the PC dimension to use...\n")
		srsc <- JackStraw(srsc, num.replicate = 100)
		srsc <- ScoreJackStraw(srsc, dims = 1:20)
		jsGG <- JackStrawPlot(srsc, dims = 1:20)
		ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6)
		cat("Jack Straw costs")
		print(Sys.time()-jst)
	}

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)

	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- FindClusters(srsc, resolution = 0.5)

	srsc <- RunUMAP(srsc, dims = 1:ndim)
	
	umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
	ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)

	if (!is.null(plotFeat)) {
		cat("Start to plot feature plots and violin plots\n")
		vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "data", log = TRUE, ncol = 3, pt.size = 0.5)
		ggsave(plot = vlnGG, filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30) 

		featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
		ggsave(plot = featGG, filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
		       dpi = res, width = 16, height = 30)

		dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
		ggsave(plot = dotGG, filename = paste(plotPf, "DotPlot_for_cell_annotation.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}

	if (!is.null(batchName)) {
		batchCheckGG <- DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, filename = paste(plotPf, "batch_effect_check_dimPlot.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}
	
	srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

	write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
	write.csv(topMarkers, file = paste(plotPf, "top10_pos_markers.csv", sep = "")) 

	topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	hm <- DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend()
	ggsave(paste(plotPf, "top5_marker_heatmap.png", sep = ""), hm, dpi = res, width = 24, height = 18)

	topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	geneMarkers <- unique(as.character(topDotMarkers$gene))

	dotMarkerGG <- DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
	ggsave(plot = dotMarkerGG, filename = paste(plotPf, "Top5Marker_DotPlot.png", sep = ""), 
	       dpi = res, width = 9, height = 16)

	
	if (rdsSave) {
		saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	}
	print(srsc)
	return(srsc)
}

resScreen <- function(srsc, plotPf, plotDir, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, jsFlag = FALSE, resolutions) {
	qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = batchName)
	ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) 

	srsc <- ScaleData(srsc)

	#############
	srsc <- RunPCA(srsc)

	elGG <- ElbowPlot(srsc, ndims = 50)
	ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)
	srsc <- FindNeighbors(srsc, dims = 1:ndim)
	srsc <- RunUMAP(srsc, dims = 1:ndim)

	if (!is.null(batchName)) {
		batchCheckGG <- DimPlot(srsc, reduction = "umap", group.by = batchName)
		ggsave(plot = batchCheckGG, filename = paste(plotPf, "batch_effect_check_dimPlot.png", sep = ""),
		       dpi = res, width = 9, height = 6)
	}
	bldGG <- FeaturePlot(srsc, features = c("CD274", "SIGLEC15"), blend = TRUE)
	ggsave(paste(plotPf, "blend_pdl1_siglec15_expr_umap.png", sep = ""), bldGG, dpi = res, width = 32, height = 6)


	for (ires in resolutions) {
		resFolder <- paste("tam_subcluster_res", ires, sep = "_")
		dir.create(file.path(plotDir, resFolder), showWarnings = FALSE)
		resDir <- paste(plotDir, resFolder, sep = "/")
		rppf <- paste(resDir, "/", resFolder, "_", sep = "")

		srsc <- FindClusters(srsc, resolution = ires)

		srsc <- BuildClusterTree(srsc)
		png(paste(rppf, "cluster_tree.png", sep = ""), res = 300, width = 12, height = 18, unit = 'in')
		PlotClusterTree(object = srsc)
		gar <- dev.off()

		umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
		ggsave(plot = umapGG, filename = paste(rppf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)


		if (!is.null(plotFeat)) {
			cat("Start to plot feature plots and violin plots\n")
			vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "data", log = FALSE, ncol = 3, pt.size = 0.5)
			ggsave(plot = vlnGG, filename = paste(rppf, "violin_plot_for_cell_annotation.png", sep = ""), 
			       dpi = res, width = 16, height = 16) 

			featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
			ggsave(plot = featGG, filename = paste(rppf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
			       dpi = res, width = 16, height = 16)

			dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()
			ggsave(plot = dotGG, filename = paste(rppf, "DotPlot_for_cell_annotation.png", sep = ""),
			       dpi = res, width = 9, height = 6)

			rdgGG <- RidgePlot(srsc, features = plotFeat, ncol = 3)
			ggsave(paste(rppf, "ridge_plot_for_cell_annotation.png", sep = ""), rdgGG, dpi = res, width = 16, height = 18)

		}

		
		srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
		topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

		write.csv(srsc.markers, file = paste(rppf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
		write.csv(topMarkers, file = paste(rppf, "top10_pos_markers.csv", sep = "")) 

		topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		hm <- DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend()
		ggsave(paste(rppf, "top5_marker_heatmap.png", sep = ""), hm, dpi = res, width = 24, height = 18)

		topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		geneMarkers <- unique(as.character(topDotMarkers$gene))

		dotMarkerGG <- DotPlot(srsc, features = geneMarkers) + RotatedAxis() + coord_flip()
		ggsave(plot = dotMarkerGG, filename = paste(rppf, "Top5Marker_DotPlot.png", sep = ""), 
		       dpi = res, width = 9, height = 16)

		
		if (rdsSave) {
			saveRDS(srsc, file = paste(rppf, "Seurat_Objects_Clustered.RDS", sep = ""))
		}
		print(srsc)
	}
	return(srsc)
}

cat("Loading packages...\n")
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(future))
suppressMessages(library(stringr))
suppressMessages(library(Seurat))
suppressMessages(library(parallel))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MAST))
suppressMessages(library(clustree))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))


workDir <- "/home/weihua/mnts/smb_plee/Group/weihua/pdl1_data_hub"
dataDir <- paste(workDir, "4_cell_azizi_scrna_seq", sep = "/") # Please change this to repeat the work
resDir <- paste(workDir, "5_cell_azizi_scrna_seq_result", sep = "/")

expID <- "cell_azizi_luminal_tumor"
readFlag <- FALSE
inteFlag <- FALSE
inteFilter <- 100 # Minimum cell number for integration
clstFlag <- FALSE
majorVisFlag <- FALSE
tamClusterFlag <- FALSE
tamVisFlag <- FALSE
tamDEFlag <- TRUE

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, sep = "/")

ppf <- paste(resDir, "/", expID, "_", sep = "")
intePf <- paste(expDir, "/", expID, "_Integrated_Results_", sep = "")

if (readFlag) {
	cat("Reading the raw data...\n")
	raw_df <- read.csv(paste(dataDir, "/GSE114725_rna_raw.csv", sep = ""), check.names = F)
	print(raw_df[1:9,1:6])
	print(dim(raw_df))
	cell_id <- str_c(raw_df$patient, "_", raw_df$tissue, "_", raw_df$replicate, "_", raw_df$cellid)
	cts_df <- t(raw_df[,6:ncol(raw_df)])
	colnames(cts_df) <- cell_id
	print(cts_df[1:9,1:6])
	meta_df <- raw_df[,1:5]
	meta_df$cell_id <- str_c(meta_df$patient, "_", meta_df$tissue, "_", meta_df$replicate, "_", meta_df$cellid)
	print(head(meta_df))
	pat_meta <- read_excel(paste(dataDir, "/data_clinical_patient_GSE114727.xlsx", sep = ""))
	meta_df <- merge(meta_df, pat_meta, by = "patient", all.x = T)
	print(head(meta_df))
	rownames(meta_df) <- meta_df$cell_id
	meta_df$sampleID <- str_c(meta_df$patient,"_", meta_df$tissue)
	meta_df <- meta_df[meta_df$er_status == "Positive",]
	meta_df <- meta_df[meta_df$tissue == "TUMOR",]
	cts_df <- cts_df[,rownames(meta_df)]
	print(head(meta_df))
#	q(save = "no")
	srsc <- CreateSeuratObject(counts = cts_df, meta.data = meta_df)
	print(srsc)
	cell_cts <- srsc@meta.data %>% group_by(patient, tissue) %>% summarise(cell_num = n())
	write.csv(cell_cts, paste(ppf, "cell_counts.csv", sep = ""))
	srsc_list <- SplitObject(srsc, split.by = "sampleID")
	saveRDS(srsc_list, paste(ppf, "raw_split_list.rds", sep = ""))
	
	saName <- names(srsc_list)
	srCleanObjs <- vector(mode = "list", length = length(saName))
	names(srCleanObjs) <- saName


	for (isa in 1:length(srsc_list)) {
		st <- Sys.time()
		cat(saName[isa], "\n")
		tmpObj <- srsc_list[[saName[isa]]]
		print(tmpObj)
		print(head(tmpObj@meta.data))
		print(unique(tmpObj@meta.data$patient))
		plotPrefix <- paste(expDir, "/", expID, "_", saName[isa], "_", sep = "")
		tmpProObj <- srQCPrePro(tmpObj, plotPf = plotPrefix)
		srCleanObjs[saName[isa]] <- tmpProObj
	}
	rawObjRDS <- paste(expDir, "/", expID, "_all_raw_seurat_objects_list.RDS", sep = "")
	cleanObjRDS <- paste(expDir, "/", expID, "_all_clean_seurat_objects_list.RDS", sep = "")
	saveRDS(srsc_list, file = rawObjRDS)
	saveRDS(srCleanObjs, file = cleanObjRDS)
}

if (inteFlag) {
	srCleanObjs <- readRDS(paste(expDir, "/", expID, "_all_clean_seurat_objects_list.RDS", sep = ""))
	st <- Sys.time()
	cat("Start to merge the samples...\n")
	integrateList <- srCleanObjs
	for (inm in names(srCleanObjs)) {
		if(length(colnames(srCleanObjs[[inm]])) <= inteFilter) {integrateList[[inm]] <- NULL}
	}
	itg_ftrs <- SelectIntegrationFeatures(object.list = integrateList)
	integrateList <- lapply(X = integrateList, FUN = function(x) {
				 x <- ScaleData(x, features = itg_ftrs, verbose = TRUE)
				 x <- RunPCA(x, features = itg_ftrs, verbose = TRUE)
			     })
	saveRDS(integrateList, paste(ppf, "used_list_after_pca.RDS", sep = ""))

	cat("Perform integration... \n")
	itg_anchors <- FindIntegrationAnchors(object.list = integrateList, reference = 1,
					      anchor.features = itg_ftrs, reduction = "rpca")
	inteObjs <- IntegrateData(anchorset = itg_anchors)

#	integrateAnchors <- FindIntegrationAnchors(object.list = integrateList, dims = 1:50, k.filter = inteFilter)
#	inteObjs <- IntegrateData(anchorset = integrateAnchors, dims = 1:50)
	DefaultAssay(inteObjs) <- "integrated"
	print(inteObjs)
	print(Sys.time() - st)

	saveRDS(inteObjs, file = paste(intePf, "unprocessed_integrated_seurat_object.RDS", sep = ""))
}

if (clstFlag) {
	cat("Analyzing integrated Seurat object...\n")
	sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "NCAM1", "KIT", "FCER1A")
	pst <- Sys.time()
	inteObjs <- readRDS(paste(intePf, "unprocessed_integrated_seurat_object.RDS", sep = ""))

	inteProObj <- inteSrPro(inteObjs, batchName = c("patient"), plotPf = intePf, rdsSave = TRUE, ndim = 20, plotFeat = sigMarkers)
	cat("Total")
	print(Sys.time()-pst)
}


if (majorVisFlag) {
	cat("Merge and visualize all the cells\n")
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))
#	raw_cts <- t(as.matrix(proObj[['RNA']]@counts))
#	print(dim(raw_cts))
#	print(raw_cts[1:9,1:6])
#	write.csv(raw_cts, paste(intePf, "raw_counts_for_celltypist.csv", sep = ""))
#	q(save = "no")
	cellAnnDf <- as.data.frame(read_excel(paste(intePf, "ALL_POS_MARKERS.xlsx", sep = "")))
	rownames(cellAnnDf) <- cellAnnDf[,1]
	cellAnnDf[,1] <- NULL
	cellMarkerDf <- cellAnnDf
	cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$cell_type),]
	print(head(cellAnnDf))
	proObj@meta.data$major_cell_type <- "NA"
	print(head(proObj@meta.data))
	for (ic in cellAnnDf$cluster) {
		proObj@meta.data$major_cell_type[proObj@meta.data$seurat_clusters == ic] <- cellAnnDf$cell_type[cellAnnDf$cluster == ic]
	}
	print(head(proObj@meta.data))

	sigMarkers <- c("CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA.DRA", "CD1C", "TPSAB1", "LILRA4")
	inteGenes <- rownames(proObj[['integrated']]@data)
	rnaGenes <- rownames(proObj[['RNA']]@data)
	olGenes <- intersect(sigMarkers, inteGenes)
	diffGenes <- setdiff(sigMarkers, inteGenes)
	diffRNAGenes <- setdiff(diffGenes, rnaGenes)
	rnaOlGenes <- intersect(diffGenes, rnaGenes)

#	print(olGenes)
#	print(diffGenes)
#	print(rnaOlGenes)
	print(diffRNAGenes)
#	q(save = "no")
	inteExpr <- proObj[['integrated']]@data[olGenes,]
	rnaExpr <- proObj[['RNA']]@data[rnaOlGenes,]

	inteExpr <- as.data.frame(inteExpr)
	rnaExpr <- as.data.frame(rnaExpr)
	sigExpr <- rbind(inteExpr, rnaExpr)
	sigExpr$gene <- rownames(sigExpr)
	print(dim(sigExpr))
	print(sigExpr[1:9,1:6])

	gathExpr <- gather(sigExpr, "cell_id", "expr", rownames(proObj@meta.data))
	gathExpr <- merge(gathExpr, proObj@meta.data, by.x = "cell_id", by.y = "row.names", all.x = TRUE)
	gathExpr$gene <- factor(gathExpr$gene, levels = sigMarkers)
#	print(head(gathExpr))

	figs1e <- ggplot(gathExpr, aes(x = seurat_clusters, y = expr, fill = seurat_clusters)) +
		geom_violin(scale = "width", trim = TRUE) +
		facet_wrap(.~gene, nrow = 1, scales = "free_x") +
		labs(x = "Cluster (r = 0.5)", y = "Normalized expression", fill = "Cluster (r = 0.5)") +
		coord_flip() +
		theme_classic()
	ggsave(plot = figs1e, filename = paste(intePf, "FigSCellE_sc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 4.5) # FigSCellx
	saveRDS(figs1e, paste(intePf, "FigSCellE_sc_violin_plot_for_cell_annotation.RDS", sep = ""))


	figs1j <- ggplot(gathExpr, aes(x = major_cell_type, y = expr, fill = major_cell_type)) +
		geom_violin(scale = "width", trim = TRUE) +
		facet_wrap(.~gene, nrow = 1, scales = "free_x") +
		scale_fill_brewer(palette = "Spectral") +
		labs(x = "General cell type in TME", y = "Normalized expression", fill = "General cell type\nin TME") +
		coord_flip() +
		theme_classic()
	ggsave(plot = figs1j, filename = paste(intePf, "FigSCellI_mc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 2.5) # FigSCellx
	saveRDS(figs1j, paste(intePf, "FigSCellI_mc_violin_plot_for_cell_annotation.RDS", sep = ""))

	cat("Heatmaps\n")
	topMarkers <- cellMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	figs1d <- DoHeatmap(object = proObj, features = topMarkers$gene) + NoLegend()
	ggsave(paste(intePf, "FigSCellD_top5_marker_heatmap.png", sep = ""), figs1d, dpi = 300, width = 24, height = 18)
	saveRDS(figs1d, paste(intePf, "FigSCellD_sc_top5_marker_heatmap.RDS", sep = ""))

	if (FALSE) { # NOTE: To save some time
	c <- 0
	for (imct in unique(proObj@meta.data$major_cell_type)) {
		cat("Caculating", imct, "\n")
		tmpMarkers <- FindMarkers(proObj, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25, group.by = "major_cell_type", ident.1 = imct)
		tmpMarkers$gene <- rownames(tmpMarkers)
		tmpMarkers$cluster <- imct
		if (c == 0) {
			gnrlMarkerDf <- tmpMarkers
		} else {
			gnrlMarkerDf <- rbind(gnrlMarkerDf, tmpMarkers)
		}
		c <- c+1
	}
	print(head(gnrlMarkerDf))
	write.csv(gnrlMarkerDf, paste(intePf, "ALL_POS_MARKERS_MAJOR_CELL_TYPE.csv", sep = ""))
	} else {
		gnrlMarkerDf <- read.csv(paste(intePf, "ALL_POS_MARKERS_MAJOR_CELL_TYPE.csv", sep = ""), row.names = 1, header = T, check.names = F)
	}

	cols <- brewer.pal(length(unique(gnrlMarkerDf$cluster)), "Spectral")

	topMarkers <- gnrlMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	topOrders <- topMarkers$gene[order(topMarkers$cluster, topMarkers$avg_log2FC, decreasing = T)]
	topOrders <- topOrders[!duplicated(topOrders)]
	figs1i <- DoHeatmap(object = proObj, features = topMarkers$gene, group.colors = cols, group.by = "major_cell_type") + NoLegend()
#	print(figs1i$data)
	figs1i$data$Feature <- factor(figs1i$data$Feature, levels = topOrders)
	ggsave(paste(intePf, "FigSCellH_top5_marker_mct_heatmap.png", sep = ""), figs1i, dpi = 300, width = 24, height = 18)
	saveRDS(figs1i, paste(intePf, "FigSCellH_mc_top10_marker_heatmap.RDS", sep = ""))

	geneExpr <- AverageExpression(proObj, features = topMarkers$gene, group.by = "major_cell_type") 
	geneExpr <- as.data.frame(geneExpr$RNA)
	cellType <- data.frame(colnames(geneExpr))
	colnames(cellType) <- 'General cell type in TME'
	ann_cols <- brewer.pal(length(unique(topMarkers$cluster)), "Spectral")
	names(ann_cols) <- unique(topMarkers$cluster)[order(unique(topMarkers$cluster))]
	print(ann_cols)
	top_anno = HeatmapAnnotation(df = cellType,
				     border = T,
				     show_annotation_name = F,
				     gp = gpar(col = 'black'),
				     col = list(`General cell type in TME` = ann_cols))
	marker_exp <- t(scale(t(geneExpr),scale = T,center = T))
	print(topMarkers$gene)
	print(rownames(marker_exp))
	print(dim(topMarkers))
	print(dim(marker_exp))
	avgHm <- Heatmap(marker_exp,
			 name = "Z-score\n(Average expression)",
		cluster_rows = F,
		cluster_columns = F,
		show_column_names = F,
		show_row_names = T,
		row_order = order(topMarkers$cluster),
		heatmap_legend_param = list(title = ""),
		col = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
		border = 'black',
		rect_gp = gpar(col = "black", lwd = 1),
		row_names_gp = gpar(fontsize = 10),
		column_names_gp = gpar(fontsize = 10),
		top_annotation = top_anno)
	png(filename = paste(intePf, "FigSCellH_mc_top5_marker_avg_heatmap.png", sep = ""),res = 600, width = 6, height = 8,units = "in")
	print(draw(avgHm, merge_legend = T))
	dev.off()
	q(save = "no")


	cat("Cell proportions\n")
	scCts <- proObj@meta.data %>%
		group_by(patient, seurat_clusters) %>%
		summarize(n = n())

	figs1c <- ggplot(scCts, aes(x = patient, y = n, fill = seurat_clusters)) +
		geom_bar(stat = "identity", position = "fill", color = "gray") +
		coord_flip() +
		labs(x = "Sample ID", y = "Cell proportion", fill = "Cluster (r = 0.5)") +
		theme_classic()
	ggsave(plot = figs1c, filename = paste(intePf, "FigSCellC_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1c, paste(intePf, "FigSCellC.RDS", sep = ""))

	mcCts <- proObj@meta.data %>%
		group_by(patient, major_cell_type) %>%
		summarize(n = n())
	print(head(mcCts))

	figs1h <- ggplot(mcCts, aes(x = patient, y = n, fill = major_cell_type)) +
		geom_bar(stat = "identity", position = "fill", color = "gray") +
		scale_fill_brewer(palette = "Spectral") +
		coord_flip() +
		labs(x = "Sample ID", y = "Cell proportion", fill = "General cell type\nin TME") +
		theme_classic()
	ggsave(plot = figs1h, filename = paste(intePf, "FigSCellG_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1h, paste(intePf, "FigSCellG.RDS", sep = ""))
	
	if (FALSE) {# FPDL1 counts here are not neccessary.
	proObj@meta.data$PDL1_expr <- proObj[["RNA"]]@data["CD274",]
	proObj@meta.data$PDL1_status <- ifelse(proObj@meta.data$PDL1_expr > 0, "PDL1+", "PDL1-")

	pdl1Cts <- proObj@meta.data %>%
		filter(major_cell_type == "Cancer cell" | major_cell_type == "TAM") %>%
		group_by(patient, major_cell_type, PDL1_status) %>%
		summarize(n = n())

	pdl1Spr <- as.data.frame(spread(pdl1Cts, "major_cell_type", "n"))
	pdl1Spr$CM_ratio <- pdl1Spr$`Cancer cell`/pdl1Spr$TAM

	cmSpr <- mcCts %>%
		filter(major_cell_type == "Cancer cell"|major_cell_type == "TAM") %>%
		spread(major_cell_type, n)
	cmSpr$CM_ratio <- cmSpr$`Cancer cell`/cmSpr$TAM
	cmSpr$PDL1_status <- "PDL1+ & PDL1-"

	cmRatio <- rbind(cmSpr, pdl1Spr)

	figs1g <- ggplot(cmRatio, aes(x = PDL1_status, y = CM_ratio, fill = patient)) +
		geom_bar(stat = "identity", color = "gray", position = "dodge") +
		scale_fill_brewer(palette = "Dark2") +
		coord_flip() +
		labs(x = "PD-L1 status", y = "Cancer cell/TAM ratio", fill = "Sample ID") +
		theme_classic()
	ggsave(plot = figs1g, filename = paste(intePf, "FigSCellG_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1g, paste(intePf, "FigSCellG.RDS", sep = ""))
	} # FPDL1

	cat("UMAPs\n")
	figs1a <- DimPlot(proObj, reduction = "umap", label = TRUE) + labs(color = "Clusters (r = 0.5)", title = "")
	ggsave(plot = figs1a, filename = paste(intePf, "FigSCellA_UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1a, paste(intePf, "FigSCellA.RDS", sep = ""))

	figs1b <- DimPlot(proObj, reduction = "umap", label = FALSE, group.by = "patient") + 
		labs(color = "Sample ID", title = "Batch effect crossing samples") +
		scale_color_brewer(palette = "Dark2")
	ggsave(plot = figs1b, filename = paste(intePf, "FigSCellB_UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs1b, paste(intePf, "FigSCellB.RDS", sep = ""))

	figs1f <- DimPlot(proObj, reduction = "umap", label = TRUE, group.by = "major_cell_type") + 
		labs(color = "General cell type\nin TME", title = "Annotated and merged cell types") +
		scale_color_brewer(palette = "Spectral")

	ggsave(plot = figs1f, filename = paste(intePf, "FigSCellF_UMAP.png", sep = ""), dpi = 600, width = 9, heigh = 6)
	saveRDS(figs1f, paste(intePf, "FigSCellF.RDS", sep = ""))

}

if (tamClusterFlag) {
	cat("Cluster TAMs...\n")
	proObj <- readRDS(paste(intePf, "Seurat_Objects_Clustered.RDS", sep = ""))

	ist <- Sys.time()
	subSrsc <- subset(proObj, idents = c(4,7,8,10,12,15,17)) # NOTE: ALL_POS_MARKER.xlsx
	colnames(subSrsc@meta.data)[colnames(subSrsc@meta.data) == "integrated_snn_res.0.5"] <- "integrated_snn_res.0.0"
#	print(head(subSrsc@meta.data))

	newSrsc <- CreateSeuratObject(counts = subSrsc[["RNA"]]@counts, meta.data = subSrsc@meta.data)
	seuratList <- SplitObject(newSrsc, split.by = "patient")
	
	subFolder <- "myeloid_sub_cluster"
	dir.create(file.path(expDir, subFolder), showWarnings = FALSE)
	subDir <- paste(expDir, subFolder, sep = "/")
	subPrefix <- paste(subDir, "/", expID, "_", subFolder, "_", sep = "")

	for (i in names(seuratList)) {
		cat("Patient", i, ":", length(colnames(seuratList[[i]])), "\n")
		seuratList[[i]] <- NormalizeData(seuratList[[i]], verbose = FALSE)
		seuratList[[i]] <- FindVariableFeatures(seuratList[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
		if (length(colnames(seuratList[[i]])) <= 30) {
			## NOTE: remove the patients with few cells (<100 cells) for integration!
			seuratList[[i]] <- NULL
			cat("\tPatient", i, "is removed due to small cell number!!!\n")
		}
	}
	integrateAnchors <- FindIntegrationAnchors(object.list = seuratList, k.filter = 30)
	inteObjs <- IntegrateData(anchorset = integrateAnchors, k.weight = 30)
	DefaultAssay(inteObjs) <- "integrated"
	saveRDS(inteObjs, paste(subPrefix, "_Integrated_results.RDS", sep = ""))
	cat("Integration time cost:")
	print(Sys.time()-ist)

	pst <- Sys.time()
	ress <- c(seq(0.1,1.2, 0.1))#, seq(1.0,2.0,0.2))
	geneOi <- c("CD14", "CD68", "HLA-DRA", "CD274", "SIGLEC15", "PDCD1LG2", "CD1C", "ITGAX", "FCER1A", "TPSAB1", "KIT", "FUT4", "LILRA4", "CLEC4C", "CD1A")
	inteProObj <- resScreen(inteObjs, plotPf = subPrefix, plotDir = subDir, rdsSave = TRUE, jsFlag = FALSE, ndim = 20,
				batchName = c("patient"), plotFeat = geneOi, resolutions = ress)
	cat("Analysis time cost:")
	print(Sys.time()-pst)

	all_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(ress)+1)
	cltr <- clustree(inteProObj, prefix = "integrated_snn_res.") + scale_color_manual(values = all_cols)
	ggsave(paste(subPrefix, "clustree.png", sep = ""), cltr, dpi = 300, width = 12, height = 12) # FigS2
	cat("++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
}

if (tamVisFlag) {
	use_res <- 0.6
	subFolder <- "myeloid_sub_cluster"
	subDir <- paste(expDir, subFolder, sep = "/")
	resFolder <- paste("tam_subcluster_res", use_res, sep = "_")
	tamPf <- paste(expDir, "/", subFolder, "/", resFolder, "/", resFolder, "_", sep = "")
	tamObj <- readRDS(paste(tamPf, "Seurat_Objects_Clustered.RDS", sep = ""))
#	print(tamObj)
	print(head(tamObj@meta.data))
	tamObj@meta.data$patient <- tamObj@meta.data$G_patient
	cellAnnDf <- as.data.frame(read_excel(paste(tamPf, "ALL_POS_MARKERS.xlsx", sep = "")))
	rownames(cellAnnDf) <- cellAnnDf[,1]
	cellAnnDf[,1] <- NULL
	cellMarkerDf <- cellAnnDf
	cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$myeloid_type),]
	print(head(cellAnnDf))

	cols <- brewer.pal(length(unique(cellAnnDf$myeloid_type)), "Set2")
	names(cols) <- unique(cellAnnDf$myeloid_type)
	cols['PD-L1- mo/macrophage'] <- "dodgerblue"
	cols['PD-L1+ mo/macrophage'] <- "firebrick"

	pdl1_cols <- c("PD-L1- mo/macrophage" = "dodgerblue", "PD-L1+ mo/macrophage" = "firebrick")

	tamObj@meta.data$myeloid_type <- "NA"
	for (ic in cellAnnDf$cluster) {
		tamObj@meta.data$myeloid_type[tamObj@meta.data$seurat_clusters == ic] <- cellAnnDf$myeloid_type[cellAnnDf$cluster == ic]
	}
#	print(head(tamObj@meta.data))
	tamObj@meta.data$PDL1_expr <- tamObj[['RNA']]@data["CD274",]
	tamObj@meta.data$SIGLEC15_expr <- tamObj[['RNA']]@data["SIGLEC15",]
	tamObj@meta.data$PDL1_SIGLEC15_PN <- ifelse(tamObj@meta.data$PDL1_expr == 0, 
						    ifelse(tamObj@meta.data$SIGLEC15_expr == 0, "PD-L1- SIGLEC15-", "SIGLEC15+"), 
						    ifelse(tamObj@meta.data$SIGLEC15_expr == 0, "PD-L1+", "PD-L1+ SIGLEC15+")
	)

	cat("Heatmaps\n")
	topMarkers <- cellMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
	figs3c <- DoHeatmap(object = tamObj, features = topMarkers$gene) + NoLegend()
	ggsave(paste(tamPf, "FigSCellTAMC_top5_marker_heatmap.png", sep = ""), figs3c, dpi = 300, width = 18, height = 12)
	saveRDS(figs3c, paste(tamPf, "FigSCellTAMC_sc_top5_marker_heatmap.RDS", sep = ""))

	true_pn <- tamObj@meta.data %>%
		group_by(seurat_clusters, PDL1_SIGLEC15_PN) %>%
		summarize(nCell = n())
	spr_pn <- spread(true_pn[,c("seurat_clusters", "PDL1_SIGLEC15_PN", "nCell")], "PDL1_SIGLEC15_PN", "nCell")
	spr_pn[is.na(spr_pn)] <- 0
	gath_pn <- gather(spr_pn, "PDL1_SIGLEC15_PN", "nCell", unique(true_pn$PDL1_SIGLEC15_PN))
	gath_pn <- merge(gath_pn, cellAnnDf[,c("cluster", "myeloid_type")], by.x = "seurat_clusters", by.y = "cluster", all.x = T)
	write.csv(gath_pn, paste(tamPf, "true_pn_counts.csv", sep = "")) #ST1

	plot_pn <- gath_pn[gath_pn$PDL1_SIGLEC15_PN %in% c("PD-L1+", "SIGLEC15+"),]
	plot_pn <- plot_pn %>% group_by("PDL1_SIGLEC15_PN") %>% mutate(pos_total_num = sum(nCell))
	plot_pn <- plot_pn[str_detect(plot_pn$myeloid_type, "macrophage"),]
	plot_pn$pos_per <- plot_pn$nCell/plot_pn$pos_total_num * 100
	pn_gg <- ggplot(plot_pn, aes(x = myeloid_type, y = pos_per)) +
		geom_boxplot(aes(fill = myeloid_type), alpha = 0.2) +
		geom_point(aes(color = seurat_clusters)) +
		stat_compare_means(aes(group = myeloid_type), label = "p.signif") +
		facet_wrap(~PDL1_SIGLEC15_PN, nrow = 1, scales = "free") +
		scale_fill_manual(values = pdl1_cols) +
		labs(x = "PD-L1 status", y = "Percentage of TP cells within each cluster to total TP cell number", color = "Cluster (r = 0.6)", fill = "PD-L1 status") +
		theme_classic() +
		theme(axis.text.x = element_blank())
	ggsave(paste(tamPf, "FigSCellTAMX_true_positive_relative_pct_boxplot.png", sep = ""), pn_gg, dpi = 300, width = 7.5, height = 6)

	if (FALSE) {#FPN
	gath_pn <- gath_pn %>%
		group_by(seurat_clusters) %>%
		mutate(nCluster = sum(nCell))
	gath_pn$PDL1_SIGLEC15_PN <- factor(gath_pn$PDL1_SIGLEC15_PN, levels = c("PD-L1+", "SIGLEC15+", "PD-L1- SIGLEC15-", "PD-L1+ SIGLEC15+"))
	figs3h <- ggplot(gath_pn, aes(x = PDL1_SIGLEC15_PN, y = nCell, fill = seurat_clusters)) +
		geom_bar(stat = "identity", position = "fill", color = "black") +
		labs(y = "Relative cell number", x = "Monocyte/macrophage type based on \nnon-zero normalized expression", fill = "Cluster (r = 0.6)") +
		coord_flip() +
		guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
		theme_classic() +
		theme(legend.position = "top")
	ggsave(plot = figs3h, filename = paste(tamPf, "FigSCellTAMH_true_pn_bar.png", sep = ""), dpi = 300, width = 9, heigh = 3)
	saveRDS(figs3h, paste(tamPf, "FigSCellTAMH_true_pn_bar.RDS", sep = ""))
	} #FPN

	figs3f <- FeaturePlot(tamObj, features = c("CD274", "SIGLEC15"), blend = T)
	ggsave(plot = figs3f, filename = paste(tamPf, "FigSCellTAMF_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 18, heigh = 4)
	saveRDS(figs3f, paste(tamPf, "FigSCellTAMF.RDS", sep = ""))

	cat("UMAPs\n")
	figs3a <- DimPlot(tamObj, reduction = "umap", label = TRUE) + labs(color = "Clusters (r = 0.6)", title = "")
	ggsave(plot = figs3a, filename = paste(tamPf, "FigSCellTAMA_UMAP_all_tam_cluster.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs3a, paste(tamPf, "FigSCellTAMA.RDS", sep = ""))

	figs3b <- DimPlot(tamObj, reduction = "umap", label = FALSE, group.by = "patient") + 
		labs(color = "Sample ID", title = "Batch effect crossing samples") +
		scale_color_brewer(palette = "Dark2")
	ggsave(plot = figs3b, filename = paste(tamPf, "FigSCellTAMB_UMAP_all_tam_batch.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(figs3b, paste(tamPf, "FigSCellTAMB.RDS", sep = ""))

	fig1a <- DimPlot(tamObj, reduction = "umap", label = TRUE, group.by = "myeloid_type") + 
		labs(color = "Myeloid cell types", title = "Annotated and merged cell types") +
#		scale_color_brewer(palette = "Set3") +
		scale_color_manual(values = cols)
	ggsave(plot = fig1a, filename = paste(tamPf, "Fig1A_UMAP_all_myeloid_pdl1.png", sep = ""), dpi = 600, width = 9, heigh = 6)
	saveRDS(fig1a, paste(tamPf, "Fig1A.RDS", sep = ""))

	sigMarkers <- c("CD14", "CD68", "HLA.DRA", "CD1C", "CST3", "GPX1", "KIT", "TPSAB1", "CLEC4C", "CEACAM8", "CTAG2", "S100A9", "MARC1")#, "CD274", "SIGLEC15", "CD3D")
	inteGenes <- rownames(tamObj[['integrated']]@data)
	rnaGenes <- rownames(tamObj[['RNA']]@data)
	olGenes <- intersect(sigMarkers, inteGenes)
	diffGenes <- setdiff(sigMarkers, inteGenes)
	diffRNAGenes <- setdiff(diffGenes, rnaGenes)
	rnaOlGenes <- intersect(diffGenes, rnaGenes)

#	print(olGenes)
#	print(diffGenes)
#	print(rnaOlGenes)
	print(diffRNAGenes)
#	q(save = "no")
	inteExpr <- tamObj[['integrated']]@data[olGenes,]
	rnaExpr <- tamObj[['RNA']]@data[rnaOlGenes,]

	inteExpr <- as.data.frame(inteExpr)
	rnaExpr <- as.data.frame(rnaExpr)
	sigExpr <- rbind(inteExpr, rnaExpr)
	sigExpr$gene <- rownames(sigExpr)
	print(dim(sigExpr))
	print(sigExpr[1:9,1:6])
#	q(save = "no")

	gathExpr <- gather(sigExpr, "cell_id", "expr", rownames(tamObj@meta.data))
	gathExpr <- merge(gathExpr, tamObj@meta.data, by.x = "cell_id", by.y = "row.names", all.x = TRUE)
	gathExpr$gene <- factor(gathExpr$gene, levels = sigMarkers)
	print(head(gathExpr))

	figs3d <- ggplot(gathExpr, aes(x = seurat_clusters, y = expr, fill = seurat_clusters)) +
		geom_violin(scale = "width", trim = TRUE) +
		coord_flip() +
		facet_grid(myeloid_type~gene, space = "free_y", scale = "free_y", switch = NULL) +
		labs(x = "Cluster (r = 0.6)", y = "Normalized expression", fill = "Cluster (r = 0.6)") +
		theme_classic() +
		theme(strip.text.y = element_text(angle = 0), axis.text.x = element_blank())
	ggsave(plot = figs3d, filename = paste(tamPf, "FigSCellTAMD_sc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 4.5)
	saveRDS(figs3d, paste(tamPf, "FigSCellTAMD_sc_violin_plot_for_cell_annotation.RDS", sep = ""))


	figs3e <- ggplot(gathExpr, aes(x = myeloid_type, y = expr, fill = myeloid_type)) +
		geom_violin(scale = "width", trim = TRUE) +
		facet_wrap(.~gene, nrow = 1, scales = "free_x") +
#		scale_fill_brewer(palette = "Spectral") +
		scale_fill_manual(values = cols) +
		labs(x = "General cell type", y = "Normalized expression", fill = "General cell type") +
		coord_flip() +
		theme_classic() +
		theme(axis.text.x = element_blank())
	ggsave(plot = figs3e, filename = paste(tamPf, "FigSCellTAME_mc_violin_plot_for_cell_annotation.png", sep = ""), 
	       dpi = 300, width = 16, height = 2.5) # FigS1x
	saveRDS(figs3e, paste(tamPf, "FigSCellTAME_mc_violin_plot_for_cell_annotation.RDS", sep = ""))
}

if (tamDEFlag) {
	use_res <- 0.6
	subFolder <- "myeloid_sub_cluster"
	subDir <- paste(expDir, subFolder, sep = "/")
	resFolder <- paste("tam_subcluster_res", use_res, sep = "_")
	tamPf <- paste(expDir, "/", subFolder, "/", resFolder, "/", resFolder, "_", sep = "")
	tamObj <- readRDS(paste(tamPf, "Seurat_Objects_Clustered.RDS", sep = ""))
	print(tamObj)
	cellAnnDf <- as.data.frame(read_excel(paste(tamPf, "ALL_POS_MARKERS.xlsx", sep = "")))
	rownames(cellAnnDf) <- cellAnnDf[,1]
	cellAnnDf[,1] <- NULL
	cellMarkerDf <- cellAnnDf
	cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$myeloid_type),]
	cellAnnCleanDf <- cellAnnDf[str_detect(cellAnnDf$myeloid_type, "macrophage"),]
	print(head(cellAnnDf))

	cols <- brewer.pal(length(unique(cellAnnDf$myeloid_type)), "Set2")
	names(cols) <- unique(cellAnnDf$myeloid_type)
	cols['PD-L1- mo/macrophage'] <- "dodgerblue"
	cols['PD-L1+ mo/macrophage'] <- "firebrick"

	pdl1_cols <- c("PD-L1- mo/macrophage" = "dodgerblue", "PD-L1+ mo/macrophage" = "firebrick")

	tamObj@meta.data$myeloid_type <- "NA"
	for (ic in cellAnnDf$cluster) {
		tamObj@meta.data$myeloid_type[tamObj@meta.data$seurat_clusters == ic] <- cellAnnDf$myeloid_type[cellAnnDf$cluster == ic]
	}
	tamObj@meta.data$tam_yn <- ifelse(str_detect(tamObj@meta.data$myeloid_type, "mo/macrophage"), "TAM", "Non-TAM")
	proObj <- subset(tamObj, subset = tam_yn == "TAM")
	print(unique(proObj@meta.data$myeloid_type))

	proObj@meta.data$seurat_clusters <- factor(proObj@meta.data$seurat_clusters, levels = cellAnnCleanDf$cluster[order(cellAnnCleanDf$myeloid_type, cellAnnCleanDf$cluster)])
	pdl1VlnGG <- VlnPlot(proObj, features = "CD274", group.by = "seurat_clusters", log = FALSE, pt.size = 0.1) + scale_y_continuous(limits = c(0.0, NA))
	ggsave(plot = pdl1VlnGG, filename = paste(tamPf, "PDL1_FigSCellTAMY_manual_selected_VlnPlot.png", sep = ""),
	       dpi = 300, width = 5, height = 5)

	cat("UMAPs\n")
	if (FALSE) { #F1a
	fig1a <- DimPlot(proObj, reduction = "umap", label = TRUE) + labs(color = "Clusters (r = 0.5)", title = "")
	ggsave(plot = fig1a, filename = paste(tamPf, "Fig1A_UMAP_clean_tam_cluster.png", sep = ""), dpi = 300, width = 9, heigh = 6)
	saveRDS(fig1a, paste(tamPf, "Fig1A.RDS", sep = ""))
	} #F1a

	fig1b <- FeaturePlot(proObj, features = c("CD274", "SIGLEC15"), blend = T)
	ggsave(plot = fig1b, filename = paste(tamPf, "Fig1B_UMAP_clean_blend_pdl1.png", sep = ""), dpi = 300, width = 18, heigh = 4)
	saveRDS(fig1b, paste(tamPf, "FigCellB.RDS", sep = ""))

	fig1d <- DimPlot(proObj, reduction = "umap", label = FALSE, group.by = "myeloid_type") + 
		labs(color = "PD-L1 status") +
		scale_color_manual(values = pdl1_cols)
	ggsave(plot = fig1d, filename = paste(tamPf, "FigCellD_UMAP_clean_tam_pdl1.png", sep = ""), dpi = 600, width = 9, heigh = 6)
	saveRDS(fig1d, paste(tamPf, "FigCellD.RDS", sep = ""))

	vlnGenes <- c("HLA-DRA", "C1QA", "IL1B", "CD74", "CD83", "FOS", "JUNB", "CEBPD", "SPP1", "FABP5", "CSTB", "IL1RN", "CD9", "CD52", "LPL", 
		      "MMP9", "COL1A1", "FN1", "CD14", "CD274", "SIGLEC15")    
	vlnGenes <- c("SPP1", "FN1", "TREM2", "CD9", "FABP5", "FABP4", "LGALS3", "IL1RN", "CSTB", "LPL", 
		     "C1QA", "CEBPD", "HLA.DQA1", "HLA.DQB1", "HLA.DQA2", "FOSB", "IL1B", "HLA.DPB1", "CD14", "CD274", "SIGLEC15")

	for (ivg in vlnGenes) {
		cat(ivg, "\n")
		pdl1VlnGG <- VlnPlot(proObj, features = ivg, group.by = "myeloid_type", log = FALSE, pt.size = 0,
				     cols = c('dodgerblue', 'firebrick')) + scale_y_continuous(limits = c(0.0, NA))
		ggsave(plot = pdl1VlnGG, filename = paste(tamPf, ivg, "_Fig1HI_manual_selected_VlnPlot.png", sep = ""),
		       dpi = 300, width = 5, height = 5)

		pdl1VlnGG <- VlnPlot(proObj, features = ivg, group.by = "myeloid_type", log = FALSE, pt.size = 0,
				     cols = c('dodgerblue', 'firebrick')) + scale_y_continuous(limits = c(0.0, NA)) + stat_compare_means(label = "p.signif")
		ggsave(plot = pdl1VlnGG, filename = paste(tamPf, ivg, "_Fig1HI_manual_selected_VlnPlot_stats.png", sep = ""),
		       dpi = 300, width = 5, height = 6)

	}

	de_res <- FindMarkers(proObj, ident.1 = "PD-L1+ mo/macrophage", ident.2 = "PD-L1- mo/macrophage", group.by = "myeloid_type", logfc.threshold = 0.5, test.use = "MAST")
	print(head(de_res))
	write.csv(de_res, paste(tamPf, "pdl1_tam_de.csv", sep = "")) # Supplementary table XX
	plot_de <- de_res[de_res$p_val_adj < 0.05 & abs(de_res$avg_log2FC)>= 1,]
	print(plot_de)
	figs4a <- DotPlot(proObj, features = rownames(plot_de)[order(plot_de$avg_log2FC)], group.by = "myeloid_type") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(paste(tamPf, "FigS4A_dot_plot_pdl1_de.png", sep = ""), figs4a, dpi = 300, width = 24, height = 4)

	labFeatures <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", 
			"C1QA", "C1QB", "C1QC", "IL1B", "CXCL2", "CXCL3", "CXCL8", "CCL3", "CCL4", "CCL18", 
			"CD74", "CD83", "FOS", "FOSB", "JUN", "JUNB", "CEBPD", "FOLR2", "MS4A6A", 
			"SPP1", "FABP4", "FABP5", "CSTB", "SPARC", "IL1RN", "CD9", "CD52", "CHI3L1", "LGALS3", 
			"PKM", "LPL", "CHIT1", "SFRP2", "LDHA", "S100A10", "DCN", "GSN", "MMP9", "COL1A1", "COL1A2", "COL3A1", "FN1", "LUM")
	labFeatures <- c("SPP1", "FN1", "TREM2", "CD9", "FABP5", "FABP4", "LGALS3", "IL1RN", "CSTB", "LPL", 
		     "C1QA", "CEBPD", "HLA.DQA1", "HLA.DQB1", "HLA.DQA2", "FOSB", "IL1B", "HLA.DPB1")

	deRes <- de_res
	keyvals <- ifelse(deRes$avg_log2FC < -0.5, 'dodgerblue', 
			  ifelse(deRes$avg_log2FC > 0.5, 'firebrick', 'grey'))
	keyvals[is.na(keyvals)] <- 'grey'
	names(keyvals)[keyvals == 'firebrick'] <- 'Upregulated\nin PD-L1+ mo/macrophage'
	names(keyvals)[keyvals == 'grey'] <- 'Insiginificant genes'
	names(keyvals)[keyvals == 'dodgerblue'] <- 'Upregulated\nin PD-L1- mo/macrophage'

	ev <- EnhancedVolcano(deRes,
			      lab = rownames(deRes),
			      x = "avg_log2FC",
			      y = "p_val_adj",
			      title = "PD-L1+ vs PD-L1- mo/macrophage",
			      subtitle = "MAST",
			      colCustom = keyvals,
			      colAlpha = 0.9,
			      labSize = 3.6,
			      pointSize = 2.0,
			      xlim = c(-1.5, 1.5),
			      pCutoff = 0.10,
			      FCcutoff = 0.5,
			      legendPosition = 'top',
			      selectLab = labFeatures,
			      ylab = bquote(~ '-' ~Log[10]~ 'Adjusted P-value'),
			      drawConnectors = TRUE
	)
	png(paste(tamPf, "Fig1G_post_MAST_pos_vs_neg_ev.png", sep = ""), 
	    res = 300, height = 7.2, width = 9, units = 'in')
	print(ev)
	gar <- dev.off()


}

q(save = "no")
sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", 
	       "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
	      "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", 
	      "FOXP3", "TBX21", "GATA3")
metaItems <- c("patient", "tissue", "dataset", "er_status", "pr_status", "her2_status")

res_dir <- paste(sc_dir, "/breast_results", sep = "")
if (any(dataset_ids == "all")) {
	cat("All the available breast scRNA-seq will be used!\n")
	#TODO: need to read or find all the datasets
} else {
	expr_name <- paste(dataset_ids, collapse="_")
	tis_name <- paste(tissues, collapse = "_")
	expr_name <- paste(expr_name, tis_name, inte_prep, subtype, expr_id, sep = "_")
}
print(expr_name)

dir.create(file.path(res_dir, expr_name), showWarnings = FALSE)
expDir <- paste(res_dir, "/", expr_name, sep = "")
ppf <- paste(res_dir, "/", expr_name, "/", expr_name, "_", sep = "")

resFolder <- paste(expr_name, "_r", glb_res, sep = "")
resFolder <- paste("r", glb_res, sep = "")
rDir <- paste(res_dir, "/", expr_name, "/", resFolder, "/", sep = "")
rplotPf <- paste(rDir, "/", resFolder, "_", expr_name, "_", sep = "")
if (!evFlag) {
srsc <- readRDS(paste(rplotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
print(srsc)
}
subExpID <- paste("tam_glbr", glb_res, sep = "") #r1.5
dir.create(file.path(expDir, subExpID), showWarnings = FALSE)
subDir <- paste(expDir, "/", subExpID, "/", sep = "")
sppf <- paste(subDir, subExpID, "_", expr_name, "_", sep = "")

if (tamInteFlag) {
	subSrsc <- subset(srsc, idents = tamCluster)
	subSrList <- SplitObject(subSrsc, split.by="batch_factor")
	for (isn in names(subSrList)) {
		tmpCts <- subSrList[[isn]]@assays$RNA@counts
		tmpMeta <- subSrList[[isn]]@meta.data
		print(dim(tmpMeta))
		if (nrow(tmpMeta) >= 30) {
			colnames(tmpMeta) <- str_c("G_",colnames(tmpMeta))
			tmpSS <- CreateSeuratObject(counts = tmpCts, project = isn, meta.data=tmpMeta)
			plotPf <- paste(sppf, isn, "_", sep = "")
			if (inte_prep == "std") {
				cat("\tUsing the standard approach for integration (log-normalization)...\n")
				tmpSS <- NormalizeData(tmpSS, assay = "RNA", 
						       normalization.method = "LogNormalize", 
						       scale.factor = 10000)
				################
				cat("\tFind most variant features...\n")
				tmpSS <- FindVariableFeatures(tmpSS, selection.method = "vst", 
							      nfeatures = 2000)
				top10Features <- head(VariableFeatures(tmpSS), 18)
				varGG <- VariableFeaturePlot(tmpSS)
				labVargg <- LabelPoints(plot = varGG, points = top10Features, 
							repel = TRUE)
				ggsave(plot = labVargg, 
				       filename = paste(plotPf, "variable_features.png", sep = ""), 
				       dpi = 300, width = 9, height = 6)
			}
			if (inte_prep == "sct") {
				cat("\tUsing SCT approach for integration prep...\n")
				plan("multiprocess", workers = detectCores()/2)
				tmpSS <- SCTransform(tmpSS, verbose = TRUE)
			}
			subSrList[[isn]] <- tmpSS
		} else {
			subSrList[[isn]] <- NULL
		}
	}
	print(sppf)
	saveRDS(subSrList, paste(sppf, "preprocessed_seurat_list.RDS", sep = ""))
	if (inte_prep == "std") {
		ist = Sys.time()
		cat("Start to merge the samples with batch correction...\n")

		cat("Standard integration process (log-normalization)\n")
		integrateAnchors <- FindIntegrationAnchors(object.list = subSrList)
		inteObjs <- IntegrateData(anchorset = integrateAnchors, k.weight = 30)
		DefaultAssay(inteObjs) = "integrated"
	}
	if (inte_prep == "sct") {
		ist = Sys.time()
		cat("SCT-based integration process...\n")
		inteFeat <- SelectIntegrationFeatures(object.list = subSrList, nfeatures = 3000)
		subSrList <- PrepSCTIntegration(object.list = subSrList, 
						anchor.features = inteFeat, verbose = TRUE)
		plan("multiprocess", workers = detectCores()/2)
		inteAnchors <- FindIntegrationAnchors(object.list = subSrList, 
						      normalization.method = "SCT", 
						      anchor.features = inteFeat, 
						      verbose = TRUE)
		inteObjs <- IntegrateData(anchorset = inteAnchors, k.weight = 30,
					  normalization.method = "SCT", verbose = TRUE)
	}
	saveRDS(inteObjs, file = paste(sppf, "integrated_seurat_object.RDS", sep = ""))
	print(inteObjs)
	cat("Total integration cost")
	print(Sys.time()-ist)
}

if (tamPostFlag) {
	st <- Sys.time()
#	if (splitFlag) {rdsFile <- paste(expDir, expID, "_", intePref, "_split.RDS", sep = "")}
	inteObjs <- readRDS(paste(sppf, "integrated_seurat_object.RDS", sep = ""))
	print(inteObjs)
	inteProObj <- srPostPro(inteObjs, plotPf = sppf, expDir = subDir, 
				rdsSave = TRUE, ndim = 30, exprName = paste(expr_name, subExpID, sep = "_"),
				jsFlag = FALSE, qcFlag = FALSE, tamFlag = TRUE,
				plotFeat = c("CD68", "HLA.DRA", "CD274", "SIGLEC15"),
				resGrid = seq(from=0.2, to=3.0, by=0.2),
				metaPlot = str_c("G_", metaItems))
	print(inteProObj)
	print(Sys.time() - st)
}

if (tamCompFlag) {
	tam_res_folder <- paste("r", tam_res, sep = "")
	tamSubPf <- paste(subDir, tam_res_folder,  "/", 
			  tam_res_folder, "_", expr_name, "_", subExpID,"_",sep = "")
	tamRDSFile <- paste(tamSubPf, "Seurat_Objects_Clustered.RDS", sep = "")
	proObj <- readRDS(tamRDSFile)
	print(proObj)

	## List for annotate PD-L1 status on TAM
	tamTypes <- c("PDL1pos", "PDL1neg")
	tamAnnList <- vector(mode = "list", length = length(tamTypes))
	names(tamAnnList) <- tamTypes
	# r0.6
	tamAnnList[["PDL1pos"]] <- c(0,2,3,4,5,12)
	tamAnnList[["PDL1neg"]] <- c(1,6,8,11)

	proObj@meta.data$PDL1 <- "NonTAM"
	for (itt in names(tamAnnList)) {
		tmpMask <- proObj@meta.data$seurat_clusters %in% tamAnnList[[itt]]
		proObj@meta.data[tmpMask, "PDL1"] <- itt
	}

	## Big group UMAP
	pdl1DimGG <- DimPlot(proObj, reduction = "umap", 
			     group.by = "PDL1", cols = c("grey", "dodgerblue", "indianred"))
	ggsave(plot = pdl1DimGG, 
	       filename = paste(tamSubPf, "post_PDL1_status_dimplot.png", sep = ""),
	       dpi = 300, width = 6, height = 4)
	proObj@meta.data$G_seurat_clusters <- as.factor(proObj@meta.data$G_seurat_clusters)
	pdl1DimGG <- DimPlot(proObj, reduction = "umap", group.by = "G_seurat_clusters")
	ggsave(plot = pdl1DimGG, 
	       filename = paste(tamSubPf, "G_seurat_cluster_dimplot.png", sep = ""),
	       dpi = 300, width = 6, height = 4)

	## Subset 
	tamObj <- proObj
	proObj <- subset(proObj, subset = PDL1 == "NonTAM", invert = TRUE)
	pdl1Col <- c("white", "dodgerblue", "indianred")
	pdl1Col <- c("dodgerblue", "indianred")

	pdl1DimGG <- DimPlot(proObj, reduction = "umap", group.by = "PDL1", cols = pdl1Col) +
		xlim(-5, 7.5)
	ggsave(plot = pdl1DimGG, 
	       filename = paste(tamSubPf, "post_PDL1_status_dimplot_clean.png", sep = ""),
	       dpi = 300, width = 6, height = 4)

	pdl1VlnGG <- VlnPlot(proObj, features = c("rna_CD274", "rna_SIGLEC15"), 
			     group.by = "PDL1", log = TRUE)
	ggsave(plot = pdl1VlnGG, 
	       filename = paste(tamSubPf, "post_PDL1_status_VlnPlot.png", sep = ""),
	       dpi = 300, width = 6, height = 4)
	pdl1_cts <- proObj@meta.data %>% 
		group_by(G_patient, PDL1) %>% 
		summarise(n = n()) %>%
		group_by(G_patient) %>%
		mutate(per = n/sum(n))
	write.csv(pdl1_cts, paste(tamSubPf, "pdl1_cell_counts.csv", sep = ""))


	cat("Start to run DE analysis -- Wilcox\n")
	deRes <- FindMarkers(proObj, ident.1 = "PDL1pos", ident.2 = "PDL1neg", 
			     group.by = "PDL1", test.use = "wilcox")
	write.csv(deRes, paste(tamSubPf, "post_wilcox_pos_vs_neg_PDL1_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_log2FC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_log2FC >= 0.5 | deRes$avg_log2FC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(tamSubPf, "post_wilcox_marker_heatmap.png", sep = ""), 
	    res = 300, width = 9, height = 12, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "PDL1", label = FALSE, 
			group.colors = pdl1Col))
	gar <- dev.off()

	cat("Start to run DE analysis -- MAST\n")
	deRes <- FindMarkers(proObj, ident.1 = "PDL1pos", ident.2 = "PDL1neg", 
			     group.by = "PDL1", test.use = "MAST")
	write.csv(deRes, paste(tamSubPf, "post_MAST_pos_vs_neg_PDL1_res.csv", sep = ""))
	deRes <- deRes[order(deRes$avg_log2FC),]
	sigMask <- deRes$p_val_adj <= 0.10 &(deRes$avg_log2FC >= 0.5 | deRes$avg_log2FC <= -0.5)
	sigFeat <- rownames(deRes)[sigMask]
	png(paste(tamSubPf, "post_MAST_marker_heatmap.png", sep = ""), 
	    res = 300, width = 9, height = 12, units = "in")
	print(DoHeatmap(object = proObj, features = sigFeat, group.by = "PDL1", label = FALSE, 
			group.colors = pdl1Col))
	gar <- dev.off()


	## Output for CIBERSORTx
	rnaExpr <- as.matrix(proObj[["RNA"]]@data)
	allGene <- rownames(rnaExpr)
#	print(rnaExpr[1:9, 1:6])
	metaData <- t(proObj@meta.data)
	cbstxOutDf <- rbind(metaData, rnaExpr)
#	print(cbstxOutDf[1:9,1:6])
#	print(dim(rnaExpr))
#	print(length(allGene))
#	cbstxOutDf <- FetchData(proObj, vars = c("PDL1", "seurat_clusters", allGene))
#	write.csv(cbstxOutDf, paste(tamSubPf, "post_PDL1_for_cbstx.csv", sep = ""))
}

if (tamVisFlag) {
	tam_res_folder <- paste("r", tam_res, sep = "")
	tamSubPf <- paste(subDir, tam_res_folder,  "/", 
			  tam_res_folder, "_", expr_name, "_", subExpID,"_",sep = "")
	tamRDSFile <- paste(tamSubPf, "Seurat_Objects_Clustered.RDS", sep = "")
	proObj <- readRDS(tamRDSFile)
	print(proObj)

	## List for annotate PD-L1 status on TAM
	tamTypes <- c("PDL1pos", "PDL1neg")
	tamAnnList <- vector(mode = "list", length = length(tamTypes))
	names(tamAnnList) <- tamTypes
	# r0.6
	tamAnnList[["PDL1pos"]] <- c(0,2,3,4,5,12)
	tamAnnList[["PDL1neg"]] <- c(1,6,8,11)
	pdl1Col <- c("dodgerblue", "indianred")

	proObj@meta.data$PDL1 <- "NonTAM"
	for (itt in names(tamAnnList)) {
		tmpMask <- proObj@meta.data$seurat_clusters %in% tamAnnList[[itt]]
		proObj@meta.data[tmpMask, "PDL1"] <- itt
	}
	proObj <- subset(proObj, subset = PDL1 == "NonTAM", invert = TRUE)
	print(proObj)

	hmFeats <- c("SPP1", "FN1", "TREM2", "CD9", "FABP5", "FABP4", "LGALS3", "IL1RN", "CSTB", "LPL", 
		     "C1QA", "CEBPD", "HLA.DQA1", "HLA.DQB1", "HLA.DQA2", "FOSB", "IL1B", "HLA.DPB1")

	proObj <- ScaleData(proObj, assay = "RNA")
	png(paste(tamSubPf, "manual_selected_marker_heatmap_v2.png", sep = ""), 
	    res = 300, width = 6, height = 6, units = "in")
	print(DoHeatmap(object = proObj, features = hmFeats, group.by = "PDL1", label = FALSE, 
			group.colors = pdl1Col, assay = "RNA"))
	gar <- dev.off()
	for (ihf in hmFeats) {
		pdl1VlnGG <- VlnPlot(proObj, features = ihf, group.by = "PDL1", log = TRUE, pt.size = 0,
				     cols = c('dodgerblue', 'firebrick'))
		ggsave(plot = pdl1VlnGG, filename = paste(tamSubPf, ihf, "_manual_selected_VlnPlot.png", sep = ""),
		       dpi = 300, width = 3, height = 4)

		pdl1Siglec15GG <- FeaturePlot(proObj, features = ihf, order = TRUE, 
					      cols = c("honeydew", "purple")) + xlim(-5, 7.5)
		ggsave(plot = pdl1Siglec15GG, filename = paste(tamSubPf, ihf, "_umap_clean.png", sep = ""),
		       dpi = 300, width = 6, height = 4)

	}

	hmDotPlot <- DotPlot(proObj, features = hmFeats, cols = c("honeydew", "purple"), group.by = "PDL1", assay = "RNA") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(plot = hmDotPlot, filename = paste(tamSubPf, "manual_selected_marker_dotplot.png", sep = ""),
	       dpi = 300, width = 12, height = 4)

}

if (evFlag) { # EnhancedVolcanoPlots
	library(EnhancedVolcano)
	## Run DE analysis between PD-L1+/PD-L1- TAMs (Highly customized!)
	cat("Visualize differential expression results between PD-L1+/- TAMs...\n")
	tam_res_folder <- paste("r", tam_res, sep = "")
	subPrefix <- paste(subDir, tam_res_folder,  "/", 
			  tam_res_folder, "_", expr_name, "_", subExpID,"_",sep = "")

	hmFeats <- c("SPP1", "FN1", "TREM2", "CD9", "FABP5", "FABP4", "LGALS3", "IL1RN", "CSTB", "LPL", 
		     "C1QA", "CEBPD", "HLA.DQA1", "HLA.DQB1", "HLA.DQA2", "FOSB", "IL1B", "HLA.DPB1")

	deRes <- read.csv(paste(subPrefix, "post_MAST_pos_vs_neg_PDL1_res.csv", sep = ""), 
			  check.names = F, header = T, row.names = 1)
	print(head(deRes))
	
	keyvals <- ifelse(deRes$avg_log2FC < -0.5, 'dodgerblue', 
			  ifelse(deRes$avg_log2FC > 0.5, 'firebrick', 'grey'))
	keyvals[is.na(keyvals)] <- 'grey'
	names(keyvals)[keyvals == 'firebrick'] <- 'Upregulated\nin PD-L1+ TAM'
	names(keyvals)[keyvals == 'grey'] <- 'Insiginificant genes'
	names(keyvals)[keyvals == 'dodgerblue'] <- 'Upregulated\nin PD-L1- TAM'

	ev <- EnhancedVolcano(deRes,
			      lab = rownames(deRes),
			      x = "avg_log2FC",
			      y = "p_val_adj",
			      title = "PD-L1+ vs PD-L1- TAM",
			      subtitle = "MAST",
			      colCustom = keyvals,
			      colAlpha = 0.9,
			      labSize = 3.6,
			      pointSize = 2.0,
			      xlim = c(-1.5, 1.5),
			      pCutoff = 0.10,
			      FCcutoff = 0.5,
			      legendPosition = 'top',
			      selectLab = hmFeats,
			      ylab = bquote(~ '-' ~Log[10]~ 'Adjusted P-value'),
			      drawConnectors = TRUE
	)
	png(paste(subPrefix, "post_MAST_pos_vs_neg_ev.png", sep = ""), 
	    res = 300, height = 7.2, width = 9, units = 'in')
	print(ev)
	gar <- dev.off()

}

