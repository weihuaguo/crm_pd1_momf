### This script is designed to analyze public scRNA-seq data Bassez et al., Nat Med 2021
### data resource http://biokey.lambrechtslab.org
### Zhikun Guo, WL laboratory, Shenzhen University
### Cleaned on 2023/11/23

rm(list = ls())
suppressMessages(library(Seurat))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MAST))
suppressMessages(library(clustree))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

set.seed(1)

workDir <- "..." # Please change this directory
dataDir <- paste(workDir, "1_Nature_scrna_seq", sep = "/") 
resDir <- paste(workDir, "2_Nature_scrna_seq_result", sep = "/")

expID <- "public_NatMed_TNBC"
readFlag <- FALSE
QCFlag <- FALSE
clstFlag <- FALSE
majorVisFlag <- FALSE
tamClusterFlag <- FALSE
tamVisFlag <- TRUE

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, sep = "/")

rnaPf <- paste(expDir, "/", expID, "_RNA_Results_", sep = "")


if (readFlag) {
  cat("Start to load Seurat object...\n")
  coh1 <- readRDS(paste(dataDir,"1863-counts_cells_cohort1.rds",sep = "/"))
  coh1 <- CreateSeuratObject(coh1,project = "bc",assay = "RNA")
  meta1 <- read.csv(paste(dataDir,"1872-BIOKEY_metaData_cohort1_web.csv",sep = "/"),row.names = 1)
  coh1 <- AddMetaData(coh1,metadata = meta1)
  
  coh2 <- readRDS(paste(dataDir,"1867-counts_cells_cohort2.rds",sep = "/"))
  coh2 <- CreateSeuratObject(coh2,project = "bc",assay = "RNA")
  meta2 <- read.csv(paste(dataDir,"1871-BIOKEY_metaData_cohort2_web.csv",sep = "/"),row.names = 1)
  coh2 <- AddMetaData(coh2,metadata = meta2)
  
  rawObjRDS <- merge(x = coh1,y = coh2)
  rawObjRDS <- subset(rawObjRDS,subset=BC_type=="TNBC")
  saveRDS(object = rawObjRDS, file = paste(expDir, "/", expID, "_unprocessed_seurat_object.RDS", sep = ""))
}

if(QCFlag){
  cat("Start to QC and pre-processing Seurat object...\n")
  rawObjRDS <- readRDS(file = paste(expDir, "/", expID, "_unprocessed_seurat_object.RDS", sep = ""))
  print(rawObjRDS)
  rawObjRDS[["percent.mt"]] <- PercentageFeatureSet(rawObjRDS, pattern = "^MT-")

  qcgg = VlnPlot(rawObjRDS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,raster=FALSE)
  ggsave(qcgg,filename = paste(expDir, "/", expID,"_QCPlot1.png",sep = ""),dpi = 300,width = 9, height = 6,units = "in")

  qcsc1 <- FeatureScatter(rawObjRDS, feature1 = "nCount_RNA", feature2 = "percent.mt")
  qcsc2 <- FeatureScatter(rawObjRDS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  cqcgg <- CombinePlots(plots = list(qcsc1, qcsc2))+theme_bw()
  ggsave(plot = cqcgg, filename = paste(expDir, "/", expID,"_QCPlot2.png", sep = ""), dpi = 300, width = 9, height = 6,units = "in")

  rawObjRDS <- subset(rawObjRDS, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA > 5 & percent.mt < 15)
  rawObjRDS <- NormalizeData(rawObjRDS, normalization.method = "LogNormalize", scale.factor = 10000)

  rawObjRDS <- FindVariableFeatures(rawObjRDS, selection.method = "vst", nfeatures = 2000)
  top10Features <- head(VariableFeatures(rawObjRDS), 10)
  varGG <- VariableFeaturePlot(rawObjRDS)
  labVargg <- LabelPoints(plot = varGG, points = top10Features, repel = TRUE)+theme_bw()
  ggsave(plot = labVargg,filename = paste(expDir, "/", expID,"_variable_features.png", sep = ""),dpi = 300,width = 9, height = 6,units = "in")

  all.genes <- rownames(rawObjRDS)
  rawObjRDS <- ScaleData(rawObjRDS, features = all.genes)
  rawObjRDS <- RunPCA(rawObjRDS, features = VariableFeatures(object = rawObjRDS))

  print(rawObjRDS[["pca"]], dims = 1:5, nfeatures = 5)
  pcaGG = DimPlot(rawObjRDS, reduction = "pca")
  ggsave(plot = pcaGG,filename = paste(expDir, "/", expID,"_pca_plot.png", sep = ""),dpi = 300,height = 5,width = 6,units = "in")
  
  png(filename = paste(expDir, "/", expID,"_PCA_HEATMAP.png", sep = ""),res = 300,width = 9, height = 16, units = "in")
  DimHeatmap(rawObjRDS, dims = 1:9, cells = 500, balanced = TRUE)
  dev.off()
  
  rawObjRDS <- JackStraw(rawObjRDS, num.replicate = 10)
  rawObjRDS <- ScoreJackStraw(rawObjRDS, dims = 1:20)
  jsGG <- JackStrawPlot(rawObjRDS, dims = 1:15)+theme_bw()
  ggsave(plot = jsGG, filename = paste(expDir, "/", expID,"_JackStraw.png", sep = ""), dpi = 300, width = 9, heigh = 6, units = "in")
  
  elGG <- ElbowPlot(rawObjRDS)+theme_bw()
  ggsave(plot = elGG, filename = paste(expDir, "/", expID,"_Elbow.png", sep = ""), dpi = 300, width = 9, heigh = 6, units = "in")

  rawObjRDS <- FindNeighbors(rawObjRDS, dims = 1:10)
  rawObjRDS <- FindClusters(rawObjRDS, resolution = 0.2)
  rawObjRDS <- RunUMAP(rawObjRDS, dims = 1:10)
  
  umapGG <- DimPlot(rawObjRDS, reduction = "umap", label = TRUE)
  ggsave(plot = umapGG, filename = paste(rnaPf, "UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6, units = "in")
  umapGG <- DimPlot(rawObjRDS, reduction = "umap", label = TRUE,group.by = "cellType")
  ggsave(plot = umapGG, filename = paste(rnaPf, "UMAP_cellType.png", sep = ""), dpi = 300, width = 9, heigh = 6, units = "in")
  
  saveRDS(rawObjRDS, file = paste(rnaPf,"Seurat_Objects_Clustered.RDS", sep = ""))
}
  
if (clstFlag) {
  cat("Analyzing Seurat object...\n")
  srObjRDS <- readRDS(file = paste(rnaPf,"Seurat_Objects_Clustered.RDS", sep = ""))
  
  Idents(srObjRDS) <- srObjRDS$cellType ####Add
  srObjRDS.markers <- FindAllMarkers(srObjRDS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  topMarkers <- srObjRDS.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  write.csv(srObjRDS.markers, file = paste(rnaPf, "ALL_POS_MARKERS.csv", sep = "")) 
  write.csv(topMarkers, file = paste(rnaPf, "top10_pos_markers.csv", sep = "")) 
  
  topMarkers <- srObjRDS.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  hm <- DoHeatmap(object = srObjRDS, features = topMarkers$gene) + NoLegend()
  ggsave(paste(rnaPf, "top5_marker_heatmap.png", sep = ""), hm, dpi = 300, width = 24, height = 18)
  
  topDotMarkers <- srObjRDS.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  geneMarkers <- unique(as.character(topDotMarkers$gene))
  
  dotMarkerGG <- DotPlot(srObjRDS, features = geneMarkers) + RotatedAxis() + coord_flip()+theme_classic()
  ggsave(plot = dotMarkerGG, filename = paste(rnaPf, "Top5Marker_DotPlot.png", sep = ""), 
         dpi = 300, width = 9, height = 16)
}

if (majorVisFlag) {
  cat("Merge and visualize all the cells\n")
  proObj <- readRDS(paste(rnaPf, "Seurat_Objects_Clustered.RDS", sep = ""))
  cellAnnDf <- as.data.frame(read_excel(paste(rnaPf, "ALL_POS_MARKERS.xlsx", sep = "")))
  
  rownames(cellAnnDf) <- cellAnnDf[,1]
  cellAnnDf[,1] <- NULL
  cellMarkerDf <- cellAnnDf
  cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$cell_type),]
  print(head(cellAnnDf))
  proObj@meta.data$major_cell_type <- "NA"
  for (ic in cellAnnDf$cluster) {
    proObj@meta.data$major_cell_type[proObj@meta.data$cellType == ic] <- cellAnnDf$cell_type[cellAnnDf$cluster == ic]
  }
  print(head(proObj@meta.data))
  
  sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "CD1C", "KIT", "ITGAM")
  sigExpr <- proObj[['RNA']]@data[sigMarkers,]
  sigExpr <- as.data.frame(sigExpr)
  sigExpr$gene <- rownames(sigExpr)
  
  gathExpr <- gather(sigExpr, "cell_id", "expr", rownames(proObj@meta.data))
  gathExpr <- merge(gathExpr, proObj@meta.data, by.x = "cell_id", by.y = "row.names", all.x = TRUE)
  gathExpr$gene <- factor(gathExpr$gene, levels = sigMarkers)
  
  ####FigS9D####
  figs9d <- ggplot(gathExpr, aes(x = major_cell_type, y = expr, fill = major_cell_type)) +
  	geom_violin(scale = "width", trim = TRUE) +
  	facet_wrap(.~gene, nrow = 1, scales = "free_x") +
  	scale_fill_brewer(palette = "Spectral") +
  	labs(x = "General cell type in TME", y = "Normalized expression", fill = "General cell type\nin TME") +
  	coord_flip() +
  	theme_classic()
  ggsave(plot = figs9d, filename = paste(rnaPf, "FigS9D_mc_violin_plot_for_cell_annotation.png", sep = ""),
         dpi = 600, width = 16, height = 2.5) 
  
  cat("Heatmaps\n")
  topMarkers <- cellMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  geneExpr <- AverageExpression(proObj,features = topMarkers$gene,group.by = "major_cell_type") 
  geneExpr <- as.data.frame(geneExpr$RNA)
  cellType <- data.frame(colnames(geneExpr))
  colnames(cellType) <- 'class'
  cols <- brewer.pal(length(unique(topMarkers$cell_type)), "Spectral")
  top_anno = HeatmapAnnotation(df = cellType,
                               border = T,
                               show_annotation_name = F,
                               gp = gpar(col = 'black'),
                               col = list(class = c("B cell" = "#D53E4F",
                                                    "Cancer cell" = "#F46D43",
                                                    "DC" = "#FDAE61",
                                                    "Endothelium" = "#FEE08B",
                                                    "Fibroblast" = "#E6F598",
                                                    "Mast cell" = "#ABDDA4",
                                                    "Monocyte/Macrophage" = "#66C2A5",
                                                    "T cell" = "#3288BD")))

  marker_exp <- t(scale(t(geneExpr),scale = T,center = T))
  ####FigS9C####
  png(filename = paste(rnaPf, "FigS9C_top5_marker_avg_heatmap.png", sep = ""),res = 600, width = 6, height = 8,units = "in")
  Heatmap(marker_exp,
          cluster_rows = F,
          cluster_columns = F,
          show_column_names = F,
          show_row_names = T,
          heatmap_legend_param = list(title = ""),
          col = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
          border = 'black',
          rect_gp = gpar(col = "black", lwd = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          top_annotation = top_anno)
  dev.off()
  
  cat("Cell proportions\n")
  mcCts <- proObj@meta.data %>%
  	group_by(patient_id, major_cell_type) %>%
  	summarize(n = n())
  print(head(mcCts))

  ####FigS9B####
  figs9b <- ggplot(mcCts, aes(x = patient_id, y = n, fill = major_cell_type)) +
  	geom_bar(stat = "identity", position = "fill", color = "gray") +
  	scale_fill_brewer(palette = "Spectral") +
  	coord_flip() +
  	labs(x = "Sample ID", y = "Cell proportion", fill = "General cell type\nin TME") +
  	theme_classic()
  ggsave(plot = figs9b, filename = paste(rnaPf, "FigS9B_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)

  ####FigS9A####
  cat("UMAPs\n")
  figs9a <- DimPlot(proObj, reduction = "umap", label = TRUE, group.by = "major_cell_type") + 
    labs(color = "General cell type\nin TME", title = "Annotated and merged cell types") +
    scale_color_brewer(palette = "Spectral")
  ggsave(plot = figs9a, filename = paste(rnaPf, "FigS9A_UMAP.png", sep = ""), dpi = 600, width = 9, heigh = 6)
  
}

if (tamClusterFlag) {
  cat("Cluster TAMs...\n")
  ist <- Sys.time()
  subSrsc <- subset(proObj, subset=major_cell_type=="Monocyte/Macrophage") 

  subFolder <- "tam_sub_cluster"
  dir.create(file.path(expDir, subFolder), showWarnings = FALSE)
  subDir <- paste(expDir, subFolder, sep = "/")
  subPrefix <- paste(subDir, "/", expID, "_", subFolder, "_", sep = "")
 
  pst <- Sys.time()
  ress <- c(seq(0.1,0.9, 0.1),seq(1.0,1.2,0.2))
  geneOi <- c("CD14", "CD68", "HLA-DRA", "CD274", "SIGLEC15", "PDCD1LG2", "CD1C", "ITGAX", "FCER1A", "TPSAB1", "KIT", "FUT4")
  
  subSrsc <- JackStraw(subSrsc, num.replicate = 20)
  subSrsc <- FindNeighbors(subSrsc, dims = 1:10)
  subSrsc <- RunUMAP(subSrsc, dims = 1:10)
  
  bldGG <- FeaturePlot(subSrsc, features = c("CD274", "SIGLEC15"), blend = TRUE)
  ggsave(paste(subPrefix, "blend_pdl1_siglec15_expr_umap.png", sep = ""), bldGG, dpi = 300, width = 25, height = 6)
  
  for (ires in ress) {
    plotFeat <- geneOi
    resFolder <- paste("tam_subcluster_res", ires, sep = "_")
    dir.create(file.path(subDir, resFolder), showWarnings = FALSE)
    resDir <- paste(subDir, resFolder, sep = "/")
    rppf <- paste(resDir, "/", resFolder, "_", sep = "")
    
    subSrsc <- FindClusters(subSrsc, resolution = ires)
    umapGG <- DimPlot(subSrsc, reduction = "umap", label = TRUE,group.by = paste0("RNA_snn_res.",ires))
    ggsave(plot = umapGG, filename = paste(rppf, "UMAP.png", sep = ""), dpi = 300, width = 9, heigh = 6)
    
    if (!is.null(plotFeat)) {
      cat("Start to plot feature plots and violin plots\n")
      vlnGG <- VlnPlot(subSrsc, features = plotFeat, slot = "data", log = FALSE, ncol = 3, pt.size = 0.5,)
      ggsave(plot = vlnGG, filename = paste(rppf, "violin_plot_for_cell_annotation.png", sep = ""), 
             dpi = 300, width = 16, height = 16) 
      
      featGG <- FeaturePlot(subSrsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
      ggsave(plot = featGG, filename = paste(rppf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
             dpi = 300, width = 16, height = 16)
      
      dotGG <- DotPlot(subSrsc, features = plotFeat) + RotatedAxis()+theme_classic()
      ggsave(plot = dotGG, filename = paste(rppf, "DotPlot_for_cell_annotation.png", sep = ""),
             dpi = 300, width = 9, height = 6)
      
      rdgGG <- RidgePlot(subSrsc, features = plotFeat, ncol = 3)
      ggsave(paste(rppf, "ridge_plot_for_cell_annotation.png", sep = ""), rdgGG, dpi = 300, width = 16, height = 18)
      
    }
  }
  saveRDS(subSrsc, paste(subPrefix, "Seurat_Objects_Clustered.RDS", sep = ""))
  
  ####FigS10E####
  all_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(ress)+1)
  figs10e_cltr <- clustree(subSrsc, prefix = c("RNA_snn_res.0"),node_text_size = 4) + scale_color_manual(values = all_cols)
  ggsave(paste(subPrefix, "FigS10E_clustree.png", sep = ""), figs10e_cltr, dpi = 600, width = 8, height = 10)
  
  cat("Analysis time cost:")
  print(Sys.time()-pst)
  cat("++++++++++++++++++++++++++++++++++++++++++++++++\n\n")
}

####Find DEG##########
resFolder <- paste("tam_subcluster_res", 0.8, sep = "_")
dir.create(file.path(subDir, resFolder), showWarnings = FALSE)
resDir <- paste(subDir, resFolder, sep = "/")
rppf <- paste(resDir, "/", resFolder, "_", sep = "")

Idents(subSrsc) <- subSrsc$RNA_snn_res.0.8
subSrsc.markers <- FindAllMarkers(subSrsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers <- subSrsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(subSrsc.markers, file = paste(rppf, "ALL_POS_MARKERS.csv", sep = "")) 
write.csv(topMarkers, file = paste(rppf, "top10_pos_markers.csv", sep = ""))

topDotMarkers <- subSrsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneMarkers <- unique(as.character(topDotMarkers$gene))

dotMarkerGG <- DotPlot(subSrsc, features = geneMarkers) + RotatedAxis() + coord_flip()+theme_classic()
ggsave(plot = dotMarkerGG, filename = paste(rppf, "Top5Marker_DotPlot.png", sep = ""),
         dpi = 300, width = 9, height = 16)


if (tamVisFlag) {
  use_res <- 0.8
  subFolder <- "tam_sub_cluster"
  subDir <- paste(expDir, subFolder, sep = "/")
  resFolder <- paste("tam_subcluster_res", use_res, sep = "_")
  tamPf <- paste(expDir, "/", subFolder, "/", resFolder, "/", resFolder, "_", sep = "")
  tamObj <- readRDS(paste(tamPf, "Seurat_Objects_Clustered.RDS", sep = ""))
  Idents(tamObj) <- tamObj$RNA_snn_res.0.8
  print(tamObj)
  cellAnnDf <- as.data.frame(read_excel(paste(tamPf, "ALL_POS_MARKERS.xlsx", sep = "")))
  rownames(cellAnnDf) <- cellAnnDf[,1]
  cellAnnDf[,1] <- NULL
  cellMarkerDf <- cellAnnDf
  cellAnnDf <- cellAnnDf[!is.na(cellAnnDf$tam_type),]
  print(head(cellAnnDf))
  
  cols <- brewer.pal(length(unique(cellAnnDf$tam_type)), "Set2")
  names(cols) <- unique(cellAnnDf$tam_type)
  cols['PD-L1- TAM'] <- "dodgerblue"
  cols['PD-L1+ TAM'] <- "firebrick"
  
  tamObj@meta.data$tam_type <- "NA"
  for (ic in cellAnnDf$cluster) {
    tamObj@meta.data$tam_type[tamObj@meta.data$RNA_snn_res.0.8 == ic] <- cellAnnDf$tam_type[cellAnnDf$cluster == ic]
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
  geneExpr <- AverageExpression(tamObj,features = topMarkers$gene,group.by = "RNA_snn_res.0.8") 
  geneExpr <- as.data.frame(geneExpr$RNA)
  cellType <- data.frame(colnames(geneExpr))
  colnames(cellType) <- 'cluster'
  top_anno = HeatmapAnnotation(df = cellType,
                               border = T,
                               show_annotation_name = F,
                               gp = gpar(col = 'black'),
                               annotation_legend_param = list(at = c(0:14)),
                               col = list(cluster = c("0" = "#F8766D",
                                                    "1" = "#E58700",
                                                    "2" = "#C99800",
                                                    "3" = "#A3A500",
                                                    "4" = "#6BB100",
                                                    "5" = "#00BA38",
                                                    "6" = "#00BF7D",
                                                    "7" = "#00C0AF",
                                                    "8" = "#00BCD8",
                                                    "9" = "#00B0F6",
                                                    "10" = "#619CFF",
                                                    "11" = "#B983FF",
                                                    "12" = "#E76BF3",
                                                    "13" = "#FD61D1",
                                                    "14" = "#FF67A4")))
  
  
  marker_exp <- t(scale(t(geneExpr),scale = T,center = T))
  ####FigS10D####
  png(filename = paste(tamPf, "FigS10D_top5_marker_avg_heatmap.png", sep = ""),res = 600, width = 6, height = 8,units = "in")
  Heatmap(marker_exp,
          cluster_rows = F,
          cluster_columns = F,
          show_column_names = F,
          show_row_names = T,
          heatmap_legend_param = list(title = ""),
          col = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
          border = 'black',
          rect_gp = gpar(col = "black", lwd = 1),
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          top_annotation = top_anno)
  dev.off()
  
  ####FigS10B####
  figs10b <- FeaturePlot(tamObj, features = c("CD274", "SIGLEC15"), blend = T)
  ggsave(plot = figs10b, filename = paste(tamPf, "FigS10B_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 25, heigh = 6)
  
  cat("UMAPs\n")
  ####FigS10A####
  figs10a <- DimPlot(tamObj, reduction = "umap", label = FALSE) + labs(color = "Clusters (r = 0.8)", title = "")
  ggsave(plot = figs10a, filename = paste(tamPf, "FigS10A_UMAP_all_tam_cluster.png", sep = ""), dpi = 300, width = 9, heigh = 6)
  
  ####Fig2G####
  fig2g <- DimPlot(tamObj, reduction = "umap", label = FALSE, group.by = "tam_type") + 
    labs(color = "TAM types", title = "Annotated and merged cell types") +
    scale_color_manual(values = cols)
  ggsave(plot = fig2g, filename = paste(tamPf, "Fig2G_UMAP_all_tam_pdl1.png", sep = ""), dpi = 600, width = 9, heigh = 6)
  
  ####Fig2H####
  fig2h <- DimPlot(tamObj, reduction = "umap", label = FALSE, group.by = "expansion",cols = c("purple","grey","grey"),pt.size = 0.1) + 
    labs(color = "Expansion", title = "Annotated and merged cell types") 
  ggsave(plot = fig2h, filename = paste(tamPf, "Fig2H_UMAP_all_tam_expansion.png", sep = ""), dpi = 600, width = 9, heigh = 6)
 
  ####FigS10F####
  figs10f <- VlnPlot(tamObj,features = "CD274",group.by = "RNA_snn_res.0.8",pt.size = 0.001,sort = T)
  ggsave(plot = figs10f, filename = paste(tamPf, "FigS10F_all_tam_PDL1_sort.png", sep = ""), dpi = 600, width = 9, heigh = 6)
  
  vlnGenes <- c("CD274")
  ####FigS10G####
  for (ivg in vlnGenes) {
    cat(ivg, "\n")
    pdl1VlnGG <- VlnPlot(tamObj, features = ivg, group.by = "tam_type", log = FALSE, pt.size = 0.001,
                         cols = c('dodgerblue', 'firebrick')) + scale_y_continuous(limits = c(0.0, NA)) + 
      stat_compare_means(comparisons = list(c("PD-L1+ TAM", "PD-L1- TAM")), label = "p.signif")
    ggsave(plot = pdl1VlnGG, filename = paste(tamPf, ivg, "_FigS10G_VlnPlot.png", sep = ""),
           dpi = 300, width = 3, height = 5)
  }
}

###finished---------------------------------------------
