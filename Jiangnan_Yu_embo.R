# This script is designed to analyze EMBO TNBC scRNA-seq data 
# Jiangnan Yu, WL Lab.
# 23/11/2023
# Data resouce: GSE161529
# RDS files resource: https://www.nature.com/articles/s41597-022-01236-2/tables/3

####begin####
rm(list = ls())
setwd("...") # Please change this directory
tst <- Sys.time()
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
options(warn = 1)

####*emboObj*####
###embo_sce.big###
embo_sce.big = readRDS("tnbc_sce.big.RDS")

###QC###
embo_sce.big[["percent.mt"]] <- PercentageFeatureSet(embo_sce.big, pattern = "^MT-")
VlnPlot(embo_sce.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(embo_sce.big, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(embo_sce.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
embo_sce.big <- subset(embo_sce.big, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 25)

###Normalization###
embo_sce.big <- NormalizeData(embo_sce.big, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
embo_sce.big <- FindVariableFeatures(embo_sce.big, selection.method = "vst", nfeatures = 2000)
top10Features <- head(VariableFeatures(embo_sce.big), 10) ###changed parameter
varGG <- VariableFeaturePlot(embo_sce.big)
ggsave(varGG,"variable_features.png",dpi = res, width = 9, height = 6)
all.genes <- rownames(embo_sce.big)
embo_sce.big <- ScaleData(embo_sce.big, features = all.genes)

###Reduction###
embo_sce.big <- RunPCA(embo_sce.big, features = VariableFeatures(object = embo_sce.big))
set.seed(1)
embo_sce.big <- JackStraw(embo_sce.big, num.replicate = 10)
embo_sce.big <- ScoreJackStraw(embo_sce.big, dims = 1:20)
embo_sce.big<- embo_sce.big %>% 
  RunUMAP(embo_sce.big, dims = 1:30) %>% 
  RunTSNE(embo_sce.big, dims = 1:30) %>% 
  FindNeighbors(embo_sce.big, dims = 1:30) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()

###Cluster markers###
res = 300
sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "NCAM1", "KIT", "FCER1A")
FeaturePlot(embo_sce.big, features = sigMarkers, reduction = "umap", ncol = 3)
ggsave("umap_featureplot.png",dpi = 400, width = 12, heigh = 18)
srsc.markers <- FindAllMarkers(embo_sce.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(srsc.markers, file = ("ALL_POS_MARKERS.csv")) 
write.csv(topMarkers, file = ("top10_pos_markers.csv")) 
topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneMarkers <- unique(as.character(topDotMarkers$gene))
dotMarkerGG <- DotPlot(embo_sce.big, features = geneMarkers) + RotatedAxis() + coord_flip()
ggsave(plot = dotMarkerGG, filename = "Top5Marker_DotPlot.png", 
       dpi = res, width = 9, height = 16)

###annotation based on top 20 marker genes and feature plots###
embo_sce.big$major_cell_type <- embo_sce.big$RNA_snn_res.0.3
cell_type <- c("0" = "T cell",
               "1" = "Cancer cell",
               "2" = "Monocyte/Macrophage",
               "3" = "Cancer cell",
               "4" = "Cancer cell",
               "5" = "NK",
               "6" = "Cancer cell",
               "7" = "Plasma cell",
               "8" = "Cancer cell",
               "9" = "B cell",
               "10" = "Fibroblast",
               "11" = "NK",
               "12" = "Cancer cell",
               "13" = "DC",
               "14" = "Endothelium",
               "15" = "SMC")
embo_sce.big[['major_cell_type']] = unname(cell_type[embo_sce.big$major_cell_type])

proObj <- embo_sce.big
sigMarkers <- c("EPCAM", "KRT19", "COL1A1", "DCN", "VWF", "PECAM1", "CD3D", "CD8A", "CD4", "MS4A1", "CD27", "JCHAIN", "CD14", "CD68", "HLA-DRA", "CD1C", "KIT", "ITGAM")
sigExpr <- proObj[['RNA']]@data[sigMarkers,] 
sigExpr <- as.data.frame(sigExpr)
sigExpr$gene <- rownames(sigExpr)
gathExpr <- gather(sigExpr, "cell_id", "expr", rownames(proObj@meta.data))
gathExpr <- merge(gathExpr, proObj@meta.data, by.x = "cell_id", by.y = "row.names", all.x = TRUE)
gathExpr$gene <- factor(gathExpr$gene, levels = sigMarkers)

#####FigS7A#####
FigS7A <- DimPlot(proObj, reduction = "umap", label = TRUE, group.by = "major_cell_type") + 
  labs(color = "General cell type\nin TME", title = "Annotated and merged cell types") +
  scale_color_brewer(palette = "Spectral")
ggsave(plot = FigS7A, filename = paste(intePf, "FigS7A_UMAP.png", sep = ""), dpi = 600, width = 9, heigh = 6)

#####FigS7B#####
mcCts <- proObj@meta.data %>%
  group_by(patient, major_cell_type) %>%
  summarize(n = n())
#print(head(mcCts))
FigS7B <- ggplot(mcCts, aes(x = patient, y = n, fill = major_cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "gray") +
  scale_fill_brewer(palette = "Spectral") +
  coord_flip() +
  labs(x = "Sample ID", y = "Cell proportion", fill = "General cell type\nin TME") +
  theme_classic()
ggsave(plot = FigS7B, filename = paste(intePf, "FigS7B_BAR.png", sep = ""), dpi = 300, width = 9, heigh = 6)

#####FigS7C#####
cat("Heatmaps\n")
cellAnnDf <- as.data.frame(read.csv("ALL_POS_MARKERS.csv"))
rownames(cellAnnDf) <- cellAnnDf[,1]
cellAnnDf[,1] <- NULL
cellMarkerDf <- cellAnnDf
topMarkers <- cellMarkerDf %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
geneExpr <- AverageExpression(proObj,features = topMarkers$gene,group.by = "major_cell_type") 
geneExpr <- as.data.frame(geneExpr$RNA)
cellType <- data.frame(colnames(geneExpr))
colnames(cellType) <- 'class'
cols <- brewer.pal(length(unique(topMarkers$cell_type)), "Spectral")
top_anno = HeatmapAnnotation(df = cellType,#cell_type name/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c("B cell" = "#9e3e4f",
                                                  "Cancer cell" = "#d53e4f",
                                                  "DC" = "#f46d43",
                                                  "Endothelium" = "#fdae61",
                                                  "Fibroblast" = "#fee08b",
                                                  "Monocyte/Macrophage" = "#e6f598",
                                                  "NK" = "#abdda4",
                                                  "Plasma cell" = "#66c2a5",
                                                  "SMC" = "#3288bd",
                                                  "T cell" = "#5e4fa2")))

marker_exp <- t(scale(t(geneExpr),scale = T,center = T))
png(filename = "FigS7C_top5_marker_avg_heatmap.png",res = 600, width = 6, height = 8,units = "in")
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        #column_title = top_anno,
        heatmap_legend_param = list(title = ""),
        col = colorRampPalette(c("#2fa1dd", "white", "#f87669"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)
dev.off()

#####FigS7D#####
FigS7D <- ggplot(gathExpr, aes(x = major_cell_type, y = expr, fill = major_cell_type)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(.~gene, nrow = 1, scales = "free_x") +
  scale_fill_brewer(palette = "Spectral") +
  labs(x = "General cell type in TME", y = "Normalized expression", fill = "General cell type\nin TME") +
  coord_flip() +
  theme_classic()
ggsave(plot = figs1j, filename = "FigS7D_mc_violin_plot_for_cell_annotation.png", 
       dpi = 300, width = 16, height = 3) # FigS1x
saveRDS(figs1j, "FigS7D_mc_violin_plot_for_cell_annotation.RDS")

####*tamObj*####
tamObj <- proObj[ ,proObj$majority_voting_cell_type %in% "Monocytes/Macrophages"]
###Normalization###
tamObj <- NormalizeData(tamObj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
tamObj <- FindVariableFeatures(tamObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(tamObj)
tamObj <- ScaleData(tamObj, features = all.genes)
###Reduction###
tamObj <- RunPCA(tamObj, features = VariableFeatures(object = tamObj))
set.seed(1)
tamObj <- JackStraw(tamObj, num.replicate = 10)
tamObj <- ScoreJackStraw(tamObj, dims = 1:20)
tamObj<- tamObj %>% 
  RunUMAP(tamObj, dims = 1:30) %>% 
  RunTSNE(tamObj, dims = 1:30) %>% 
  FindNeighbors(tamObj, dims = 1:30) %>% 
  FindClusters(resolution = 0.3) %>% 
  identity()
tamPf <- NULL

#####FigS8A#####
FigS8A <- DimPlot(tamObj, reduction = "umap", label = F) + labs(color = "Clusters (r = 0.8)", title = "")
ggsave(plot = FigS8A, filename = paste("FigS8A_UMAP_all_tam_cluster.png", sep = ""), dpi = 300, width = 9, heigh = 6)

#####FigS8B#####
FigS8B <- FeaturePlot(tamObj, features = c("rna_CD274", "rna_SIGLEC15"), blend = T)
ggsave(plot = FigS8B, filename = paste(tamPf, "FigS8B_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 18, heigh = 4)

#####FigS8C#####
FigS8C <- VlnPlot(tamObj, features = "rna_CD274", group.by = "RNA_snn_res.0.8")
ggsave(plot = FigS8C, filename = paste(tamPf, "FigS8C_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 18, heigh = 4)

#####FigS8D#####
geneExpr <- AverageExpression(tamObj,features = topMarkers$gene,group.by = "RNA_snn_res.0.8") 
geneExpr <- as.data.frame(geneExpr$RNA)
cellType <- data.frame(colnames(geneExpr))
colnames(cellType) <- 'class'
cols <- brewer.pal(length(unique(topMarkers$cluster)), "Spectral")
top_anno = HeatmapAnnotation(df = cellType,#cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             annotation_legend_param = list(at = c(0:14)),
                             col = list(class = c("0" = "#F8766D",
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
colnames(marker_exp) = factor(colnames(marker_exp), levels = c(0:14))

png(filename = "FigS8D_top5_tam_marker_avg_heatmap.png",res = 600, width = 8, height = 10,units = "in")
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

#####FigS8E#####
ress <- c(seq(0.1,0.9, 0.1))
tamObj <- FindClusters(tamObj, resolution = ress)
ress = 0.8
all_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(ress)+1)
FigS8E <- clustree(tamObj, prefix = "RNA_snn_res.", node_text_size = 4) + scale_color_manual(values = all_cols)
ggsave("FigS8E_clustree.png", FigS8E, dpi = 600, width = 8, height = 8)###res.0.8

#####Fig2A#####
Fig2A <- DimPlot(tamObj, reduction = "umap", label = F) + labs(color = "Clusters (r = 0.8)", title = "")
ggsave(plot = Fig2A, filename = paste("Fig2A_UMAP_all_tam_cluster.png", sep = ""), dpi = 300, width = 9, heigh = 6)

#####Fig2B#####
Fig2B <- FeaturePlot(tamObj, features = c("rna_CD274", "rna_SIGLEC15"), blend = T)
ggsave(plot = Fig2B, filename = paste(tamPf, "Fig2B_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 18, heigh = 4)

#####Fig2C#####
tamObj@meta.data$tam_type = ifelse(tamObj$RNA_snn_res.0.8%in%c("2","3","4","6",'7',"14"), "PD-L1- TAM","PD-L1+ TAM")%>%
                            factor(tamObj@meta.data$tam_type, levels = c("PD-L1- TAM","PD-L1+ TAM"))
Fig2C <- DimPlot(tamObj, group.by = "tam_type", reduction = "umap", cols = c("#1E90FF","#CD5C5C"),pt.size = 1)
ggsave(plot = Fig2C, filename = paste(tamPf, "Fig2C_UMAP_blend_pdl1.png", sep = ""), dpi = 300, width = 9, heigh = 6)

#####Fig2D#####
Fig2D <- VlnPlot(tamObj, features = "rna_CD274", group.by = "tam_type", y.max = 10,cols = c("#1E90FF","#CD5C5C"))+
  stat_compare_means(comparisons = list(c("PD-L1- TAM","PD-L1+ TAM")), label = "p.signif")
ggsave(Fig2D, "Fig2D_rna_cd274_vlnplot.png", width = 4, height = 6, dpi = 600, units = "in")

#####Fig2E#####
de_res <- FindMarkers(tamObj, ident.1 = "PD-L1+ TAM", ident.2 = "PD-L1- TAM", group.by = "tam_type", logfc.threshold = 0.5, test.use = "wilcox")
print(head(de_res))
write.csv(de_res, paste(tamPf, "pdl1_tam_de.csv", sep = "")) 
plot_de <- de_res[de_res$p_val_adj < 0.05 & abs(de_res$avg_log2FC)>= 1,]
print(plot_de)
labFeatures <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", 
                 "C1QA", "C1QB", "C1QC", "IL1B", "CXCL2", "CXCL3", "CXCL8", "CCL3", "CCL4", "CCL18", 
                 "CD74", "CD83", "FOS", "FOSB", "JUN", "JUNB", "CEBPD", "FOLR2", "MS4A6A", 
                 "SPP1", "FABP4", "FABP5", "CSTB", "SPARC", "IL1RN", "CD9", "CD52", "CHI3L1", "LGALS3", 
                 "PKM", "LPL", "CHIT1", "SFRP2", "LDHA", "S100A10", "DCN", "GSN", "MMP9", "COL1A1", "COL1A2", "COL3A1", "FN1", "LUM")

deRes <- de_res
keyvals <- ifelse(deRes$avg_log2FC < -0.5, 'dodgerblue', 
                  ifelse(deRes$avg_log2FC > 0.5, 'firebrick', 'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'firebrick'] <- 'Upregulated\nin PD-L1+ TAM'
names(keyvals)[keyvals == 'grey'] <- 'Insiginificant genes'
names(keyvals)[keyvals == 'dodgerblue'] <- 'Upregulated\nin PD-L1- TAM'

Fig2E <- EnhancedVolcano(deRes,
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
                      selectLab = labFeatures,
                      ylab = bquote(~ '-' ~Log[10]~ 'Adjusted P-value'),
                      drawConnectors = TRUE
)
png(paste(tamPf, "Fig2E_post_MAST_pos_vs_neg_ev.png", sep = ""), 
    res = 300, height = 7.2, width = 9, units = 'in')
print(Fig2E)
gar <- dev.off()

#####Fig2F#####
Maturation_sig <- c("rna_CD74", "rna_HLA-DRA", "rna_HLA-DQA1", "rna_C1QA")
Pro_inflammatory_sig <- c("rna_C1QA", "rna_CXCL9","rna_CXCL10")
Anti_inflammatory_sig <- c("rna_MARCO", "rna_CD52", "rna_IL1RN", "rna_CSTB")
Pro_tumor_sig <- c("rna_SPP1", "rna_VEGFA","rna_CXCL8", "rna_S100A8")
for (sig_gene in Maturation_sig) #change sig_genelist here
  { 
  cat(sig_gene, "\n")
  pdl1VlnGG <- VlnPlot(tamObj, features = sig_gene, group.by = "tam_type", log = FALSE, pt.size = 0,
                       cols = c('dodgerblue', 'firebrick')) + scale_y_continuous(limits = c(0.0, NA))
  ggsave(plot = pdl1VlnGG, filename = paste(tamPf, sig_gene, "_Fig2F_manual_selected_VlnPlot.png", sep = ""),
         dpi = 300, width = 3, height = 4)
}

####finished####
