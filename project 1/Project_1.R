library(Seurat)
library(dplyr)
library(SeuratObject)
library(ggplot2)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)
library(tidyr)
library(celldex)
setwd("~/Third Semester/Single Cell/")


BMMC_T1 <- readRDS("GSM4138872_scRNA_BMMC_D1T1.rds")
BMMC_T2 <- readRDS("GSM4138873_scRNA_BMMC_D1T2.rds")
CD34_D2 <- readRDS("GSM4138874_scRNA_CD34_D2T1.rds")
CD34_D3 <- readRDS("GSM4138875_scRNA_CD34_D3T1.rds")


BMMC_T1_obj <- CreateSeuratObject(counts = BMMC_T1, project = "BMMC_T1", min.cells = 3, min.features = 200)
BMMC_T2_obj <- CreateSeuratObject(counts = BMMC_T2, project = "BMMC_T2", min.cells = 3, min.features = 200)
CD34_D2_obj <- CreateSeuratObject(counts = CD34_D2, project = "CD34_D2", min.cells = 3, min.features = 200)
CD34_D3_obj <- CreateSeuratObject(counts = CD34_D3, project = "CD34_D3", min.cells = 3, min.features = 200)

# BMMC_T1_obj
# BMMC_T2_obj
# CD34_D2_obj
# CD34_D3_obj

BMMC_T1_obj$Donor <- "D1"
BMMC_T1_obj$Replicate <- "T1"
BMMC_T1_obj$Sex <- "F"

BMMC_T2_obj$Donor <- "D1"
BMMC_T2_obj$Replicate <- "T2"
BMMC_T2_obj$Sex <- "F"

CD34_D2_obj$Donor <- "D2"
CD34_D2_obj$Replicate <- "T1"
CD34_D2_obj$Sex <- "M"

CD34_D3_obj$Donor <- "D3"
CD34_D3_obj$Replicate <- "T1"
CD34_D3_obj$Sex <- "F"

#head(BMMC_T1_obj, n = 5)
#BMMC_T1_obj

BMMC_T1_obj[["percent.mt"]] <- PercentageFeatureSet(BMMC_T1_obj, pattern = "^MT-")
VlnPlot(BMMC_T1_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot_1 <- FeatureScatter(BMMC_T1_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot_2 <- FeatureScatter(BMMC_T1_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_1 + plot_2
BMMC_T1_obj <- subset(BMMC_T1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
BMMC_T1_obj

BMMC_T1_obj <- NormalizeData(BMMC_T1_obj) 
BMMC_T1_obj <- FindVariableFeatures(BMMC_T1_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(BMMC_T1_obj), 10)
plot1 <- VariableFeaturePlot(BMMC_T1_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2

BMMC_T1_obj_fil <- ScaleData(BMMC_T1_obj)
BMMC_T1_obj_fil <- RunPCA(BMMC_T1_obj_fil)
#ElbowPlot(BMMC_T1_obj_fil)
BMMC_T1_obj_fil <- FindNeighbors(BMMC_T1_obj_fil, dim = 1:20)
BMMC_T1_obj_fil <- FindClusters(BMMC_T1_obj_fil)
BMMC_T1_obj_fil <- RunUMAP(BMMC_T1_obj_fil, dim = 1:20)


sweep.res <- paramSweep_v3(BMMC_T1_obj_fil, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn_BMMCT1 <- find.pK(sweep.stats)
#ggplot(bcmvn_BMMCT1, aes(pK, BCmetric, group = 1))+
#  geom_point() + 
 # geom_line()

pK <- bcmvn_BMMCT1 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK <- as.numeric(as.character(pK[[1]]))

colnames(BMMC_T1_obj_fil@meta.data)

annotations <- BMMC_T1_obj_fil@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.75*nrow(BMMC_T1_obj_fil@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

BMMC_T1_obj_fil <- doubletFinder_v3(BMMC_T1_obj_fil, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
colnames(BMMC_T1_obj_fil@meta.data)
#DimPlot(BMMC_T1_obj_fil, reduction = "umap", group.by = "DF.classifications_0.25_0.005_264")


###BMMC_T2----------------------------------------------------------------------

BMMC_T2_obj[["percent.mt"]] <- PercentageFeatureSet(BMMC_T2_obj, pattern = "^MT-")
VlnPlot(BMMC_T2_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BMMC_T2_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BMMC_T2_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
BMMC_T2_obj <- subset(BMMC_T2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
#BMMC_T2_obj

BMMC_T2_obj <- NormalizeData(BMMC_T2_obj)
BMMC_T2_obj <- FindVariableFeatures(BMMC_T2_obj, selection.method = "vst", nfeatures = 2000)
BMMCT2_top10 <- head(VariableFeatures(BMMC_T2_obj), 10)
FT_plot1 <- VariableFeaturePlot(BMMC_T2_obj)
FT_plot2 <- LabelPoints(plot = FT_plot1, points = BMMCT2_top10, repel = TRUE)
FT_plot1 + FT_plot2

BMMC_T2_obj_fil <- ScaleData(BMMC_T2_obj)
BMMC_T2_obj_fil <- RunPCA(BMMC_T2_obj_fil)
ElbowPlot(BMMC_T2_obj_fil)
BMMC_T2_obj_fil <- FindNeighbors(BMMC_T2_obj_fil, dim = 1:20)
BMMC_T2_obj_fil <- FindClusters(BMMC_T2_obj_fil)
BMMC_T2_obj_fil <- RunUMAP(BMMC_T2_obj_fil, dim = 1:20)


sweep.res_2 <- paramSweep_v3(BMMC_T2_obj_fil, PCs = 1:20, sct = FALSE)
sweep.stats_2 <- summarizeSweep(sweep.res_2, GT = FALSE)
bcmvn_BMMCT2 <- find.pK(sweep.stats_2)

ggplot(bcmvn_BMMCT2, aes(pK, BCmetric, group = 1))+
  geom_point() + 
  geom_line()

pK_2 <- bcmvn_BMMCT2 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK_2 <- as.numeric(as.character(pK_2[[1]]))

colnames(BMMC_T2_obj_fil@meta.data)

annotations <- BMMC_T2_obj_fil@meta.data$seurat_clusters
homotypic.prop_2 <- modelHomotypic(annotations)
nExp_poi_2 <- round(0.75*nrow(BMMC_T2_obj_fil@meta.data))
nExp_poi.adj_2 <- round(nExp_poi_2*(1-homotypic.prop_2))

BMMC_T2_obj_fil <- doubletFinder_v3(BMMC_T2_obj_fil, PCs = 1:20, pN = 0.25, pK = pK_2, nExp = nExp_poi_2, reuse.pANN = FALSE, sct = FALSE)
colnames(BMMC_T2_obj_fil@meta.data)
#DimPlot(BMMC_T2_obj_fil, reduction = "umap", group.by = "DF.classifications_0.25_0.005_475")


###CD34_D2----------------------------------------------------------------------

CD34_D2_obj[["percent.mt"]] <- PercentageFeatureSet(CD34_D2_obj, pattern = "^MT-")
VlnPlot(CD34_D2_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pl1 <- FeatureScatter(CD34_D2_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
pl2 <- FeatureScatter(CD34_D2_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pl1 + pl2
CD34_D2_obj <- subset(CD34_D2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)


CD34_D2_obj <- NormalizeData(CD34_D2_obj)
CD34_D2_obj <- FindVariableFeatures(CD34_D2_obj, selection.method = "vst", nfeatures = 2000)
CD34_d2_top10 <- head(VariableFeatures(CD34_D2_obj), 10)
ft_plot1 <- VariableFeaturePlot(CD34_D2_obj)
ft_plot2 <- LabelPoints(plot = ft_plot1, points = CD34_d2_top10, repel = TRUE)
ft_plot1 + ft_plot2

CD34_D2_obj_fil <- ScaleData(CD34_D2_obj)
CD34_D2_obj_fil <- RunPCA(CD34_D2_obj_fil)
ElbowPlot(CD34_D2_obj_fil)
CD34_D2_obj_fil <- FindNeighbors(CD34_D2_obj_fil, dim = 1:20)
CD34_D2_obj_fil <- FindClusters(CD34_D2_obj_fil)
CD34_D2_obj_fil <- RunUMAP(CD34_D2_obj_fil, dim = 1:20)


sweep.res_cd <- paramSweep_v3(CD34_D2_obj_fil, PCs = 1:20, sct = FALSE)
sweep.stats_cd <- summarizeSweep(sweep.res_cd, GT = FALSE)
bcmvn_CD34_d2 <- find.pK(sweep.stats_cd)

ggplot(bcmvn_CD34_d2, aes(pK, BCmetric, group = 1))+
  geom_point() + 
  geom_line()

pK_cd <- bcmvn_CD34_d2 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK_cd <- as.numeric(as.character(pK_cd[[1]]))

colnames(CD34_D2_obj_fil@meta.data)

annotations_cd <- CD34_D2_obj_fil@meta.data$seurat_clusters
homotypic.prop_cd <- modelHomotypic(annotations_cd)
nExp_poi_cd <- round(0.75*nrow(CD34_D2_obj_fil@meta.data))
nExp_poi.adj_cd <- round(nExp_poi_cd*(1-homotypic.prop_cd))

CD34_D2_obj_fil <- doubletFinder_v3(CD34_D2_obj_fil, PCs = 1:20, pN = 0.25, pK = pK_cd, nExp = nExp_poi_cd, reuse.pANN = FALSE, sct = FALSE)
colnames(CD34_D2_obj_fil@meta.data)
#DimPlot(CD34_D2_obj_fil, reduction = "umap", group.by = "DF.classifications_0.25_0.3_116")


#CD34_D3------------------------------------------------------------------------

CD34_D3_obj[["percent.mt"]] <- PercentageFeatureSet(CD34_D3_obj, pattern = "^MT-")
VlnPlot(CD34_D3_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plt1 <- FeatureScatter(CD34_D3_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plt2 <- FeatureScatter(CD34_D3_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plt1 + plt2
CD34_D3_obj <- subset(CD34_D3_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)


CD34_D3_obj <- NormalizeData(CD34_D3_obj)

CD34_D3_obj <- FindVariableFeatures(CD34_D3_obj, selection.method = "vst", nfeatures = 2000)
CD34_d3_top10 <- head(VariableFeatures(CD34_D3_obj), 10)
f_plot1 <- VariableFeaturePlot(CD34_D3_obj)
f_plot2 <- LabelPoints(plot = f_plot1, points = CD34_d3_top10, repel = TRUE)
f_plot1 + f_plot2

CD34_D3_obj_fil <- ScaleData(CD34_D3_obj)
CD34_D3_obj_fil <- RunPCA(CD34_D3_obj_fil)
ElbowPlot(CD34_D3_obj_fil)
CD34_D3_obj_fil <- FindNeighbors(CD34_D3_obj_fil, dim = 1:20)
CD34_D3_obj_fil <- FindClusters(CD34_D3_obj_fil)
CD34_D3_obj_fil <- RunUMAP(CD34_D3_obj_fil, dim = 1:20)


sweep.res_cd2 <- paramSweep_v3(CD34_D3_obj_fil, PCs = 1:20, sct = FALSE)
sweep.stats_cd2 <- summarizeSweep(sweep.res_cd2, GT = FALSE)
bcmvn_CD34_d3 <- find.pK(sweep.stats_cd2)

ggplot(bcmvn_CD34_d3, aes(pK, BCmetric, group = 1))+
  geom_point() + 
  geom_line()

pK_cd2 <- bcmvn_CD34_d3 %>%
  filter(BCmetric == max(BCmetric)) %>%
  select(pK)
pK_cd2 <- as.numeric(as.character(pK_cd2[[1]]))

colnames(CD34_D3_obj_fil@meta.data)

annotations_cd2 <- CD34_D3_obj_fil@meta.data$seurat_clusters
homotypic.prop_cd2 <- modelHomotypic(annotations_cd2)
nExp_poi_cd2 <- round(0.75*nrow(CD34_D3_obj_fil@meta.data))
nExp_poi.adj_cd2 <- round(nExp_poi_cd2*(1-homotypic.prop_cd2))

CD34_D3_obj_fil <- doubletFinder_v3(CD34_D3_obj_fil, PCs = 1:20, pN = 0.25, pK = pK_cd2, nExp = nExp_poi_cd2, reuse.pANN = FALSE, sct = FALSE)
colnames(CD34_D3_obj_fil@meta.data)
#DimPlot(CD34_D3_obj_fil, reduction = "umap", group.by = "DF.classifications_0.25_0.005_403")



#Merging the Seurat Objects and Finding Batch Effects in the Data---------------

allData <- merge(x = BMMC_T1_obj_fil, y = list(BMMC_T2_obj_fil, CD34_D2_obj_fil, CD34_D3_obj_fil))
#allData$Sample <- rownames(allData@meta.data)
#allData@meta.data <- separate(allData@meta.data, col = "Sample", into = c("Sample"), sep = ":")

saveRDS(allData, file = "merged_data.rds")
allData <- readRDS("merged_data.rds")

allData_fil <- NormalizeData(allData)
allData_fil <- FindVariableFeatures(allData_fil)
allData_fil <- ScaleData(allData_fil, verbose = FALSE)
allData_fil <- RunPCA(allData_fil, npcs = 30, verbose = FALSE)
allData_fil <- RunUMAP(allData_fil, reduction = "pca", dims = 1:20)
ElbowPlot(allData_fil)
allData_fil <- FindNeighbors(allData_fil, reduction="pca", dims = 1:20)
allData_fil <- FindClusters(allData_fil, resolution = 0.26)

DimPlot(allData_fil, reduction = "umap", group.by = "orig.ident")
#DimPlot(allData_fil, reduction = "umap", group.by = "Sample")
DimPlot(allData_fil, reduction = "umap", group.by = "Sex")
DimPlot(allData_fil, reduction = "umap", group.by = "Replicate")
DimPlot(allData_fil, reduction = "umap", group.by = "Donor")


obj.list <- SplitObject(allData_fil, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(obj.list)
anchors <- FindIntegrationAnchors(obj.list, anchor.features = features)

saveRDS(anchors, file = "anchors.rds")
anchors <- readRDS("anchors.rds")

allData_fil <- IntegrateData(anchorset = anchors)
#allData_fil <- FindVariableFeatures(allData_fil)
allData_fil <- ScaleData(allData_fil)
allData_fil <- RunPCA(allData_fil, dims =1:20)
#ElbowPlot(allData_fil)
allData_fil <- FindNeighbors(allData_fil, dims = 1:20)
allData_fil <- FindClusters(allData_fil, resolution = 0.26)
allData_fil <- RunUMAP(allData_fil, dims = 1:20)
DimPlot(allData_fil, reduction = "umap", group.by = "orig.ident")
#DimPlot(allData_fil, reduction = "umap", group.by = "Sample")
DimPlot(allData_fil, reduction = "umap", group.by = "Sex")
DimPlot(allData_fil, reduction = "umap", group.by = "Replicate")
DimPlot(allData_fil, reduction = "umap", group.by = "Donor")

#Dimensionality Reduction with PCA and UMAP-------------------------------------

allData_dim <- RunPCA(allData_fil, features = VariableFeatures(allData_fil))
allData_dim <- JackStraw(allData_dim, num.replicate = 100)
allData_dim <- ScoreJackStraw(allData_dim, dims = 1:20)
ElbowPlot(allData_dim)
JackStrawPlot(allData_dim, dims = 1:15)
allData_dim <- RunUMAP(allData_dim, dims = 1:20) #Non-linear Reduction with UMAP
DimPlot(allData_dim, reduction = "umap")


saveRDS(allData_dim, file = "allData_dime.rds")


#CLustering---------------------------------------------------------------------
allData_dim <- readRDS("allData_dime.rds")
allData_dim <- FindNeighbors(allData_dim, dims = 1:10)
allData_dim <- FindClusters(allData_dim, resolution = 0.26)
allData_dim <- RunUMAP(allData_dim, dims = 1:20)

DimPlot(allData_dim, reduction = "umap")

saveRDS(allData_dim, file = "allData_clus.rds")


#Cell Annotation __ Automatic Annotation----------------------------------------
allData_dim <- readRDS("allData_clus.rds")
ref <- celldex :: HumanPrimaryCellAtlasData()
allData_annot <- as.SingleCellExperiment(DietSeurat(allData_dim))
annot <- SingleR(test = allData_annot, ref = ref, assay.type.test=1, labels = ref$label.main)
allData_dim@meta.data$Auto_Annotate <- annot$pruned.labels
DimPlot(allData_dim, reduction = "umap", group.by = "Auto_Annotate" ,label = T, repel = T) + NoLegend()
saveRDS(allData_dim, file = "allData_annotation_auto.rds")


#Manual Cell Annotation---------------------------------------------------------

allData_dim <- readRDS("allData_clus.rds")
#cluster2.markers <- FindMarkers(allData_dim, ident.1 = 2, min.pct = 0.25)
#head(cluster2.markers, n = 5)
#cluster5.markers <- FindMarkers(allData_dim, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
#head(cluster5.markers, n = 5)
DimPlot(allData_dim, reduction = "umap", group.by = "seurat_clusters", label = T)

markers <- FindAllMarkers(allData_dim, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

feature_HSC <- c("CD34", "CD38", "Sca1", "Kit") #6
feature_LMPP <- c("CD38", "CD52", "CSF3R", "ca1", "Kit", "CD34", "Flk2")#1
feature_CLP <- c("IL7R")#2
feature_GMPN <- c("ELANE")#6
feature_CMP <- c("IL3", "GM-CSF", "M-CSF")
feature_B <- c("CD19", "CD20", "CD38")#8
feature_PreB <- c("CD19", "CD34")#10
feature_Plasma <- c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN")#11
feature_Cd8 <- c("CD3D", "CD3E", "CD8A", "CD8B")#0,5
feature_Cd4 <- c("CD3D", "CD3E", "CD4")#0
feature_NK <- c("FCGR3A", "NCAM1", "NKG7", "KLRB1")#4
feature_Erythrocytes <- c("GATA1", "HBB", "HBA1", "HBA2")#5
feature_pDC <- c("IRF8", "IRF4", "IRF7")#9
feature_cDC <- c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA")#1
feature_CD14 <- c("CD14", "CCL3", "CCL4", "IL1B")#1
feature_CD16 <- c("FCGR3A", "CD68", "S100A12")#1
feature_Basophils <- c("GATA2")#2

#FeaturePlot(allData_dim, features = c("CD14"), min.cutoff = "q10")

#new.cluster.ids <- markers %>% group_by(cluster) %>% top_n(1, avg_log2FC) %>% pull(gene)
#names(new.cluster.ids) <- levels(allData_dim)
#allData_dim <- RenameIdents(allData_dim, new.cluster.ids)

DimPlot(allData_dim, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


FeaturePlot(allData_dim, features = feature_pDC, label = T)

cluster.ids <- c("CD4/CD8", "Myleoid-Cells", "Basophils", "LMPP", "Erthrocytes", "NK",
                 "HSC", "B", "pDC", "PreB", "Plasma")

names(cluster.ids) <- levels(allData_dim)
allData_dim <- RenameIdents(allData_dim, cluster.ids)
DimPlot(allData_dim, reduction = "umap", label = TRUE)

VlnPlot(allData_dim, features = c("CD1C", "CD19", "ELANE"))
#DimPlot()


#Differential Analysis----------------------------------------------------------
data <- readRDS("allData_clus.rds")
De.markers <- FindMarkers(data, ident.1 = 0 , ident.2 = 7, min.pct = 0.25) # T-cells vs B-cells
head(De.markers, 5)
De.markers$expression <- "No"
De.markers$expression[De.markers$avg_log2FC > 0.25 & De.markers$p_val_adj  < 0.01] <- "UP" 
De.markers$expression[De.markers$avg_log2FC < 0.25 & De.markers$p_val_adj < 0.01] <- "Down"

De.markers$labels[De.markers$p_val<= thresh] <- (De.markers$gene[De.markers$p_val <= thresh])
head(arrange(De.markers,p_val),10)

ggplot(De.markers, aes(x=avg_log2FC, y=-log10(p_val), col=expression)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c('blue', 'red', 'black')) +
  theme(text=element_text(size=20))



De.markers_cd <- FindMarkers(data, ident.1 = 0 , ident.2 = 1, min.pct = 0.25) # Cd4 vs Cd14
De.markers_cd$expression <- "No"
De.markers_cd$expression[De.markers_cd$avg_log2FC > 0.25 & De.markers_cd$p_val_adj  < 0.05] <- "UP" 
De.markers_cd$expression[De.markers_cd$avg_log2FC < 0.25 & De.markers_cd$p_val_adj < 0.05] <- "Down"
ggplot(De.markers_cd, aes(x=avg_log2FC, y=-log10(p_val), col=expression)) + 
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c('blue', 'black', 'red')) +
  theme(text=element_text(size=20))

head(FindMarkers(allData_dim, ident.1 = "CD4/CD8", ident.2 = "B", logfc.threshold = log(2)))

