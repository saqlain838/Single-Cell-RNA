setwd("C:/Users/Public/Single-cell")
library(ArchR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
addArchRGenome("hg19")
set.seed(1)
InputFiles <- list.files(pattern = ".tsv")

names(InputFiles) <- c("scATAC_BMMC_D5T1", "scATAC_CD34_D7T1", "scATAC_CD34_D8T1")

Arrow <- createArrowFiles(
    inputFiles = InputFiles,
    sampleNames = names(InputFiles),
    filterTSS = 5,
    filterFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
    )

doubScores <- addDoubletScores(
    input = Arrow,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
)

projR <- ArchRProject(
  ArrowFiles = Arrow, 
  outputDirectory = "Project2",
  copyArrows = TRUE
)

saveArchRProject(ArchRProj = projR, outputDirectory = "ArchR_project", load = FALSE)


proj1 <-  loadArchRProject(path = "ArchR_project", force = FALSE, showLogo = TRUE)
colnames(proj1@cellColData)

# metadata <- read.csv("metadata.csv")
# names(metadata)
# for(meta in names(metadata)){
#     proj1@cellColData[[meta]] <- metadata[match(as.character(proj1@cellColData$Sample), metadata$Sample), meta]
# }
# #head(proj1@cellColData)
#proj1
#proj1[proj1$cellNames[1:100], ]


##### Adding Metadata according to the Sample in the ArchR Project ######
samples <- gsub("scATAC_","",proj1$Sample)
proj1$samples <- samples
Type <- gsub("scATAC_|_D5T1|_D8T1|_D7T1_","",proj1$Sample)
proj1$Type <- Type
Donor <-  gsub('BMMC_D5T1', 'D5',
               gsub('CD34_D8T1', 'D7',
                    gsub('CD34_D7T1', 'D6', proj1$samples)))
proj1$Donor <- Donor
Sex <- gsub('BMMC_D5T1', 'F',
            gsub('CD34_D8T1', 'M',
                 gsub('CD34_D7T1', 'M', proj1$samples)))          
proj1$Sex <- Sex

##### Cell Counts in Each Sample #####

idxSample <- BiocGenerics::which(proj1$Sample %in% "scATAC_BMMC_D5T1") ### Cell count in BMMC_D5T1 Sample
cellsSample <- proj1$cellNames[idxSample]
proj1[cellsSample, ]

idxSample <- BiocGenerics::which(proj1$Sample %in% "scATAC_CD34_D7T1") ### Cell count in CD34_D7T1 Sample
cellsSample <- proj1$cellNames[idxSample]
proj1[cellsSample, ]

idxSample <- BiocGenerics::which(proj1$Sample %in% "scATAC_CD34_D8T1") ### Cell count in CD34_D7T1 Sample
cellsSample <- proj1$cellNames[idxSample]
proj1[cellsSample, ]

##### Plot for TSS Enrichment distribution of each Sample #####

plot_TSS <- plotGroups(
    ArchRProj = proj1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

plotPDF(plot_TSS, name = "TSSEnrichment-plot-ridges-2.pdf", ArchRProj = proj1, addDOC = FALSE, width = 4, height = 4)


##### Plot for Fragment Size Distribution

plot_frag <- plotFragmentSizes(ArchRProj = proj1)
#plot_frag
plotPDF(plot_frag, name = "FragmentDistribution-plot.pdf", ArchRProj = proj1, addDOC = FALSE, width = 4, height = 4)


## Task 1.6: Filtering
dc <- as.matrix(proj1)
proj2 <- filterDoublets(proj1)
proj2 <- proj2[which(proj2$TSSEnrichment > 7 & proj2$nFrags > 2000 & proj2$nFrags < 30000)]

idxSample <- BiocGenerics::which(proj2$Sample %in% "scATAC_BMMC_D5T1") ### Cell count in BMMC_D5T1 Sample
cellsSample <- proj2$cellNames[idxSample]
proj2[cellsSample, ]

idxSample <- BiocGenerics::which(proj2$Sample %in% "scATAC_CD34_D7T1") ### Cell count in CD34_D7T1 Sample
cellsSample <- proj2$cellNames[idxSample]
proj2[cellsSample, ]

idxSample <- BiocGenerics::which(proj2$Sample %in% "scATAC_CD34_D8T1") ### Cell count in CD34_D7T1 Sample
cellsSample <- proj2$cellNames[idxSample]
proj2[cellsSample, ]

#############################################  WEEK 02 #####################################################

##### Dimensionality Reduction #####

proj2 <- addIterativeLSI(ArchRProj = proj2, 
    useMatrix = "TileMatrix", 
    name = "IterativeLSI",
    iterations = 2, 
    clusterParams = list(
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30)


proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

plot_dimensionality <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = c("Sample", "TSSEnrichment", "nFrags"), embedding = "UMAP")
#plot_dimensionality
plotPDF(plot_dimensionality, name = "Dimensionality.pdf", ArchRProj = proj2, addDOC = FALSE, width = 4, height = 4)


###### Batch Effect using Harmony #######

proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

plot_Batch <- plotEmbedding(ArchRProj = proj2, colorBy = "cellColData", name = c("Sample", "TSSEnrichment", "nFrags"), embedding = "UMAPHarmony")
#p3
plotPDF(plot_Batch, name = "BatchCorrection-Harmony.pdf", ArchRProj = proj2, addDOC = FALSE, width = 7, height = 6)



####### CLUSTERING ###########

proj3 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2
)


head(proj3$Clusters)
table(proj3$Clusters)
cM <- confusionMatrix(paste0(proj3$Clusters), paste0(proj3$Sample))
cM
cluster_plot <- plotEmbedding(ArchRProj = proj3, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
#cluster_plot
plotPDF(cluster_plot, name = "Clusters-plot.pdf", ArchRProj = proj3, addDOC = FALSE, width = 7, height = 7)

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
plotPDF(p, name = "Clusters-plot-heatmap.pdf", ArchRProj = proj3, addDOC = FALSE, width = 7, height = 7)

########### PEAK CALLING #############
proj3 <- addGroupCoverages(ArchRProj = proj3, groupBy = "Clusters", force = TRUE)


proj3 <- addReproduciblePeakSet(
    ArchRProj = proj3, 
    groupBy = "Clusters",
    peakMethod = "Tiles",
    method = "p", 
    threads = 1
)
getPeakSet(proj3)
proj4 <- addPeakMatrix(proj3)
getAvailableMatrices(proj4)


markersPeaks <- getMarkerFeatures(
    ArchRProj = proj4, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  threads = 1
)

markersPeaks

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList$C3

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Heatmap-PEAKCalling.pdf", ArchRProj = proj4, addDOC = FALSE, width = 20, height = 10)

gene_browser <- c("CD19", "CD34", "CD3D", "CD69")
plot_Browser <- plotBrowserTrack(
    ArchRProj = proj4, 
    groupBy = "Clusters", 
    geneSymbol = gene_browser,
    upstream = 50000,
    downstream = 50000
)
plotPDF(plot_Browser, name = "BrowserTrack.pdf", ArchRProj = proj4, addDOC = FALSE, width = 8, height = 8)

grid::grid.draw(p$CD19)


################ WEEK 03 ###########################
###### GENE ACTIVITY ########
markersGS <- getMarkerFeatures(
    ArchRProj = proj4, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    threads = 1
)

markerListGS <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
#markerListGS$C3

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  nLabel = 5,
  nPrint = 5,
  transpose = TRUE
)
plotPDF(heatmapGS, name = "Gene-Activity-Heatmap.pdf", ArchRProj = proj4, addDOC = FALSE, width = 15, height = 10)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

marker_CLust3 <- c("NPPA", "LCE3A", "GLMP", "MIR3620", "OR2T6")

plot_marker <- plotEmbedding(
    ArchRProj = proj4, 
    colorBy = "GeneScoreMatrix", 
    name = marker_CLust3, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

plotPDF(plot_marker, name = "Marker-Plot-WithoutMAGIC.pdf", ArchRProj = proj4, addDOC = FALSE, width = 4, height = 4)
saveArchRProject(ArchRProj = proj4, outputDirectory = "ArchRProject", load = FALSE)


proj4 <-  loadArchRProject(path = "ArchRProject", force = FALSE, showLogo = TRUE)
proj4 <- addImputeWeights(proj4)

plot_marker_MAGIC <- plotEmbedding(
    ArchRProj = proj4, 
    colorBy = "GeneScoreMatrix", 
    name = marker_CLust3, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(proj4)
)

plotPDF(plot_marker_MAGIC, name = "MAGIC-Plot.pdf", ArchRProj = proj4, addDOC = FALSE, width = 4, height = 4)


######### TF MOTIF ACTIVITY ###########
proj4 <- addMotifAnnotations(ArchRProj = proj4, motifSet = "cisbp", name = "Motif")

motifsUp <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj4,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

proj4 <- addBgdPeaks(proj4)
proj4 <- addDeviationsMatrix(
  ArchRProj = proj4, 
  peakAnnotation = "Motif",
  force = TRUE
)
getAvailableMatrices(proj4)
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df, 5)

motifs <- c("GATA2", "GATA1")
motifM <- getFeatures(proj4, select = paste(motifs, collapse = "|"), useMatrix = "MotifMatrix")
motifM <- grep("z:", motifM, value = TRUE)

plot_marker_TF <- plotEmbedding(
    ArchRProj = proj4, 
    colorBy = "MotifMatrix", 
    name = sort(motifM), 
    embedding = "UMAP"
)
plotPDF(plot_marker_TF, name = "TF-Motif-Activity.pdf", ArchRProj = proj4, addDOC = FALSE, width = 4, height = 4)
saveArchRProject(ArchRProj = proj4, outputDirectory = "ArchRProject", load = FALSE)


plot_TF <- plotGroups(
    ArchRProj = proj4, 
    groupBy = "Clusters", 
    colorBy = "MotifMatrix",
    name = motifM,
    imputeWeights = getImputeWeights(proj4)
   )

plotPDF(plot_TF, name = "TF-Motif-Activity-Distribution.pdf", ArchRProj = proj4, addDOC = FALSE, width = 6, height = 6)

########### INTEGRATION WITH scRNA  ###############
###TASK 8.1
proj4 <-  loadArchRProject(path = "ArchRProject", force = FALSE, showLogo = TRUE)
sRna <- readRDS("allData_annotation_auto.rds")

proj4 <- addGeneIntegrationMatrix(
    ArchRProj = proj4, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sRna,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "Auto_Annotate",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
getAvailableMatrices(proj4)
proj4 <- addImputeWeights(proj4)
saveArchRProject(ArchRProj = proj4, outputDirectory = "ArchRProject", load = FALSE)




markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "CD19", #B-Cell Trajectory
    "CD3D", #TCells
    "CD14", #Monocytes
    "TBX21", #NK
    "CD38" #LMPP
    #"AP351" #GMP
  )

p1 <- plotEmbedding(
    ArchRProj = proj4, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj4)
)


p2 <- plotEmbedding(
    ArchRProj = proj4, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj4)
)

plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation-GeneInte.pdf", 
    ArchRProj = proj4, 
    addDOC = FALSE, width = 5, height = 5)

plotPDF(plotList = p2, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation-GeneScore.pdf", 
    ArchRProj = proj4, 
    addDOC = FALSE, width = 5, height = 5)


####TASK 8.2
corGSM_MM <- correlateMatrices(
    ArchRProj = proj4,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI", 
    threads = 1
)
s <- sort(corGSM_MM)
head(s, 5)
tail(s, 5)

####TASK 8.3
cM <- confusionMatrix(proj4$Clusters, proj4$predictedGroup)
labelOld <- rownames(cM)
labelOld
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

#labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
#labelNew2

proj4$Clusters <- mapLabels(proj4$Clusters, newLabels = labelNew, oldLabels = labelOld)
#head(proj4$Clusters, 5)
plot_clust <- plotEmbedding(proj4, colorBy = "cellColData", name = "Clusters")
#plot_clust
plotPDF(plot_clust, name = "label-transfer.pdf", ArchRProj = proj4, addDOC = FALSE, width = 4, height = 4)
saveArchRProject(ArchRProj = proj4, outputDirectory = "ArchRProject", load = FALSE)


##################### PEAK GENE CALLING #####################
####TASK 8.4.1
proj4 <-  loadArchRProject(path = "ArchRProject", force = FALSE, showLogo = TRUE)

proj4 <- addPeak2GeneLinks(
    ArchRProj = proj4,
    reducedDims = "IterativeLSI",
    threads = 1
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj4,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = TRUE
)

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "CD19", #B-Cell Trajectory
    "CD3D", #TCells
    "CD14", #Monocytes
    "TBX21", #NK
    "CD38" #LMPP
    #"AP351" #GMP
  )

p <- plotBrowserTrack(
    ArchRProj = proj4, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(proj4)
)
plotPDF(p, name = "Peak2Gene-BrowserTrack.pdf", ArchRProj = proj4, addDOC = FALSE, width = 10, height = 10)

p2g_heatmap <- plotPeak2GeneHeatmap(ArchRProj = proj4, groupBy = "Clusters")
#p2g_heatmap
plotPDF(p2g_heatmap, name = "Peak2Gene-HeatMap.pdf", ArchRProj = proj4, addDOC = FALSE, width = 10, height = 10)
saveArchRProject(ArchRProj = proj4, outputDirectory = "ArchRProject", load = FALSE)

############################## WEEK 04 #######################################
###### DIFFERENTIAL ACCESSIBILITY ########
proj4 <-  loadArchRProject(path = "ArchRProject", force = FALSE, showLogo = TRUE)
proj5 <- addMotifAnnotations(ArchRProj = proj4, motifSet = "cisbp", name = "Motif", force = TRUE)

markerTestT_M <- getMarkerFeatures(
  ArchRProj = proj4, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "T_cells",
  bgdGroups = "Monocyte",
  threads = 1
)

pma <- markerPlot(seMarker = markerTestT_M, name = "T_cells", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markerTestT_M, name = "T_cells", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")

plotPDF(pma, pv, name = "TCells-vs-Monocyte-Markers-MA-Volcano--222", width = 8, height = 5, ArchRProj = proj4, addDOC = FALSE)

markersPeaks <- getMarkerFeatures(
    ArchRProj = proj5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  threads = 1
)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE)
plotPDF(heatmapEM, name = "TF-Motif-Enrichment-Heatmap", width = 20, height = 15, ArchRProj = proj5, addDOC = FALSE)


motifsUp <- peakAnnoEnrichment(
    seMarker = markerTestT_M,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
)
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))


motifsDo <- peakAnnoEnrichment(
    seMarker = markerTestT_M,
    ArchRProj = proj5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5",
  )

df2 <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))
head(df2)

ggDo <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))




plotPDF(ggUp, ggDo, name = "T-cells-Monocytes-Enriched", width = 5, height = 5, ArchRProj = proj5, addDOC = FALSE)


########### MOTIF FOOTPRINTING #################

getAvailableMatrices(proj5)
motifPositions <- getPositions(proj5)

motifs <- c("ZBTB7A_258", "WT1_266", "ZNF148_222", "CBFB_801", "RUNX2_732")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))


seFoot <- getFootprints(
  ArchRProj = proj5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters", 
  threads = 1
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)



######### CO-ACCESSIBILITY ##########
proj5 <- addCoAccessibility(
    ArchRProj = proj5,
    reducedDims = "IterativeLSI",
    threads = 1
)

cA <- getCoAccessibility(
    ArchRProj = proj5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "CD19", #B-Cell Trajectory
    "CD3D", #TCells
    "CD14", #Monocytes
    "TBX21", #NK
    "CD38" #LMPP
    #"AP351" #GMP
  )


plot_cA <- plotBrowserTrack(
    ArchRProj = proj5, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(proj5)
)

plotPDF(plotList = plot_cA, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = proj5, 
    addDOC = FALSE, width = 15, height = 8)
