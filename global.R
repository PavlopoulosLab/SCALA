library(shinydashboard)
library(DT)
library(shiny)
library(shinyjs)
library(SeuratObject)
library(Seurat)
library(plotly)
library(igraph)
library(rgl)
library(RColorBrewer)
library(dplyr)
library(visNetwork)
library(heatmaply)
library(gprofiler2)
library(ggplot2)
library(ggpubr)
library(CIPR) # devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = F)
library(dittoSeq) # BiocManager::install("dittoSeq")
library(slingshot) # BiocManager::install("slingshot")
library(nichenetr) # devtools::install_github("saeyslab/nichenetr") # BiocManager::install("limma")
library(tidyverse)
library(destiny) #remotes::install_github("theislab/destiny")
library(UCell) #remotes::install_github("carmonalab/UCell")
library(colorspace)
library(missMDA)
library(phateR) #pip install phate //\\ #install.packages("phateR") //\\ *devtools::install_github("scottgigante/seurat", ref="patch/add-PHATE-again") //\\ #reticulate::py_install("phate", pip=TRUE)
#ATAC libraries
library(ArchR)
library(pheatmap)
#library(DelayedMatrixStats) #BiocManager::install("DelayedMatrixStats")

#Global variables

#tab Upload
#objectInputType <- "Input10x"
seurat_object <- NULL # readRDS("seurat_processed.RDS") # NULL
init_seurat_object <- NULL
#my_metadata <- NULL
minimum_cells <- 3
minimum_features <- 200
organism <- "mouse" #or human

#tab Quality Control
qc_minFeatures <- 500
qc_maxFeatures <- 6000
qc_maxMtPercent <- 10

#tab Normalization
normalize_normMethod <- "LogNormalize"
normalize_normScaleFactor <- 10000
normalize_hvgMethod <- "vst"
normalize_hvgNGenes <- 2000
normalize_scaleRegressOut <- NULL

#tab Clustering
snn_dims <- 15
snn_k <- 20
cluster_res <- 0.6
cluster_dims <- 15

#tab DEA
markers_test <- "wilcox"
markers_minPct <- "0.1"
markers_minLogfc <- "0.25"
markers_minPval <- "0.01"
markers_logFCBase <- "avg_logFC"

#tabs Umap/tsne, DEA, Cell cycle, Trajectory
reductions_choices <- c("-")

#ATAC variables
ArrowFiles <- NULL
proj_default <- NULL

js.enrich <- "
  shinyjs.Enrich = function(url) {
    window.open(url[0]);
  }
"



# proj_default <<- loadArchRProject(path = "default/")
# 
# proj_default <<- addIterativeLSI(ArchRProj = proj_default, useMatrix = "TileMatrix", name = "IterativeLSI",
#                                  iterations = 1, varFeatures = 25000,
#                                  clusterParams = list( resolution = 0.8, sampleCells = 10000, n.start = 10),dimsToUse=1:30, force = T)
# 
# proj_default <<- addClusters(input = proj_default, reducedDims = "IterativeLSI", method = "Seurat", neme = "Clusters", #name = paste0("Clusters_res_", input$clusterResATAC),
#                              resolution = 0.8, dimsToUse = 1:30, force = T)
# 
# proj_default <<- addUMAP(ArchRProj = proj_default, reducedDims = "IterativeLSI", name = "umap", nNeighbors = 30, minDist = 0.5, metric = "cosine",
#                          force = T, n_components = 3, dimsToUse = 1:30)
# 
# proj_default <<- addUMAP(ArchRProj = proj_default, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine",
#                          force = T, n_components = 2, dimsToUse = 1:30)
# 
# cluster_names <- unique(proj_default$Clusters)
# 
# rD <- getEmbedding(ArchRProj = proj_default, embedding = "umap")
# groups <- getCellColData(ArchRProj = proj_default, select = "Clusters")
# sds <- slingshot(
#   data = rD[, 1:3], 
#   clusterLabels = groups[rownames(rD), ], 
#   start.clus = "C1", end.clus = "C7"
# )
# 
# names(sds@lineages)

# #-----Trajectory
# proj_default <<- addSlingShotTrajectories(
#   ArchRProj = proj_default,
#   groupBy = "Clusters",
#   name = "Slingshot",
#   useGroups = c("C1", "C2", "C3"),#cluster_names,
#   principalGroup = "C1",
#   embedding = "UMAP",
#   force=TRUE
# )
# p <- plotTrajectory(proj_default, trajectory = "Slingshot.Curve1", colorBy = "cellColData", name = "Slingshot.Curve1", embedding = "UMAP")
# plot(p[[1]])
# 
# if(length(colnames(proj_default@cellColData)[grep(pattern = "Slingshot", colnames(proj_default@cellColData))]) != 0)
# {
#   #delete
#   proj_default@cellColData <- proj_default@cellColData[, -grep(pattern = "Slingshot", colnames(proj_default@cellColData))]
# }
# 
# plotEmbedding(proj_default, embedding = "UMAP", colorBy = "cellColData", name = "Clusters")
# 
# meta <- as.data.frame(getCellColData(proj_default))
# meta$Cell_id <- rownames(meta)
# reduc_data <- data.frame()
# 
# #prepare colors
# cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, "Clusters"])))
# 
# #for all reductions
# archr_object_reduc <- proj_default@embeddings[["UMAP"]]$df #as.data.frame(seurat_object@reductions[[input$umapType]]@cell.embeddings)
# archr_object_reduc <- archr_object_reduc[, c(1:ncol(archr_object_reduc))]
# archr_object_reduc$Cell_id <- rownames(archr_object_reduc)
# reduc_data <- left_join(archr_object_reduc, meta)
# print(head(reduc_data))
# reduc_data$Clusters <- factor(reduc_data$Clusters, levels = sort(unique(reduc_data$Clusters))) 
# 
# fibro6a_Co_NEW_tg <- addTrajectory(
#   ArchRProj = proj_default,
#   name = "Archr_backbone",
#   groupBy = "Clusters",
#   trajectory = c("C2", "C1","C5","C6","C4"),
#   embedding = "UMAP",
#   force = TRUE
# )
# p <- plotTrajectory(fibro6a_Co_NEW_tg, trajectory = "Archr_backbone", colorBy = "cellColData", name = "Archr_backbone", embedding = "UMAP")
# plot(p[[1]])

# 
# proj_default <- addGroupCoverages(ArchRProj = proj_default, groupBy = "Clusters")
# 
# proj_default <- addSlingShotTrajectories(
#   ArchRProj = proj_default,
#   groupBy = "Clusters",
#   name = "SlingS5",
#   useGroups = c("C9","C5","C6","C10","C7","C8","C4","C2","C1","C3" ),
#   principalGroup = "C1",
#   embedding = "umap",
#   force=TRUE
# )
# p <- plotTrajectory(proj_default, trajectory = "SlingS5.Curve1", colorBy = "cellColData", name = "SlingS5.Curve1", embedding = "umap")
# df_traj <- proj_default@projectMetadata@metadata
# 
# fibro6a_Co_NEW_tg <- addTrajectory(
#   ArchRProj = proj_default,
#   name = "Archr_backbone",
#   groupBy = "Clusters",
#   trajectory = c("C1", "C2","C5","C7","C9"),
#   embedding = "umap",
#   force = TRUE
# )
# 
# plotTrajectory(fibro6a_Co_NEW_tg, trajectory = "Archr_backbone", colorBy = "cellColData", name = "Archr_backbone", embedding = "umap")
# 
# addMonocleTrajectory(
#   ArchRProj = NULL,
#   name = "Trajectory",
#   useGroups = NULL,
#   groupBy = "Clusters",
#   monocleCDS = NULL,
#   force = FALSE
# )
# 
# p <- plotBrowserTrack(
#   ArchRProj = proj_default,
#   groupBy = "Clusters",
#   geneSymbol = "Thy1",
#   upstream = 1000,
#   downstream = 1000,
#   baseSize = 15,
#   facetbaseSize = 10,
#   sizes = c(10, 4, 3, 4)
# )
# plot(p[["Thy1"]])
# 
# 
# grid::grid.newpage()
# grid::grid.draw(p$Prg4)

# 
# markers_cluster <- getMarkerFeatures(
#   ArchRProj = proj_default,
#   useMatrix = "GeneScoreMatrix",
#   groupBy = "Clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# 
# heatmap_matrix <- plotMarkerHeatmap(
#   seMarker = markers_cluster,
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.36",
#   labelMarkers = NULL,
#   transpose = TRUE,
#   returnMatrix = TRUE
# )
# 
# markers_clusterList <- unlist(getMarkers(markers_cluster, cutOff = paste0("FDR <= ",0.01," & Log2FC >= ",0.35), n = 5))
# markers_clusterList$name
# heatmaply(t(heatmap_matrix[1:10, markers_clusterList$name]))

# 
#  proj_default <<- addImputeWeights(proj_default,sampleCells = 5000)
#  
#  p <- plotEmbedding(
#    ArchRProj = proj_default, 
#    colorBy = "GeneScoreMatrix", 
#    name = "Prg4", 
#    embedding = "umap",
#    imputeWeights = getImputeWeights(proj_default)
#  )
#  p
#  
#  p$data
# # ggplotly(p)
# 
# p_mat <- getMatrixFromProject(
#   ArchRProj = proj_default,
#   useMatrix = "GeneScoreMatrix",
#   useSeqnames = NULL,
#   verbose = TRUE,
#   binarize = FALSE,
#   threads = 1,
#   logFile = createLogFile("getMatrixFromProject")
# )
# perCellGeneScore<-as.data.frame(assays(p_mat)$GeneScoreMatrix)
# perCellGeneScore$Gene_name <- p_mat@elementMetadata@listData$name
# 
# # mat <- p_mat@assays$GeneScoreMatrix
#  
#  imp_matrix
# 

# #---
# meta <- as.data.frame(getCellColData(proj_default))
# meta$Cell_id <- rownames(meta)
# reduc_data <- data.frame()
# 
# #prepare colors
# cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorByATAC])))
# 
# #for all reductions
# archr_object_reduc <- proj_default@embeddings[["umap"]]$df #as.data.frame(seurat_object@reductions[[input$umapType]]@cell.embeddings)
# archr_object_reduc <- archr_object_reduc[, c(1:ncol(archr_object_reduc))]
# archr_object_reduc$Cell_id <- rownames(archr_object_reduc)
# reduc_data <- left_join(archr_object_reduc, meta)
# print(head(reduc_data))

#gProfiler
#set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")

#TO DO
# Rename + delete buttons gia ta objects
# Selectbox gia kathe tab gia to arxeio poy ginetai processed
# input color file by user
# bug barplot in clustering doesn't appear when loading a new object