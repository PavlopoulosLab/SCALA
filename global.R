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


# #
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