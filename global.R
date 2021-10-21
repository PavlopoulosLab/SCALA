library(shinydashboard)
library(DT)
library(shiny)
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

#Global variables

#tab Upload
objectInputType <- "Input10x"
seurat_object <- NULL #readRDS("seurat_processed.RDS")#NULL
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

#gProfiler
#set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")

#TO DO
# Rename + delete buttons gia ta objects
# Selectbox gia kathe tab gia to arxeio poy ginetai processed
# input color file by user
# bug barplot in clustering doesn't appear when loading a new object