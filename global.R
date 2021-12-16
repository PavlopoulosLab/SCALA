library(shinydashboard)
library(DT)
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinyalert)
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
library(SCENIC)
library(SCopeLoomR)
library(AUCell)
library(GSEABase)
library(RcisTarget)
library(stringr)
library(readr)
library(parallel)
library(chromVAR)
library(chromVARmotifs)
library(reticulate)
library(JASPAR2020)
library(JASPAR2018)
library(JASPAR2016) #BiocManager::install("JASPAR2020"), BiocManager::install("JASPAR2018"), BiocManager::install("JASPAR2016")

#Global variables
user_dir <- "" #user's folder in temp

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

#export tables RNA
export_metadata_RNA <- ""
export_clustertable_RNA <- ""
export_markerGenes_RNA <- ""
export_enrichedTerms_RNA <- ""
export_annotation_RNA <- ""
export_ligandReceptor_full_RNA <- ""
export_ligandReceptor_short_RNA <- ""


#ATAC variables
ArrowFiles <- NULL
proj_default <- NULL
#export tables
export_metadata_ATAC <- ""
export_clustertable_ATAC <- ""
export_markerGenes_ATAC <- ""
export_markerPeaks_ATAC <- ""
export_motifs_ATAC <- ""
export_positiveRegulators_ATAC <- ""
export_peakToGenelinks_ATAC <- ""

userMode <- FALSE

if(userMode == F)
{
  #to improve speed
}

js.enrich <- "
  shinyjs.Enrich = function(url) {
    window.open(url[0]);
  }
"

# Functions ####

# This is a void function that hides all shiny css loaders
hideAllLoaders <- function(){
  shinyjs::hide("hvgScatter_loader")
  shinyjs::hide("nFeatureViolin_loader")
  shinyjs::hide("totalCountsViolin_loader")
  shinyjs::hide("mitoViolin_loader")
  shinyjs::hide("genesCounts_loader")
  shinyjs::hide("mtCounts_loader")
  shinyjs::hide("filteredNFeatureViolin_loader")
  shinyjs::hide("filteredTotalCountsViolin_loader")
  shinyjs::hide("filteredMitoViolin_loader")
  shinyjs::hide("filteredGenesCounts_loader")
  shinyjs::hide("filteredMtCounts_loader")
  shinyjs::hide("TSS_plot_loader")
  shinyjs::hide("nFrag_plot_loader")
  shinyjs::hide("TSS_nFrag_plot_loader")
  shinyjs::hide("elbowPlotPCA_loader")
  shinyjs::hide("PCAscatter_loader")
  shinyjs::hide("PCAloadings_loader")
  shinyjs::hide("PCAheatmap_loader")
  shinyjs::hide("clusterBarplot_loader")
  shinyjs::hide("clusterBarplotATAC_loader")
  shinyjs::hide("umapPlot_loader")
  shinyjs::hide("umapPlotATAC_loader")
  shinyjs::hide("findMarkersHeatmap_loader")
  shinyjs::hide("findMarkersDotplot_loader")
  shinyjs::hide("findMarkersFeaturePlot_loader")
  shinyjs::hide("findMarkersFPfeature1_loader")
  shinyjs::hide("findMarkersFPfeature2_loader")
  shinyjs::hide("findMarkersFPfeature1_2_loader")
  shinyjs::hide("findMarkersFPcolorbox_loader")
  shinyjs::hide("findMarkersViolinPlot_loader")
  shinyjs::hide("findMarkersVolcanoPlot_loader")
  shinyjs::hide("findMarkersFeaturePlotATAC_loader")
  shinyjs::hide("snnSNN_loader")
  shinyjs::hide("findMarkersGenesHeatmapATAC_loader")
  shinyjs::hide("findMarkersGenesATACTable_loader")
  shinyjs::hide("findMarkersPeaksATACTable_loader")
  shinyjs::hide("findMarkersPeaksHeatmapATAC_loader")
  shinyjs::hide("cellCyclePCA_loader")
  shinyjs::hide("cellCycleBarplot_loader")
  shinyjs::hide("gProfilerManhatan_loader")
  shinyjs::hide("findMotifsHeatmapATAC_loader")
  shinyjs::hide("findMotifsATACTable_loader")
  shinyjs::hide("annotateClustersCIPRDotplot_loader")
  shinyjs::hide("ligandReceptorFullHeatmap_loader")
  shinyjs::hide("ligandReceptorCuratedHeatmap_loader")
  shinyjs::hide("trajectoryPlot_loader")
  shinyjs::hide("trajectoryPseudotimePlot_loader")
  shinyjs::hide("trajectoryPseudotimePlotATAC_loader")
  shinyjs::hide("grnHeatmapRNA_loader")
  shinyjs::hide("grnHeatmapATAC_loader")
  shinyjs::hide("grnATACTable_loader")
  shinyjs::hide("grnATACTable2_loader")
  shinyjs::hide("visualizeTracksOutput_loader")
}

#gProfiler
#set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")

#TO DO
# bug barplot in clustering doesn't appear when loading a new object