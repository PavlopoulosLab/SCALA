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

#gProfiler
#set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e102_eg49_p15")


#**** summary(df) check dataframes columns before and after binding

# seurat_object <- readRDS("Tg4_clustered.RDS")
# reduction <- Embeddings(seurat_object, "umap")
# sds <- slingshot(reduction[, 1:3], clusterLabels = seurat_object$seurat_clusters, 
#                  start.clus = 0, end.clus = 4, stretch = 0)
# ptime <- as.data.frame(slingPseudotime(sds, na=F))
# ptime$Cell_id <- rownames(ptime)
# seurat_object@meta.data$Cell_id <- rownames(seurat_object@meta.data)
# seurat_object@meta.data <- left_join(seurat_object@meta.data, ptime)
# 
# dittoHeatmap(seurat_object, genes=c("Prg4", "Cd44", "Timp1", "Smoc2", "Thy1", "Ccl11"),
#              annot.by = c("curve1"),
#              order.by = "orig.ident", cluster_cols = F, annot.colors = "orange")
             
#dev.off()
#picker input
# pickerInput("datasources", "Select datasources", 
#             choices=list('Gene Ontology'=list("Gene Ontology-Molecular Function (GO:MF)"="GO:MF", "Gene Ontology-Cellular Component (GO:CC)"= "GO:CC", "Gene Ontology-Biological Process (GO:BP)"="GO:BP"),
#                          'Biological Pathways'= list("KEGG PATHWAY"="KEGG", "Reactome"="REAC", "WikiPathways"="WP"),
#                          
#                          options = list('actions-box' = TRUE), multiple = TRUE,
#                          selected = c("Gene Ontology-Molecular Function (GO:MF)"="GO:MF","Gene Ontology-Cellular Component (GO:CC)"= "GO:CC",
#                                       "Gene Ontology-Biological Process (GO:BP)"="GO:BP","KEGG PATHWAY"="KEGG")

#
# my_seurat <- readRDS("Tg4_degs.RDS")
# all_clusters <- as.character(unique(my_seurat@misc$markers$cluster))
# 
# gene_lists <- list()
# for(i in 1:length(all_clusters))
# {
#   gene_lists[[i]] <- my_seurat@misc$markers[which(my_seurat@misc$markers$cluster == all_clusters[i] & my_seurat@misc$markers$avg_logFC > 0.25), 'gene']
# }
# 
# names(gene_lists) <- all_clusters
# gostres <- gost(query = gene_lists,
#                 organism = "mmusculus", ordered_query = FALSE,
#                 multi_query = F, significant = TRUE, exclude_iea = F,
#                 measure_underrepresentation = FALSE, evcodes = TRUE,
#                 user_threshold = 0.05, correction_method = "g_SCS",
#                 domain_scope = "annotated", custom_bg = NULL,
#                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
# plot <- gostplot(gostres, capped = T, interactive = F)
# #ggplotly(plot, tooltip = c("x", "y")) 
# gostplot(gostres, capped = T, interactive = T) %>% 
#   layout(width=8*500, height=8*500)

# This function calculates the height of the barplot plot. The height is variable and depends on the value of slider
# @param num_entries: integer value of slider, with number of entries to print
# @return height: total calculated pixels to be assigned to div height
# height_barplot<-function(num_entries){
#   height <- paste( ((num_entries*20) + 100), "px", sep="")
#   return(height)
# }
# 
# output$barplot <- renderUI({plotOutput("barplot1", height = height_barplot(sliderBarplot))})

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#output$barplot <-renderUI({plotOutput("barplot1", height = height_barplot(input$slider_barplot))})
#height_barplot<-function(max_height){
  
#  height<- paste( ((max_height*20)+ 100), "px", sep="")
#}

# seurat_object <- readRDS("Tg4_degs.RDS")
# colnames(seurat_object@meta.data)
# metaD <- seurat_object@meta.data
# f <- sapply(metaD, is.factor)
# factors <- colnames(metaD[, f])

# set.seed(9)
# mygraph <- as.matrix(seurat_object@graphs$RNA_snn)
# graphOut <- graph_from_adjacency_matrix(mygraph, mode = "undirected", weighted = T)
# graphSimple <- simplify(graphOut, remove.loops=T)
# weights <- E(graphSimple)$weight
# sub_nodes <- V(graphSimple)$name
# 
# tableCl <- seurat_object@meta.data[, ]
# tableCl$Cell_id <- rownames(tableCl)
# tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
# tableCl <- tableCl[order(tableCl$seurat_clusters), ]
# tableCl <- tableCl[sub_nodes, ]
# colors_cl <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
# tableCl$color <- colors_cl[as.numeric(tableCl$seurat_clusters)]
# 
# V(graphSimple)$color <- tableCl$color
# 
# visIgraph(graphSimple, layout = "layout_with_lgl") %>% 
#   visInteraction(navigationButtons = TRUE, hover = TRUE)


# seurat_object <- readRDS("Tg4_degs.RDS")
# top10 <- seurat_object@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# scaled_tabe <- as.data.frame(seurat_object@assays$RNA@scale.data)
# scaled_tabe$gene <- rownames(scaled_tabe)
# scaled_tabe_order <- as.data.frame(top10$gene)
# colnames(scaled_tabe_order)[1] <- "gene"
# scaled_tabe_final <- left_join(scaled_tabe_order, scaled_tabe)
# scaled_tabe_final <- na.omit(scaled_tabe_final)
# #scaled_tabe_final <- scaled_tabe[which(rownames(scaled_tabe) %in% top10$gene), ]
# 
# tableCl <- seurat_object@meta.data[, ]
# tableCl$Cell_id <- rownames(tableCl)
# tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
# tableCl <- tableCl[order(tableCl$seurat_clusters), ]
# 
# clip<-function(x, min=-2, max=2) {
#   x[x<min]<-min; 
#   x[x>max]<-max; 
#   x
# }
# 
# final_mat <- scaled_tabe_final[, -1]
# final_mat <- final_mat[, tableCl$Cell_id]
# 
# cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
# names(cols) <- unique(tableCl$seurat_clusters)
# heatmaply(clip(final_mat), Rowv = F, Colv = F, colors = rev(RdBu(256)), showticklabels = c(F, T), labRow  = scaled_tabe_final$gene, 
#           col_side_colors = tableCl$seurat_clusters, col_side_palette =  cols)


# p <- FeaturePlot(seurat_object, features = "Prg4", label = T, label.size = 5, pt.size = 2)
# pdata <- p$data
# pdata$Cell_id <- rownames(pdata)
# seurat_object_reduc <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
# seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
# final_df <- left_join(seurat_object_reduc, pdata)
# #plotly 2D
# Noax <- list(
#   title = "",
#   zeroline = FALSE,
#   showline = T,
#   showticklabels = T,
#   showgrid = T
# )
# 
# p <- plot_ly(final_df, x= ~UMAP_1, y= ~UMAP_2
#         marker = list(size = 10) 
#         ) %>%
#   layout(xaxis = Noax, yaxis = Noax) %>% add_trace(color = ~Prg4, type="scatter", name="Prg4")
# p
# ggplotly(p)
# 
# 
# d <- diamonds[sample(nrow(diamonds), 1000), ]
# 
# fig <- plot_ly(
#   d, x = ~carat, y = ~price,
#   color = ~carat, size = ~carat
# )
# 
# fig
# 
# geneS <- "Actb"
# plot <- FeaturePlot(seurat_object, features = geneS, pt.size = 1.5, label = T, label.size = 5, cols = c("lightgrey", "red"), order = T, reduction = "umap") + 
#      theme_bw() +
#      theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
#            axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
#            axis.title.y = element_text(face = "bold", color = "black", size = 25),
#            axis.title.x = element_text(face = "bold", color = "black", size = 25),
#            legend.text = element_text(face = "bold", color = "black", size = 9),
#            legend.title = element_text(face = "bold", color = "black", size = 9),
#            #legend.position="right",
#            title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
#      labs(x="UMAP 1", y="UMAP 2", title = geneS, fill="Expression")
# plot
# gp <- ggplotly(plot, tooltip = c("x", "y", geneS))
# gp
# 
# hl <- HoverLocator(plot = plot, information = FetchData(seurat_object, vars = c("seurat_clusters", "Prg4")))
# ggplotly(plot)
# 
# plot_ly(final_df, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", alpha = 1, size = 10, color =~Prg4)

#seurat version
# seurat_object <- readRDS("Tg4_clustered.RDS")
# p <- DimPlot(seurat_object, label = T, pt.size = 1, label.size = 5, group.by = "seurat_clusters", cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data$seurat_clusters))))
# p
# ggplotly(p)
# 
# #ggplot2
# meta <- seurat_object@meta.data
# meta$Cell_id <- rownames(meta)
# meta <- meta[, c('Cell_id', 'seurat_clusters', 'orig.ident')]
# 
# seurat_object_reduc <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
# seurat_object_reduc <- seurat_object_reduc[, c(1:3)]
# seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
# reduc_data <- left_join(seurat_object_reduc, meta)
# 
# cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(reduc_data$seurat_clusters)))
# level_order3 <- factor(reduc_data$seurat_clusters)
# 
# p <- ggplot(data=reduc_data, aes(x=UMAP_1, y=UMAP_2, fill=seurat_clusters)) +
#   geom_point(size=4, shape=21)+
#   scale_fill_manual(values = cols)+
#   scale_size()+
#   theme_bw() +
#   theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
#         axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
#         axis.title.y = element_text(face = "bold", color = "black", size = 25),
#         axis.title.x = element_text(face = "bold", color = "black", size = 25),
#         legend.text = element_text(face = "bold", color = "black", size = 9),
#         legend.title = element_text(face = "bold", color = "black", size = 9),
#         legend.position="right",
#         title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
#   labs(x="UMAP 1", y="UMAP 2", color="Cell type", title = "", fill="Cluster")
# p
# 
# ggplotly(p) %>%
#   layout(legend = list(orientation = "h", x = 0.4, y = -0.2))
# 
# #plotly 3D
# plot_ly(reduc_data, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", alpha = 1, size = 10, color =~seurat_clusters, 
#         colors =colorRampPalette(brewer.pal(12, "Paired"))(length(unique(reduc_data$seurat_clusters))) )
# 
# #plotly 2D
# Noax <- list(
#   title = "",
#   zeroline = FALSE,
#   showline = T,
#   showticklabels = T,
#   showgrid = T
# )
# 
# plot_ly(reduc_data, x=~UMAP_1, y=~UMAP_2, alpha = 1, color =~seurat_clusters, 
#         colors =colorRampPalette(brewer.pal(12, "Paired"))(length(unique(reduc_data$seurat_clusters))),
#         marker = list(size = 10,
#                       line = list(color = 'black',
#                                   width = 2))
#         ) %>% 
#   layout(xaxis = Noax, yaxis = Noax) 
# 
# #seurat_object <- RunUMAP(object=seurat_object, dims = 1:15, n.components = 3, reduction = "pca")
# #seurat_object <- RunTSNE(seurat_object, reduction = "pca", dims = 1:15, dim.embed = 3)
# #p <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:15, do.plot = T)
# # 
# mygraph <- as.matrix(seurat_object@graphs$RNA_snn)
# graphOut <- graph_from_adjacency_matrix(mygraph, mode = "undirected", weighted = T)
# weights <- E(graphOut)$weight
# sub_nodes <- V(graphOut)$name
# mylayout <- layout_with_kk(graphOut, dim = 3, weights = weights)
# plot(graphOut, vertex.label=NA, vertex.size=2, layout=mylayout)
# graphSimple <- simplify(graphOut, remove.loops=T)
# rglplot(graphSimple, layout=mylayout,vertex.size=5,vertex.label=NA)
# 
# rgl.open()
# rglplot(graphOut, layout=mylayout,vertex.size=5,vertex.label=NA)
# rgl.close()
# 
# netPlot <- visIgraph(graphSimple, layout = "layout_with_lgl")

# graphOut
# head(layout)
# #DimPlot(seurat_object, reduction = "umap")
# #mygraphs <- as.matrix(seurat_object@graphs$RNA_snn)
# d3_tsne <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
# d3_tsne$Cell_id <- rownames(d3_tsne)
# plot_ly(d3_tsne, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, type="scatter3d", mode="markers", alpha = 1, size = 2)
# plot_ly(as.data.frame(layout), x=~V1, y=~V2, z=~V3, type="scatter3d", mode="markers", alpha = 1, size = 2)
# #p <- ElbowPlot(seurat_object, ndims = 50)
#plot1 <- p$data
#p <- DimPlot(seurat_object, reduction = "pca")
#plot1 <- p$data
#ggplotly(p)

#p <- DimPlot(seurat_object, reduction = "tsne", group.by = "RNA_snn_res.0.4", pt.size = 2)
#ggplotly(p)
#seurat_object <- FindClusters(seurat_object, resolution = 1.0)

# p <- ggplot(init_seurat_object@meta.data, aes(y=percent.mt, x=orig.ident, fill=orig.ident)) +
#   geom_violin() +
#   geom_hline(yintercept=0, alpha=0.5) +
#   labs(title = "",
#        x = "",
#        y = "") +
#   theme_bw()
# 
# ggplotly(p)

#TO DO
# Rename + delete buttons gia ta objects
# Selectbox gia kathe tab gia to arxeio poy ginetai processed
# Check oti exei ginei to proigoumeno bhma
# Weights + edges of k-NN
# clusters + cells belonging there
# input color file by user
# update Power slider sto 0.5