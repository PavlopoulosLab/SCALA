server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) #increase upload limit
  source("global.R", local=TRUE)
  
  metaD <- reactiveValues(my_project_name="-", my_seurat="-", my_activePC=1, all_lin="0")
  
  #------------------Upload tab--------------------------------
  observeEvent(input$uploadConfirm, {
    session$sendCustomMessage("handler_startLoader", c("input_loader", 10))
    session$sendCustomMessage("handler_disableButton", "uploadConfirm")
    tryCatch({
      # Create the user directory for the input and output of the analysis
      metaD$my_project_name <- input$projectID
      minimum_cells <<- input$minCells
      minimum_features <<- input$minFeatures
      organism <<- input$radioSpecies
      objectInputType <<- input$inputType
      #  userId <- "User1_"
      #  user_dir <- paste0("D:\\BSRC_Fleming\\aaa_PhD_template\\", userId, metaD$my_project_name)
      #  dir.create(user_dir)
      #  file.copy(from = input$matrix$datapath, to = paste0(user_dir, "\\", input$matrix$name), overwrite = TRUE)
      #  file.copy(from = input$barcodes$datapath, to = paste0(user_dir, "\\", input$barcodes$name), overwrite = TRUE)
      #  file.copy(from = input$genes$datapath, to = paste0(user_dir, "\\", input$genes$name), overwrite = TRUE)
      #
      #   seurat_data <- Read10X(user_dir)
      if(objectInputType == "Input10x")
      {
        seurat_data <- Read10X("User1_Tg4w/")#"hg19/"
        seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
        
        init_seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      }
      else
      {
        testMatrix <- read.table("testMatrix.txt")
        seurat_object <<- CreateSeuratObject(counts = testMatrix, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
        
        init_seurat_object <<- CreateSeuratObject(counts = testMatrix, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      }
      
      session$sendCustomMessage("handler_startLoader", c("input_loader", 25))
      if(organism == "mouse")
      {
        seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^mt-")
        init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^mt-")
      }
      else
      {
        seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^MT-")
        init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^MT-")
      }
      
      session$sendCustomMessage("handler_startLoader", c("input_loader", 50))
      #seurat_object <<- readRDS("Large_Pooled_example.RDS")#("pbmc.RDS")#("Tg4_degs.RDS") #../ALEX_Project/Fibroblast_Only_Analysis/imcs.combined_res0.5.RDS
      metaD$my_seurat <- seurat_object
      output$metadataTable <- renderDataTable(metaD$my_seurat@meta.data, options = list(pageLength = 20))
      session$sendCustomMessage("handler_startLoader", c("input_loader", 75))
      updateSelInpColor()
      updateInputGeneList()
      updateGeneSearchFP()
      # updateInputLRclusters()
      # updateInpuTrajectoryClusters()
      # print(organism)
    }, warning = function(w) {
      print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
    }, finally = { # with or without error
      session$sendCustomMessage("handler_startLoader", c("input_loader", 100))
      Sys.sleep(2) # giving some time for renderer for smoother transition
      session$sendCustomMessage("handler_finishLoader", "input_loader")
      session$sendCustomMessage("handler_enableButton", "uploadConfirm")
    })
  })
   
  #------------------Quality Control tab--------------------------------
  observeEvent(input$qcDisplay, {
    session$sendCustomMessage("handler_startLoader", c("qc_loader", 10))
    session$sendCustomMessage("handler_disableButton", "qcDisplay")
    output$nFeatureViolin <- renderPlotly(
      {
        p <- VlnPlot(seurat_object, features = c("nFeature_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
                     cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
          theme_bw() + 
          theme(
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            #axis.title.y = element_blank(),
            legend.position = "none") +
          labs(title = "", y="Genes detected/cell")
        plotly::ggplotly(p, tooltip = c("x", "y")) 
      }
    )
    output$totalCountsViolin <- renderPlotly(
      {
        p <- VlnPlot(seurat_object, features = c("nCount_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
                     cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
          theme_bw() + 
          theme(
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            #axis.title.y = element_blank(),
            legend.position = "none")+
          labs(title = "", y="Total counts/cell")
        plotly::ggplotly(p, tooltip = c("x", "y")) 
      }
    )
    output$mitoViolin <- renderPlotly(
      {
        p <- VlnPlot(seurat_object, features = c("percent.mt"), pt.size = 0.5, group.by = input$qcColorBy,
                     cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
          theme_bw() + 
          theme(
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            #axis.title.y = element_blank(),
            legend.position = "none")+
          labs(title = "", y="% of reads mapped to mitochonrial genome/cell")
        plotly::ggplotly(p, tooltip = c("x", "y"))
      }
    )
    output$mtCounts <- renderPlotly(
      {
        # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=percent.mt, color=input$qcColorBy)) + 
        #    geom_point() +
        #    theme_bw()+
        #    labs(x="Total reads/ cell", y="% of reads mapped to mitochondrial genome/ cell") +
        #    theme(legend.position = "none")
        p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = input$qcColorBy, 
                            cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
        gp <- plotly::ggplotly(p)
        print(gp)
      }
    )
    output$genesCounts <- renderPlotly(
      {
        # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=input$qcColorBy)) + 
        #    geom_point() +
        #    theme_bw()+
        #    labs(x="Total reads/ cell", y="Genes detected/ cell") +
        #    theme(legend.position = "none")
        p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = input$qcColorBy, 
                            cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
        gp <- plotly::ggplotly(p)
        print(gp)
      }
    )
    output$cellStats <- renderPrint(
      {
        print(paste0("Total number of cells: ", nrow(seurat_object@meta.data)))
      }
    )
    
    session$sendCustomMessage("handler_startLoader", c("qc_loader", 100))
    Sys.sleep(2)
    session$sendCustomMessage("handler_finishLoader", "qc_loader")
    session$sendCustomMessage("handler_enableButton", "qcDisplay")
  })
  
   observeEvent(input$qcConfirm, {
     session$sendCustomMessage("handler_startLoader", c("qc_loader", 10))
     session$sendCustomMessage("handler_disableButton", "qcConfirm")
     tryCatch({
       qc_minFeatures <<- input$minUniqueGenes
       qc_maxFeatures <<- input$maxUniqueGenes
       qc_maxMtPercent <<- input$maxMtReads
       
       session$sendCustomMessage("handler_startLoader", c("qc_loader", 50))
       seurat_object <<- subset(init_seurat_object, subset = nFeature_RNA > as.numeric(qc_minFeatures) & nFeature_RNA < as.numeric(qc_maxFeatures) & percent.mt < as.double(qc_maxMtPercent)) #filter object
       
       session$sendCustomMessage("handler_startLoader", c("qc_loader", 75))
       output$nFeatureViolin <- renderPlotly(
         {
           p <- VlnPlot(seurat_object, features = c("nFeature_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
                        cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
             theme_bw() + 
             theme(
               plot.title = element_blank(),
               axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               legend.position = "none") +
             labs(title = "", y="Genes detected/cell")
           plotly::ggplotly(p, tooltip = c("x", "y")) 
         }
       )
       output$totalCountsViolin <- renderPlotly(
         {
           p <- VlnPlot(seurat_object, features = c("nCount_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
                        cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
             theme_bw() + 
             theme(
               plot.title = element_blank(),
               axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               legend.position = "none")+
             labs(title = "", y="Total counts/cell")
           plotly::ggplotly(p, tooltip = c("x", "y")) 
         }
       )
       output$mitoViolin <- renderPlotly(
         {
           p <- VlnPlot(seurat_object, features = c("percent.mt"), pt.size = 0.5, group.by = input$qcColorBy,
                        cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy])))) + 
             theme_bw() + 
             theme(
               plot.title = element_blank(),
               axis.title.x = element_blank(),
               #axis.title.y = element_blank(),
               legend.position = "none")+
             labs(title = "", y="% of reads mapped to mitochonrial genome/cell")
           plotly::ggplotly(p, tooltip = c("x", "y"))
         }
       )
       output$mtCounts <- renderPlotly(
         {
           # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=percent.mt, color=input$qcColorBy)) + 
           #    geom_point() +
           #    theme_bw()+
           #    labs(x="Total reads/ cell", y="% of reads mapped to mitochondrial genome/ cell") +
           #    theme(legend.position = "none")
           p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = input$qcColorBy, 
                               cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
           gp <- plotly::ggplotly(p)
           print(gp)
         }
       )
       output$genesCounts <- renderPlotly(
         {
           # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=input$qcColorBy)) + 
           #    geom_point() +
           #    theme_bw()+
           #    labs(x="Total reads/ cell", y="Genes detected/ cell") +
           #    theme(legend.position = "none")
           p <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = input$qcColorBy, 
                               cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, input$qcColorBy]))))
           gp <- plotly::ggplotly(p)
           print(gp)
         }
       )
       output$cellStats <- renderPrint(
         {
           print(paste0("Total number of cells: ", nrow(seurat_object@meta.data)))
         }
       )
       
       #metaD$my_seurat <- seurat_object
       #output$metadataTable <- renderDataTable(metaD$my_seurat@meta.data, options = list(pageLength = 20))
     }, warning = function(w) {
       print(paste("Warning:  ", w))
     }, error = function(e) {
       print(paste("Error :  ", e))
       session$sendCustomMessage("handler_alert", "The selected Quality Control arguments cannot produce meaningful visualizations.")
     }, finally = {
       session$sendCustomMessage("handler_startLoader", c("qc_loader", 100))
       Sys.sleep(2)
       session$sendCustomMessage("handler_finishLoader", "qc_loader")
       session$sendCustomMessage("handler_enableButton", "qcConfirm")
     })
   })
  
  #------------------Normalization tab--------------------------------
  observeEvent(input$normalizeConfirm, {
    session$sendCustomMessage("handler_log", " ### Starting normalization procedure ###")
    session$sendCustomMessage("handler_startLoader", c("normalize_loader", 10))
    session$sendCustomMessage("handler_disableButton", "normalizeConfirm")
    tryCatch({
      normalize_normMethod <- input$radioNormalize
      normalize_normScaleFactor <- input$normScaleFactor
      seurat_object <<- NormalizeData(seurat_object, normalization.method = normalize_normMethod, scale.factor = as.numeric(normalize_normScaleFactor))
      session$sendCustomMessage("handler_log", "Finished performing log-normalization.")
      session$sendCustomMessage("handler_startLoader", c("normalize_loader", 25))
      
      normalize_hvgMethod <<- input$radioHVG
      normalize_hvgNGenes <<- input$nHVGs
      seurat_object <<- FindVariableFeatures(seurat_object, selection.method = normalize_hvgMethod, nfeatures = as.numeric(normalize_hvgNGenes))
      session$sendCustomMessage("handler_log", "Finished calculating gene and feature variances of standardized and clipped values.")
      session$sendCustomMessage("handler_startLoader", c("normalize_loader", 50))
      
      normalize_scaleRegressOut <- input$radioScaling
      # all.genes <- rownames(seurat_object) # TODO use below
      if(normalize_scaleRegressOut == "Y") seurat_object <<- ScaleData(seurat_object, vars.to.regress = "percent.mt")
      else seurat_object <<- ScaleData(seurat_object)
      session$sendCustomMessage("handler_log", "Finished centering and scaling data matrix.")
      session$sendCustomMessage("handler_startLoader", c("normalize_loader", 75))
      
      #seurat_object <<- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) # TODO move RunPCA to its own tab
      metaD$my_seurat <- seurat_object # TODO possibly remove metaD$my_seurat across script
    }, warning = function(w) {
      print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "The selected Normalization arguments cannot produce meaningful visualizations.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("normalize_loader", 100))
      Sys.sleep(2)
      session$sendCustomMessage("handler_finishLoader", "normalize_loader")
      session$sendCustomMessage("handler_enableButton", "normalizeConfirm")
      session$sendCustomMessage("handler_log", " ### Finished normalization procedure ###")
    })
  })
  
  #------------------PCA tab------------------------------------------
  observeEvent(input$PCrunPCA, {
    seurat_object <<- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) # TODO move RunPCA to its own tab
    
    output$elbowPlotPCA <- renderPlotly(
      {
        plot1 <- ElbowPlot(seurat_object, ndims = 50)
        plot1_data <- plot1$data
        colnames(plot1_data)[1] <- "PC"
        colnames(plot1_data)[2] <- "SD"
        
        p <- ggplot(plot1_data, aes(x=PC, y=SD))+
          geom_point() +
          theme_bw() +
          labs(x="PC", y="Standard Deviation")
        
        gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
        print(gp)
      }
    )
    
    output$PCAscatter <- renderPlotly(
      {
        plot1 <- DimPlot(seurat_object, reduction = "pca")
        plot1_data <- plot1$data
        
        p <- ggplot(plot1_data, aes(x=PC_1, y=PC_2, color=ident))+
          geom_point() +
          theme_bw() +
          labs(x="PC1", y="PC2")+
          theme(legend.position = "none")
        
        gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
        print(gp)
      }
    )
  })
    
  observeEvent(input$PCconfirm, {
    #metaD$my_activePC <- as.numeric(input$PCin)
    activePC <- as.numeric(input$PCin)
    
    output$PCAloadings <- renderPlotly(
      {
        plot1 <- VizDimLoadings(seurat_object, dims = activePC, reduction = "pca", balanced = TRUE)
        plot1_data <- plot1$data
        plot1_data$orig.ident <- unique(seurat_object@meta.data$orig.ident)
        colnames(plot1_data)[1] <- "PC"
        activePC <- paste0("PC", "_", activePC)
        
        p <- ggplot(plot1_data, aes(x=PC, y=feature, color=orig.ident))+
          geom_point() +
          theme_bw() +
          scale_color_manual(values="blue")+
          theme(legend.position = "none")+
          labs(x=activePC, y="")
        
        gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
        print(gp)
      }
    )
    
    output$PCAheatmap <- renderPlotly(
      {
        p <- DimHeatmap(seurat_object, dims = activePC, cells = 500, balanced = TRUE, fast = F)
        p
        plotly::ggplotly(p) 
      }
    )
    #updateSelInpColor()
  })
  #------------------Clustering tab------------------------------------------
  observeEvent(input$snnConfirm, {
    snn_dims <<- input$snnPCs
    snn_k <<- input$snnK
    cluster_res <<- input$clusterRes
    cluster_dims <<- input$clusterPCs
    
    seurat_object <<- FindNeighbors(seurat_object, k.param = as.numeric(snn_k), dims = 1:as.numeric(snn_dims), reduction = "pca")
    seurat_object <<- FindClusters(seurat_object, resolution = as.numeric(cluster_res), dims = 1:as.numeric(cluster_dims))
    seurat_object <<- RunTSNE(seurat_object, dims = 1:as.numeric(cluster_dims), seed.use = 42, dim.embed = 3, reduction = "pca")
    seurat_object <<- RunUMAP(seurat_object, dims = 1:as.numeric(cluster_dims), seed.use = 42, n.components = 3, reduction = "pca")
    
    metaD$my_seurat <- seurat_object
    
    cluster_df <- as.data.frame(table(Idents(metaD$my_seurat)))
    colnames(cluster_df)[1] <- "Cluster"
    colnames(cluster_df)[2] <- "Number of cells"
    cluster_df$`% of cells per cluster` <- cluster_df$`Number of cells`/length(metaD$my_seurat@meta.data$orig.ident)*100
    
    output$clusterTable <- renderDataTable(cluster_df, options = list(pageLength = 10), rownames = F)
    updateSelInpColor()
    updateInputLRclusters()
    updateInpuTrajectoryClusters()
    
    #******cell cycle analysis
    if(organism == "mouse")
    {
      s.genes_ver2 <- readRDS("S.phase.human.mouse.conversion.rds")
      s.genes_ver2 <- s.genes_ver2$MGI.symbol
      
      g2m.genes_ver2 <- readRDS("G2M.phase.human.mouse.conversion.rds")
      g2m.genes_ver2 <- g2m.genes_ver2$MGI.symbol
      s.genes <- s.genes_ver2
      g2m.genes <- g2m.genes_ver2
    }
    else
    {
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
    }
    
    seurat_object <<- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
    seurat_object@meta.data$CC.Difference <<- seurat_object$S.Score - seurat_object$G2M.Score
    seurat_object@meta.data$Phase <<- factor(seurat_object@meta.data$Phase, levels = c("G1", "S", "G2M"))
    metaD$my_seurat <- seurat_object
    
    print("Cell_cycle_executed")
  })
  
  #------------------DEA tab-----------------------------------------------
  observeEvent(input$findMarkersConfirm, {
    seurat_object@misc$markers <<- FindAllMarkers(seurat_object, test.use = input$findMarkersTest, min.pct = as.numeric(input$findMarkersMinPct), logfc.threshold = as.numeric(input$findMarkersLogFC), 
                                                  return.thresh = as.numeric(input$findMarkersPval), base = exp(1))
    metaD$my_seurat <- seurat_object
    updateInputGeneList()
  })
  
  #------------------gProfiler tab-----------------------------------------------
  observeEvent(input$gProfilerConfirm, {
    cluster_temp <- input$gProfilerList
    temp_df <- data.frame()
    if(cluster_temp == "all_clusters")
    {
      #all clusters used
      all_clusters <- as.character(unique(metaD$my_seurat@misc$markers$cluster))
      
      gene_lists <- list()
      for(i in 1:length(all_clusters))
      {
        #gene_lists[[i]] <- metaD$my_seurat@misc$markers[which(metaD$my_seurat@misc$markers$cluster == all_clusters[i] & metaD$my_seurat@misc$markers$avg_logFC > 0.25), 'gene']
        if(input$gProfilerLFCRadio == "Up")#UP regulated
        {
          gene_lists[[i]] <- metaD$my_seurat@misc$markers[which(metaD$my_seurat@misc$markers$cluster == all_clusters[i] & 
                                                                  metaD$my_seurat@misc$markers$avg_logFC >= as.numeric(input$gProfilerSliderLogFC) &
                                                                  metaD$my_seurat@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
        }
        else #down
        {
          gene_lists[[i]] <- metaD$my_seurat@misc$markers[which(metaD$my_seurat@misc$markers$cluster == all_clusters[i] & 
                                                                  metaD$my_seurat@misc$markers$avg_logFC <= (as.numeric(input$gProfilerSliderLogFC)*(-1)) &
                                                                  metaD$my_seurat@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
        }
      }
      
      names(gene_lists) <- all_clusters
      
      # gostres <- gost(query = gene_lists, 
      #                 organism = "mmusculus", ordered_query = FALSE, 
      #                 multi_query = F, significant = TRUE, exclude_iea = T, 
      #                 measure_underrepresentation = FALSE, evcodes = TRUE, 
      #                 user_threshold = 0.05, correction_method = "g_SCS", 
      #                 domain_scope = "annotated", custom_bg = NULL, 
      #                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
      
      gostres <- gost(query = gene_lists, 
                      organism = input$gProfilerOrganism, ordered_query = FALSE, 
                      multi_query = F, significant = TRUE, exclude_iea = F, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = as.numeric(input$gProfilerSliderSignificanceTerms), 
                      correction_method = input$gprofilerRadioCorrection, 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = input$gProfilerDatasources, as_short_link = FALSE)
      
      temp_df <- gostres$result
      
      output$gProfilerManhatan <- renderPlotly({ gostplot(gostres, capped = T, interactive = T) %>% 
      layout(width=1500, height=length(all_clusters)*400) })
      #output$gProfilerManhatan <- renderPlotly({ print(gostplot(gostres, capped = TRUE, interactive = TRUE)) })
      #output$gProfilerManhatan <- renderPlotly({ plotly::ggplotly(gostplot(gostres, capped = TRUE, interactive = TRUE)) })
      #output$gProfilerManhatan <- renderPlot({gostplot(gostres, capped = TRUE, interactive = FALSE)}, height = 4000)
    }
    else
    {
      gene_lists <- list()
      
      if(input$gProfilerLFCRadio == "Up")#UP regulated
      {
        gene_lists[[1]] <- metaD$my_seurat@misc$markers[which(metaD$my_seurat@misc$markers$cluster == cluster_temp & 
                                                                metaD$my_seurat@misc$markers$avg_logFC >= as.numeric(input$gProfilerSliderLogFC) &
                                                                metaD$my_seurat@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
      }
      else #down
      {
        gene_lists[[1]] <- metaD$my_seurat@misc$markers[which(metaD$my_seurat@misc$markers$cluster == cluster_temp & 
                                                                metaD$my_seurat@misc$markers$avg_logFC <= (as.numeric(input$gProfilerSliderLogFC)*(-1)) &
                                                                metaD$my_seurat@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
      }
      
      names(gene_lists) <- cluster_temp
      
      gostres <- gost(query = gene_lists, 
                      organism = input$gProfilerOrganism, ordered_query = FALSE, 
                      multi_query = F, significant = TRUE, exclude_iea = F, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = as.numeric(input$gProfilerSliderSignificanceTerms), 
                      correction_method = input$gprofilerRadioCorrection, 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = input$gProfilerDatasources, as_short_link = FALSE)
      
      temp_df <- gostres$result
      output$gProfilerManhatan <- renderPlotly({ print(gostplot(gostres, capped = TRUE, interactive = TRUE)) })
    }
    
    output$gProfilerTable <- renderDataTable(temp_df[, c(1, 3:6, 9:11, 16)], options = list(pageLength = 10), rownames = F)
  })
  #------------------CIPR tab-----------------------------------------------
  observeEvent(input$annotateClustersConfirm, {
    marker_genes <- metaD$my_seurat@misc$markers
    
    avgexp <- AverageExpression(metaD$my_seurat)
    avgexp <- as.data.frame(avgexp$RNA)
    avgexp$gene <- rownames(avgexp)
    input_CIPR <- marker_genes
    
    #Average expresssion methods
    if(input$annotateClustersMethod %in% c("all_genes_spearman", "all_genes_pearson"))
    {
      input_CIPR <- avgexp
    }
    
    CIPR(input_dat = input_CIPR,
         comp_method = input$annotateClustersMethod, 
         reference = input$annotateClustersReference, 
         keep_top_var = as.numeric(input$annotateClustersSlider),
         plot_ind = F,
         plot_top = F)
    
    CIPR_top_results$index <- as.character(CIPR_top_results$index)
    CIPR_top_results$index <- factor(CIPR_top_results$index, levels = as.character(seq(1:length(metaD$my_seurat$seurat_clusters))))#1:45
    cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(CIPR_top_results$cluster)))
    
    p <- ggplot(CIPR_top_results, aes(x=index, y=identity_score, fill = cluster, size=7)) +
      theme_bw() +
      geom_point(alpha=1, shape=21) +
      scale_size(range = c(7, 7), guide = F) + 
      scale_x_discrete(labels=CIPR_top_results$long_name) +
      scale_fill_manual(values = cols) +
      labs(x="") + 
      theme(legend.position="top",
            axis.text.x = element_text(vjust=0.5, hjust=1, angle = 90))
    
    output$annotateClustersCIPRTable <- renderDataTable(CIPR_top_results[], options = list(pageLength = 20), rownames = F) #remove Description
    output$annotateClustersCIPRDotplot <- renderPlotly({ print(p)})
  })
  
  #--------------------trajectory tab----------------------------------------
  observeEvent(input$trajectoryConfirm, {
    
    #delete previous lineages columns
    for_delete <- grep("Lineage", colnames(metaD$my_seurat@meta.data))
    if(length(for_delete) != 0)
    {
      metaD$my_seurat@meta.data <- metaD$my_seurat@meta.data[, -for_delete]
    }
    
    reduction <- Embeddings(metaD$my_seurat, input$trajectoryReduction)
    
    sds <- slingshot(reduction[, 1:as.numeric(input$trajectorySliderDimensions)], clusterLabels = metaD$my_seurat$seurat_clusters, 
                     start.clus = as.numeric(input$trajectoryStart), end.clus = as.numeric(input$trajectoryEnd), stretch = 0)
    
    metaD$all_lin <- slingLineages(sds)
    pt <- slingPseudotime(sds)
    
    for(i in 1:ncol(pt))
    {
      temp_lin <- pt[, i]
      metaD$my_seurat <- AddMetaData(
        object = metaD$my_seurat,
        metadata = temp_lin,
        col.name = paste0("Lineage", i)
      )
    }

    updateInputLineageList(names(metaD$all_lin))
    
    #print(paste0("after update:", names(metaD$all_lin)))
    
    plot_D <- dittoDimPlot(metaD$my_seurat, "seurat_clusters",
                           do.label = TRUE,
                           labels.repel = TRUE, 
                           reduction.use = "umap",
                           add.trajectory.lineages = slingLineages(sds),
                           trajectory.cluster.meta = "seurat_clusters",
                           color.panel = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, 'seurat_clusters']))),
                           data.out = F)
    
    output$trajectoryPlot <- renderPlot({ print(plot_D)})
    #***delete lineage columns before second run
    plot_P <- dittoDimPlot(metaD$my_seurat, var = "Lineage1",  
                           do.label = F,
                           labels.repel = F, 
                           reduction.use = "umap",
                           add.trajectory.lineages = list((metaD$all_lin[["Lineage1"]])),
                           trajectory.cluster.meta = "seurat_clusters",
                           size = 2,
                           data.out = F) + scale_colour_gradientn(colours = plasma(100), na.value = "grey90") +
                                           labs(colour = "Pseudotime") #new column paste(cluster, annotation)
    output$trajectoryPseudotimePlot <- renderPlot({ print(plot_P) })
    output$trajectoryText <- renderPrint({print(metaD$all_lin)})
  })
  
  observeEvent(input$trajectoryConfirmLineage, {
    #Lineage view
    plot_P <- dittoDimPlot(metaD$my_seurat, var = input$trajectoryLineageSelect,
                           do.label = F,
                           labels.repel = F,
                           reduction.use = "umap",
                           add.trajectory.lineages = list(metaD$all_lin[[input$trajectoryLineageSelect]]),
                           trajectory.cluster.meta = "seurat_clusters",
                           size = 2,
                           rename.var.groups = "Pseudotime",
                           data.out = F) + scale_colour_gradientn(colours = plasma(100), na.value = "grey90") +
                                           labs(colour = "Pseudotime")
    output$trajectoryPseudotimePlot <- renderPlot({ print(plot_P) })
  })
  
  #--------------------Ligand Receptor tab---------------------------
  observeEvent(input$ligandReceptorConfirm, {
    #load interactions
    ligand_target_matrix = readRDS("ligand_target_matrix.rds")
    lr_network = readRDS("lr_network.rds")
    weighted_networks = readRDS("weighted_networks.rds")
    
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
    
    if(organism == "mouse")
    {
      lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
      colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
      rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
      ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
      weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
    }
    
    #expressed genes
    ## receiver
    receiver = input$ligandReceptorSender
    expressed_genes_receiver = get_expressed_genes(receiver, metaD$my_seurat, pct = 0.10, assay_oi = "RNA")
    
    ## sender
    sender = input$ligandReceptorSender
    expressed_genes_sender = get_expressed_genes(sender, metaD$my_seurat, pct = 0.10, assay_oi = "RNA")
    
    #active ligands-receptors
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    
    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
    
    ## L-R analysis
    lr_network_top = lr_network %>% filter(from %in% potential_ligands & to %in% expressed_receptors) %>% distinct(from,to)
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    
    lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% potential_ligands & to %in% best_upstream_receptors)
    
    lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
    
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
    
    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
    p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
    #p_ligand_receptor_network
    
    output$ligandReceptorFullHeatmap <- renderPlotly({ plotly::ggplotly(p_ligand_receptor_network) })
    
    #only curated interactions
    lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
    ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
    receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
    
    lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
    lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
    
    lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
    lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
    
    dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
    dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
    
    vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
    p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "brown", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
    #p_ligand_receptor_network_strict
    
    output$ligandReceptorCuratedHeatmap <- renderPlotly({ plotly::ggplotly(p_ligand_receptor_network_strict) })
  })
  
  #------------------------------------------------------------------Output rendering
  #  output$nFeatureViolin <- renderPlotly(
  #    {
  #       p <- VlnPlot(metaD$my_seurat, features = c("nFeature_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
  #                    cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, input$qcColorBy])))) + 
  #         theme_bw() + 
  #         theme(
  #          plot.title = element_blank(),
  #          axis.title.x = element_blank(),
  #          #axis.title.y = element_blank(),
  #          legend.position = "none") +
  #          labs(title = "", y="Genes detected/cell")
  #       plotly::ggplotly(p, tooltip = c("x", "y")) 
  #    }
  #  )
  #  output$totalCountsViolin <- renderPlotly(
  #    {
  #       p <- VlnPlot(metaD$my_seurat, features = c("nCount_RNA"), pt.size = 0.5, group.by = input$qcColorBy,
  #                    cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, input$qcColorBy])))) + 
  #         theme_bw() + 
  #         theme(
  #          plot.title = element_blank(),
  #          axis.title.x = element_blank(),
  #          #axis.title.y = element_blank(),
  #          legend.position = "none")+
  #          labs(title = "", y="Total counts/cell")
  #       plotly::ggplotly(p, tooltip = c("x", "y")) 
  #    }
  #  )
  #  output$mitoViolin <- renderPlotly(
  #    {
  #      p <- VlnPlot(metaD$my_seurat, features = c("percent.mt"), pt.size = 0.5, group.by = input$qcColorBy,
  #                   cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, input$qcColorBy])))) + 
  #        theme_bw() + 
  #        theme(
  #          plot.title = element_blank(),
  #          axis.title.x = element_blank(),
  #          #axis.title.y = element_blank(),
  #          legend.position = "none")+
  #         labs(title = "", y="% of reads mapped to mitochonrial genome/cell")
  #      plotly::ggplotly(p, tooltip = c("x", "y"))
  #    }
  #  )
  #  output$mtCounts <- renderPlotly(
  #     {
  #        # p <- ggplot(metaD$my_seurat@meta.data, aes(x=nCount_RNA, y=percent.mt, color=input$qcColorBy)) + 
  #        #    geom_point() +
  #        #    theme_bw()+
  #        #    labs(x="Total reads/ cell", y="% of reads mapped to mitochondrial genome/ cell") +
  #        #    theme(legend.position = "none")
  #        p <- FeatureScatter(metaD$my_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = input$qcColorBy, 
  #                            cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, input$qcColorBy]))))
  #        gp <- plotly::ggplotly(p)
  #        print(gp)
  #     }
  #  )
  #  output$genesCounts <- renderPlotly(
  #     {
  #        # p <- ggplot(metaD$my_seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=input$qcColorBy)) + 
  #        #    geom_point() +
  #        #    theme_bw()+
  #        #    labs(x="Total reads/ cell", y="Genes detected/ cell") +
  #        #    theme(legend.position = "none")
  #       p <- FeatureScatter(metaD$my_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = input$qcColorBy, 
  #                           cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, input$qcColorBy]))))
  #        gp <- plotly::ggplotly(p)
  #        print(gp)
  #     }
  #  )
  # output$nFeatureStats <- renderPrint(
  #   {
  #      summary(metaD$my_seurat@meta.data$nFeature_RNA)
  #   }
  #   
  # )
  # output$totalCountsStats <- renderPrint(
  #    {
  #       summary(metaD$my_seurat@meta.data$nCount_RNA)
  #    }
  #    
  # )
  # output$mitoStats <- renderPrint(
  #    {
  #       if(!identical(seurat_object, NULL))
  #       {
  #         summary(metaD$my_seurat@meta.data$percent.mt)
  #       }
  #    }
  # )
  
  output$hvgScatter <- renderPlotly(
    {
      if(length(VariableFeatures(metaD$my_seurat)) != 0)
      {
        plot1 <- VariableFeaturePlot(metaD$my_seurat)
        varplot <- plot1$data
        varplot$gene <- rownames(varplot)
        varplot$colors[varplot$colors == "yes"] <- paste0("Highly Variable genes(", length(VariableFeatures(metaD$my_seurat)), ")")
        varplot$colors[varplot$colors == "no"] <- paste0("Not Variable genes(", length(rownames(metaD$my_seurat)) - length(VariableFeatures(metaD$my_seurat)), ")")
        
        if(normalize_hvgMethod == "vst")
        {
          p <- ggplot(varplot, aes(x=log10(mean), y=variance.standardized, color=colors, label=gene)) + 
            geom_point() +
            theme_bw() +
            #scale_color_manual(values = c("black", "red")) +
            scale_color_manual(
              #labels = paste(c('Non-variable', 'Variable'), 'count:', table(varplot$colors)),
              values = c("red", "black")
            )+
            labs(x="Average Expression", y="Standardized Variance", color="")  
        }
        else
        {
          p <- ggplot(varplot, aes(x=mean, y=dispersion.scaled, color=colors, label=gene)) + 
            geom_point() +
            theme_bw() +
            scale_color_manual(
              #labels = paste(c('Non-variable', 'Variable'), 'count:', table(varplot$colors)),
              values = c("red", "black")
            )+
            labs(x="Average Expression", y="Scaled Dispersion", color="")
        }
        
        gp <- plotly::ggplotly(p, tooltip = c("label", "x", "y"))
        print(gp)  
      }
    }
  )
  
  # output$elbowPlotPCA <- renderPlotly(
  #   {
  #     plot1 <- ElbowPlot(metaD$my_seurat, ndims = 50)
  #     plot1_data <- plot1$data
  #     colnames(plot1_data)[1] <- "PC"
  #     colnames(plot1_data)[2] <- "SD"
  #     
  #     p <- ggplot(plot1_data, aes(x=PC, y=SD))+
  #       geom_point() +
  #       theme_bw() +
  #       labs(x="PC", y="Standard Deviation")
  #     
  #     gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
  #     print(gp)
  #   }
  # )
  # 
  # output$PCAscatter <- renderPlotly(
  #   {
  #     plot1 <- DimPlot(metaD$my_seurat, reduction = "pca")
  #     plot1_data <- plot1$data
  #     
  #     p <- ggplot(plot1_data, aes(x=PC_1, y=PC_2, color=ident))+
  #       geom_point() +
  #       theme_bw() +
  #       labs(x="PC1", y="PC2")+
  #       theme(legend.position = "none")
  #     
  #     gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
  #     print(gp)
  #   }
  # )
  
  # output$PCAloadings <- renderPlotly(
  #   {
  #     plot1 <- VizDimLoadings(metaD$my_seurat, dims = metaD$my_activePC, reduction = "pca", balanced = TRUE)
  #     plot1_data <- plot1$data
  #     plot1_data$orig.ident <- unique(metaD$my_seurat@meta.data$orig.ident)
  #     colnames(plot1_data)[1] <- "PC"
  #     activePC <- paste0("PC", "_", metaD$my_activePC)
  #     
  #     p <- ggplot(plot1_data, aes(x=PC, y=feature, color=orig.ident))+
  #       geom_point() +
  #       theme_bw() +
  #       scale_color_manual(values="blue")+
  #       theme(legend.position = "none")+
  #       labs(x=activePC, y="")
  #     
  #     gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
  #     print(gp)
  #   }
  # )
  # 
  # output$PCAheatmap <- renderPlotly(
  #   {
  #     p <- DimHeatmap(metaD$my_seurat, dims = metaD$my_activePC, cells = 500, balanced = TRUE, fast = F)
  #     p
  #     plotly::ggplotly(p) 
  #   }
  # )
  
  output$cellCyclePCA <- renderPlotly(
  {
    if(!is.null(metaD$my_seurat))
    {
      meta <- metaD$my_seurat@meta.data
      meta$Cell_id <- rownames(meta)
      label_x <- "PC_1"
      label_y <- "PC_2"
      selected_Reduction <- input$cellCycleReduction
      #
      if(selected_Reduction == "umap")
      {
        label_x <- "UMAP_1"
        label_y <- "UMAP_2"
      }
      else if(selected_Reduction == "tsne")
      {
        label_x <- "tSNE_1"
        label_y <- "tSNE_2"
      }
      
      plot1 <- DimPlot(metaD$my_seurat, reduction = selected_Reduction)
      plot1_data <- plot1$data
      plot1_data$Cell_id <- rownames(plot1_data)
      plot1_data <- left_join(plot1_data, meta)

      p <- ggplot(plot1_data, aes_string(x=label_x, y=label_y, color="Phase"))+
        geom_point() +
        theme_bw() +
        labs(x=label_x, y=label_y)+
        theme(legend.position = "right")

      gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
      print(gp)
    }
   }
  )
  
  output$cellCycleBarplot <- renderPlotly(
    {
      #barplot
      fsize <- 24
      
      obj_meta <- metaD$my_seurat@meta.data
      obj_meta$cell_id <- row.names(obj_meta)
      obj_meta$Cluster_name <- Idents(metaD$my_seurat)
      obj_meta <- obj_meta[, c('cell_id', 'Cluster_name', 'Phase')]
      final_df_seurat <- obj_meta %>% dplyr::group_by(Cluster_name) %>% dplyr::count(Phase)
      colnames(final_df_seurat)[3] <- "Cells"
      
      total_cells_per_cluster <- as.data.frame(table(Idents(metaD$my_seurat)))
      colnames(total_cells_per_cluster)[1] <- "Cluster_name"
      colnames(total_cells_per_cluster)[2] <- "Total"
      
      final_df_seurat <- left_join(final_df_seurat, total_cells_per_cluster)
      final_df_seurat$Percentage <- (final_df_seurat$Cells/final_df_seurat$Total)*100
      
      p <- ggplot(final_df_seurat, aes(x=Phase, y=Percentage, fill=Phase))+
        geom_bar(stat='identity')+
        facet_wrap(~Cluster_name) +
        theme_bw() +
        #scale_fill_manual(values = c("red", "green", "blue")) +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = fsize, angle = 90, vjust = 0.5),
              axis.text.y = element_text(face = "bold", color = "black", size = 18, angle = 0),
              axis.title.y = element_text(face = "bold", color = "black", size = fsize),
              axis.title.x = element_text(face = "bold", color = "black", size = fsize),
              legend.text = element_text(size = 18, color = "black", face = "bold.italic"),
              legend.title = element_text(size = 18, color = "black", face = "bold.italic"),
              strip.text.x = element_text(size = 18, color = "black", face = "bold.italic"),
              strip.background = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid")) +
        labs(y="% of cells", x="Phase", color="Phase")
      p
      plotly::ggplotly(p)
    }
  )
  
  output$snnSNN <- renderVisNetwork(
    {
      if(!is.null(metaD$my_seurat@graphs$RNA_snn))
      {
        #set.seed(9)
        mygraph <- as.matrix(metaD$my_seurat@graphs$RNA_snn)
        graphOut <- graph_from_adjacency_matrix(mygraph, mode = "undirected", weighted = T, diag = F)
        graphSimple <- simplify(graphOut) #, remove.loops=T)
        weights <- E(graphSimple)$weight
        sub_nodes <- V(graphSimple)$name
        
        tableCl <- metaD$my_seurat@meta.data[, ]
        tableCl$Cell_id <- rownames(tableCl)
        tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
        tableCl <- tableCl[order(tableCl$seurat_clusters), ]
        tableCl <- tableCl[sub_nodes, ]
        colors_cl <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
        tableCl$color <- colors_cl[as.numeric(tableCl$seurat_clusters)]
        
        V(graphSimple)$color <- tableCl$color
        
        visIgraph(graphSimple, layout = "layout_with_lgl", randomSeed = 9) %>% 
        visInteraction(navigationButtons = TRUE, hover = TRUE)    
      }
    }
  )
  
  output$clusterTable <- renderDataTable(
    {
      if(!is.null(metaD$my_seurat@meta.data))
      {
        cluster_df <- as.data.frame(table(Idents(metaD$my_seurat)))
        colnames(cluster_df)[1] <- "Cluster"
        colnames(cluster_df)[2] <- "Number of cells"
        cluster_df$`% of cells per cluster` <- cluster_df$`Number of cells`/length(metaD$my_seurat@meta.data$orig.ident)*100
        
        output$clusterTable <- renderDataTable(cluster_df, options = list(pageLength = 20))
      }
    }
  )
  
  output$clusterBarplot <- renderPlotly(
    {
      #barplot for cell distribution per cluster
      
      
      clusterTable <- as.data.frame(table(Idents(metaD$my_seurat)))
      totalCells <- sum(clusterTable$Freq)
      
      clusterTable$Perc <- (clusterTable$Freq*100)/totalCells
      colnames(clusterTable)[1] <- "Cluster"
      colnames(clusterTable)[2] <-  "Cells"
      colnames(clusterTable)[3] <- "Percentage"
      
      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(clusterTable$Cluster)))
      
      p <- ggplot(clusterTable) + theme_bw() +
        geom_bar( mapping = aes(x = Cluster, y = Percentage, fill=Cluster), stat = "identity" ) +
        scale_fill_manual(values = cols ) +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
              axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
              axis.title.y = element_text(face = "bold", color = "black", size = 12),
              axis.title.x = element_text(face = "bold", color = "black", size = 12),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        labs(x="", y="Percentage % of cells", fill="Clusters")
      gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
      print(gp)
    }
  )
  
  output$umapPlot <- renderPlotly(
    {
      #get input
      dims <- as.numeric(input$umapDimensions)
      type <- input$umapType
      
      #prepare metadata
      meta <- metaD$my_seurat@meta.data
      meta$Cell_id <- rownames(meta)
      meta <- meta[, ]#meta[, c('Cell_id', 'seurat_clusters', 'orig.ident')]
      reduc_data <- data.frame()
      
      #prepare colors
      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy])))
      
      #umap data frame
      seurat_object_reduc <- as.data.frame(metaD$my_seurat@reductions$umap@cell.embeddings)
      seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
      seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
      reduc_data <- left_join(seurat_object_reduc, meta)
      
      #tsne data frame
      seurat_object_reduc2 <- as.data.frame(metaD$my_seurat@reductions$tsne@cell.embeddings)
      seurat_object_reduc2 <- seurat_object_reduc2[, c(1:ncol(seurat_object_reduc2))]
      seurat_object_reduc2$Cell_id <- rownames(seurat_object_reduc2)
      reduc_data2 <- left_join(seurat_object_reduc2, meta)
      
      if(type == "umap" & dims == 2)
      {
        p <- ggplot(data=reduc_data, aes_string(x="UMAP_1", y="UMAP_2", fill=input$umapColorBy)) +
          geom_point(size=4, shape=21)+
          scale_fill_manual(values = cols)+
          scale_size()+
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.title.y = element_text(face = "bold", color = "black", size = 25),
                axis.title.x = element_text(face = "bold", color = "black", size = 25),
                legend.text = element_text(face = "bold", color = "black", size = 9),
                legend.title = element_text(face = "bold", color = "black", size = 9),
                legend.position="right",
                title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
          labs(x="UMAP 1", y="UMAP 2", color="Cell type", title = "", fill="Color")
        plotly::ggplotly(p)  
      }
      else if(type == "umap" & dims == 3)
      {
        plot_ly(reduc_data, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, type="scatter3d", alpha = 1, mode="markers", color=as.formula(paste0('~', input$umapColorBy)),
                marker = list(size = 6, 
                              line = list(color = 'black', width = 1)
                              ),
                colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
      }
      else if(type == "tsne" & dims == 2)
      {
        p <- ggplot(data=reduc_data2, aes_string(x="tSNE_1", y="tSNE_2", fill=input$umapColorBy)) +
          geom_point(size=4, shape=21)+
          scale_fill_manual(values = cols)+
          scale_size()+
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.title.y = element_text(face = "bold", color = "black", size = 25),
                axis.title.x = element_text(face = "bold", color = "black", size = 25),
                legend.text = element_text(face = "bold", color = "black", size = 9),
                legend.title = element_text(face = "bold", color = "black", size = 9),
                legend.position="right",
                title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
          labs(x="tSNE 1", y="tSNE 2", color="Cell type", title = "", fill="Color")
        plotly::ggplotly(p)  
      }
      else if(type == "tsne" & dims == 3)
      {
        plot_ly(reduc_data2, x=~tSNE_1, y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", alpha = 1, size = 40, color=as.formula(paste0('~', input$umapColorBy)), 
                marker = list(size = 6, 
                              line = list(color = 'black', width = 1)
                ),
                colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
      }
    }
  )
  
  output$findMarkersTable <- renderDataTable(
    {
      if (!is.null(metaD$my_seurat@misc$markers))
      {
        output$findMarkersTable <- renderDataTable(metaD$my_seurat@misc$markers, options = list(pageLength = 20), filter = 'top', rownames = FALSE)  
      }
    }
  )
  
  output$findMarkersHeatmap <- renderPlotly(
    {
      top10 <- metaD$my_seurat@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
      set.seed(9)
      downsampled <- subset(metaD$my_seurat, cells = sample(Cells(metaD$my_seurat), 1500))
      
      scaled_tabe <- as.data.frame(downsampled@assays$RNA@scale.data)
      scaled_tabe$gene <- rownames(scaled_tabe)
      scaled_tabe_order <- as.data.frame(top10$gene)
      colnames(scaled_tabe_order)[1] <- "gene"
      scaled_tabe_final <- left_join(scaled_tabe_order, scaled_tabe)
      scaled_tabe_final <- na.omit(scaled_tabe_final)
      
      tableCl <- downsampled@meta.data[, ]
      tableCl$Cell_id <- rownames(tableCl)
      tableCl <- tableCl[, c('Cell_id', 'seurat_clusters')]
      tableCl <- tableCl[order(tableCl$seurat_clusters), ]
      
      clip<-function(x, min=-2, max=2) {
        x[x<min]<-min; 
        x[x>max]<-max; 
        x
      }
      
      final_mat <- scaled_tabe_final[, -1]
      final_mat <- final_mat[, tableCl$Cell_id]
      
      cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(tableCl$seurat_clusters)))
      names(cols) <- unique(tableCl$seurat_clusters)
      heatmaply(clip(final_mat), Rowv = F, Colv = F, colors = rev(RdBu(256)), showticklabels = c(F, T), labRow  = scaled_tabe_final$gene, 
                col_side_colors = tableCl$seurat_clusters, col_side_palette =  cols)
    }
  )
  
  output$findMarkersDotplot <- renderPlotly(
    {
      top10 <- metaD$my_seurat@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
      p <- DotPlot(metaD$my_seurat, features = rev(unique(top10$gene)), dot.scale = 6, cols = c("grey", "red")) + RotatedAxis()
      plotly::ggplotly(p)
    }
  )
  
  output$findMarkersFeaturePlot <- renderPlotly(
    {
      geneS <- input$findMarkersGeneSelect
      plot <- FeaturePlot(metaD$my_seurat, features = geneS, pt.size = 1.5, label = T, label.size = 5, cols = c("lightgrey", "red"), order = T, reduction = "umap") + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.title.y = element_text(face = "bold", color = "black", size = 25),
              axis.title.x = element_text(face = "bold", color = "black", size = 25),
              legend.text = element_text(face = "bold", color = "black", size = 9),
              legend.title = element_text(face = "bold", color = "black", size = 9),
              #legend.position="right",
              title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
        labs(x="UMAP 1", y="UMAP 2", title = geneS, fill="Expression")
      gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
      gp
    }
  )
  
  output$findMarkersViolinPlot <- renderPlotly(
    {
      geneS <- input$findMarkersGeneSelect2
      plot <- VlnPlot(metaD$my_seurat, features = geneS, pt.size = 0, 
                      cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(metaD$my_seurat@meta.data[, 'seurat_clusters'])))) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.title.y = element_text(face = "bold", color = "black", size = 25),
              axis.title.x = element_text(face = "bold", color = "black", size = 25),
              legend.text = element_text(face = "bold", color = "black", size = 9),
              legend.title = element_text(face = "bold", color = "black", size = 9),
              #legend.position="right",
              title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
        labs(x="Cluster", y="Expression", title = geneS, fill="Cluster")
      gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
      gp
    }
  )
  
  output$findMarkersVolcanoPlot <- renderPlotly(
    {
      diff_exp_genes <- metaD$my_seurat@misc$markers
      cluster_degs <- diff_exp_genes[which(diff_exp_genes$cluster == input$findMarkersClusterSelect), ]
      cluster_degs$status <- "Down regulated"
      for(i in 1:length(cluster_degs$gene))
      {
        if(cluster_degs$avg_logFC[i] > 0)
        {
          cluster_degs$status[i] <- "Up regulated"
        }
      }
      cluster_degs$log10Pval <- -log10(cluster_degs$p_val)
      
      p <- ggplot(data=cluster_degs, aes_string(x="avg_logFC", y="log10Pval", fill="status", label="gene", color="status")) +
        geom_point(size=1, shape=16)+
        scale_fill_manual(values = c("cyan3", "orange"))+
        scale_color_manual(values = c("cyan3", "orange"))+
        scale_size()+
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
              axis.title.y = element_text(face = "bold", color = "black", size = 25),
              axis.title.x = element_text(face = "bold", color = "black", size = 25),
              legend.text = element_text(face = "bold", color = "black", size = 9),
              legend.title = element_text(face = "bold", color = "black", size = 9),
              legend.position="right",
              title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
        labs(x="average LogFC", y="-log10(Pvalue)", fill="Color", color="")
      plotly::ggplotly(p, tooltip = c("x", "y", "label"))
    }
  )
  
  
  #function update selectInput
  updateSelInpColor <- function()
  {
    colnames(metaD$my_seurat@meta.data)
    tableMeta <- metaD$my_seurat@meta.data
    f <- sapply(tableMeta, is.factor)
    factors <- colnames(tableMeta[, f])
    # print(factors)
    updateSelectInput(session, "umapColorBy", choices = factors)
    updateSelectInput(session, "qcColorBy", choices = factors)
  }
  
  updateGeneSearchFP <- function()
  {
    total_genes <- rownames(metaD$my_seurat@assays$RNA@counts)
    updateSelectizeInput(session, 'findMarkersGeneSelect', choices = total_genes, server = TRUE) # server-side selectize drastically improves performance
    updateSelectizeInput(session, 'findMarkersGeneSelect2', choices = total_genes, server = TRUE)
  }
  
  updateInputGeneList <- function()
  {
    all_cluster_names <- as.character(unique(metaD$my_seurat@misc$markers$cluster))
    gene_list_choises <- c(all_cluster_names, "all_clusters")
    updateSelectInput(session, "gProfilerList", choices = gene_list_choises)
  }
  
  updateInputLineageList <- function(lin_names)
  {
    updateSelectInput(session, "trajectoryLineageSelect", choices = lin_names)
  }
  
  updateInputLRclusters <- function()
  {
    all_cluster_names <- (levels(metaD$my_seurat@meta.data[, 'seurat_clusters']))
    updateSelectInput(session, "ligandReceptorSender", choices = all_cluster_names)
    updateSelectInput(session, "ligandReceptorReciever", choices = all_cluster_names)
  }
  
  updateInpuTrajectoryClusters <- function()
  {
    all_cluster_names <- (levels(metaD$my_seurat@meta.data[, 'seurat_clusters']))
    updateSelectInput(session, "trajectoryStart", choices = all_cluster_names)
    updateSelectInput(session, "trajectoryEnd", choices = all_cluster_names)
  }
}



