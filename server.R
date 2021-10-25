server <- function(input, output, session) {
  options(shiny.maxRequestSize=30*1024^2) #increase upload limit
  source("global.R", local=TRUE)
  
  #session$sendCustomMessage("handler_disableTabs", "sidebarMenu") # disable all tab panels (except Data Input) until files are uploaded
  metaD <- reactiveValues(my_project_name="-", all_lin="0")
  
  #------------------Upload tab--------------------------------
  observeEvent(input$uploadCountMatrixConfirm, {
    #session$sendCustomMessage("handler_disableTabs", "sidebarMenu") # disable all tab panels (except Data Input) until files are uploaded
    session$sendCustomMessage("handler_startLoader", c("input_loader", 10))
    session$sendCustomMessage("handler_disableButton", "uploadCountMatrixConfirm")
    tryCatch({
      # Create the user directory for the input and output of the analysis
      metaD$my_project_name <- input$uploadCountMatrixprojectID
      minimum_cells <<- input$uploadCountMatrixminCells
      minimum_features <<- input$uploadCountMatrixminFeatures
      organism <<- input$uploadCountMatrixRadioSpecies
      #  userId <- "User1_"
      #  user_dir <- paste0("D:\\BSRC_Fleming\\aaa_PhD_template\\", userId, metaD$my_project_name) #TODO remove 2 random barcodes, does it crash?
      #  dir.create(user_dir)
      #  file.copy(from = input$matrix$datapath, to = paste0(user_dir, "\\", input$matrix$name), overwrite = TRUE)
      #  file.copy(from = input$barcodes$datapath, to = paste0(user_dir, "\\", input$barcodes$name), overwrite = TRUE)
      #  file.copy(from = input$genes$datapath, to = paste0(user_dir, "\\", input$genes$name), overwrite = TRUE)
      #
      #   seurat_data <- Read10X(user_dir)
      
      testMatrix <- read.table("testMatrix.txt")
      seurat_object <<- CreateSeuratObject(counts = testMatrix, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
        
      init_seurat_object <<- CreateSeuratObject(counts = testMatrix, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      
      
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
      
      output$metadataTable <- renderDataTable(seurat_object@meta.data, options = list(pageLength = 20))
      session$sendCustomMessage("handler_startLoader", c("input_loader", 75))
      updateSelInpColor()
      updateInputGeneList()
      updateGeneSearchFP()
      updateQC_choices()
      cleanAllPlots(T) # fromDataInput -> TRUE
      # updateInputLRclusters()
      # updateInpuTrajectoryClusters()
      # print(organism)
      # saveRDS(seurat_object, "seurat_object.RDS")
      session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " QUALITY CONTROL", " DATA NORMALIZATION\n& SCALING", "UTILITY OPTIONS"))
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
    }, finally = { # with or without error
      session$sendCustomMessage("handler_startLoader", c("input_loader", 100))
      Sys.sleep(1) # giving some time for renderer for smoother transition
      session$sendCustomMessage("handler_finishLoader", "input_loader")
      session$sendCustomMessage("handler_enableButton", "uploadCountMatrixConfirm")
    })
  })
  
  observeEvent(input$upload10xRNAConfirm, {
    #session$sendCustomMessage("handler_disableTabs", "sidebarMenu") # disable all tab panels (except Data Input) until files are uploaded
    session$sendCustomMessage("handler_startLoader", c("input_loader", 10))
    session$sendCustomMessage("handler_disableButton", "upload10xRNAConfirm")
    tryCatch({
      # Create the user directory for the input and output of the analysis
      metaD$my_project_name <- input$upload10xRNAprojectID
      minimum_cells <<- input$upload10xRNAminCells
      minimum_features <<- input$upload10xRNAminFeatures
      organism <<- input$upload10xRNARadioSpecies
      #  userId <- "User1_"
      #  user_dir <- paste0("D:\\BSRC_Fleming\\aaa_PhD_template\\", userId, metaD$my_project_name) #TODO remove 2 random barcodes, does it crash?
      #  dir.create(user_dir)
      #  file.copy(from = input$matrix$datapath, to = paste0(user_dir, "\\", input$matrix$name), overwrite = TRUE)
      #  file.copy(from = input$barcodes$datapath, to = paste0(user_dir, "\\", input$barcodes$name), overwrite = TRUE)
      #  file.copy(from = input$genes$datapath, to = paste0(user_dir, "\\", input$genes$name), overwrite = TRUE)
      #
      #   seurat_data <- Read10X(user_dir)
      seurat_data <- Read10X("User1_Tg4w/")#"hg19/"
      seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
        
      init_seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      
      
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
      
      output$metadataTable <- renderDataTable(seurat_object@meta.data, options = list(pageLength = 20))
      session$sendCustomMessage("handler_startLoader", c("input_loader", 75))
      updateSelInpColor()
      updateInputGeneList()
      updateGeneSearchFP()
      updateQC_choices()
      cleanAllPlots(T) # fromDataInput -> TRUE
      # updateInputLRclusters()
      # updateInpuTrajectoryClusters()
      # print(organism)
      # saveRDS(seurat_object, "seurat_object.RDS")
      session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " QUALITY CONTROL", " DATA NORMALIZATION\n& SCALING", "UTILITY OPTIONS"))
      # }, warning = function(w) {
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
    }, finally = { # with or without error
      session$sendCustomMessage("handler_startLoader", c("input_loader", 100))
      Sys.sleep(1) # giving some time for renderer for smoother transition
      session$sendCustomMessage("handler_finishLoader", "input_loader")
      session$sendCustomMessage("handler_enableButton", "upload10xRNAConfirm")
    })
  })
  
  observeEvent(input$upload10xExampleRNAConfirm, {
    #session$sendCustomMessage("handler_disableTabs", "sidebarMenu") # disable all tab panels (except Data Input) until files are uploaded
    session$sendCustomMessage("handler_startLoader", c("input_loader", 10))
    session$sendCustomMessage("handler_disableButton", "upload10xExampleRNAConfirm")
    tryCatch({
      # Create the user directory for the input and output of the analysis
      seurat_data <- Read10X("hg19/")#"hg19/"
      metaD$my_project_name <- "PBMC"
      organism <<- "human"
      
      seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      
      init_seurat_object <<- CreateSeuratObject(counts = seurat_data, project = metaD$my_project_name, min.cells = as.numeric(minimum_cells), min.features = as.numeric(minimum_features))
      
      seurat_object[["percent.mt"]] <<- PercentageFeatureSet(seurat_object, pattern = "^MT-")
      init_seurat_object[["percent.mt"]] <<- PercentageFeatureSet(init_seurat_object, pattern = "^MT-")
      
      
      session$sendCustomMessage("handler_startLoader", c("input_loader", 50))
      
      output$metadataTable <- renderDataTable(seurat_object@meta.data, options = list(pageLength = 20))
      session$sendCustomMessage("handler_startLoader", c("input_loader", 75))
      updateSelInpColor()
      updateInputGeneList()
      updateGeneSearchFP()
      updateQC_choices()
      cleanAllPlots(T) # fromDataInput -> TRUE
      # updateInputLRclusters()
      # updateInpuTrajectoryClusters()
      # print(organism)
      # saveRDS(seurat_object, "seurat_object.RDS")
      session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", "UTILITY OPTIONS", " QUALITY CONTROL", " DATA NORMALIZATION\n& SCALING"))
      # }, warning = function(w) {
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "Data Input error. Please, refer to the help pages for input format.")
    }, finally = { # with or without error
      session$sendCustomMessage("handler_startLoader", c("input_loader", 100))
      Sys.sleep(1) # giving some time for renderer for smoother transition
      session$sendCustomMessage("handler_finishLoader", "input_loader")
      session$sendCustomMessage("handler_enableButton", "upload10xExampleRNAConfirm")
    })
  })
  
  #------------------Utilities------------------------------------------
  # observeEvent(input$utilitiesColorPicker, {
  #   user_palette <- choose_palette()
  #   print(class(user_palette))
  #   print("\n")
  #   print(user_palette)
  # })
  
  observeEvent(input$utilitiesConfirmRename, {
    if(!identical(seurat_object, NULL) & input$utilitiesRenameOldName != "-" & input$utilitiesRenameNewName != "")
    {
      old_idents <- unique(seurat_object$seurat_clusters)
      #print(old_idents)
      new_idents <- gsub(input$utilitiesRenameOldName, input$utilitiesRenameNewName, old_idents)
      #print(new_idents)
      names(new_idents) <- levels(old_idents)
      #print(new_idents)
      seurat_object <<- RenameIdents(seurat_object, new_idents)
      seurat_object$seurat_clusters <<- Idents(seurat_object)
      updateInputLRclusters()
    }
  })
  
  observeEvent(input$utilitiesConfirmDelete, {
    if(!identical(seurat_object, NULL) & input$utilitiesDeleteCluster != "-")
    {
      ident_to_remove <- input$utilitiesDeleteCluster
      seurat_object <<- subset(seurat_object, idents = ident_to_remove, invert = TRUE)
      seurat_object$seurat_clusters <<- Idents(seurat_object)
      updateInputLRclusters()
    }
  })
  
  output$utilitiesConfirmExport <- downloadHandler(
    filename = function() { 
      paste("processed_seurat_object-", Sys.Date(), ".RDS", sep="")
    },
    content = function(file) {
      saveRDS(seurat_object, file)
    })
  
  #------------------Quality Control tab--------------------------------
  observeEvent(input$qcDisplay, {
    session$sendCustomMessage("handler_startLoader", c("qc_loader", 10))
    session$sendCustomMessage("handler_disableButton", "qcDisplay")
    tryCatch({
      if (identical(init_seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else{
        output$nFeatureViolin <- renderPlotly(
          {
            p <- VlnPlot(init_seurat_object, features = c("nFeature_RNA"), pt.size = 0.5, group.by = "orig.ident",
                         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) + 
              theme_bw() + 
              geom_hline(yintercept=c(as.numeric(input$minUniqueGenes), as.numeric(input$maxUniqueGenes)), linetype="dashed", color = "red", size=1) + 
              theme(
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                #axis.title.y = element_blank(),
                legend.position = "none") +
              labs(title = "", y="Genes detected/cell")
            plotly::ggplotly(p, tooltip = c("x", "y")) 
          }
        )
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 25))
        
        output$totalCountsViolin <- renderPlotly(
          {
            p <- VlnPlot(init_seurat_object, features = c("nCount_RNA"), pt.size = 0.5, group.by = "orig.ident",
                         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) + 
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
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 35))
        
        output$mitoViolin <- renderPlotly(
          {
            p <- VlnPlot(init_seurat_object, features = c("percent.mt"), pt.size = 0.5, group.by = "orig.ident",
                         cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident'])))) + 
              theme_bw() + 
              geom_hline(yintercept= as.numeric(input$maxMtReads), linetype="dashed", color = "red", size=1) +
              theme(
                plot.title = element_blank(),
                axis.title.x = element_blank(),
                #axis.title.y = element_blank(),
                legend.position = "none")+
              labs(title = "", y="% of reads mapped to mitochonrial genome/cell")
            plotly::ggplotly(p, tooltip = c("x", "y"))
          }
        )
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 50))
        
        output$mtCounts <- renderPlotly(
          {
            # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=percent.mt, color=input$qcColorBy)) + 
            #    geom_point() +
            #    theme_bw()+
            #    labs(x="Total reads/ cell", y="% of reads mapped to mitochondrial genome/ cell") +
            #    theme(legend.position = "none")
            p <- FeatureScatter(init_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident", 
                                cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident']))))
            gp <- plotly::ggplotly(p)
            print(gp)
          }
        )
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 60))
        
        output$genesCounts <- renderPlotly(
          {
            # p <- ggplot(seurat_object@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=input$qcColorBy)) + 
            #    geom_point() +
            #    theme_bw()+
            #    labs(x="Total reads/ cell", y="Genes detected/ cell") +
            #    theme(legend.position = "none")
            p <- FeatureScatter(init_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident", 
                                cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(init_seurat_object@meta.data[, 'orig.ident']))))
            gp <- plotly::ggplotly(p)
            print(gp)
          }
        )
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 75))
        
        output$cellStats <- renderPrint(
          {
            cat(paste0("Total number of cells: ", nrow(init_seurat_object@meta.data)))
          }
        )
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "The selected Quality Control arguments cannot produce meaningful visualizations.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("qc_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "qc_loader")
      session$sendCustomMessage("handler_enableButton", "qcDisplay")
    })
  })
  
  observeEvent(input$qcConfirm, {
    session$sendCustomMessage("handler_startLoader", c("qc_loader", 10))
    session$sendCustomMessage("handler_disableButton", "qcConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else{
        
        qc_minFeatures <<- input$minUniqueGenes
        qc_maxFeatures <<- input$maxUniqueGenes
        qc_maxMtPercent <<- input$maxMtReads
        
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 50))
        seurat_object <<- subset(init_seurat_object, subset = nFeature_RNA > as.numeric(qc_minFeatures) & nFeature_RNA < as.numeric(qc_maxFeatures) & percent.mt < as.double(qc_maxMtPercent)) #filter object
        
        session$sendCustomMessage("handler_startLoader", c("qc_loader", 75))
        output$filteredNFeatureViolin <- renderPlotly(
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
        output$filteredTotalCountsViolin <- renderPlotly(
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
        output$filteredMitoViolin <- renderPlotly(
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
        output$filteredMtCounts <- renderPlotly(
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
        output$filteredGenesCounts <- renderPlotly(
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
        output$filteredCellStats <- renderPrint(
          {
            cat(paste0("\nTotal number of cells after filtering: ", nrow(seurat_object@meta.data)))
          }
        )
        updateRegressOut()
        updateSelInpColor()
        #session$sendCustomMessage("handler_disableTabs", "sidebarMenu")
        cleanAllPlots(F) # fromDataInput -> FALSE
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " QUALITY CONTROL", " DATA NORMALIZATION\n& SCALING", "UTILITY OPTIONS"))
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "The selected Quality Control arguments cannot produce meaningful visualizations.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("qc_loader", 100))
      Sys.sleep(1)
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
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else{
        normalize_normMethod <- normalize_normMethod
        normalize_normScaleFactor <- input$normScaleFactor
        seurat_object <<- NormalizeData(seurat_object, normalization.method = normalize_normMethod, scale.factor = as.numeric(normalize_normScaleFactor))
        session$sendCustomMessage("handler_log", "Finished performing log-normalization.")
        session$sendCustomMessage("handler_startLoader", c("normalize_loader", 25))
        
        normalize_hvgMethod <<- input$radioHVG
        normalize_hvgNGenes <<- input$nHVGs
        seurat_object <<- FindVariableFeatures(seurat_object, selection.method = normalize_hvgMethod, nfeatures = as.numeric(normalize_hvgNGenes))
        session$sendCustomMessage("handler_log", "Finished calculating gene and feature variances of standardized and clipped values.")
        session$sendCustomMessage("handler_startLoader", c("normalize_loader", 50))
        
        normalize_scaleRegressOut <- input$normalizeRegressColumns
        # all.genes <- rownames(seurat_object) # TODO use below
        if(is.null(normalize_scaleRegressOut)) seurat_object <<- ScaleData(seurat_object) 
        else seurat_object <<- ScaleData(seurat_object, vars.to.regress=normalize_scaleRegressOut)
        session$sendCustomMessage("handler_log", "Finished centering and scaling data matrix.")
        session$sendCustomMessage("handler_startLoader", c("normalize_loader", 75))
        updateSignatures()
        
        output$hvgScatter <- renderPlotly(#tooltip example TODO the same in all plots
          {
            if(length(VariableFeatures(seurat_object)) != 0)
            {
              plot1 <- VariableFeaturePlot(seurat_object)
              varplot <- plot1$data
              varplot$gene <- rownames(varplot)
              varplot$colors[varplot$colors == "yes"] <- paste0("Highly Variable genes(", length(VariableFeatures(seurat_object)), ")")
              varplot$colors[varplot$colors == "no"] <- paste0("Not Variable genes(", length(rownames(seurat_object)) - length(VariableFeatures(seurat_object)), ")")
              
              if(normalize_hvgMethod == "vst")
              {
                #p <- ggplot(varplot, aes(x=log10(mean), y=variance.standardized, color=colors, label=gene)) + 
                p <- ggplot(varplot, aes(x=log10(mean), y=variance.standardized, color=colors, text = paste0("log10(mean expression): ", log10(mean),
                                                                                                             "\nstandardized variance: ", variance.standardized,
                                                                                                             "\ngene: ", gene)
                                         )) +
                  geom_point() +
                  theme_bw() +
                  #scale_color_manual(values = c("black", "red")) +
                  scale_color_manual(
                    values = c("red", "black")
                  )+
                  labs(x="Average Expression", y="Standardized Variance", color="")  
              }
              else
              {
                p <- ggplot(varplot, aes(x=mean, y=dispersion.scaled, color=colors, text = paste0("mean expression: ", mean,
                                                                                                  "\nscaled dispersion: ", dispersion.scaled,
                                                                                                  "\ngene: ", gene)
                )) +
                  geom_point() +
                  theme_bw() +
                  scale_color_manual(
                    values = c("red", "black")
                  )+ 
                  labs(x="Average Expression", y="Scaled Dispersion", color="")
              }
              
              gp <- plotly::ggplotly(p, tooltip = "text")#c("label", "x", "y"))
              print(gp)  
            }
          }
        )
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " PRINCIPAL COMPONENT\nANALYSIS"))
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "The selected Normalization arguments cannot produce meaningful visualizations.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("normalize_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "normalize_loader")
      session$sendCustomMessage("handler_enableButton", "normalizeConfirm")
      session$sendCustomMessage("handler_log", " ### Finished normalization procedure ###")
    })
  })
  
  #------------------PCA tab------------------------------------------
  observeEvent(input$PCrunPCA, {
    session$sendCustomMessage("handler_log", " ### Starting PCA Analysis ###")
    session$sendCustomMessage("handler_startLoader", c("PCA1_loader", 10))
    session$sendCustomMessage("handler_disableButton", "PCrunPCA")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"ScaleData.RNA" %in% names(seurat_object@commands)) session$sendCustomMessage("handler_alert", "Data need to be scaled first. Please, execute the previous step in DATA NORMALIZATION & SCALING.")
      else {
        seurat_object <<- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
        updateUmapTypeChoices("pca")
        
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
        
        session$sendCustomMessage("handler_startLoader", c("PCA1_loader", 50))
        output$PCAscatter <- renderPlotly( #TODO fix error
          {
            #prepare metadata
            meta <- seurat_object@meta.data
            meta$Cell_id <- rownames(meta)
            meta <- meta[, ]#meta[, c('Cell_id', 'seurat_clusters', 'orig.ident')]
            reduc_data <- data.frame()
            
            #prepare colors
            cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$pcaColorBy])))
            
            #umap data frame
            seurat_object_reduc <- as.data.frame(seurat_object@reductions$pca@cell.embeddings)
            seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
            seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
            reduc_data <- left_join(seurat_object_reduc, meta)
            
            p <- ggplot(data=reduc_data, aes_string(x="PC_1", y="PC_2", fill=input$pcaColorBy)) +
              geom_point(size=4, shape=19)+
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
              labs(x="PC 1", y="PC 2", color="Cell type", title = "", fill="Color")
            print(plotly::ggplotly(p))
          }
        )
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " CLUSTERING", " CELL CYCLE PHASE ANALYSIS"))
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the PCA analysis.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("PCA1_loader", 75))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "PCA1_loader")
      session$sendCustomMessage("handler_enableButton", "PCrunPCA")
      session$sendCustomMessage("handler_log", " ### Finished PCA Analysis ###")
    })
  })
  
  observeEvent(input$PCconfirm, {
    session$sendCustomMessage("handler_log", " ### Starting PC Exploration ###")
    session$sendCustomMessage("handler_startLoader", c("PCA2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "PCconfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first Run PCA above.")
      else {
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
        
        session$sendCustomMessage("handler_startLoader", c("PCA2_loader", 75))
        output$PCAheatmap <- renderPlotly(
          {
            p <- DimHeatmap(seurat_object, dims = activePC, cells = 500, balanced = TRUE, fast = F)
            p
            plotly::ggplotly(p) 
          }
        )
        #updateSelInpColor()
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the PCA analysis.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("PCA2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "PCA2_loader")
      session$sendCustomMessage("handler_enableButton", "PCconfirm")
      session$sendCustomMessage("handler_log", " ### Finished PC Exploration ###")
    })
  })
  
  #------------------Clustering tab------------------------------------------
  observeEvent(input$snnConfirm, {
    session$sendCustomMessage("handler_startLoader", c("clust1_loader", 10))
    session$sendCustomMessage("handler_disableButton", "snnConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first execute PRINCIPAL COMPONENT ANALYSIS.")
      else {
        snn_dims <<- input$snnPCs
        snn_k <<- input$snnK
        cluster_res <<- input$clusterRes
        cluster_dims <<- input$clusterPCs
        
        seurat_object <<- FindNeighbors(seurat_object, k.param = as.numeric(snn_k), dims = 1:as.numeric(snn_dims), reduction = "pca")
        seurat_object <<- FindClusters(seurat_object, resolution = as.numeric(cluster_res), dims = 1:as.numeric(cluster_dims))
        
        session$sendCustomMessage("handler_startLoader", c("clust1_loader", 25))
        cluster_df <- as.data.frame(table(Idents(seurat_object)))
        colnames(cluster_df)[1] <- "Cluster"
        colnames(cluster_df)[2] <- "Number of cells"
        cluster_df$`% of cells per cluster` <- cluster_df$`Number of cells`/length(seurat_object@meta.data$orig.ident)*100
        
        session$sendCustomMessage("handler_startLoader", c("clust1_loader", 50))
        output$clusterTable <- renderDataTable(cluster_df, options = list(pageLength = 10), rownames = F)
        
        session$sendCustomMessage("handler_startLoader", c("clust1_loader", 70))
        output$snnSNN <- renderVisNetwork(
          {
            if(!is.null(seurat_object@graphs$RNA_snn))
            {
              #set.seed(9)
              mygraph <- as.matrix(seurat_object@graphs$RNA_snn)
              graphOut <- graph_from_adjacency_matrix(mygraph, mode = "undirected", weighted = T, diag = F)
              graphSimple <- simplify(graphOut) #, remove.loops=T)
              weights <- E(graphSimple)$weight
              sub_nodes <- V(graphSimple)$name
              
              tableCl <- seurat_object@meta.data[, ]
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
        session$sendCustomMessage("handler_startLoader", c("clust1_loader", 80))
        updateSelInpColor()
        updateInputLRclusters()
        updateInpuTrajectoryClusters()
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " ADDITIONAL DIMENSIONALITY\nREDUCTION METHODS", " TRAJECTORY ANALYSIS",
                                                          " MARKERS' IDENTIFICATION", " LIGAND - RECEPTOR\nANALYSIS"))
      }
      # }, warning = function(w) { # if this is not commented, the table does not render
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the Clustering procedure.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("clust1_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "clust1_loader")
      session$sendCustomMessage("handler_enableButton", "snnConfirm")
    })
  })
  
  observeEvent(input$clusterBarplotConfirm, {
    session$sendCustomMessage("handler_startLoader", c("clust2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "clusterBarplotConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first.")
      else {
        if(input$clusterGroupBy == "seurat_clusters")
        {
          #barplot for cell distribution per cluster
          clusterTable <- as.data.frame(table(Idents(seurat_object)))
          totalCells <- sum(clusterTable$Freq)
          
          #clusterTable$Perc <- (clusterTable$Freq*100)/totalCells
          clusterTable$Perc <- (clusterTable$Freq)/totalCells
          colnames(clusterTable)[1] <- "Cluster"
          colnames(clusterTable)[2] <-  "Cells"
          colnames(clusterTable)[3] <- "Percentage"
          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(clusterTable$Cluster)))
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 25))
          
          p <- ggplot(clusterTable) + theme_bw() +
            geom_bar( mapping = aes(x = Cluster, y = Percentage, fill=Cluster), stat = "identity" ) +
            scale_fill_manual(values = cols ) +
            theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
                  axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
                  axis.title.y = element_text(face = "bold", color = "black", size = 12),
                  axis.title.x = element_text(face = "bold", color = "black", size = 12),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) +
            labs(x="", y="Percent of cells", fill="Clusters")
          gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 50))
          
          output$clusterBarplot <- renderPlotly({print(gp)})
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 75))
        }
        else
        {
          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object$seurat_clusters)))
          
          p <- dittoBarPlot(seurat_object, var = "seurat_clusters",group.by = input$clusterGroupBy, scale = "percent") + theme_bw() +
            scale_fill_manual(values = cols ) +
            theme(axis.text.x = element_text(face = "bold", color = "black", size = 12, angle = 0),
                  axis.text.y = element_text(face = "bold", color = "black", size = 12, angle = 0),
                  axis.title.y = element_text(face = "bold", color = "black", size = 12),
                  axis.title.x = element_text(face = "bold", color = "black", size = 12),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) #+
          #labs(x="", y="Percentage % of cells", fill="Clusters")
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 25))
          
          gp <- plotly::ggplotly(p, tooltip = c("x", "y"))
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 50))
          
          output$clusterBarplot <- renderPlotly({print(gp)})
          session$sendCustomMessage("handler_startLoader", c("clust2_loader", 75))
        }
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the Clustering procedure.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("clust2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "clust2_loader")
      session$sendCustomMessage("handler_enableButton", "clusterBarplotConfirm")
    })
  })
  
  #------------------Umap/tSNE/DFM tab---------------------------------------
  observeEvent(input$umapRunUmap, {
    session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 10))
    session$sendCustomMessage("handler_disableButton", "umapRunUmap")
    session$sendCustomMessage("handler_disableButton", "umapRunTsne")
    session$sendCustomMessage("handler_disableButton", "umapRunDFM")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first execute PRINCIPAL COMPONENT ANALYSIS.")
      else {
        session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 25))
        seurat_object <<- RunUMAP(seurat_object, dims = 1:as.numeric(input$umapPCs), seed.use = 42, n.components = as.numeric(input$umapOutComponents), reduction = "pca") #TODO add diffusion map, addition of extra dimensions UMAP, select dimensions to plot, alpha and dot size
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " TRAJECTORY ANALYSIS"))
        updateUmapTypeChoices("umap")
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the UMAP procedure.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "dim_red1_loader")
      session$sendCustomMessage("handler_enableButton", "umapRunUmap")
      session$sendCustomMessage("handler_enableButton", "umapRunTsne")
      session$sendCustomMessage("handler_enableButton", "umapRunDFM")
    })
  })
  
  observeEvent(input$umapRunTsne, {
    session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 10))
    session$sendCustomMessage("handler_disableButton", "umapRunUmap")
    session$sendCustomMessage("handler_disableButton", "umapRunTsne")
    session$sendCustomMessage("handler_disableButton", "umapRunDFM")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first execute PRINCIPAL COMPONENT ANALYSIS.")
      else {
        session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 25))
        seurat_object <<- RunTSNE(seurat_object, dims = 1:as.numeric(input$umapPCs), seed.use = 42, dim.embed = 3, reduction = "pca", verbose = T)
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " TRAJECTORY ANALYSIS"))
        updateUmapTypeChoices("tsne")
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the tSNE procedure.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "dim_red1_loader")
      session$sendCustomMessage("handler_enableButton", "umapRunUmap")
      session$sendCustomMessage("handler_enableButton", "umapRunTsne")
      session$sendCustomMessage("handler_enableButton", "umapRunDFM")
    })
  })
  
  observeEvent(input$umapRunDFM, {
    session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 10))
    session$sendCustomMessage("handler_disableButton", "umapRunUmap")
    session$sendCustomMessage("handler_disableButton", "umapRunTsne")
    session$sendCustomMessage("handler_disableButton", "umapRunDFM")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first execute PRINCIPAL COMPONENT ANALYSIS.")
      else {
        session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 25))
        #prepare input
        dfm_in <- as.data.frame(seurat_object@reductions[['pca']]@cell.embeddings)
        dfm_in$Cell_id <- rownames(dfm_in)
        dfm_in <- dfm_in[, c(51, 1:as.numeric(input$umapPCs))]
        #print(head(dfm_in))
    
        #run DFM default
        dfm <- DiffusionMap(dfm_in, n_eigs = as.numeric(input$umapOutComponents))
        dfm_out <- dfm@eigenvectors
        colnames(dfm_out) <- gsub("DC", "DC_", colnames(dfm_out))
        #print(head(dfm_out))
    
        #add new reduction in seurat_object
        seurat_object[["dfm"]] <<- CreateDimReducObject(embeddings = dfm_out, key = "DC_", assay = DefaultAssay(seurat_object), global = T)
        # print("finished DFM")
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " TRAJECTORY ANALYSIS"))
        updateUmapTypeChoices("dfm")
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the Diffusion Map procedure.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("dim_red1_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "dim_red1_loader")
      session$sendCustomMessage("handler_enableButton", "umapRunUmap")
      session$sendCustomMessage("handler_enableButton", "umapRunTsne")
      session$sendCustomMessage("handler_enableButton", "umapRunDFM")
    })
  })
  
  # observeEvent(input$umapType, {
  #   if(input$umapType != "-")
  #     updateReduction()
  # })
  # 
  # observeEvent(input$umapColorBy, { 
  #   if(input$umapType != "-")
  #   updateReduction()
  # })
  # 
  # observeEvent(input$umapDimensions, { 
  #   if(input$umapType != "-")
  #   updateReduction()
  # })
  # 
  # observeEvent(input$umapDotSize, { 
  #   if(input$umapType != "-")
  #   updateReduction()
  # })
  # 
  # observeEvent(input$umapDotOpacity, { 
  #   if(input$umapType != "-")
  #   updateReduction()
  # })
  # 
  # observeEvent(input$umapDotBorder, { 
  #   if(input$umapType != "-")
  #   updateReduction()
  # })
  
  observeEvent(input$umapConfirm, { 
    if(input$umapType != "-")
      updateReduction()
  })
  
  #------------------DEA tab-----------------------------------------------
  observeEvent(input$findMarkersConfirm, { #TODO selectinput gia clustering column + check for errors in extra tabs with readRDS("seurat_processed.RDS")
    session$sendCustomMessage("handler_startLoader", c("DEA1_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("DEA2_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("DEA3_loader", 10))
    #session$sendCustomMessage("handler_startLoader", c("DEA4_loader", 10))
    #session$sendCustomMessage("handler_startLoader", c("DEA5_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("DEA6_loader", 10))
    session$sendCustomMessage("handler_disableButton", "umapConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first and then re-run UMAP or tSNE or Diffusion Map above.")
      else {
        markers_logFCBase <<- input$findMarkersLogBase
        
        if(markers_logFCBase == "avg_logFC")
        {
          seurat_object@misc$markers <<- FindAllMarkers(seurat_object, test.use = input$findMarkersTest, min.pct = as.numeric(input$findMarkersMinPct), logfc.threshold = as.numeric(input$findMarkersLogFC), 
                                                        return.thresh = as.numeric(input$findMarkersPval), base = exp(1))
        }
        else
        {
          seurat_object@misc$markers <<- FindAllMarkers(seurat_object, test.use = input$findMarkersTest, min.pct = as.numeric(input$findMarkersMinPct), logfc.threshold = as.numeric(input$findMarkersLogFC), 
                                                        return.thresh = as.numeric(input$findMarkersPval), base = 2)
        }
        updateInputGeneList()
        
        output$findMarkersTable <- renderDataTable(
          {
            if (!is.null(seurat_object@misc$markers))
            {
              output$findMarkersTable <- renderDataTable(seurat_object@misc$markers, options = list(pageLength = 20), filter = 'top', rownames = FALSE)  
            }
          }
        )
        session$sendCustomMessage("handler_startLoader", c("DEA1_loader", 75))
        
        #DEA output rendering
        output$findMarkersHeatmap <- renderPlotly(
          {
            top10 <- seurat_object@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = eval(parse(text=markers_logFCBase))) #wt = avg_logFC)
            set.seed(9)
            downsampled <- subset(seurat_object, cells = sample(Cells(seurat_object), 1500))
            
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
        session$sendCustomMessage("handler_startLoader", c("DEA2_loader", 75))
        
        output$findMarkersDotplot <- renderPlotly(
          {
            top10 <- seurat_object@misc$markers %>% group_by(cluster) %>% top_n(n = 10, wt = eval(parse(text=markers_logFCBase)))#wt = avg_logFC)
            p <- DotPlot(seurat_object, features = rev(unique(top10$gene)), dot.scale = 6, cols = c("grey", "red")) + RotatedAxis() + labs(fill="Average\nexpression")
            plotly::ggplotly(p)
          }
        )
        
        #session$sendCustomMessage("handler_startLoader", c("DEA3_loader", 75))
        
        #session$sendCustomMessage("handler_startLoader", c("DEA4_loader", 75))
        
        #session$sendCustomMessage("handler_startLoader", c("DEA5_loader", 75))
        
        output$findMarkersVolcanoPlot <- renderPlotly(
          {
            diff_exp_genes <- seurat_object@misc$markers
            cluster_degs <- diff_exp_genes[which(diff_exp_genes$cluster == input$findMarkersClusterSelect), ]
            cluster_degs$status <- "Down regulated"
            for(i in 1:length(cluster_degs$gene))
            {
              if(cluster_degs[i, markers_logFCBase] > 0) #(cluster_degs$avg_logFC[i] > 0)
              {
                cluster_degs$status[i] <- "Up regulated"
              }
            }
            cluster_degs$log10Pval <- -log10(cluster_degs$p_val)
            
            p <- ggplot(data=cluster_degs, aes_string(x=markers_logFCBase, y="log10Pval", fill="status", label="gene", color="status")) + #x="avg_logFC"
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
              labs(x=paste0("", markers_logFCBase), y="-log10(Pvalue)", fill="Color", color="")
            plotly::ggplotly(p, tooltip = c("x", "y", "label"))
          }
        )
        session$sendCustomMessage("handler_startLoader", c("DEA6_loader", 75))
        session$sendCustomMessage("handler_enableTabs", c("sidebarMenu", " FUNCTIONAL ENRICHMENT\nANALYSIS", " CLUSTERS' ANNOTATION"))
      }
      # }, warning = function(w) {
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the DE Analysis")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("DEA1_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("DEA2_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("DEA3_loader", 100))
      #session$sendCustomMessage("handler_startLoader", c("DEA4_loader", 100))
      #session$sendCustomMessage("handler_startLoader", c("DEA5_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("DEA6_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", c("DEA1_loader"))
      session$sendCustomMessage("handler_finishLoader", c("DEA2_loader"))
      session$sendCustomMessage("handler_finishLoader", c("DEA3_loader"))
      #session$sendCustomMessage("handler_finishLoader", c("DEA4_loader"))
      #session$sendCustomMessage("handler_finishLoader", c("DEA5_loader"))
      session$sendCustomMessage("handler_finishLoader", c("DEA6_loader"))
      session$sendCustomMessage("handler_enableButton", "umapConfirm")
    })
  })

  observeEvent(input$findMarkersSignatureAdd, {
    markers <- list()
    varTextarea <- input$findMarkersSignatureMembers
    markers[[1]] <- unlist(strsplit(varTextarea, "\\n")) #c("Prg4", "Tspan14", "Clic5", "Htra4")
    sig_name <- input$findMarkersSignatureName
    print(sig_name)
    names(markers)[1] <- sig_name
    print(markers)
    
    seurat_object <<- AddModuleScore_UCell(seurat_object, features = markers)
    updateSignatures()
  })
    
observeEvent(input$findMarkersFPConfirm, {
  if(input$findMarkersReductionType != "-")
  {
    output$findMarkersFeaturePlot <- renderPlotly(
      {
        geneS <- ""
        if(input$findMarkersFeatureSignature == "signature")
        {
          geneS <- input$findMarkersSignatureSelect
        }
        else
        {
          geneS <- input$findMarkersGeneSelect
        }
        label_x <- ""
        label_y <- ""
        show_label <- as.logical(input$findMarkersLabels)
        order_exp <- as.logical(input$findMarkersOrder)
        minq <- paste0("q", input$findMarkersMinCutoff)
        maxq <- paste0("q", input$findMarkersMaxCutoff)
        
        if(input$findMarkersReductionType == "umap")
        {
          label_x <- "UMAP_1"
          label_y <- "UMAP_2"
        }
        else if(input$findMarkersReductionType == "tsne")
        {
          label_x <- "tSNE_1"
          label_y <- "tSNE_2"
        }
        else if(input$findMarkersReductionType == "dfm")
        {
          label_x <- "DC_1"
          label_y <- "DC_2"
        }
        else if(input$findMarkersReductionType == "pca")
        {
          label_x <- "PC_1"
          label_y <- "PC_2"
        }
        
        plot <- FeaturePlot(seurat_object, features = geneS, pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red"), 
                            order = order_exp, reduction = input$findMarkersReductionType, max.cutoff = maxq, min.cutoff = minq) +
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.title.y = element_text(face = "bold", color = "black", size = 25),
                axis.title.x = element_text(face = "bold", color = "black", size = 25),
                legend.text = element_text(face = "bold", color = "black", size = 9),
                legend.title = element_text(face = "bold", color = "black", size = 9),
                legend.position="right",
                title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
          labs(x=label_x, y=label_y, title = geneS, color="")
        gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
        gp
      }
    )
  }
})

observeEvent(input$findMarkersFeaturePairConfirm, {
   updateFeaturePair()
})

# observeEvent(input$findMarkersBlendThreshold, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePairReductionType, {
#     updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePairMaxCutoff, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePairMinCutoff, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePair1, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePair2, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePairLabels, {
#   updateFeaturePair()
# })
# 
# observeEvent(input$findMarkersFeaturePairOrder, {
#   updateFeaturePair()
# })

observeEvent(input$findMarkersViolinConfirm, {
  if(!identical(seurat_object, NULL)){
    if ((input$findMarkersViolinFeaturesSignature == "gene" & !identical(input$findMarkersGeneSelect2, NULL))
        | (input$findMarkersViolinFeaturesSignature == "signature" & input$findMarkersViolinSignatureSelect != "-")){
      output$findMarkersViolinPlot <- renderPlotly(
        {
          geneS <- ""
          if(input$findMarkersViolinFeaturesSignature == "gene")
          {
            geneS <- input$findMarkersGeneSelect2
          }
          else
          {
            geneS <- input$findMarkersViolinSignatureSelect
          }
          print(geneS)
          plot <- VlnPlot(seurat_object, features = geneS, pt.size = 0,
                          cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'seurat_clusters'])))) +
            theme_bw() +
            theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                  axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                  axis.title.y = element_text(face = "bold", color = "black", size = 25),
                  axis.title.x = element_text(face = "bold", color = "black", size = 25),
                  legend.text = element_text(face = "bold", color = "black", size = 9),
                  legend.title = element_text(face = "bold", color = "black", size = 9),
                  #legend.position="right",
                  title = element_text(face = "bold", color = "black", size = 25, angle = 0)) +
            labs(x="Cluster", y="", title = geneS, fill="Cluster")
          gp <- plotly::ggplotly(plot, tooltip = c("x", "y", geneS))
          gp
        }
      )
    } else session$sendCustomMessage("handler_alert", "Please select signature.")
  } else session$sendCustomMessage("handler_alert", "Please upload some data first.")
})
  
  #------------------Cell cycle tab------------------------------------------
  observeEvent(input$cellCycleRun, { # observe selectInput cellCycleReduction instead of cellCycleRun actionButton
    session$sendCustomMessage("handler_startLoader", c("CC1_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("CC2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "cellCycleRun")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (!"pca" %in% names(seurat_object)) session$sendCustomMessage("handler_alert", "Please, first execute PRINCIPAL COMPONENT ANALYSIS.") # TODO add if conditions for UMAP, tSNE etc
      else {
        #******cell cycle analysis #TODO before and after regression plot
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
        
        output$cellCyclePCA <- renderPlotly(
          {
            if(!is.null(seurat_object) & input$cellCycleReduction != "-")
            {
              meta <- seurat_object@meta.data
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
              else if(selected_Reduction == "dfm")
              {
                label_x <- "DC_1"
                label_y <- "DC_2"
              }
              
              plot1 <- DimPlot(seurat_object, reduction = selected_Reduction)
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
        session$sendCustomMessage("handler_startLoader", c("CC1_loader", 80))
        
        output$cellCycleBarplot <- renderPlotly(
          {
            #barplot
            fsize <- 24
            
            obj_meta <- seurat_object@meta.data
            obj_meta$cell_id <- row.names(obj_meta)
            obj_meta$Cluster_name <- Idents(seurat_object)
            obj_meta <- obj_meta[, c('cell_id', 'Cluster_name', 'Phase')]
            final_df_seurat <- obj_meta %>% dplyr::group_by(Cluster_name) %>% dplyr::count(Phase)
            colnames(final_df_seurat)[3] <- "Cells"
            
            total_cells_per_cluster <- as.data.frame(table(Idents(seurat_object)))
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
        session$sendCustomMessage("handler_startLoader", c("CC2_loader", 80))
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with drawing the resutls.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("CC1_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("CC2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "CC1_loader")
      session$sendCustomMessage("handler_finishLoader", "CC2_loader")
      session$sendCustomMessage("handler_enableButton", "cellCycleRun")
    })
  })
  
  #------------------gProfiler tab-----------------------------------------------
  observeEvent(input$gProfilerConfirm, {
    session$sendCustomMessage("handler_startLoader", c("gprof1_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("gprof2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "gProfilerConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@misc$markers, NULL)) session$sendCustomMessage("handler_alert", "Please, execute the gene differential analysis at the MARKERS' IDENTIFICATION tab first.")
      else {
        cluster_temp <- input$gProfilerList
        temp_df <- data.frame()
        # if(cluster_temp == "all_clusters")
        # {
        #   #all clusters used
        #   all_clusters <- as.character(unique(seurat_object@misc$markers$cluster))
        #   
        #   gene_lists <- list()
        #   for(i in 1:length(all_clusters))
        #   {
        #     if(input$gProfilerLFCRadio == "Up")#UP regulated
        #     {
        #       gene_lists[[i]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == all_clusters[i] & 
        #                                                             seurat_object@misc$markers[, markers_logFCBase] >= as.numeric(input$gProfilerSliderLogFC) &
        #                                                             seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
        #     }
        #     else #down
        #     {
        #       gene_lists[[i]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == all_clusters[i] & 
        #                                                             seurat_object@misc$markers[, markers_logFCBase] <= (as.numeric(input$gProfilerSliderLogFC)*(-1)) &
        #                                                             seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
        #     }
        #   }
        #   
        #   names(gene_lists) <- all_clusters
        #   session$sendCustomMessage("handler_startLoader", c("gprof1_loader", 60))
        #   
        #   # gostres <- gost(query = gene_lists, 
        #   #                 organism = "mmusculus", ordered_query = FALSE, 
        #   #                 multi_query = F, significant = TRUE, exclude_iea = T, 
        #   #                 measure_underrepresentation = FALSE, evcodes = TRUE, 
        #   #                 user_threshold = 0.05, correction_method = "g_SCS", 
        #   #                 domain_scope = "annotated", custom_bg = NULL, 
        #   #                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
        #   
        #   gostres <- gost(query = gene_lists, 
        #                   organism = input$gProfilerOrganism, ordered_query = FALSE, 
        #                   multi_query = F, significant = TRUE, exclude_iea = F, 
        #                   measure_underrepresentation = FALSE, evcodes = TRUE, 
        #                   user_threshold = as.numeric(input$gProfilerSliderSignificanceTerms), 
        #                   correction_method = input$gprofilerRadioCorrection, 
        #                   domain_scope = "annotated", custom_bg = NULL, 
        #                   numeric_ns = "", sources = input$gProfilerDatasources, as_short_link = FALSE)
        #   
        #   temp_df <- gostres$result
        #   
        #   # output$gProfilerManhatan <- renderPlotly({ gostplot(gostres, capped = T, interactive = T) %>% 
        #   # layout(width=1500, height=length(all_clusters)*400) })
        #   if (!is.na(all_clusters)){ # TODO check if null or sth else
        #     session$sendCustomMessage("handler_fixHeight", c("gProfilerManhatan", length(all_clusters)))
        #     output$gProfilerManhatan <- renderPlotly({ print(gostplot(gostres, capped = TRUE, interactive = TRUE)) })
        #     session$sendCustomMessage("handler_startLoader", c("gprof2_loader", 80))
        #   } else session$sendCustomMessage("handler_alert", "Please, first calculate clusters") # TODO desicde sentence
        #   #output$gProfilerManhatan <- renderPlotly({ plotly::ggplotly(gostplot(gostres, capped = TRUE, interactive = TRUE)) })
        #   #output$gProfilerManhatan <- renderPlot({gostplot(gostres, capped = TRUE, interactive = FALSE)}, height = 4000)
        # }
        # else
        # {
          gene_lists <- list()
          
          if(input$gProfilerLFCRadio == "Up")#UP regulated
          {
            gene_lists[[1]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == cluster_temp & 
                                                                  #seurat_object@misc$markers$avg_logFC >= as.numeric(input$gProfilerSliderLogFC) & 
                                                                  seurat_object@misc$markers[, markers_logFCBase] >= as.numeric(input$gProfilerSliderLogFC) &
                                                                  seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
          }
          else #down
          {
            gene_lists[[1]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == cluster_temp & 
                                                                  seurat_object@misc$markers[, markers_logFCBase] <= (as.numeric(input$gProfilerSliderLogFC)*(-1)) &
                                                                  seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
          }
          session$sendCustomMessage("handler_startLoader", c("gprof1_loader", 60))
          
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
          session$sendCustomMessage("handler_startLoader", c("gprof2_loader", 80))
        
        session$sendCustomMessage("handler_startLoader", c("gprof1_loader", 80))
        
        output$gProfilerTable <- renderDataTable(temp_df[, c(1, 3:6, 9:11, 16)], options = list(pageLength = 10), rownames = F)
      }
      # }, warning = function(w) {
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with the Enrichment Analysis.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("gprof1_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("gprof2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "gprof1_loader")
      session$sendCustomMessage("handler_finishLoader", "gprof2_loader")
      session$sendCustomMessage("handler_enableButton", "gProfilerConfirm")
    })
  })

observeEvent(input$sendToFlame, {
  tryCatch({
    if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
    else if (identical(seurat_object@misc$markers, NULL)) session$sendCustomMessage("handler_alert", "Please, execute the gene differential analysis at the MARKERS' IDENTIFICATION tab first.")
    else {
      cluster_temp <- input$gProfilerList
      temp_df <- data.frame()
      gene_lists <- list()
      
      if(input$gProfilerLFCRadio == "Up")#UP regulated
      {
        gene_lists[[1]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == cluster_temp & 
                                                              #seurat_object@misc$markers$avg_logFC >= as.numeric(input$gProfilerSliderLogFC) & 
                                                              seurat_object@misc$markers[, markers_logFCBase] >= as.numeric(input$gProfilerSliderLogFC) &
                                                              seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
      }
      else #down
      {
        gene_lists[[1]] <- seurat_object@misc$markers[which(seurat_object@misc$markers$cluster == cluster_temp & 
                                                              seurat_object@misc$markers[, markers_logFCBase] <= (as.numeric(input$gProfilerSliderLogFC)*(-1)) &
                                                              seurat_object@misc$markers[, input$gprofilerRadio] < as.numeric(input$gProfilerSliderSignificance)), 'gene']
      }
      js$Enrich(paste0("http://bib.fleming.gr:3838/Flame/?url_genes=", paste(gene_lists[[1]], collapse = ",")))
    }
  }, error = function(e) {
    print(paste("Error :  ", e))
    session$sendCustomMessage("handler_alert", "There was an error with the Enrichment Analysis.")
  })
})
  
  #------------------CIPR tab-----------------------------------------------
  observeEvent(input$annotateClustersConfirm, {
    session$sendCustomMessage("handler_startLoader", c("annot1_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("annot2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "annotateClustersConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@misc$markers, NULL)) session$sendCustomMessage("handler_alert", "Please, execute the gene differential analysis at the MARKERS' IDENTIFICATION tab first.")
      else {
        marker_genes <- seurat_object@misc$markers
        
        avgexp <- AverageExpression(seurat_object)
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
        CIPR_top_results$index <- factor(CIPR_top_results$index, levels = as.character(seq(1:length(seurat_object$seurat_clusters))))#1:45
        cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(CIPR_top_results$cluster)))
        session$sendCustomMessage("handler_startLoader", c("annot1_loader", 50))
        session$sendCustomMessage("handler_startLoader", c("annot2_loader", 50))
        
        p <- ggplot(CIPR_top_results, aes(x=index, y=identity_score, fill = cluster, size=7)) +
          theme_bw() +
          geom_point(alpha=1, shape=21) +
          scale_size(range = c(7, 7), guide = F) + 
          scale_x_discrete(labels=CIPR_top_results$long_name) +
          scale_fill_manual(values = cols) +
          labs(x="") + 
          theme(legend.position="top",
                axis.text.x = element_text(vjust=0.5, hjust=1, angle = 90))
        
        session$sendCustomMessage("handler_startLoader", c("annot1_loader", 80))
        output$annotateClustersCIPRTable <- renderDataTable(CIPR_top_results[], options = list(pageLength = 20), rownames = F) #remove Description
        
        session$sendCustomMessage("handler_startLoader", c("annot2_loader", 80))
        output$annotateClustersCIPRDotplot <- renderPlotly({ print(p)})
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with Cluster Annotation.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("annot1_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("annot2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "annot1_loader")
      session$sendCustomMessage("handler_finishLoader", "annot2_loader")
      session$sendCustomMessage("handler_enableButton", "annotateClustersConfirm")
    })
  })
  
  #--------------------trajectory tab----------------------------------------
  observeEvent(input$trajectoryConfirm, {

    session$sendCustomMessage("handler_startLoader", c("traj1_loader", 10))
    session$sendCustomMessage("handler_startLoader", c("traj2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "trajectoryConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first and then run UMAP, in the NON LINEAR DIMENSIONALITY REDUCTION tab.")
      else if (identical(seurat_object$umap, NULL)) session$sendCustomMessage("handler_alert", "Please, calculate UMAP first, in the NON LINEAR DIMENSIONALITY REDUCTION tab.")
      else {
        #delete previous lineages columns
        for_delete <- grep("Lineage", colnames(seurat_object@meta.data))
        if(length(for_delete) != 0)
        {
          seurat_object@meta.data <<- seurat_object@meta.data[, -for_delete]
        }
        
        reduction <- Embeddings(seurat_object, input$trajectoryReduction)
        
        sds <- slingshot(reduction[, 1:as.numeric(input$trajectorySliderDimensions)], clusterLabels = seurat_object$seurat_clusters, 
                         start.clus = as.numeric(input$trajectoryStart), end.clus = as.numeric(input$trajectoryEnd), stretch = 0)
        
        metaD$all_lin <- slingLineages(sds)
        pt <- slingPseudotime(sds)
        
        for(i in 1:ncol(pt))
        {
          temp_lin <- pt[, i]
          seurat_object <<- AddMetaData(
            object = seurat_object,
            metadata = temp_lin,
            col.name = paste0("Lineage", i)
          )
        }
        
        updateInputLineageList(names(metaD$all_lin))
        session$sendCustomMessage("handler_startLoader", c("traj1_loader", 50))
        session$sendCustomMessage("handler_startLoader", c("traj2_loader", 50))
        #print(paste0("after update:", names(metaD$all_lin)))
        
        plot_D <- dittoDimPlot(seurat_object, "seurat_clusters",
                               do.label = TRUE,
                               labels.repel = TRUE, 
                               reduction.use = "umap", # TODO, why always umap here? can't this be pca?
                               add.trajectory.lineages = slingLineages(sds),
                               trajectory.cluster.meta = "seurat_clusters",
                               color.panel = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(seurat_object@meta.data[, 'seurat_clusters']))),
                               data.out = F)
        
        session$sendCustomMessage("handler_startLoader", c("traj1_loader", 80))
        output$trajectoryPlot <- renderPlot({ print(plot_D)})
        #***delete lineage columns before second run
        plot_P <- dittoDimPlot(seurat_object, var = "Lineage1",  
                               do.label = F,
                               labels.repel = F, 
                               reduction.use = "umap",
                               add.trajectory.lineages = list((metaD$all_lin[["Lineage1"]])),
                               trajectory.cluster.meta = "seurat_clusters",
                               size = 2,
                               data.out = F) + scale_colour_gradientn(colours = plasma(100), na.value = "grey90") +
          labs(colour = "Pseudotime") #new column paste(cluster, annotation)
        
        session$sendCustomMessage("handler_startLoader", c("traj2_loader", 80))
        output$trajectoryPseudotimePlot <- renderPlot({ print(plot_P) })
        output$trajectoryText <- renderPrint({print(metaD$all_lin)})
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with Trajectory Analysis.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("traj1_loader", 100))
      session$sendCustomMessage("handler_startLoader", c("traj2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "traj1_loader")
      session$sendCustomMessage("handler_finishLoader", "traj2_loader")
      session$sendCustomMessage("handler_enableButton", "trajectoryConfirm")
    })
  })
  
  observeEvent(input$trajectoryConfirmLineage, {
    session$sendCustomMessage("handler_startLoader", c("traj2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "trajectoryConfirmLineage")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else {
        #Lineage view
        plot_P <- dittoDimPlot(seurat_object, var = input$trajectoryLineageSelect,
                               do.label = F,
                               labels.repel = F,
                               reduction.use = "umap",
                               add.trajectory.lineages = list(metaD$all_lin[[input$trajectoryLineageSelect]]),
                               trajectory.cluster.meta = "seurat_clusters",
                               size = 2,
                               rename.var.groups = "Pseudotime",
                               data.out = F) + scale_colour_gradientn(colours = plasma(100), na.value = "grey90") +
          labs(colour = "Pseudotime")
        session$sendCustomMessage("handler_startLoader", c("traj2_loader", 80))
        
        output$trajectoryPseudotimePlot <- renderPlot({ print(plot_P) })
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with viewing Trajectory Lineage.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("traj2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "traj2_loader")
      session$sendCustomMessage("handler_enableButton", "trajectoryConfirmLineage")
    })
  })
  
  #--------------------Ligand Receptor tab---------------------------
  observeEvent(input$ligandReceptorConfirm, {
    session$sendCustomMessage("handler_startLoader", c("lr_loader", 10))
    session$sendCustomMessage("handler_disableButton", "ligandReceptorConfirm")
    tryCatch({
      # if (!"assays" %in% slotNames(seurat_object)) session$sendCustomMessage("handler_alert", "Slot assays has not been calculated yet. Execute Normalization procedure first.")
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first.")
      else {
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
        
        session$sendCustomMessage("handler_startLoader", c("lr_loader", 40))
        #expressed genes
        ## receiver
        receiver = input$ligandReceptorSender
        
        expressed_genes_receiver = get_expressed_genes(receiver, seurat_object, pct = 0.10, assay_oi = "RNA")
        ## sender
        sender = input$ligandReceptorSender
        expressed_genes_sender = get_expressed_genes(sender, seurat_object, pct = 0.10, assay_oi = "RNA")
        
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
        
        session$sendCustomMessage("handler_startLoader", c("lr_loader", 60))
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
        # p_ligand_receptor_network_strict
        session$sendCustomMessage("handler_startLoader", c("lr_loader", 80))
        
        output$ligandReceptorCuratedHeatmap <- renderPlotly({ plotly::ggplotly(p_ligand_receptor_network_strict) })
      }
    # }, warning = function(w) {
    #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with Ligand-Receptor Analysis.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("lr_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "lr_loader")
      session$sendCustomMessage("handler_enableButton", "ligandReceptorConfirm")
    })
  })
  
  #updates on UI elements ----------
  updateReduction <- function()
  {
    session$sendCustomMessage("handler_startLoader", c("dim_red2_loader", 10))
    session$sendCustomMessage("handler_disableButton", "umapConfirm")
    tryCatch({
      if (identical(seurat_object, NULL)) session$sendCustomMessage("handler_alert", "Please, upload some data via the DATA INPUT tab first.")
      else if (identical(seurat_object@meta.data$seurat_clusters, NULL)) session$sendCustomMessage("handler_alert", "Please, execute CLUSTERING first and then re-run UMAP or tSNE or Diffusion Map above.")
      else {
        #get input
        dims <- as.numeric(input$umapDimensions)
        type <- input$umapType
        
        #prepare metadata
        meta <- seurat_object@meta.data
        meta$Cell_id <- rownames(meta)
        meta <- meta[, ]#meta[, c('Cell_id', 'seurat_clusters', 'orig.ident')]
        reduc_data <- data.frame()
        
        #prepare colors
        cols = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy])))
        
        #for all reductions
        seurat_object_reduc <- as.data.frame(seurat_object@reductions[[input$umapType]]@cell.embeddings)
        seurat_object_reduc <- seurat_object_reduc[, c(1:ncol(seurat_object_reduc))]
        seurat_object_reduc$Cell_id <- rownames(seurat_object_reduc)
        reduc_data <- left_join(seurat_object_reduc, meta)
        
        session$sendCustomMessage("handler_startLoader", c("dim_red2_loader", 50))
        
        if(type == "umap" & dims == 2)
        {
          p <- ggplot(data=reduc_data, aes_string(x="UMAP_1", y="UMAP_2", fill=input$umapColorBy)) +
            geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder))+
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
          output$umapPlot <- renderPlotly({plotly::ggplotly(p)})
        }
        else if(type == "umap" & dims == 3)
        {
          p <- plot_ly(reduc_data, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, type="scatter3d", alpha = as.numeric(input$umapDotOpacity), mode="markers", color=as.formula(paste0('~', input$umapColorBy)),
                       marker = list(size = as.numeric(input$umapDotSize), 
                                     line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                       ),
                       colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) )
          
          output$umapPlot <- renderPlotly({print(p)})
        }
        else if(type == "tsne" & dims == 2)
        {
          p <- ggplot(data=reduc_data, aes_string(x="tSNE_1", y="tSNE_2", fill=input$umapColorBy)) +
            geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
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
          output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
        }
        else if(type == "tsne" & dims == 3)
        {
          p <- plot_ly(reduc_data, x=~tSNE_1, y=~tSNE_2, z=~tSNE_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                       marker = list(size = as.numeric(input$umapDotSize), 
                                     line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                       ),
                       colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
          output$umapPlot <- renderPlotly({print(p)})
        }
        else if(type == "dfm" & dims == 2)
        {
          p <- ggplot(data=reduc_data, aes_string(x="DC_1", y="DC_2", fill=input$umapColorBy)) +
            geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
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
            labs(x="DC 1", y="DC 2", color="Cell type", title = "", fill="Color")
          output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
        }
        else if(type == "dfm" & dims == 3)
        {
          p <- plot_ly(reduc_data, x=~DC_1, y=~DC_2, z=~DC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                       marker = list(size = as.numeric(input$umapDotSize), 
                                     line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                       ),
                       colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
          output$umapPlot <- renderPlotly({print(p)})
        }
        else if(type == "pca" & dims == 2)
        {
          p <- ggplot(data=reduc_data, aes_string(x="PC_1", y="PC_2", fill=input$umapColorBy)) +
            geom_point(size= as.numeric(input$umapDotSize), shape=21, alpha= as.numeric(input$umapDotOpacity), stroke=as.numeric(input$umapDotBorder)) +
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
            labs(x="PC 1", y="PC 2", color="Cell type", title = "", fill="Color")
          output$umapPlot <- renderPlotly({plotly::ggplotly(p)})  
        }
        else if(type == "pca" & dims == 3)
        {
          p <- plot_ly(reduc_data, x=~PC_1, y=~PC_2, z=~PC_3, type="scatter3d", mode="markers", alpha = as.numeric(input$umapDotOpacity), color=as.formula(paste0('~', input$umapColorBy)), 
                       marker = list(size = as.numeric(input$umapDotSize), 
                                     line = list(color = 'black', width = as.numeric(input$umapDotBorder))
                       ),
                       colors = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(meta[, input$umapColorBy]))) ) 
          output$umapPlot <- renderPlotly({print(p)})
        }
      }
      # }, warning = function(w) {
      #   print(paste("Warning:  ", w))
    }, error = function(e) {
      print(paste("Error :  ", e))
      session$sendCustomMessage("handler_alert", "There was an error with drawing the resutls.")
    }, finally = {
      session$sendCustomMessage("handler_startLoader", c("dim_red2_loader", 100))
      Sys.sleep(1)
      session$sendCustomMessage("handler_finishLoader", "dim_red2_loader")
      session$sendCustomMessage("handler_enableButton", "umapConfirm")
    })
  }
  
  updateFeaturePair <- function()
  {
      if (!identical(seurat_object, NULL) & input$findMarkersFeaturePairReductionType != "-" & input$findMarkersFeaturePair1 %in% rownames(seurat_object) & input$findMarkersFeaturePair2 %in% rownames(seurat_object))
      {
        geneS1 <- input$findMarkersFeaturePair1
        geneS2 <- input$findMarkersFeaturePair2
        label_x <- ""
        label_y <- ""
        show_label <- as.logical(input$findMarkersFeaturePairLabels)
        order_exp <- as.logical(input$findMarkersFeaturePairOrder)
        minq <- paste0("q", input$findMarkersFeaturePairMinCutoff)
        maxq <- paste0("q", input$findMarkersFeaturePairMaxCutoff)
        blendThr <- as.numeric(input$findMarkersBlendThreshold)
        
        if(input$findMarkersFeaturePairReductionType == "umap")
        {
          label_x <- "UMAP_1"
          label_y <- "UMAP_2"
        }
        else if(input$findMarkersFeaturePairReductionType == "tsne")
        {
          label_x <- "tSNE_1"
          label_y <- "tSNE_2"
        }
        else if(input$findMarkersFeaturePairReductionType == "dfm")
        {
          label_x <- "DC_1"
          label_y <- "DC_2"
        }
        else if(input$findMarkersFeaturePairReductionType == "pca")
        {
          label_x <- "PC_1"
          label_y <- "PC_2"
        }
        
        plot <- FeaturePlot(seurat_object, features = c(geneS1, geneS2), blend.threshold = blendThr, 
                            pt.size = 1.5, label = show_label, label.size = 5, cols = c("lightgrey", "red", "dodgerblue4"), 
                            order = order_exp, reduction = input$findMarkersFeaturePairReductionType, blend = TRUE, max.cutoff = maxq, min.cutoff = minq) +
          theme_bw() +
          theme(axis.text.x = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.text.y = element_text(face = "bold", color = "black", size = 25, angle = 0),
                axis.title.y = element_text(face = "bold", color = "black", size = 25),
                axis.title.x = element_text(face = "bold", color = "black", size = 25),
                legend.text = element_text(face = "bold", color = "black", size = 9),
                legend.title = element_text(face = "bold", color = "black", size = 9),
                legend.position="none",
                title = element_text(face = "bold", color = "black", size = 25, angle = 0)) #+
        labs(x=label_x, y=label_y, color="Normalized\nexpression")
        gp1 <- plotly::ggplotly(plot[[1]]) #, tooltip = c("x", "y", geneS))
        gp2 <- plotly::ggplotly(plot[[2]])
        gp3 <- plotly::ggplotly(plot[[3]])
        gp4 <- plotly::ggplotly(plot[[4]]+theme(title = element_text(face = "bold", color = "black", size = 15)))
        
        output$findMarkersFPfeature1 <- renderPlotly({ gp1 })
        output$findMarkersFPfeature2 <- renderPlotly({ gp2 })
        output$findMarkersFPfeature1_2 <- renderPlotly({ gp3 })
        output$findMarkersFPcolorbox <- renderPlotly({ gp4 })  
      }
  }
  
  updateQC_choices <- function()
  {
    updateSliderInput(session, "minUniqueGenes", min = min(init_seurat_object$nFeature_RNA), max = max(init_seurat_object$nFeature_RNA)-2)
    updateSliderInput(session, "maxUniqueGenes", min = min(init_seurat_object$nFeature_RNA)+2, max = max(init_seurat_object$nFeature_RNA))
  }
  
  
  
  updateUmapTypeChoices <- function(type)
  {
    reductions_choices <<- c(reductions_choices, type)
    updateSelectInput(session, "umapType", choices = reductions_choices)  #umapType, findMarkersReductionType, cellCycleReduction, trajectoryReduction
    updateSelectInput(session, "findMarkersReductionType", choices = reductions_choices)
    updateSelectInput(session, "findMarkersFeaturePairReductionType", choices = reductions_choices)
    updateSelectInput(session, "cellCycleReduction", choices = reductions_choices)
    updateSelectInput(session, "trajectoryReduction", choices = reductions_choices)
  }
  
  updateSignatures <- function()
  {
    #sig_names <- grep(pattern = "_UCell", x = colnames(seurat_object@meta.data))
    colnames(seurat_object@meta.data)
    tableMeta <- seurat_object@meta.data
    f <- sapply(tableMeta, is.numeric)
    sig_names <- colnames(tableMeta[, f])
    updateSelectInput(session, "findMarkersSignatureSelect", choices = sig_names)#colnames(seurat_object@meta.data)[sig_names])
    updateSelectInput(session, "findMarkersViolinSignatureSelect", choices = sig_names) #colnames(seurat_object@meta.data)[sig_names])
  }
  
  #function update selectInput
  updateSelInpColor <- function()
  {
    colnames(seurat_object@meta.data)
    tableMeta <- seurat_object@meta.data
    f <- sapply(tableMeta, is.factor)
    factors <- colnames(tableMeta[, f])
    # print(factors)
    updateSelectInput(session, "umapColorBy", choices = factors)
    updateSelectInput(session, "qcColorBy", choices = factors)
    updateSelectInput(session, "pcaColorBy", choices = factors)
    updateSelectInput(session, "clusterGroupBy", choices = factors)
  }
  
  updateRegressOut <- function()
  {
    variables <- colnames(seurat_object@meta.data)
    updateSelectInput(session, "normalizeRegressColumns", choices = variables)
  }
  
  updateGeneSearchFP <- function()
  {
    total_genes <- rownames(seurat_object@assays$RNA@counts)
    updateSelectizeInput(session, 'findMarkersGeneSelect', choices = total_genes, server = TRUE) # server-side selectize drastically improves performance
    updateSelectizeInput(session, 'findMarkersGeneSelect2', choices = total_genes, server = TRUE)
    updateSelectizeInput(session, 'findMarkersFeaturePair1', choices = total_genes, server = TRUE)
    updateSelectizeInput(session, 'findMarkersFeaturePair2', choices = total_genes, server = TRUE)
  }
  
  updateInputGeneList <- function()
  {
    all_cluster_names <- as.character(unique(seurat_object@misc$markers$cluster))
    updateSelectInput(session, "gProfilerList", choices = all_cluster_names)
  }
  
  updateInputLineageList <- function(lin_names)
  {
    updateSelectInput(session, "trajectoryLineageSelect", choices = lin_names)
  }
  
  updateInputLRclusters <- function()
  {
    all_cluster_names <- (levels(seurat_object@meta.data[, 'seurat_clusters']))
    updateSelectInput(session, "ligandReceptorSender", choices = all_cluster_names)
    updateSelectInput(session, "ligandReceptorReciever", choices = all_cluster_names)
    updateSelectInput(session, "utilitiesRenameOldName", choices = all_cluster_names)
    updateSelectInput(session, "utilitiesDeleteCluster", choices = all_cluster_names)
  }
  
  updateInpuTrajectoryClusters <- function()
  {
    all_cluster_names <- (levels(seurat_object@meta.data[, 'seurat_clusters']))
    updateSelectInput(session, "trajectoryStart", choices = all_cluster_names)
    updateSelectInput(session, "trajectoryEnd", choices = all_cluster_names)
  }
  
  # Helper Functions ####
  
  # This function is called after a new input file has been uploaded
  # and is responsible for clearing all generated plots across all tabs
  # @param fromDataInput: If TRUE, clears all, including QC, else skips clearing QC
  cleanAllPlots <- function(fromDataInput){
    # renderPlotly
    if (fromDataInput){
      output$nFeatureViolin <- NULL
      output$totalCountsViolin <- NULL
      output$mitoViolin <- NULL
      output$mtCounts <- NULL
      output$genesCounts <- NULL
      output$cellStats <- NULL
      output$filteredNFeatureViolin <- NULL
      output$filteredTotalCountsViolin <- NULL
      output$filteredMitoViolin <- NULL
      output$filteredMtCounts <- NULL
      output$filteredGenesCounts <- NULL
      output$filteredCellStats <- NULL
    }
    output$elbowPlotPCA <- NULL
    output$PCAscatter <- NULL
    output$PCAloadings <- NULL
    output$PCAheatmap <- NULL
    output$gProfilerManhatan <- NULL
    output$annotateClustersCIPRDotplot <- NULL
    output$ligandReceptorFullHeatmap <- NULL
    output$ligandReceptorCuratedHeatmap <- NULL
    output$hvgScatter <- NULL
    output$cellCyclePCA <- NULL
    output$cellCycleBarplot <- NULL
    output$clusterBarplot <- NULL
    output$umapPlot <- NULL
    output$findMarkersHeatmap <- NULL
    output$findMarkersDotplot <- NULL
    output$findMarkersFeaturePlot <- NULL
    output$findMarkersFPcolorbox <- NULL
    output$findMarkersFPfeature1_2 <- NULL
    output$findMarkersFPfeature1 <- NULL
    output$findMarkersFPfeature2 <- NULL
    output$findMarkersViolinPlot <- NULL
    output$findMarkersVolcanoPlot <- NULL
    # renderPlot
    output$trajectoryPlot <- NULL
    output$trajectoryPseudotimePlot <- NULL
    
    # renderDataTable
    output$clusterTable <- NULL
    output$gProfilerTable <- NULL
    output$annotateClustersCIPRTable <- NULL
    output$findMarkersTable <- NULL
    
    # renderPrint
    if (fromDataInput) output$cellStats <- NULL
    output$trajectoryText <- NULL
    
    # dittoDimPlot
    plot_D <- NULL
    plot_P <- NULL
    
    reductions_choices <<- c("-")
    updateSelectInput(session, "umapType", choices = reductions_choices)  #umapType, findMarkersReductionType, cellCycleReduction, trajectoryReduction
    updateSelectInput(session, "findMarkersReductionType", choices = reductions_choices)
    updateSelectInput(session, "findMarkersFeaturePairReductionType", choices = reductions_choices)
    updateSelectInput(session, "cellCycleReduction", choices = reductions_choices)
    updateSelectInput(session, "trajectoryReduction", choices = reductions_choices)
  }
}



