source("help.R", local=TRUE)
library(bsplus)

ui <- dashboardPage(
  title="SCALA",
  skin = "black",
  #------------------------------------------------------------Header
  dashboardHeader(
    titleWidth = "285px",
    title = tags$img(src='logo.png')
  ),
  
  #------------------------------------------------------------Sidebar
  dashboardSidebar(
    width = "285px",
    sidebarMenu(id = "sidebarMenu",
      menuItem(text = "HOME", tabName = "home", icon = icon("home")),
      tags$hr(),
      menuItem(text = "DATA INPUT", tabName = "upload", icon = icon("upload")),
      menuItem(text = "QUALITY CONTROL", tabName = "qc", icon = icon("check-circle")),
      menuItem(tags$div("DATA NORMALIZATION",
                        tags$br(),
                        "& SCALING", class = "menu_item_div"), tabName = "normalize", icon = icon("balance-scale")),
      tags$hr(),
      menuItem(text = "PCA/LSI", tabName = "pca", icon = icon("chart-line")),
      menuItem(text = "CLUSTERING", tabName = "clustering", icon = icon("project-diagram")),
      menuItem(text = "UTILITY OPTIONS", tabName = "utilities", icon = icon("edit")),
      menuItem(tags$div("ADDITIONAL DIMENSIONALITY",
                        tags$br(),
                        "REDUCTION METHODS", class = "menu_item_div"), tabName = "umap", icon = icon("draw-polygon")),
      menuItem(text = "MARKERS' IDENTIFICATION", tabName = "findMarkers", icon = icon("map-marker-alt")),
      menuItem(text = "FEATURE INSPECTION", tabName = "features", icon = icon("braille")),
      menuItem(text = "DOUBLETS' DETECTION", tabName = "doubletDetection", icon = icon("check-double")),
      menuItem(text = "CELL CYCLE PHASE ANALYSIS", tabName = "cellCycle", icon = icon("circle-notch")),
      tags$hr(),
      menuItem(tags$div("FUNCTIONAL/MOTIF",
                        tags$br(),
                        "ENRICHMENT ANALYSIS", class = "menu_item_div"), tabName = "gProfiler", icon=icon("chart-bar")),
      menuItem(text = "CLUSTERS' ANNOTATION", tabName = "annotateClusters", icon = icon("id-card")),
      menuItem(text = "TRAJECTORY ANALYSIS", tabName = "trajectory", icon = icon("route")),
      menuItem(tags$div("LIGAND - RECEPTOR",
                        tags$br(),
                        "ANALYSIS", class = "menu_item_div"), tabName = "ligandReceptor", icon = icon("satellite-dish")), #icon("satellite-dish")),
      menuItem(tags$div("GENE REGULATORY NETWORK",
                        tags$br(),
                        "ANALYSIS", class = "menu_item_div"), tabName = "grn", icon = icon("network-wired")),
      menuItem(text = "TRACKS", tabName = "visualizeTracks", icon = icon("compact-disc")),
      tags$hr(),
      menuItem(text = "Help", tabName = "help", icon = icon("question")),
      menuItem(text = "About", tabName = "about", icon = icon("info"))
    )
  ),
  #------------------------------------------------------------Body
  dashboardBody(
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "main.css")),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href="loading-bar.css")), # loading bar CSS
    tags$head(tags$script(src = "rshiny_handlers.js")), # R to JS
    tags$head(tags$script(src = "loading-bar.js")), # loading bar JS
    #tags$head(tags$script(src = "sliderfix.js")),
    useShinyjs(),
    extendShinyjs(text = js.enrich, functions = c("Enrich")),
    tabItems(
      #home tab
      tabItem(tabName = "home", 
              div(id = "home_div", class = "div_container",
                  h1(class = "container_title", "Welcome to SCALA"),
                  HTML("<p class=container_text> SCALA is a web application and stand-alone toolkit, that handles the analysis of scRNA-seq and scATAC-seq datasets, 
                  from quality control and normalization, to dimensionality reduction, differential expression/accessibility analysis, cell clustering, functional enrichment analysis, 
                  trajectory inference, ligand – receptor analysis, gene regulatory network inference, and visualization. Try out our sample data and visit the Help pages for guidance. </p>"
                  ),
              )
      ),
      
      #Upload tab
      tabItem(tabName = "upload",
              bsCollapse(id = 'countMatrixRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with the upload of count matrix?', ih_inputCountMatrix_rna, style = 'warning')
              ),
              bsCollapse(id = '10xRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with the upload of 10x files?', ih_input10x_rna, style = 'warning')
              ),
              bsCollapse(id = 'rdsRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with the upload of an RDS file?', ih_inputRDS_rna, style = 'warning')
              ),
              bsCollapse(id = 'arrowATAC_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with the upload of an arrow file?', ih_inputArrow_atac, style = 'warning')
              ),
              fluidRow(
                box(
                  width = 3, status = "info", solidHeader = TRUE,
                  title = "Upload your data",
                  tabsetPanel(type = "tabs",
                              tabPanel("Gene-count matrix (scRNA-seq)",
                                       tags$h3("Load PBMC 10x dataset (example scRNA-seq)", class="h3-example"),
                                       tags$hr(class="hr-example"),
                                       actionButton(inputId = "upload10xExampleRNACountMatrixConfirm", label = "Load example", class="btn-example"),
                                       tags$br(),
                                       tags$h3("OR"),
                                       tags$h3("Upload your file"),
                                       tags$hr(),
                                       textInput(inputId = "uploadCountMatrixprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "countMatrix", label = "1. Genes-Cells count matrix", accept = ".txt"),
                                       sliderInput(inputId = "uploadCountMatrixminCells", label = "Include features detected in at least this many cells :", min = 1, max = 20, value = 3, step = 1),
                                       sliderInput(inputId = "uploadCountMatrixminFeatures", label = "Include cells where at least this many features are detected :", min = 1, max = 1000, value = 200, step = 1),
                                       radioButtons("uploadCountMatrixRadioSpecies", label = h3("Select organism : "),
                                                    choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                   "Homo sapiens (Human)" = "human"
                                                    ), 
                                                    selected = "mouse"),
                                       actionButton(inputId = "uploadCountMatrixConfirm", label = "Submit", class="btn-run", icon = icon("check-circle")),
                                       tags$h3("Export working object as .RDS file"),
                                       tags$hr(),
                                       downloadButton(outputId = "utilitiesConfirmExport1", label = "Export .RDS"),
                                       ),
                              tabPanel("10x input files (scRNA-seq)", 
                                       tags$h3("Load PBMC 10x dataset (example scRNA-seq)", class="h3-example"),
                                       tags$hr(class="hr-example"),
                                       actionButton(inputId = "upload10xExampleRNA10xFilesConfirm", label = "Load example", class="btn-example"),
                                       tags$br(),
                                       tags$h3("OR"),
                                       tags$h3("Upload your files"),
                                       tags$hr(),
                                       textInput(inputId = "upload10xRNAprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "barcodes", label = "1. Choose barcodes.tsv.gz file", accept = ".gz"),
                                       fileInput(inputId = "genes", label = "2. Choose features.tsv.gz file", accept = ".gz"),
                                       fileInput(inputId = "matrix", label = "3. Choose matrix.mtx.gz file", accept = ".gz"),
                                       sliderInput(inputId = "upload10xRNAminCells", label = "Include features detected in at least this many cells :", min = 0, max = 20, value = 3, step = 1),
                                       sliderInput(inputId = "upload10xRNAminFeatures", label = "Include cells where at least this many features are detected :", min = 0, max = 1000, value = 200, step = 1),
                                                    radioButtons("upload10xRNARadioSpecies", label = h3("Select organism : "),
                                                                 choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                                "Homo sapiens (Human)" = "human"
                                                                 ), 
                                                                 selected = "mouse"),
                                       actionButton(inputId = "upload10xRNAConfirm", label = "Submit", class="btn-run", icon = icon("check-circle")),
                                       tags$h3("Export working object as .RDS file"),
                                       tags$hr(),
                                       downloadButton(outputId = "utilitiesConfirmExport2", label = "Export .RDS"),
                                       ),
                              tabPanel("RDS Seurat object input (scRNA-seq)",
                                       fileInput(inputId = "uploadRdsFile", label = "Choose a Seurat object saved in .RDS format", accept = ".RDS"),
                                       radioButtons("uploadRdsRadioSpecies", label = h3("Select organism : "),
                                                    choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                   "Homo sapiens (Human)" = "human"
                                                    ), 
                                                    selected = "mouse"),
                                       actionButton(inputId = "uploadSeuratRdsConfirm", label = "Load Seurat object", class="btn-run", icon = icon("check-circle")),
                                       tags$br(),
                                       tags$br(),
                                       tags$hr(),
                                       tags$br(),
                                       selectInput("utilitiesActiveAssay", "(Optional) Change active assay:",
                                                   c("Assay" = "RNA")),
                                       actionButton(inputId = "utilitiesConfirmChangeAssay", label = "Change assay"),
                                       tags$br(),
                                       tags$h3("Export working object as .RDS file"),
                                       tags$hr(),
                                       downloadButton(outputId = "utilitiesConfirmExport", label = "Export .RDS")
                                       ),
                              tabPanel("Arrow input files (scATAC-seq)",
                                       tags$h3("PBMC 10x dataset (example scATAC-seq)", class="h3-example"),
                                       tags$hr(class="hr-example"),
                                       actionButton(inputId = "upload10xExampleATACConfirm", label = "Load example", class="btn-example"),
                                       tags$br(),
                                       tags$h3("OR"),
                                       tags$h3("Upload your file"),
                                       tags$hr(),
                                       textInput(inputId = "uploadATACprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "uploadATACArrow", label = "Please upload an .arrow file", accept = ".arrow"),
                                       radioButtons("upload10xATACRadioSpecies", label = h3("Select organism and genome version: "),
                                                    choices = list("Mus musculus (Mouse) - mm10" = "mm10", 
                                                                   "Homo sapiens (Human) - hg19" = "hg19",
                                                                   "Homo sapiens (Human) - hg38" = "hg38"
                                                    ), 
                                                    selected = "mm10"),
                                       sliderInput(inputId = "upload10xATACThreads", label = "Threads to be used:", min = 1, max = 100, value = 2, step = 1), #max=2 in the online version
                                       actionButton(inputId = "upload10xATACConfirm", label = "Submit", class="btn-run", icon = icon("check-circle"))
                              )
                  )
                ),
                box(
                  width = 8, solidHeader = TRUE, status = "info",
                  title = "Metadata table",
                  
                  tabsetPanel(type = "tabs", id = "uploadTabPanel",
                              tabPanel("scRNA-seq",
                                       dataTableOutput("metadataTable"),
                                       downloadButton(outputId = "uploadMetadataExportRNA", label = "Save table")
                                       ),
                              tabPanel("scATAC-seq",
                                       dataTableOutput("metadataTableATAC"),
                                       downloadButton(outputId = "uploadMetadataExport", label = "Save table")
                                       )
                              )
                )
              ),
      ),
      
      #utilities tab
      tabItem(tabName = "utilities",
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Edit/export working object",
                  tags$h3("Rename cluster"),
                  tags$hr(),
                  selectInput(inputId = "utilitiesRenameOldName", label = "Cluster to be renamed (old name):", choices = "-", multiple = F),
                  textInput(inputId = "utilitiesRenameNewName", label = "New name of the cluster:", value = "New_name_1"),
                  actionButton(inputId = "utilitiesConfirmRename", label = "Rename", class="btn-run", icon = icon("check-circle")),
                  tags$h3("Delete cluster"),
                  tags$hr(),
                  selectInput(inputId = "utilitiesDeleteCluster", label = "Cluster to be deleted:", choices = "-", multiple = F),
                  actionButton(inputId = "utilitiesConfirmDelete", label = "Delete", class="btn-run", icon = icon("check-circle")),
                  
                  tags$h3("Active clusters"),
                  tags$hr(),
                  selectInput("utilitiesActiveClusters", "Set active clustering column:",
                              c("Cluster" = "seurat_clusters")),
                  actionButton(inputId = "utilitiesConfirmChangeCluster", label = "Change clustering column!", class="btn-run", icon = icon("check-circle")),
                  tags$h3("Command history"),
                  tags$hr(),
                  actionButton(inputId = "utilitiesPlotCommands", label = "View command history!", class="btn-run", icon = icon("check-circle")),
                  verbatimTextOutput(outputId = "history")
                )
              )
      ),
      
      #QC tab
      tabItem(tabName = "qc",
              bsCollapse(id = 'qcRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with filtering your data?', ih_qc_rna, style = 'warning')
              ),
              tabsetPanel(type = "tabs", id = "qcTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Quality control",
                                       tags$h3("1. Display quality control plots before filtering"),
                                       actionButton(inputId = "qcDisplay", label = "Display plots", class="btn-run", icon = icon("check-circle")),
                                       tags$hr(),
                                       tags$h3("2. Filter out low quality cells"),
                                       tags$hr(),
                                       
                                       sliderInput(inputId = "minUniqueGenes", label = "Minimum features detected", min = 200, max = 2000, value = 500, step = 1)%>%
                                         shinyInput_label_embed(
                                           shiny_iconlink() %>%
                                             bs_embed_popover(
                                               title = "Filter out cells that have unique feature counts less than:", placement = "left"
                                             )
                                         ), 
                                       sliderInput(inputId = "maxUniqueGenes", label = "Maximum features detected", min = 2001, max = 7000, value = 4500, step = 1)%>%
                                         shinyInput_label_embed(
                                           shiny_iconlink() %>%
                                             bs_embed_popover(
                                               title = "Filter out cells that have unique feature counts over than:", placement = "left"
                                             )
                                         ),
                                       
                                       sliderInput(inputId = "maxMtReads", label = "Mitochondrial %", min = 1, max = 100, value = 10, step = 1)%>%
                                         shinyInput_label_embed(
                                           shiny_iconlink() %>%
                                             bs_embed_popover(
                                               title = "Filter out cells that their percentage of genes mapped to mitochondrial genome exceeds:", placement = "left"
                                             )
                                         ), 
                                       
                                       selectInput("qcColorBy", "Color by:",
                                                   c("orig.ident" = "orig.ident")),
                                       actionButton(inputId = "qcConfirm", label = "Perform filtering", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE,
                                       title = "Quality control plots",

                                       tabsetPanel(type="tabs", id = "qc_tabs_rna",
                                                   tabPanel("Pre-filtering plots",
                                                            column(
                                                              div(id="nFeatureViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "nFeatureViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="totalCountsViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "totalCountsViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="mitoViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "mitoViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="genesCounts_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "genesCounts", height= "100%")
                                                                  )
                                                              ), width = 6),
                                                            column(
                                                              div(id="mtCounts_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "mtCounts", height= "100%")
                                                                  )
                                                              ), width = 6),
                                                            column(verbatimTextOutput(outputId = "cellStats"), width = 4)
                                                   ),
                                                   tabPanel("Post-filtering plots",
                                                            column(
                                                              div(id="filteredNFeatureViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "filteredNFeatureViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="filteredTotalCountsViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "filteredTotalCountsViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="filteredMitoViolin_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "filteredMitoViolin", height = "100%")
                                                                  )
                                                              ), width = 4),
                                                            column(
                                                              div(id="filteredGenesCounts_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "filteredGenesCounts", height= "100%")
                                                                  )
                                                              ), width = 6),
                                                            column(
                                                              div(id="filteredMtCounts_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "filteredMtCounts", height= "100%")
                                                                  )
                                                              ), width = 6),
                                                            column(verbatimTextOutput(outputId = "filteredCellStats"), width = 4)
                                                   )
                                       )
                                   )
                                  )
                          ), 
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Quality control",
                                       tags$h3("Display soft filtered quality control plots"),
                                       actionButton(inputId = "qcDisplayATAC", label = "Display plot!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE,
                                       title = "Quality control plots",
                                       div(
                                         column(
                                           div(id="TSS_plot_loader",
                                               shinycssloaders::withSpinner(
                                                 plotlyOutput(outputId = "TSS_plot", height = "100%")
                                                 )
                                               ), width = 4),
                                         column(
                                           div(id="nFrag_plot_loader",
                                               shinycssloaders::withSpinner(
                                                 plotOutput(outputId = "nFrag_plot", height = "100%")
                                                 )
                                               ), width = 4),
                                         column(
                                           div(id="TSS_nFrag_plot_loader",
                                               shinycssloaders::withSpinner(
                                                 plotlyOutput(outputId = "TSS_nFrag_plot", height = "100%")
                                                 )
                                               ), width = 4),
                                         column(verbatimTextOutput(outputId = "CellStatsATAC"), width = 5)
                                       )
                                     )
                                   )
                          )
              )
      ),
      
      #Normalization tab
      tabItem(tabName = "normalize",
              bsCollapse(id = 'scalingRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with the regress out functionality?', ih_scaling_rna, style = 'warning')
              ),
              tags$div("Normalization and scaling: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 50sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              fluidRow(
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "Normalize and scale the data",
                  tags$h3("1. Log-normalization"),
                  tags$hr(),
                  sliderInput(inputId = "normScaleFactor", label = "Scale factor :", min = 1000, max = 1000000, value = 10000, step = 1000)%>%
                    shinyInput_label_embed(
                      shiny_iconlink() %>%
                        bs_embed_popover(
                          title = "It normalizes the count data per cell and transforms the result to log scale", placement = "left"
                        )
                    ), 
                  tags$h3("2. Identification of highly variable features"),
                  tags$hr(),
                  radioButtons("radioHVG", label = h3("Select one of the following methods : "),
                               choices = list("Variance Stabilizing Transformation method" = "vst", 
                                              
                                              "Mean-Variance method" = "mvp", 
                                              
                                              "Dispersion method" = "disp"), 
                               selected = "vst")%>%
                    shinyInput_label_embed(
                      shiny_iconlink() %>%
                        bs_embed_popover(
                          title = paste0("- vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).\n\n",
                                         "- mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.\n\n",
                                         "- dispersion (disp): selects the genes with the highest dispersion values"), placement = "left"
                        )
                    ), 
                  sliderInput(inputId = "nHVGs", label = "Number of genes to select as top variable genes (applicable only to the first and third option) :", min = 200, max = 8000, value = 2000, step = 100),
                  tags$h3("3. Scaling the data"),
                  tags$hr(),
                  radioButtons("normalizeScaleGenes", label = h3("Genes to be scaled : "),
                               choices = list("All genes" = "all_genes", 
                                              "Only most variable genes" = "mv_genes"), 
                               selected = "mv_genes"),
                  tags$br(),
                  selectInput("normalizeRegressColumns", "Select variables to regress out", list(), selected = NULL, multiple = TRUE, selectize = TRUE, width = NULL, size = NULL)%>%
                    shinyInput_label_embed(
                      shiny_iconlink() %>%
                        bs_embed_popover(
                          title = "Scales and centers features in the dataset. If variables are provided in this text input, they are individually regressed against each feature, and the resulting residuals are then scaled and centered. Variables stored in metadata are valid options.\nThis operation is slow when variables are provided in combination with the \"All genes\" option above.", placement = "left"
                        )
                    ),
                  
                  actionButton(inputId = "normalizeConfirm", label = "Run Normalization and Scaling process!", class="btn-run", icon = icon("check-circle")),
                ),
                box(
                  width = 8, status = "info", solidHeader = TRUE,
                  title = "Highly variable genes",
                  div(id="hvgScatter_loader",
                      shinycssloaders::withSpinner(
                        plotlyOutput(outputId = "hvgScatter", height = "800px")
                      )
                  ),
                  column(verbatimTextOutput(outputId = "hvgTop10Stats"), width = 8)
                ),  
              )
      ),
      
      #PCA tab
      tabItem(tabName = "pca", 
              tags$br(),
              tags$div("PCA estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 8sec (quick version), ~25min (slow version)", tags$br(),
                       "LSI estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 38 sec", tags$br(),
                       "*For datasets containing more than 10,000 cells the slow version of PCA is not suggested", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "pcaTabPanel",
                          tabPanel("scRNA-seq: PCA",
                                   fluidRow(
                                     box(
                                       width = 12, status = "info", solidHeader = TRUE,
                                       title = "PCA results", height = "1200px",
                                       tabsetPanel(type = "tabs", 
                                                   tabPanel("PCA run",
                                                            column(radioButtons("pcaRadio", label = h3("Suggest optimal number of PCs Using 10-fold SVA-CV: "),
                                                                                choices = list("Yes (slow operation)" = "yes", 
                                                                                               "No" = "no"), 
                                                                                selected = "no"), width = 12),
                                                            column(actionButton(inputId = "PCrunPCA", label = "Run PCA!", class="btn-run", icon = icon("check-circle")), width = 12),
                                                            div(
                                                              column(
                                                                div(id="elbowPlotPCA_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      plotlyOutput(outputId = "elbowPlotPCA", height = "750px")
                                                                    )
                                                                ), width = 6),
                                                              column(
                                                                div(id="PCAscatter_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      plotlyOutput(outputId = "PCAscatter", height = "750px")
                                                                    )
                                                                ), width = 6)
                                                            )
                                                   ),
                                                   tabPanel("PCA exploration",
                                                            selectInput("PCin", "Select a principal component :", choices=1:100, selected = 1, multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                                                            column(actionButton(inputId = "PCconfirm", label = "Explore principal component!", class="btn-run", icon = icon("check-circle")),
                                                                   width = 12),
                                                            div(
                                                              column(
                                                                tags$h3("PCA loading scores (top-30 genes for this PC)"),
                                                                div(id="PCAloadings_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      plotlyOutput(outputId = "PCAloadings", height = "700px")
                                                                    )
                                                                ), width = 6),
                                                              column(
                                                                tags$h3("Heatmap of scaled expression (top-30 genes for this PC)"),
                                                                div(id="PCAheatmap_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      plotlyOutput(outputId = "PCAheatmap", height = "700px")
                                                                    )
                                                                ), width = 6)
                                                            ),
                                                            downloadButton(outputId = "pcaRNAExport", label = "Save table")
                                                   )
                                       )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq: LSI",
                                   fluidRow(
                                     box(
                                       width = 6, status = "info", solidHeader = TRUE,
                                       title = "Latent Semantic Indexing", height = "1290px",
                                       tags$h3("Input parameters"),
                                       tags$hr(),
                                       sliderInput(inputId = "lsiVarFeatures", label = "Number of variable feures: ", min = 5000, max = 100000, value = 25000, step = 1000),#varFeatures
                                       sliderInput(inputId = "lsiDmensions", label = "Number of dimensions to use: ", min = 1, max = 100, value = 30, step = 1),#dimensions
                                       sliderInput(inputId = "lsiResolution", label = "Resolution :", min = 0.1, max = 10, value = 1, step = 0.1),#resolution
                                       sliderInput(inputId = "lsiIterations", label = "Number of iterations: ", min = 1, max = 10, value = 1, step = 1),#iterations
                                       actionButton(inputId = "lsiConfirm", label = "Run LSI!", class="btn-run", icon = icon("check-circle")),
                                       tags$hr(),
                                       verbatimTextOutput(outputId = "lsiOutput")
                                     )
                                   )
                          )
              )
      ),
      
      #Clustering tab
      tabItem(tabName = "clustering", 
              tags$br(),
              tags$div("Clustering: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 8sec", tags$br(),
                       "Clustering: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 38sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "clusteringTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(
                                       width = 4, status = "info", solidHeader = TRUE,
                                       title = "Clustering options",
                                       tags$h3("1. Construction of the shared nearest neighbour (SNN) graph"),
                                       tags$hr(),
                                       sliderInput(inputId = "snnK", label = "Number of neighbours for each cell [k]:", min = 1, max = 200, value = 20, step = 1),
                                       sliderInput(inputId = "snnPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                       tags$h3("2. Communities' detection (Louvain algorithm)"),
                                       tags$hr(),
                                       sliderInput(inputId = "clusterRes", label = "Clustering resolution :", min = 0.1, max = 5, value = 0.5, step = 0.1),
                                       actionButton(inputId = "snnConfirm", label = "Run clustering", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Clustering output",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Clustering results",
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Cluster table",
                                                                                 dataTableOutput(outputId="clusterTable"),
                                                                                 downloadButton(outputId = "clusterTableRNAExport", label = "Save table")
                                                                        ),
                                                                        tabPanel("Cluster barplot",
                                                                                 selectInput("clusterGroupBy", "Grouping variable:",
                                                                                             c("orig.ident" = "orig.ident")),
                                                                                 actionButton(inputId = "clusterBarplotConfirm", label = "Display barchart!", class="btn-run", icon = icon("check-circle")),
                                                                                 div(id="clusterBarplot_loader",
                                                                                     shinycssloaders::withSpinner(
                                                                                       plotlyOutput(outputId = "clusterBarplot", height = "700px")
                                                                                     )
                                                                                 )
                                                                        )
                                                            )															
                                                   ),
                                                   tabPanel("Shared Nearest Neighbour (SNN) Graph", 
                                                            actionButton(inputId = "snnDisplayConfirm", label = "Display SNN graph!", class="btn-run", icon = icon("check-circle")),
                                                            div(id="snnSNN_loader",
                                                                shinycssloaders::withSpinner(
                                                                  visNetworkOutput(outputId="snnSNN", height = "1300px")
                                                                )
                                                            )
                                                   )
                                       ),
                                     ),
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 4, status = "info", solidHeader = TRUE,
                                       title = "Clustering options",
                                       sliderInput(inputId = "clusterDimensionsATAC", label = "Number of dimensions to use: ", min = 1, max = 100, value = 30, step = 1),
                                       sliderInput(inputId = "clusterResATAC", label = "Clustering resolution :", min = 0.1, max = 60, value = 0.6, step = 0.1),
                                       actionButton(inputId = "clusterConfirmATAC", label = "Perform clustering!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Clustering output",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Clustering results", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Cluster table",
                                                                                 dataTableOutput(outputId="clusterTableATAC"),
                                                                                 downloadButton(outputId = "clusterTableExportATAC", label = "Save table")
                                                                        ),
                                                                        tabPanel("Cluster barplot",
                                                                                 div(id="clusterBarplotATAC_loader",
                                                                                     shinycssloaders::withSpinner(
                                                                                       plotlyOutput(outputId = "clusterBarplotATAC", height = "700")
                                                                                     )
                                                                                  )
                                                                        )
                                                                      )
                                                   )
                                       ),
                                     ),
                                   )
                          )
              )
      ),
      
      #UMAP tab
      tabItem(tabName = "umap", 
              tags$br(),
              tags$div("UMAP and tSNE: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 47sec", tags$br(),
                       "UMAP and tSNE: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 1min 20sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "umapTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Cells visualization options in reduced space",
                                         sliderInput(inputId = "umapSeed", label = "Set seed :", min = 1, max = 500, value = 42, step = 1),
                                         sliderInput(inputId = "umapPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                         sliderInput(inputId = "umapOutComponents", label = "Number of dimensions to fit output:", min = 2, max = 100, value = 3, step = 1)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "If PHATE is selected, the runtime increases when a value > 3 is used.\nPlease note that tSNE doesn't return more than 3 dimensions.", placement = "bottom"
                                               )
                                           ),
                                         actionButton(inputId = "umapRunUmap", label = "Run UMAP!", class="btn-run", icon = icon("check-circle")),
                                         actionButton(inputId = "umapRunTsne", label = "Run tSNE!", class="btn-run", icon = icon("check-circle")),
                                         actionButton(inputId = "umapRunDFM", label = "Run Diffusion Map!", class="btn-run", icon = icon("check-circle")),
                                         actionButton(inputId = "umapRunPhate", label = "Run PHATE!", class="btn-run", icon = icon("check-circle")),
                                         tags$h3("Display settings"),
                                         tags$hr(),
                                         selectInput("umapType", "Plot type:",
                                                     c("-" = "-")
                                         ),
                                         selectInput("umapDimensions", "Dimensions:",
                                                     c("2D" = "2",
                                                       "3D" = "3")),
                                         selectInput("umapColorBy", "Color by:",
                                                     c("Cluster" = "seurat_clusters")),
                                         
                                         sliderInput("umapDotSize", "Size:", min = 1, max = 20, value = 5, step = 0.5),
                                         sliderInput("umapDotOpacity", "Opacity:", min = 0, max = 1, value = 1, step = 0.1),
                                         sliderInput("umapDotBorder", "Border width:", min = 0, max = 10, value = 0.5, step = 0.1),
                                         actionButton(inputId = "umapConfirm", label = "Update plot!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                                         div(id="umapPlot_loader",
                                             shinycssloaders::withSpinner(
                                               plotlyOutput(outputId = "umapPlot", height = "1100px")
                                             )
                                         )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Cells visualization options in reduced space",
                                         sliderInput(inputId = "umapDimensionsATAC", label = "Number of input dimensions to use :", min = 1, max = 100, value = 30, step = 1),
                                         sliderInput(inputId = "umapOutComponentsATAC", label = "Number of dimensions to fit output:", min = 2, max = 100, value = 3, step = 1)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "Please note that tSNE doesn't return more than 2 dimensions.", placement = "bottom"
                                               )
                                           ),
                                         actionButton(inputId = "umapRunUmapTsneATAC", label = "Run UMAP and tSNE!", class="btn-run", icon = icon("check-circle")),
                                         tags$h3("Display settings"),
                                         tags$hr(),
                                         selectInput("umapTypeATAC", "Plot type:",
                                                     c("UMAP" = "umap",
                                                       "tSNE" = "tsne")
                                         ),
                                         selectInput("umapDimensionsPlotATAC", "Dimensions:",
                                                     c("2D" = "2",
                                                       "3D" = "3")),
                                         selectInput("umapColorByATAC", "Color by:",
                                                     c("Clusters" = "Clusters")),
                                         
                                         sliderInput("umapDotSizeATAC", "Size:", min = 1, max = 20, value = 5, step = 0.5),
                                         sliderInput("umapDotOpacityATAC", "Opacity:", min = 0, max = 1, value = 1, step = 0.1),
                                         sliderInput("umapDotBorderATAC", "Border width:", min = 0, max = 10, value = 0.5, step = 0.1),
                                         actionButton(inputId = "umapConfirmATAC", label = "Display plot!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                                         div(id="umapPlotATAC_loader",
                                             shinycssloaders::withSpinner(
                                               plotlyOutput(outputId = "umapPlotATAC", height = "1100px")
                                             )
                                         )
                                     )
                                   )
                          )
              )
      ),
      
      #Feature inspection
      tabItem(tabName = "features",
              bsCollapse(id = 'signature_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with adding a new signature?', ih_signatureFP_rna, style = 'warning')
              ),
              tags$div("Signature scoring: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 33sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be 
                                                          achieved by using the stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "featuresTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Feature plot", fluidRow(
                                                   box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                       selectizeInput(inputId = 'findMarkersGeneSelect',
                                                                      label = 'Type the name of a gene, gene signature or numeric metadata column: 
                                                                      (e.g. "Ccl2" or "Signature1_Ucell", or "nCount_RNA")',
                                                                      choices = NULL,
                                                                      selected = NULL,
                                                                      multiple = FALSE),
                                                       selectInput("findMarkersReductionType", "Plot type:",
                                                                   c("-" = "-")
                                                       ),
                                                       radioButtons("findMarkersLabels", label = "Show cluster labels: ",
                                                                    choices = list("Yes" = TRUE, 
                                                                                   "No" = FALSE)
                                                       ),
                                                       radioButtons("findMarkersOrder", label = "Prioritize expressing cells: ",
                                                                    choices = list("Yes" = TRUE, 
                                                                                   "No" = FALSE)
                                                       ),
                                                       sliderInput("findMarkersMaxCutoff", "Set max expression value: (quantile)", min = 0, max = 99, value = 99, step = 1),
                                                       sliderInput("findMarkersMinCutoff", "Set minimum expression value: (quantile)", min = 0, max = 99, value = 0, step = 1),
                                                       actionButton(inputId = "findMarkersFPConfirm", label = "Display plot!", class="btn-run", icon = icon("check-circle")),
                                                       tags$hr(),
                                                       tags$h3("Add a new signature"),
                                                       textInput(inputId = "findMarkersSignatureName", label = "Gene signature name :", value = "Signature1"),
                                                       textAreaInput(inputId = "findMarkersSignatureMembers", label = "Paste a list of genes", cols = 80, rows = 15, placeholder = "Prg4\nTspan15\nCol22a1\nHtra4"),
                                                       actionButton(inputId = "findMarkersSignatureAdd", label = "Calculate signature score!", class="btn-run", icon = icon("check-circle"))
                                                   ),
                                                   box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                       div(id="findMarkersFeaturePlot_loader",
                                                           shinycssloaders::withSpinner(
                                                             plotlyOutput(outputId = "findMarkersFeaturePlot", height = "1300px")
                                                           )
                                                       )
                                                   )
                                                 )),
                                                 tabPanel("Multi-feature vizualization", fluidRow(
                                                   box(width=3, status="info", solidHeader=T, title="Options",
                                                       selectizeInput(inputId = 'findMarkersFeaturePair1',
                                                                      label = 'Select 1st feature:',
                                                                      choices = NULL,
                                                                      selected = NULL,
                                                                      multiple = FALSE),
                                                       selectizeInput(inputId = 'findMarkersFeaturePair2',
                                                                      label = 'Select 2nd Feature:',
                                                                      choices = NULL,
                                                                      selected = NULL,
                                                                      multiple = FALSE),
                                                       sliderInput("findMarkersBlendThreshold", "Select threshold for blending:", min = 0, max = 1, value = 0.5, step = 0.1),
                                                       selectInput("findMarkersFeaturePairReductionType", "Plot type:",
                                                                   c("-" = "-")
                                                       ),
                                                       radioButtons("findMarkersFeaturePairLabels", label = "Show cluster labels: ",
                                                                    choices = list("Yes" = TRUE, 
                                                                                   "No" = FALSE)
                                                       ),
                                                       radioButtons("findMarkersFeaturePairOrder", label = "Prioritize expressing cells: ",
                                                                    choices = list("Yes" = TRUE, 
                                                                                   "No" = FALSE)
                                                       ),
                                                       sliderInput("findMarkersFeaturePairMaxCutoff", "Set max expression value: (quantile)", min = 0, max = 99, value = 99, step = 1),
                                                       sliderInput("findMarkersFeaturePairMinCutoff", "Set minimum expression value: (quantile)", min = 0, max = 99, value = 0, step = 1),
                                                       actionButton(inputId = "findMarkersFeaturePairConfirm", label = "Display plot!", class="btn-run", icon = icon("check-circle"))
                                                   ),
                                                   box(width=9, status="info", solidHeader=TRUE, title="Plot",
                                                       div(
                                                         column(
                                                           div(id="findMarkersFPfeature1_loader",
                                                               shinycssloaders::withSpinner(
                                                                 plotlyOutput(outputId="findMarkersFPfeature1", height = "650px")
                                                               )
                                                           ), width = 6),
                                                         column(
                                                           div(id="findMarkersFPfeature2_loader",
                                                               shinycssloaders::withSpinner(
                                                                 plotlyOutput(outputId="findMarkersFPfeature2", height = "650px")
                                                               )
                                                           ), width = 6),
                                                         column(
                                                           div(id="findMarkersFPfeature1_2_loader",
                                                               shinycssloaders::withSpinner(
                                                                 plotlyOutput(outputId="findMarkersFPfeature1_2", height = "650px")
                                                               )
                                                           ), width = 6),
                                                         column(
                                                           div(id="findMarkersFPcolorbox_loader",
                                                               shinycssloaders::withSpinner(
                                                                 plotlyOutput(outputId="findMarkersFPcolorbox", height = "650px")
                                                               )
                                                           ), width = 6),
                                                       )
                                                   )
                                                 )
                                                 ),
                                                 tabPanel("Violin plot", fluidRow(
                                                   box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                       selectizeInput(inputId = 'findMarkersGeneSelect2',
                                                                      label = 'Type the name of a gene, gene signature or numeric metadata column: 
                                                                      (e.g. "Ccl2" or "Signature1_Ucell", or "nCount_RNA")',
                                                                      choices = NULL,
                                                                      selected = NULL,
                                                                      multiple = FALSE), # allow for multiple inputs
                                                       actionButton(inputId = "findMarkersViolinConfirm", label = "Display plot!", class="btn-run", icon = icon("check-circle"))
                                                   ),
                                                   box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                       div(id="findMarkersViolinPlot_loader",
                                                           shinycssloaders::withSpinner(
                                                             plotlyOutput(outputId = "findMarkersViolinPlot", height = "800px")
                                                           )
                                                       )
                                                   )
                                                 )
                                                 )
                                     )
                                   )
                                  ),
                          tabPanel("scATAC-seq",
                                    fluidRow(
                                      box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                          selectizeInput(inputId = 'findMarkersGeneSelectATAC',
                                                         label = 'Select a gene:',
                                                         choices = NULL,
                                                         selected = NULL,
                                                         multiple = FALSE),
                                          selectInput("findMarkersReductionTypeATAC", "Plot type:",
                                                      c("UMAP" = "umap",
                                                        "tSNE" = "tsne")
                                          ),
                                          actionButton(inputId = "findMarkersFPConfirmATAC", label = "Display plot!", class="btn-run", icon = icon("check-circle")),
                                      ),
                                      box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                          div(id="findMarkersFeaturePlotATAC_loader",
                                              shinycssloaders::withSpinner(
                                                plotOutput(outputId = "findMarkersFeaturePlotATAC", height = "1100px")
                                              )
                                          )
                                      )
                                    ) 
                                  )
                          )
      ),
      
      #ATAC
      #DEA tab
      tabItem(tabName = "findMarkers", 
              tags$br(),
              tags$div("Marker genes: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 1min 36sec", tags$br(),
                       "Marker genes: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 1min 15sec", tags$br(),
                       "Marker peaks: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 8min 30sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "findMarkersTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Differential Expression Analysis options", 
                                         selectInput("findMarkersTest", "Test used:",
                                                     c("Wilcoxon rank sum test" = "wilcox",
                                                       "Likelihood-ratio test for single cell feature expression" = "bimod",
                                                       "Standard AUC classifier" = "roc",
                                                       "Student's t-test" = "t",
                                                       "MAST" = "MAST",
                                                       "DESeq2" = "DESeq2"
                                                     )),
                                         radioButtons("findMarkersLogBase", label = "Base used for average logFC calculation: ",
                                                      choices = list("log(e)" = "avg_logFC", 
                                                                     "log(2)" = "avg_log2FC"
                                                      ), 
                                                      selected = "avg_logFC"),
                                         
                                         sliderInput(inputId = "findMarkersMinPct", label = "Minimum % of expression", min = 0, max = 1, value = 0.25, step = 0.05)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "Only test genes that are detected in a minimum fraction of cells in either of the two populations:", placement = "bottom"
                                               )
                                           ),
                                         
                                         sliderInput(inputId = "findMarkersLogFC", label = "Avg Log FC threshold", min = 0, max = 3, value = 0.25, step = 0.05)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells:", placement = "bottom"
                                               )
                                           ),
                                         
                                         sliderInput(inputId = "findMarkersPval", label = "P-value threshold", min = 0, max = 1, value = 0.01, step = 0.01)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "Only return markers that have a p-value < slected threshold, or a power > selected threshold (if the test is ROC) :", placement = "bottom"
                                               )
                                           ),
                                         actionButton(inputId = "findMarkersConfirm", label = "Find Marker genes!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE, title = "DEA results",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Marker genes", 
                                                            dataTableOutput(outputId="findMarkersTable"),
                                                            downloadButton(outputId = "findMarkersRNAExport", label = "Save table")),
                                                   tabPanel("Heatmap", 
                                                            actionButton(inputId = "findMarkersTop10HeatmapConfirm", label = "Display top-10 marker genes heatmap!", class="btn-run", icon = icon("check-circle")),
                                                            div(id="findMarkersHeatmap_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "findMarkersHeatmap", height = "1300px")
                                                                  )
                                                                )
                                                            ),
                                                   
                                                   tabPanel("Dotplot", 
                                                            actionButton(inputId = "findMarkersTop10DotplotConfirm", label = "Display top-10 marker genes dotplot!", class="btn-run", icon = icon("check-circle")),
                                                            div(id="findMarkersDotplot_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "findMarkersDotplot", height = "1300px")
                                                                )
                                                            )
                                                   ),
                                                   tabPanel("VolcanoPlot", fluidRow(
                                                     box(width = 3, status = "info", solidHeader = TRUE, title = "Cluster selection",
                                                         selectInput("findMarkersClusterSelect", "Cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                                                         actionButton(inputId = "findMarkersVolcanoConfirm", "Display volcano plot!", class="btn-run", icon = icon("check-circle"))
                                                         ),
                                                     
                                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Volcano plot",
                                                         div(id="findMarkersVolcanoPlot_loader",
                                                             shinycssloaders::withSpinner(
                                                               plotlyOutput(outputId = "findMarkersVolcanoPlot", height = "800px")
                                                               )
                                                         )
                                                     )
                                                  )
                                                )
                                       )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Marker genes/peaks detection options (scATAC-seq)", 
                                         tags$h3("Marker genes"),
                                         tags$hr(),
                                         selectInput("findMarkersTestATAC", "Test used:",
                                                     c("Wilcoxon" = "wilcoxon",
                                                       "Binomial" = "binomial",
                                                       "T-test" = "ttest"
                                                     )),
                                         selectInput("findMarkersGroupByATAC", "Cells group by:",
                                                     c("Clusters" = "Clusters",
                                                       "Integration predicted clusters" = "predictedGroup_Co"
                                                     )),
                                         sliderInput(inputId = "findMarkersLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.01),
                                         
                                         sliderInput(inputId = "findMarkersFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         actionButton(inputId = "findMarkersConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle")),
                                         
                                         tags$h3("Marker peaks"),
                                         tags$hr(),
                                         selectInput("findMarkersPeaksTestATAC", "Test used:",
                                                     c("Wilcoxon" = "wilcoxon",
                                                       "Binomial" = "binomial",
                                                       "T-test" = "ttest"
                                                     )),
                                         selectInput("findMarkersPeaksGroupByATAC", "Cells group by:",
                                                     c("Clusters" = "Clusters",
                                                       "Integration predicted clusters" = "predictedGroup_Co"
                                                     )),
                                         sliderInput(inputId = "findMarkersPeaksLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.01),
                                         
                                         sliderInput(inputId = "findMarkersPeaksFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         fileInput(inputId = "findMarkersPeaksCustomPeaks", label = "Please upload a .bed file (If you are using the example dataset you can upload the peakset file (.bed) provided in Help > Examples)", accept = ".bed"),
                                         #--activation in local version--
                                         # textInput(inputId = "pathToMacs2", label = "Absolute path to MACS2")%>%
                                         #   shinyInput_label_embed(
                                         #     shiny_iconlink() %>%
                                         #       bs_embed_popover(
                                         #         title = "Absolute path to MACS2 installation folder:
                                         #         Windows OS: the path will be detected automatically.
                                         #         Linux OS: provide the path for the MACS2, e.g. /home/user/anaconda3/bin/macs2", placement = "left"
                                         #       )
                                         #   ),
                                         
                                         actionButton(inputId = "findMarkersPeaksConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE,
                                       tabsetPanel(type = "tabs", id = "ATAC_markers_tabs",
                                                   tabPanel("Marker genes (ATAC)", fluidRow(
                                                     tabsetPanel(type = "tabs", id = "marker_genes_tab_id",
                                                                 tabPanel("Marker genes table",
                                                                          
                                                                          div(id="findMarkersGenesATACTable_loader",
                                                                              shinycssloaders::withSpinner(
                                                                          dataTableOutput(outputId="findMarkersGenesTableATAC")
                                                                            )
                                                                          ),
                                                                          downloadButton(outputId = "findMarkersGenesATACExport", label = "Save table"),
                                                                 ),
                                                                 tabPanel("Marker genes heatmap (top-10)",
                                                                          div(id="findMarkersGenesHeatmapATAC_loader",
                                                                              shinycssloaders::withSpinner(
                                                                                plotlyOutput(outputId = "findMarkersGenesHeatmapATAC", height = "700px")
                                                                              )
                                                                          )
                                                                 )
                                                     )
                                                    )
                                                   ),
                                                   tabPanel("Marker peaks (ATAC)", fluidRow(
                                                     tabsetPanel(type = "tabs", id = "marker_peaks_tab_id",
                                                       tabPanel("Marker peaks table",
                                                                div(id="findMarkersPeaksATACTable_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      dataTableOutput(outputId="findMarkersPeaksTableATAC"),
                                                                    )
                                                                ),
                                                                downloadButton(outputId = "findMarkersPeaksATACExport", label = "Save table")
                                                       ),
                                                       tabPanel("Marker peaks heatmap (top-10)",
                                                                div(id="findMarkersPeaksHeatmapATAC_loader",
                                                                    shinycssloaders::withSpinner(
                                                                      plotlyOutput(outputId = "findMarkersPeaksHeatmapATAC", height = "700px")
                                                                    )
                                                                )
                                                       )
                                                     )
                                                    )
                                                   )
                                       )
                                     )
                                   )
                          )
              )	
      ),
      
      #Doublets' detection
      tabItem(tabName = "doubletDetection", 
              tags$br(),
              tags$div("Doublet detection: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 2min 16sec", tags$br(),
                       "Doublet detection: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 7min 50sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "doubletDetectionTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(
                                       width = 4, status = "info", solidHeader = TRUE,
                                       title = "Doublet detection parameters",
                                       tags$h3("1. Options for doublet dection"),
                                       tags$hr(),
                                       sliderInput(inputId = "doubletsPN", label = "Artificial doublet's rate :", min = 0.01, max = 1, value = 0.25, step = 0.01),
                                       sliderInput(inputId = "doubletsPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                       radioButtons("doubletsPKRadio", label = "PC neighborhood size estimation: ",
                                                    choices = list("Automatic" = "auto", 
                                                                   "Manual" = "manual"
                                                    ), 
                                                    selected = "auto"),
                                       sliderInput(inputId = "doubletsPK", label = "PC neighborhood size (used only when manual is selected):", min = 0, max = 1, value = 0.09, step = 0.01),
                                       sliderInput(inputId = "doubletsNExp", label = "Percentage of doubles expected :", min = 0.01, max = 1, value = 0.03, step = 0.01),
                                       selectInput("doubletsReduction", "Plot type:",
                                                   c("-" = "-")),
                                       actionButton(inputId = "doubletsConfirm", label = "Perform doublets' detection!", class="btn-run", icon = icon("check-circle")),
                                       tags$h3("2. Options for doublet removal (optional)"),
                                       tags$hr(),
                                       radioButtons("doubletsCellRemoval", label = "Remove doublet cells: ",
                                                    choices = list("All predicted doublets" = "all", 
                                                                   "Only heterotypic doublets" = "heterotypic"
                                                    ), 
                                                    selected = "all"),
                                       actionButton(inputId = "doubletsRemove", label = "Delete doublets!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Doublet detection output",
                                       verbatimTextOutput(outputId = "doubletsInfo"),
                                       plotlyOutput(outputId = "doubletPCAplot", height = "700px")
                                       )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 4, status = "info", solidHeader = TRUE,
                                       title = "Doublet detection parameters",
                                       tags$h3("1. Options for doublet dection"),
                                       tags$hr(),
                                       sliderInput(inputId = "doubletsATACk", label = "The number of cells neighboring a simulated doublet to be considered as putative doublets :", min = 5, max = 100, value = 10, step = 1),
                                       radioButtons("doubletsATACLSI", label = "Order of operations in the TF-IDF normalization: ",
                                                    choices = list("tf-logidf" = "1", 
                                                                   "log(tf-idf)" = "2",
                                                                   "logtf-logidf" = "3"
                                                    ), 
                                                    selected = "1"),
                                       actionButton(inputId = "doubletsATACConfirm", label = "Perform doublets' detection!", class="btn-run", icon = icon("check-circle")), 
                                       tags$h3("2. Options for doublet removal (optional)"),
                                       tags$hr(),
                                       sliderInput(inputId = "doubletsATACfilterRatio", label = "The maximum ratio of predicted doublets to filter 
                                                   based on the number of pass-filter cells:", min = 0.1, max = 5, value = 1, step = 0.1),
                                       actionButton(inputId = "doubletsATACDelete", label = "Delete doublets", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Doublet detection output",
                                       div(id="doubletATAC_loader3",
                                           shinycssloaders::withSpinner(
                                            plotOutput(outputId = "doubletsScoreATAC")
                                            ), width = 6
                                           ),
                                       div(id="doubletATAC_loader4",
                                           shinycssloaders::withSpinner(
                                            plotOutput(outputId = "doubletEnrichmentATAC")
                                           ), width = 6
                                       )
                                     ),
                                   )
                          )
              )
      ),
      
      #Cell cycle phase analysis
      tabItem(tabName = "cellCycle",
              tags$br(),
              tags$div("Cell cycle phase analyis: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 4sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = T,
                  title = "Cell cycle phase analysis",
                  tabsetPanel(type = "tabs",
                              tabPanel("Dimensionality reduction plot", 
                                       tags$br(),
                                       tags$div("After succesfully running cell cycle 
                                       phase analysis the columns S.Score, G2M.Score and CC.Difference are stored in the metadata
                                       table. You can regress out the cell cycle effect by returning to the tab DATA NORMALIZATION
                                       & SCALING and repeating step3, after selecting the preferred metadata variables.",id="cellCycleMessage"),
                                       tags$br(),
                                       tags$hr(),
                                       selectInput("cellCycleReduction", "Plot type:",
                                                   c("-" = "-")
                                       ),
                                       actionButton(inputId = "cellCycleRun", label = "Run cell cycle analysis!", class="btn-run", icon = icon("check-circle")),
                                       div(id="cellCyclePCA_loader",
                                           shinycssloaders::withSpinner(
                                             plotlyOutput(outputId = "cellCyclePCA", height = "700px")
                                           )
                                       ),
                              ),
                              tabPanel("Barplot",
                                       div(id="cellCycleBarplot_loader",
                                           shinycssloaders::withSpinner(
                                             plotlyOutput(outputId = "cellCycleBarplot", height = "1100px")
                                           )
                                       )
                              )
                  )
                )
              )
      ),
      
      #Enrichment analysis -gProfiler
      tabItem(tabName = "gProfiler", 
              bsCollapse(id = 'flameRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with enrichment analysis of multiple lists using Flame?', ih_flame_rna, style = 'warning')
              ),
              tags$div("Functional enrichment analysis: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 3sec (per genelist)", tags$br(),
                       "Motif enrichment analysis: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 7min 05sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type = "tabs", id = "gProfilerTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(width = 2, status = "info", solidHeader = TRUE,
                                         title = "Enrichment analysis options",
                                         tags$h3("1. Options for input list"),
                                         tags$hr(),
                                         selectInput("gProfilerList", "Input list of genes:",
                                                     c("-" = "-")),
                                         radioButtons("gprofilerRadio", label = "Sigificance threshold : ",
                                                      choices = list("P-value" = "p_val", 
                                                                     "Adjusted P-value" = "p_val_adj",
                                                                     "Power" = "power"
                                                      ), 
                                                      selected = "p_val"),
                                         sliderInput("gProfilerSliderSignificance", "", min = 0, max = 1, value = 0.01, step = 0.01),
                                         radioButtons("gProfilerLFCRadio", label = "Direction of deregulation : ",
                                                      choices = list("Up-regulated" = "Up",
                                                                     "Down-regulated" = "Down"),
                                                      selected = "Up"
                                         ),
                                         sliderInput("gProfilerSliderLogFC", "Log FC threshold:", min = 0, max = 3, value = 0.25, step = 0.01),
                                         tags$h3("2. Options for enrichment analysis"),
                                         tags$hr(),
                                         selectInput("gProfilerDatasources", "Select datasources", list('Gene Ontology'=list("Gene Ontology-Molecular Function (GO:MF)"="GO:MF", "Gene Ontology-Cellular Component (GO:CC)"= "GO:CC", "Gene Ontology-Biological Process (GO:BP)"="GO:BP"),
                                                                                                        'Biological Pathways'= list("KEGG PATHWAY"="KEGG", "Reactome"="REAC", "WikiPathways"="WP"),
                                                                                                        'Regulatory motifs in DNA'= list("TRANSFAC"= "TF","miRTarBase"= "MIRNA"),
                                                                                                        'Protein Databases'=list("CORUM"= "CORUM", "Human Protein Atlas (HPA)"="HPA"),
                                                                                                        'Human Phenotype Ontology' =list("Human Phenotype Ontology"= "HP")),
                                                     selected = list("GO:MF", "GO:CC", "GO:BP", "KEGG"),
                                                     multiple = TRUE, selectize = TRUE, width = NULL, size = NULL),
                                         selectInput("gProfilerOrganism", "Select organism", choices=c("Homo sapiens (Human)"="hsapiens","Mus musculus (Mouse)"="mmusculus"), selected = NULL, multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                                         radioButtons("gprofilerRadioCorrection", label = "Correction method for multiple testing : ",
                                                      choices = list("Analytical(g:SCS)" = "gSCS", 
                                                                     "Benjamini-Hochberg false discovery rate" = "fdr",
                                                                     "Bonferroni correction" = "bonferroni"
                                                      ),
                                                      selected = "bonferroni"
                                         ), 
                                         sliderInput("gProfilerSliderSignificanceTerms", "Significance for enriched terms :", min = 0, max = 1, value = 0.05, step = 0.01),
                                         actionButton(inputId = "gProfilerConfirm", label = "Run enrichment analysis!", class="btn-run", icon = icon("check-circle")),
                                         tags$h3("3. Multiple cluster enrichment analysis"),
                                         tags$hr(),
                                         selectizeInput(
                                           "gProfilerFlameSelection",
                                           label = "Select up to 10 clusters",
                                           choices = c("-"="-"),
                                           multiple = TRUE,
                                           options = list(maxItems = 10)
                                         ),
                                         actionButton(inputId = "sendToFlame", label = "Send to Flame!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     box(
                                       width = 10, status = "info", solidHeader = TRUE, title = "Enrichment analysis results",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Table of functional terms", 
                                                            dataTableOutput(outputId = "gProfilerTable"),
                                                            downloadButton(outputId = "gProfilerRNAExport", label = "Save table")),
                                                   tabPanel("Manhattan plot", 
                                                            div(id="gProfilerManhattan_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "gProfilerManhattan")
                                                                )
                                                            )
                                                   )
                                       )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Motif enrichment analysis (scATAC-seq)", 
                                         selectInput("findMotifsSetATAC", "Motif set:",
                                                     c("Cisb" = "cisbp",
                                                       "ENCODE" = "encode",
                                                       "Homer" = "homer",
                                                       "JASPAR 2016" = "JASPAR2016",
                                                       "JASPAR 2018" = "JASPAR2018",
                                                       "JASPAR 2020" = "JASPAR2020"
                                                     )),
                                         selectInput("findMotifsGroupByATAC", "Cells group by:",
                                                     c("Clusters" = "Clusters",
                                                       "Integration predicted clusters" = "predictedGroup_Co"
                                                     )),
                                         sliderInput(inputId = "findMotifsLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.01),
                                         sliderInput(inputId = "findMotifsFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         actionButton(inputId = "findMotifsConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Motif enrichment analysis results",
                                         tabsetPanel(type = "tabs",
                                                     tabPanel("Table of enriched motifs", 
                                                                div(id="findMotifsATACTable_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    dataTableOutput(outputId="findMotifsTableATAC")
                                                                  )
                                                                ),
                                                                downloadButton(outputId = "findMotifsATACExport", label = "Save table")
                                                              ),
                                                     tabPanel("Heatmap of enriched motifs (top-10)", 
                                                              div(id="findMotifsHeatmapATAC_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "findMotifsHeatmapATAC", height = "800px")
                                                                  )
                                                                )
                                                              )
                                                     )
                                     )
                                   )
                          )
              )  
      ),
      
      #Clusters' annotation
      tabItem(tabName = "annotateClusters",
              tags$br(),
              tags$div("Cell type annotation: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 22sec", tags$br(),
                       "scATAC-scRNA integration: estimated time in web server for two datasets of 12,000 cells ~ 11min 12sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type="tabs", id = "annotateClustersTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Annotation parameters",
                                       radioButtons("annotateClustersReference", label = "Reference dataset : ",
                                                    choices = list("ImmGen (mouse)" = "immgen", 
                                                                   "Presorted RNAseq (mouse)" = "mmrnaseq",
                                                                   "Blueprint-Encode (human)" = "blueprint",
                                                                   "Primary Cell Atlas (human)" = "hpca",
                                                                   "DICE (human)" = "dice",
                                                                   "Hematopoietic diff (human)" = "hema",
                                                                   "Presorted RNA seq (human)" = "hsrnaseq"
                                                    ),
                                                    selected = "mmrnaseq"
                                       ),
                                       tags$hr(),
                                       sliderInput("annotateClustersSlider", "Keep top Nth % of variable genes in reference :", min = 0, max = 100, value = 100, step = 1),
                                       tags$hr(),
                                       radioButtons("annotateClustersMethod", label = "Select method for comparisons : ",
                                                    choices = list("logFC dot product" = "logfc_dot_product", 
                                                                   "logFC Spearman" = "logfc_spearman",
                                                                   "logFC Pearson" = "logfc_pearson",
                                                                   "Spearman (all genes)" = "all_genes_spearman",
                                                                   "Pearson (all genes)" = "all_genes_pearson"
                                                    ),
                                                    selected = "all_genes_pearson"
                                       ),
                                       tags$hr(),
                                       actionButton(inputId = "annotateClustersConfirm", label = "Run cluster annotation analysis!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE, title = "Cell type annotation",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Top-5 hits table", 
                                                            dataTableOutput(outputId="annotateClustersCIPRTable"),
                                                            downloadButton(outputId = "annotationRNAExport", label = "Save table")
                                                   ),
                                                   tabPanel("Top-5 hits dotplot", 
                                                            div(id="annotateClustersCIPRDotplot_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId="annotateClustersCIPRDotplot", height = "1100px")
                                                                )
                                                            )
                                                   )
                                       )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE, title = "Annotation options",
                                       fileInput(inputId = "annotateClustersRDSInput", label = "Upload an .RDS file", accept = ".RDS"),
                                       actionButton(inputId = "annotateClustersConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle")),
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE, title = "Cell type annotation from scRNA-seq",
                                       div(id="annotateClustersUMAP_loader",
                                           shinycssloaders::withSpinner(
                                             plotlyOutput(outputId="annotateClustersUMAPplot", height = "1100px")
                                           )
                                       )
                                     )
                                   )
                          )
              ) 
      ),
      
      #Trajectory analysis
      tabItem(tabName = "trajectory",
              tags$br(),
              tags$div("Trajectory inference: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 14sec", tags$br(),
                       "Trajectory inference: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 1min 09sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type="tabs", id = "trajectoryTabPanel",
                          tabPanel("scRNA-seq", 
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Trajectory parameters",
                                       selectInput("trajectoryReduction", "Dimensionality reduction method:", 
                                                   choices=c("PCA"="pca","UMAP"="umap", "tSNE"="tsne", "Diffusion Map"="dfm"), selected = "PCA",
                                                   multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                                       sliderInput("trajectorySliderDimensions", "Number of dimensions to use :", min = 0, max = 100, value = 3, step = 1),
                                       selectInput("trajectoryStart", "Initial state:", choices=c("-"="-"), selected = "-", multiple = F, selectize = F),
                                       selectInput("trajectoryEnd", "Final state:", choices=c("-"="-"), selected = "0", multiple = F, selectize = F),
                                       actionButton(inputId = "trajectoryConfirm", label = "Run trajectory analysis!", class="btn-run", icon = icon("check-circle"))
                                     ),

                                   box(
                                     width = 9, status = "info", solidHeader = TRUE, title = "Trajectory analysis results",
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Structure overview", 
                                                          div(id="trajectoryPlot_loader",
                                                              shinycssloaders::withSpinner(
                                                                plotOutput(outputId="trajectoryPlot", height = "1100px")
                                                                )
                                                              ),
                                                          verbatimTextOutput(outputId="trajectoryText")),
                                                 tabPanel("Lineage-Pseudotime view", fluidRow(
                                                   box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                       selectInput(inputId = 'trajectoryLineageSelect',
                                                                   label = 'Select lineage:',
                                                                   choices = c("Lineage1"),
                                                                   selected = "Lineage1",
                                                                   multiple = FALSE),
                                                       actionButton(inputId = "trajectoryConfirmLineage", label = "Select lineage!", class="btn-run", icon = icon("check-circle"))
                                                   ),
                                                   box(width = 9, status = "info", solidHeader = TRUE, title = "Pseudotime plot",
                                                       div(id="trajectoryPseudotimePlot_loader",
                                                           shinycssloaders::withSpinner(
                                                             plotOutput(outputId = "trajectoryPseudotimePlot", height = "1100px")
                                                             )
                                                           )
                                                       )
                                                 )
                                                )
                                   )
                          )
                  )
                ),
                tabPanel("scATAC-seq",
                         fluidRow(
                           box(
                             width = 3, status = "info", solidHeader = TRUE,
                             title = "Trajectory parameters",
                             sliderInput("trajectorySliderDimensionsATAC", "Number of UMAP dimensions to use :", min = 0, max = 100, value = 3, step = 1),
                             selectInput(inputId = "trajectoryGroupByATAC", 
                                         label = "Cells group by: ", 
                                         choices = c("-"="-", 
                                                     "Clusters"="Clusters", 
                                                     "Integration predicted clusters"="predictedGroup_Co"), 
                                         selected = "-"),
                             selectInput("trajectoryStartATAC", "Initial state:", choices=c("-"="-"), selected = "-", multiple = F, selectize = F),
                             selectInput("trajectoryEndATAC", "Final state:", choices=c("-"="-"), selected = "-", multiple = F, selectize = F),
                             actionButton(inputId = "trajectoryConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle"))
                           ),
                           box(
                             width = 9, status = "info", solidHeader = TRUE, 
                             title = "Pseudotime plot",
                             selectInput(inputId = 'trajectoryLineageSelectATAC',
                                         label = 'Select lineage:',
                                         choices = c("Lineage1"),
                                         selected = "Lineage1",
                                         multiple = FALSE),
                             actionButton(inputId = "trajectoryConfirmLineageATAC", label = "Display pseudotime ranking!", class="btn-run", icon = icon("check-circle")),
                             div(id="trajectoryPseudotimePlotATAC_loader",
                                 shinycssloaders::withSpinner(
                                   plotOutput(outputId = "trajectoryPseudotimePlotATAC", height = "1100px")
                                   )
                                 ),
                             verbatimTextOutput(outputId="trajectoryTextATAC")
                           )
                         )
                )
              )
      ),
      
      #L-R analysis
      tabItem(tabName = "ligandReceptor",
              tags$br(),
              tags$div("L-R analysis: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 15sec (per cluster pair)", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
               fluidRow(
                 box(
                   width = 3, status = "info", solidHeader = TRUE,
                   title = "L-R analysis parameters",
                   selectInput("ligandReceptorSender", "Ligand expressing cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                   selectInput("ligandReceptorReciever", "Receptor expressing cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                   actionButton(inputId = "ligandReceptorConfirm", label = "Run ligand-receptor analysis!", class="btn-run", icon = icon("check-circle"))
                 ),
                 box(
                   width = 9, status = "info", solidHeader = TRUE, title = "L-R analysis results",
                   div(
                     tabsetPanel(type = "tabs",
                                 tabPanel("All interactions",
                                          div(id="ligandReceptorFullHeatmap_loader",
                                              shinycssloaders::withSpinner(
                                                plotlyOutput(outputId="ligandReceptorFullHeatmap", height = "1100px")
                                                )
                                              ),
                                          downloadButton(outputId = "ligandReceptorFullExport", label = "Save table")),
                                 tabPanel("Curated interactions (documented in literature and publicly available databases)",
                                          div(id="ligandReceptorCuratedHeatmap_loader",
                                              shinycssloaders::withSpinner(
                                                plotlyOutput(outputId="ligandReceptorCuratedHeatmap", height = "1100px")
                                                )
                                              ),
                                          downloadButton(outputId = "ligandReceptorShortExport", label = "Save table"))
                                 )
                   )
                 )
               )
      ),
      
      #GRN analysis
      tabItem(tabName = "grn",
              bsCollapse(id = 'scenicRNA_collapse', multiple = T,
                         bsCollapsePanel('Do you need help with SCENIC analysis?', ih_scenic_rna, style = 'warning')
              ),
              tags$div("TF activity inference: estimated time in web server for a scRNA-seq dataset of 6,000 cells ~ 4min 05sec", tags$br(),
                       "Detection of positive regulators: estimated time in web server for a scATAC-seq dataset of 6,000 cells ~ 36min 50sec", tags$br(),
                       "(The execution times were measured in the web version of the tool. However, improved performance can be achieved by using the 
                       stand-alone version on PCs with appropriate CPU and RAM specifications.)",
                       class="execTimeMessage"),
              tags$br(),
              tabsetPanel(type="tabs", id = "grnTabPanel",
                          tabPanel("scRNA-seq",
                                   tabsetPanel(type="tabs", 
                                               tabPanel("Transcription factor activity inference",
                                                        fluidRow(
                                                          box(width = 3, status = "info", solidHeader = TRUE, title = "Analysis parameters",
                                                              radioButtons("grnComplexes", label = h3("Complexes:"),
                                                                           choices = list("Keep complexes together (suggested)" = "keep",
                                                                             "Split complexes into subunits" = "split"
                                                                           ),
                                                                           selected = "keep"
                                                                           ),
                                                              radioButtons("grnAnalysisMethod", label = h3("Model used in the analysis:"),
                                                                           choices = list("Multivariate linear model" = "mlm",
                                                                                          "Univariate linear model" = "ulm",
                                                                                          "Weighted sum" = "wsum"
                                                                           ),
                                                                           selected = "ulm",
                                                              ),
                                                              actionButton(inputId = "grnRunDecoupler", label = "Run analysis!", class="btn-run", icon = icon("check-circle"))
                                                              ),
                                                          box(width = 9, status = "info", solidHeader = TRUE, title = "TF activity analysis output",
                                                              dataTableOutput(outputId= "grnMatrixRNA_DecoupleR"),
                                                              downloadButton(outputId = "grnMatrixDecoupleRRNAExport", label = "Save table"),
                                                              tags$hr(),
                                                              plotlyOutput(outputId = "grnHeatmapRNA_DecoupleR", height = "800px")
                                                          )
                                                          )
                                               ),
                                               tabPanel("GRN analysis", 
                                                        fluidRow(
                                                          box(width = 3, status = "info", solidHeader = TRUE, title = "GRN input parameters",
                                                              tags$h3("Prepare files for pyscenic"),
                                                              tags$hr(),
                                                              radioButtons("grnGenomeBuild", label = h3("Select genome build : "),
                                                                           choices = list("Mus musculus (Mouse) - mm10" = "mm10", 
                                                                                          "Homo sapiens (Human) - hg19" = "hg19",
                                                                                          "Homo sapiens (Human) - hg38" = "hg38"
                                                                           ), 
                                                                           selected = "mm10"),
                                                              
                                                              actionButton(inputId = "grnProduceLoom", label = "Prepare .RDS and .loom files for download!", class="btn-run", icon = icon("check-circle")),
                                                              tags$br(),
                                                              div(
                                                                tags$br(),
                                                                downloadButton(outputId = "grnDownloadLoom", label = "Export loom file"),
                                                                downloadButton(outputId = "grnDownloadRDS", label = "Export RDS file"),
                                                              ),
                                                              
                                                              tags$h3("Analyze pyscenic output"),
                                                              tags$hr(),
                                                              fileInput(inputId = "grnLoomInput", label = "Upload a .loom file", accept = ".loom"),
                                                              fileInput(inputId = "grnRDSInput", label = "Upload an .RDS file", accept = ".RDS"),
                                                              actionButton(inputId = "grnLoomAnalysis", label = "Analyze files", class="btn-run", icon = icon("check-circle")),
                                                              
                                                              tags$h3("Visualization options"),
                                                              tags$hr(),
                                                              selectInput(inputId = "grnMatrixSelectionRNA", label = "Regulons - display:", choices = c("Matrix of AUC values"="auc",
                                                                                                                                                        "Matrix of RSS scores"="rss")),
                                                              sliderInput(inputId = "grnTopRegulonsRNA", label = "Display top regulons:", min = 5, max = 100, value = 10, step = 1),
                                                              actionButton(inputId = "grnConfirmVisualizationRNA", label = "Display Plot!", class="btn-run", icon = icon("check-circle"))
                                                          ),
                                                          box(width = 9, status = "info", solidHeader = TRUE, title = "GRN output",
                                                              div( dataTableOutput(outputId="grnMatrixRNA") ),
                                                              downloadButton(outputId = "grnMatrixScenicRNAExport", label = "Save table"),
                                                              tags$br(),
                                                              tags$hr(),
                                                              div(id="grnHeatmapRNA_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "grnHeatmapRNA", height = "800px")
                                                                  )
                                                              )
                                                          )
                                                        )
                                               )
                                   )
                                   
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE, title = "Analysis options",
                                       selectInput("grnGroupByATAC", "Cells group by:",
                                                   c("Clusters" = "Clusters",
                                                     "Integration predicted clusters" = "predictedGroup_Co"
                                                   )),
                                       selectInput("grnMatrixATAC", "Matrix:",
                                                   c("Gene activity score matrix" = "GeneScoreMatrix",
                                                     "Gene Integration matrix" = "GeneIntegrationMatrix"
                                                   )),
                                       sliderInput(inputId = "grnFdrATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.1, step = 0.01),
                                       sliderInput(inputId = "grnCorrlationATTAC", label = "Correllation threshold:", min = 0, max = 1, value = 0.7, step = 0.1),
                                       actionButton(inputId = "grnConfirmATAC", label = "Run analysis!", class="btn-run", icon = icon("check-circle"))
                                     ),
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Gene regulatory networks results",
                                         tabsetPanel(type = "tabs",
                                                     tabPanel("Positive regulators table",
                                                              div(id="grnATACTable_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    dataTableOutput(outputId="grnMatrixATAC")
                                                                  )
                                                              ),
                                                              downloadButton(outputId = "grnPositiveRegulatorsATACExport", label = "Save table"),
                                                     ),
                                                     tabPanel("Positive regulators heatmap (top-10)",
                                                              div(id="grnHeatmapATAC_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    plotlyOutput(outputId = "grnHeatmapATAC", height = "800px")
                                                                  )
                                                              )
                                                     ),
                                                     tabPanel("Peak to gene links",
                                                              div(id="grnATACTable2_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    dataTableOutput(outputId="grnP2GlinksTable"),
                                                                    )
                                                                  ),
                                                              downloadButton(outputId = "grnPeakToGeneLinksATACExport", label = "Save table")
                                                              ),
                                                     tabPanel("Peak x motif occurence matrix",
                                                              
                                                              div(id="grnATACTable3_loader",
                                                                  shinycssloaders::withSpinner(
                                                                    dataTableOutput(outputId="grnMotifTable"),
                                                                  )
                                                              ),
                                                              downloadButton(outputId = "grnPeakMotifTableATACExport", label = "Save table")
                                                     )
                                         )
                                     )
                                   )
                          )
              )
      ),
      
      #ATAC-seq tracks
      tabItem(tabName = "visualizeTracks", 
              fluidRow(
                box(width = 3, status = "info", solidHeader = TRUE,
                    title = "scATAC-seq tracks options",
                    selectInput("visualizeTracksGroupByATAC", "Cells group by:", 
                                c("Clusters" = "Clusters",
                                  "Integration predicted clusters" = "predictedGroup_Co"
                                )),
                    selectizeInput(inputId = 'visualizeTracksGene',
                                   label = 'Select a gene:',
                                   choices = NULL,
                                   selected = NULL,
                                   multiple = FALSE),
                    sliderInput("visualizeTracksBPupstream", "BP upstream :", min = 100, max = 100000, value = 50000, step = 1000),
                    sliderInput("visualizeTracksBPdownstream", "BP downstream :", min = 100, max = 100000, value = 50000, step = 1000),
                    actionButton(inputId = "visualizeTracksConfirm", label = "Visualize tracks!", class="btn-run", icon = icon("check-circle"))
                    ),
                box(width = 9, status = "info", solidHeader = TRUE,
                    title = "scATAC-seq tracks",
                    div(id="visualizeTracksOutput_loader",
                        shinycssloaders::withSpinner(
                          plotOutput(outputId="visualizeTracksOutput", height = "1100px")
                        )
                    )
                )
              )
            ),
      
      #---------------------------------HELP------------------------------------
      tabItem(tabName = "help",
              fluidRow(
                column(12, 
                       tabsetPanel(
                         tabPanel("Examples",
                                  div(class = "div_container",
                                      examples_help
                                  ), 
                         ),
                         tabPanel("Data Input",
                                  tabsetPanel(type = "tabs",
                                    tabPanel("scRNA-seq: count matrix input",
                                      br(),
                                      file_upload_tab_intro,
                                      file_upload_txt,
                                      br(),
                                      file_upload_tab_new_project
                                    ),
                                    tabPanel("scRNA-seq: 10x files input",
                                      br(),
                                      file_upload_tab_intro,
                                      file_upload_10x,
                                      br(),
                                      file_upload_tab_new_project
                                    ),
                                    tabPanel("scRNA-seq: Seurat RDS file input",
                                      br(),
                                      file_upload_tab_intro,
                                      file_upload_RDS,
                                      br(),
                                      file_upload_tab_new_project
                                    ),
                                    tabPanel("scATAC-seq: arrow file input",
                                      br(),
                                      file_upload_tab_intro,
                                      file_upload_arrow,
                                      br(),
                                      file_upload_tab_new_project
                                    ),
                                    tabPanel("Metadata output",
                                              tabsetPanel(type = "tabs",
                                                            tabPanel("Metadata RNA",
                                                              br(),
                                                              file_upload_metadata_RNA
                                                            ),
                                                            tabPanel("Metadata ATAC",
                                                              br(),
                                                              file_upload_metadata_ATAC
                                                            )
                                                          )
                                             )
                                  )
                         ),
                         tabPanel("Quality Control",
                                  tabsetPanel(type = "tabs",
                                    tabPanel("scRNA-seq QC: prior-filtering",
                                             br(),
                                             qc_tab_intro,
                                             rna_qc
                                    ),
                                    tabPanel("scRNA-seq QC: post-filtering",
                                             br(),
                                             qc_tab_intro,
                                             rna_qc_pf
                                    ),
                                    tabPanel("scATAC-seq QC: soft-filtering",
                                      br(),
                                      qc_tab_intro,
                                      atac_qc
                                    )
                                  )
                         ),
                         tabPanel("Normalization",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("Normalization and scaling options",
                                                       br(),
                                                       norm_tab_intro,
                                                       rna_normalization_param
                                              ),
                                              tabPanel("Most variable genes",
                                                       br(),
                                                       norm_tab_intro,
                                                       rna_normalization_output
                                              )
                                  )     
                         ),
                         tabPanel("PCA/LSI",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Optimal number of PCs",
                                                       br(),
                                                       pca_tab_intro,
                                                       br(),
                                                       pca_optimal_pcs
                                              ),
                                              tabPanel("scRNA-seq: Exploration of PCs",
                                                       br(),
                                                       pca_tab_intro,
                                                       br(),
                                                       pca_explore_pcs
                                              ),
                                              tabPanel("scATAC-seq: LSI",
                                                       br(),
                                                       pca_tab_intro,
                                                       br(),
                                                       pca_lsi
                                              )
                                  )
                         ),
                         tabPanel("Clustering",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Clustering parameters",
                                                       br(),
                                                       clustering_tab_intro,
                                                       br(),
                                                       clustering_rna_input
                                              ),
                                              tabPanel("scRNA-seq: Clustering output",
                                                       br(),
                                                       clustering_tab_intro,
                                                       br(),
                                                       clustering_rna_output
                                              ),
                                              tabPanel("scATAC-seq: Clustering parameters",
                                                       br(),
                                                       clustering_tab_intro,
                                                       br(),
                                                       clustering_atac_input
                                              ),
                                              tabPanel("scATAC-seq: Clustering output",
                                                       br(),
                                                       clustering_tab_intro,
                                                       br(),
                                                       clustering_atac_output
                                              )
                                  )
                                  
                         ),
                         tabPanel("Additional dimensionality reduction methods",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Input parameters",
                                                       br(),
                                                       umap_tab_intro,
                                                       br(),
                                                       umap_rna_input
                                              ),
                                              tabPanel("scRNA-seq: Visualization",
                                                       br(),
                                                       umap_tab_intro,
                                                       br(),
                                                       umap_rna_output
                                              ),
                                              tabPanel("scATAC-seq: Input parameters",
                                                       br(),
                                                       umap_tab_intro,
                                                       br(),
                                                       umap_atac_input
                                              ),
                                              tabPanel("scATAC-seq: Visualization",
                                                       br(),
                                                       umap_tab_intro,
                                                       br(),
                                                       umap_atac_output
                                              )
                                  )    
                                  
                         ),
                         tabPanel("Markers identification analysis",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Marker genes",
                                                       br(),
                                                       dea_tab_intro,
                                                       br(),
                                                       dea_rna_input
                                              ),
                                              tabPanel("scRNA-seq: Feature and signature visualization",
                                                       br(),
                                                       dea_tab_intro,
                                                       br(),
                                                       dea_rna_signature
                                              ),
                                              tabPanel("scATAC-seq: Marker genes",
                                                       br(),
                                                       dea_tab_intro,
                                                       br(),
                                                       dea_atac_genes
                                              ),
                                              tabPanel("scATAC-seq: Marker peaks",
                                                       br(),
                                                       dea_tab_intro,
                                                       br(),
                                                       dea_atac_peaks
                                              ),
                                              tabPanel("scATAC-seq: Gene activity score",
                                                       br(),
                                                       dea_tab_intro,
                                                       br(),
                                                       dea_atac_activity
                                              )
                                  )    
                                  
                         ),
                         tabPanel("Doublets' detection",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Doublets' detection",
                                                       br(),
                                                       doublets_tab_intro,
                                                       br(),
                                                       doublets_tab_rna
                                              ),
                                              tabPanel("scATAC-seq: Doublets' detection",
                                                       br(),
                                                       doublets_tab_intro,
                                                       br(),
                                                       doublets_tab_atac       
                                              )
                                  )
                         ),
                         tabPanel("Cell Cycle phase",
                                  br(),
                                  cellCycle_tab_intro,
                                  br(),
                                  cell_cycle_rna
                         ),
                         tabPanel("Functional/Motif Enrichment",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Functional enrichment analysis",
                                                br(),
                                                functional_tab_intro,
                                                br(),
                                                grpofiler_tab_rna
                                              ),
                                              tabPanel("scATAC-seq: Motif enrichment analysis",
                                                br(),
                                                functional_tab_intro,
                                                br(),
                                                motif_tab_atac       
                                              )
                                            )
                         ),
                         tabPanel("Cluster Annotation",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Cluster annotation",
                                                       br(),
                                                       annot_tab_intro,
                                                       br(),
                                                       annot_cipr_rna
                                              ),
                                              tabPanel("scATAC-seq: Integration with scRNA-seq",
                                                       br(),
                                                       annot_tab_intro,
                                                       br(),
                                                       annot_atac
                                              )
                                  )
                         ),
                         tabPanel("Trajectory Inference",
                                  tabsetPanel(type = "tabs",
                                              tabPanel("scRNA-seq: Trajectory inference analysis",
                                                       br(),
                                                       traj_tab_intro,
                                                       br(),
                                                       traj_rna_slingshot
                                              ),
                                              tabPanel("scATAC-seq: Trajectory inference analysis",
                                                       br(),
                                                       traj_tab_intro,
                                                       br(),
                                                       traj_atac_slingshot
                                              )
                                  )
                         ),
                         tabPanel("Ligand-Receptor Analysis",
                                 br(),  
                                 lr_tab_intro,
                                 br(),
                                 lr_rna_nichnet
                         ),
                         tabPanel("Gene regulatory networks analysis",
                                  tabsetPanel(type="tabs",
                                    tabPanel("scRNA-seq: GRN inference analysis",
                                             br(),
                                             grn_tab_intro,
                                             br(),
                                             grn_tab_rna
                                    ),
                                    tabPanel("scATAC-seq: GRN inference analysis",
                                             br(),
                                             grn_tab_intro,
                                             br(),
                                             grn_tab_atac
                                    )
                                  )
                         ),
                         tabPanel("Tracks",
                                  br(),
                                  tracks_tab_intro,
                                  br(),
                                  tracks_tab_atac
                                  ),
                         tabPanel("Utilities",
                                  br(),
                                  utility_tab_intro,
                                  br(),
                                  utility_tab_rna
                         )
                       )
                )
              ) #fluidRow end
      ),
      
      tabItem (tabName = "about",
               div(id = "about_div", class = "div_container",
                   h1(class = "container_title", "About SCALA"),
                   HTML(" 
                              <hr>
                              <h2 class=sub_title> Research team </h2>
                              <ul>
                              <li> Christos Tzaferis, tzaferis[at]fleming[dot]com
                              <li> Evangelos Karatzas, karatzas[at]fleming[dot]gr
                              <li> Fotis Baltoumas, baltoumas[at]fleming[dot]gr
                              <li> Georgios A. Pavlopoulos, pavlopoulos[at]fleming[dot]gr 
                              <li> George Kollias, kollias[at]fleming[dot]gr
                              <li> Dimitris Konstantopoulos, konstantopoulos[at]fleming[dot]gr
                              </ul>
                              <footer>
                              
                              <h2 class=sub_title> Developers </h2>
                              <ul>
                                <li> Christos Tzaferis, tzaferis[at]fleming[dot]com
                                <li> Dimitris Konstantopoulos, konstantopoulos[at]fleming[dot]gr
                              </ul>
                              
                              <h3 class=sub_title>Code Availability</h3>
                              <p>The source code for SCALA can be found in <a href='https://github.com/PavlopoulosLab/SCALA/' target='_blank'>this</a> repository.</p>
                              
                              <h3 class=sub_title> Cite SCALA </h3>
                              <p style='font-size:15px'>If you find SCALA useful in your work please cite: </br>Tzaferis C., Karatzas E., Baltoumas F.A., Pavlopoulos G.A., 
                              Kollias G., Konstantopoulos D. (2022) <b>SCALA: A web application for multimodal analysis of single cell next generation sequencing data.</b> 
                              <i>bioRxiv 2022.11.24.517826; doi: <a href='https://doi.org/10.1101/2022.11.24.517826' target='_blank'>https://doi.org/10.1101/2022.11.24.517826</a></i></p>
                              
                              &copy;"), 
                              sprintf("%s", YEAR), 
                              HTML("
                              <a href=\"https://fleming.gr/kollias-lab-single-cell-analysis-unit\" target=\"_blank\">Single Cell Analysis Unit</a> | 
                              <a href=\"https://sites.google.com/site/pavlopoulossite\" target=\"_blank\">Bioinformatics and Integrative Biology Lab</a> | 
                              <a href=\"https://www.fleming.gr\" target=\"_blank\">Biomedical Sciences Research Center \"Alexander Fleming\"</a>
                              </footer>"
                                   )     
                   )
               )
      )
      
    )# tab item list
  )
