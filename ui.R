# TODO add loading bar divs and source js file
# TODO version2 add all additional parameters per function
source("help.R", local=TRUE)
library(bsplus)

ui <- dashboardPage(
  title="scAnner",
  skin = "black",
  #------------------------------------------------------------Header
  dashboardHeader(
    titleWidth = "280px",
    title = tags$img(src='logo.png')
  ),
  
  #------------------------------------------------------------Sidebar
  dashboardSidebar(
    width = "280px",
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
      menuItem(tags$div("ADDITIONAL DIMENSIONALITY",
                        tags$br(),
                        "REDUCTION METHODS", class = "menu_item_div"), tabName = "umap", icon = icon("draw-polygon")),
      menuItem(text = "MARKERS' IDENTIFICATION", tabName = "findMarkers", icon = icon("map-marker-alt")),
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
                  h1(class = "container_title", "Welcome to scAnner"),
                  HTML("<p class=container_text> scAnner is a web application and stand-alone toolkit, that handles the analysis of scRNA-seq and scATAC-seq datasets, 
                  from quality control and normalization, to dimensionality reduction, differential expression/accessibility analysis, cell clustering, functional enrichment analysis, 
                  trajectory inference, ligand – receptor analysis, gene regulatory network inference, and visualization. Try out our sample data and visit the Help pages for guidance. </p>"
                  ),
                  actionButton(inputId = "debugRNA", label = "Fast debug RNA"),
                  actionButton(inputId = "debugATAC", label = "Fast debug ATAC")
              )
      ),
      
      #Upload tab
      tabItem(tabName = "upload", #TODO 3 tabs
              # two boxes inside upload tab
              fluidRow(
                #useShinyalert(),
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "Upload your data",
                  tabsetPanel(type = "tabs",
                              tabPanel("Gene-count matrix (scRNA-seq)",
                                       tags$h3("Upload your file"),
                                       tags$hr(),
                                       textInput(inputId = "uploadCountMatrixprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "countMatrix", label = "1. Genes-Cells count matrix", accept = ".txt"),
                                       sliderInput(inputId = "uploadCountMatrixminCells", label = "Include features detected in at least this many cells :", min = 1, max = 20, value = 3, step = 1),
                                       sliderInput(inputId = "uploadCountMatrixminFeatures", label = "Include cells where at least this many features are detected :", min = 1, max = 200, value = 200, step = 1),
                                       radioButtons("uploadCountMatrixRadioSpecies", label = h3("Select organism : "),
                                                    choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                   "Homo sapiens (Human)" = "human"
                                                    ), 
                                                    selected = "mouse"),
                                       actionButton(inputId = "uploadCountMatrixConfirm", label = "Submit"),
                                       tags$h3("Load PBMC 10x dataset (example scRNA-seq)"),
                                       tags$hr(),
                                       actionButton(inputId = "upload10xExampleRNACountMatrixConfirm", label = "Load example"),
                                       tags$h3("Export working object as .RDS file"),
                                       tags$hr(),
                                       downloadButton(outputId = "utilitiesConfirmExport1", label = "Export .RDS"),
                                       ),
                              tabPanel("10x input files (scRNA-seq)", 
                                       tags$h3("Upload your files"),
                                       tags$hr(),
                                       textInput(inputId = "upload10xRNAprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "barcodes", label = "1. Choose barcodes.csv.gz file", accept = ".gz"),
                                       fileInput(inputId = "genes", label = "2. Choose features.csv.gz file", accept = ".gz"),
                                       fileInput(inputId = "matrix", label = "3. Choose matrix.mtx.gz file", accept = ".gz"),
                                       sliderInput(inputId = "upload10xRNAminCells", label = "Include features detected in at least this many cells :", min = 1, max = 20, value = 3, step = 1),
                                       sliderInput(inputId = "upload10xRNAminFeatures", label = "Include cells where at least this many features are detected :", min = 1, max = 200, value = 200, step = 1),
                                                    radioButtons("upload10xRNARadioSpecies", label = h3("Select organism : "),
                                                                 choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                                "Homo sapiens (Human)" = "human"
                                                                 ), 
                                                                 selected = "mouse"),
                                       actionButton(inputId = "upload10xRNAConfirm", label = "Submit"),
                                       tags$h3("Load PBMC 10x dataset (example scRNA-seq)"),
                                       tags$hr(),
                                       actionButton(inputId = "upload10xExampleRNA10xFilesConfirm", label = "Load example"),
                                       tags$h3("Export working object as .RDS file"),
                                       tags$hr(),
                                       downloadButton(outputId = "utilitiesConfirmExport2", label = "Export .RDS"),
                                       ),
                              tabPanel("Arrow input files (scATAC-seq)", 
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
                                       sliderInput(inputId = "upload10xATACThreads", label = "Threads to be used:", min = 1, max = 200, value = 2, step = 1), #TODO floor(parallel::detectCores()/2)
                                       actionButton(inputId = "upload10xATACConfirm", label = "Submit"),
                                       tags$h3("Load an example dataset"),
                                       tags$hr(),
                                       actionButton(inputId = "upload10xExampleATACConfirm", label = "Load PBMC 10x dataset (example scATAC-seq)")
                              )
                  )
                ),
                box(
                  width = 8, solidHeader = TRUE, status = "info",
                  title = "Metadata table",
                  div(class="ldBar", id="input_loader", "data-preset"="circle"),
                  
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
      
      #QC tab
      tabItem(tabName = "qc",
              tabsetPanel(type = "tabs", id = "qcTabPanel",
                          tabPanel("scRNA-seq",
                                   #two boxes inside QC tab
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Quality control",
                                       tags$h3("1. Display quality control plots before filtering"),
                                       actionButton(inputId = "qcDisplay", label = "Display plots"),
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
                                       actionButton(inputId = "qcConfirm", label = "Perform filtering"),
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE,
                                       title = "Quality control plots",

                                       div(class="ldBar", id="qc_loader", "data-preset"="circle"),
                                       tabsetPanel(type="tabs", id = "qc_tabs_rna",
                                                   tabPanel("Pre-filtering plots",
                                                            #column(tags$h3("Pre-filtering plots"), width=12),
                                                            #column(tags$hr(), width = 12),
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
                                                            #column(tags$h3("Post-filtering plots"), width=12),
                                                            #column(tags$hr(), width = 12),
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
                                   #two boxes inside QC tab
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Quality control",
                                       tags$h3("Display soft filtered quality control plots"),
                                       actionButton(inputId = "qcDisplayATAC", label = "OK"),
                                       
                                       #TODO list of chromosomes to be analysed
                                       #TODO mark doublets (ATAC, RNA) incorporate after tSNE, UMAP generation and visualize on UMAP, output .csv - exclude? [optional]
                                     ),
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE,
                                       title = "Quality control plots",
                                       div(class="ldBar", id="qc_loader3", "data-preset"="circle"),
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
              fluidRow(
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "Normalize and scale the data",
                  tags$h3("1. Log-normalization"),
                  tags$hr(),
                  #textInput(inputId = "normScaleFactor", label = "Scale factor :", value = "10000"), 
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
                  #textInput(inputId = "nHVGs", label = "Number of genes to select as top variable genes (applicable only to the first and third option) :", value = "2000"),
                  sliderInput(inputId = "nHVGs", label = "Number of genes to select as top variable genes (applicable only to the first and third option) :", min = 200, max = 4000, value = 2000, step = 100),
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
                  
                  actionButton(inputId = "normalizeConfirm", label = "OK"),
                ),
                box(
                  width = 8, status = "info", solidHeader = TRUE,
                  title = "Highly variable genes",
                  div(class="ldBar", id="normalize_loader", "data-preset"="circle"),
                  div(id="hvgScatter_loader",
                      shinycssloaders::withSpinner(
                        plotlyOutput(outputId = "hvgScatter", height = "800px")
                      )
                  )
                ),  
              )
      ),
      
      #PCA tab
      tabItem(tabName = "pca", 
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
                                                            #column(sliderInput(inputId = "pcaStepBy", label = "Resolution/step-by: (applicable only in SVA-CV)", min = 1, max = 5, value = 3, step = 1), width = 12),
                                                            column(actionButton(inputId = "PCrunPCA", label = "Run PCA"), width = 12),
                                                            div(class="ldBar", id="PCA1_loader", "data-preset"="circle"),
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
                                                            column(actionButton(inputId = "PCconfirm", label = "OK"), width = 12),
                                                            div(class="ldBar", id="PCA2_loader", "data-preset"="circle"),
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
                                                            )
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
                                       actionButton(inputId = "lsiConfirm", label = "Run LSI"),
                                       tags$hr(),
                                       div(class="ldBar", id="lsi_loader", "data-preset"="circle"),
                                       verbatimTextOutput(outputId = "lsiOutput")
                                     )
                                   )
                          )
              )
      ),
      #ATAC 
      #input: varFeatures[5000, 100000] step=1000, default=25000, dimsToUse[1, 50] default=30, resolution[0, 10] default=2, iterations [1, 10] default=2
      #visualize LSI?
      
      #Clustering tab
      tabItem(tabName = "clustering", 
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
                                       sliderInput(inputId = "clusterRes", label = "Clustering resolution :", min = 0.1, max = 60, value = 0.5, step = 0.1),
                                       sliderInput(inputId = "clusterPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                       actionButton(inputId = "snnConfirm", label = "Perform clustering"),
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Clustering output", #height = "1500px",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Clustering results",
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Cluster table",
                                                                                 div(class="ldBar", id="clust1_loader", "data-preset"="circle"),
                                                                                 dataTableOutput(outputId="clusterTable"),
                                                                                 downloadButton(outputId = "clusterTableRNAExport", label = "Save table")
                                                                        ),
                                                                        tabPanel("Cluster barplot",
                                                                                 selectInput("clusterGroupBy", "Grouping variable:",
                                                                                             c("orig.ident" = "orig.ident")),
                                                                                 actionButton(inputId = "clusterBarplotConfirm", label = "Display barchart"),
                                                                                 div(class="ldBar", id="clust2_loader", "data-preset"="circle"),
                                                                                 div(id="clusterBarplot_loader",
                                                                                     shinycssloaders::withSpinner(
                                                                                       plotlyOutput(outputId = "clusterBarplot", height = "700px")
                                                                                     )
                                                                                 )
                                                                        )
                                                            )															
                                                   ),
                                                   tabPanel("Shared Nearest Neighbour (SNN) Graph", 
                                                            div(class="ldBar", id="clust3_loader", "data-preset"="circle"),
                                                            actionButton(inputId = "snnDisplayConfirm", label = "Display SNN graph"),
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
                                       actionButton(inputId = "clusterConfirmATAC", label = "Perform clustering"),
                                     ),
                                     box(
                                       width = 8, status = "info", solidHeader = TRUE, title = "Clustering output",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Clustering results", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Cluster table",
                                                                                 div(class="ldBar", id="clust4_loader", "data-preset"="circle"),
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
      #ATAC 
      # resolution, dimsToUse
      # correlation heatmap (RNA, ATAC) using most variable features/or markers [optional]

      #UMAP tab
      tabItem(tabName = "umap", 
              tabsetPanel(type = "tabs", id = "umapTabPanel",
                          tabPanel("scRNA-seq",
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE,
                                         title = "Cells visualization options in reduced space",
                                         sliderInput(inputId = "umapPCs", label = "Number of principal components to use :", min = 1, max = 100, value = 10, step = 1),
                                         sliderInput(inputId = "umapOutComponents", label = "Number of dimensions to fit output:", min = 2, max = 100, value = 3, step = 1)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "If PHATE is selected, the runtime increases when a value > 3 is used.\nPlease note that tSNE doesn't return more than 3 dimensions.", placement = "bottom"
                                               )
                                           ),
                                         actionButton(inputId = "umapRunUmap", label = "Run UMAP"),
                                         actionButton(inputId = "umapRunTsne", label = "Run tSNE"),
                                         actionButton(inputId = "umapRunDFM", label = "Run Diffusion Map"),
                                         actionButton(inputId = "umapRunPhate", label = "Run PHATE"),
                                         tags$h3("Display settings"),
                                         tags$hr(),
                                         div(class="ldBar", id="dim_red1_loader", "data-preset"="circle"),
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
                                         actionButton(inputId = "umapConfirm", label = "Update plot")
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                                         div(class="ldBar", id="dim_red2_loader", "data-preset"="circle"),
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
                                         sliderInput(inputId = "umapDimensionsATAC", label = "Number of input dimensions to use :", min = 1, max = 100, value = 15, step = 1),
                                         sliderInput(inputId = "umapOutComponentsATAC", label = "Number of dimensions to fit output:", min = 2, max = 100, value = 3, step = 1)%>%
                                           shinyInput_label_embed(
                                             shiny_iconlink() %>%
                                               bs_embed_popover(
                                                 title = "Please note that tSNE doesn't return more than 2 dimensions.", placement = "bottom"
                                               )
                                           ),
                                         actionButton(inputId = "umapRunUmapTsneATAC", label = "Run UMAP and tSNE"),
                                         tags$h3("Display settings"),
                                         tags$hr(),
                                         div(class="ldBar", id="dim_red3_loader", "data-preset"="circle"),
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
                                         actionButton(inputId = "umapConfirmATAC", label = "Display plot")
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                                         div(class="ldBar", id="dim_red4_loader", "data-preset"="circle"),
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
      #ATAC
      #UMAP: n_components, ndimsToUse, distance=cosine, neighbors=30, minDist=0.5
      #tSNE: n_components, ndimsToUse, perplexity =30, check for 3D
      
      #DEA tab
      tabItem(tabName = "findMarkers", 
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
                                         actionButton(inputId = "findMarkersConfirm", label = "OK")
                                     ),
                                     
                                     box(
                                       width = 9, status = "info", solidHeader = TRUE, title = "DEA results", #height = "1500px",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Marker genes", 
                                                            div(class="ldBar", id="DEA1_loader", "data-preset"="circle"),
                                                            dataTableOutput(outputId="findMarkersTable"),
                                                            downloadButton(outputId = "findMarkersRNAExport", label = "Save table")),
                                                   tabPanel("Heatmap", 
                                                            div(class="ldBar", id="DEA2_loader", "data-preset"="circle"),
                                                            actionButton(inputId = "findMarkersTop10HeatmapConfirm", label = "Display top-10 marker genes heatmap"),
                                                            div(id="findMarkersHeatmap_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "findMarkersHeatmap", height = "1300px")
                                                                  )
                                                                )
                                                            ),
                                                   
                                                   tabPanel("Dotplot", 
                                                            div(class="ldBar", id="DEA3_loader", "data-preset"="circle"),
                                                            actionButton(inputId = "findMarkersTop10DotplotConfirm", label = "Display top-10 marker genes dotplot"),
                                                            div(id="findMarkersDotplot_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "findMarkersDotplot", height = "1300px")
                                                                )
                                                            )
                                                   ),
                                                   
                                                   tabPanel("Feature plot", fluidRow(
                                                     box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                         radioButtons("findMarkersFeatureSignature", label = "Select between gene or signature to plot: ",
                                                                      choices = list("Gene" = "gene", 
                                                                                     "Gene signature" = "signature"
                                                                      ), 
                                                                      selected = "gene"),
                                                         selectizeInput(inputId = 'findMarkersGeneSelect',
                                                                        label = 'Select a gene:',
                                                                        choices = NULL,
                                                                        selected = NULL,
                                                                        multiple = FALSE),
                                                         selectizeInput(inputId = 'findMarkersSignatureSelect',
                                                                        label = 'Select signature/numeric variable:',
                                                                        choices = "-",
                                                                        selected = "-",
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
                                                         actionButton(inputId = "findMarkersFPConfirm", label = "Display plot"),
                                                         tags$hr(),
                                                         tags$h3("Add a new signature"),
                                                         textInput(inputId = "findMarkersSignatureName", label = "Gene signature name :", value = "Signature1"),
                                                         textAreaInput(inputId = "findMarkersSignatureMembers", label = "Paste a list of genes", cols = 80, rows = 15, placeholder = "Prg4\nTspan15\nCol22a1\nHtra4"),
                                                         actionButton(inputId = "findMarkersSignatureAdd", label = "Calculate signature score")
                                                     ),
                                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                         div(class="ldBar", id="DEA4_loader", "data-preset"="circle"),
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
                                                         actionButton(inputId = "findMarkersFeaturePairConfirm", label = "Display plot")
                                                     ),
                                                     box(width=9, status="info", solidHeader=TRUE, title="Plot",
                                                         div(class="ldBar", id="DEA5_loader", "data-preset"="circle"),
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
                                                         #textInput(inputId = "findMarkersGeneSelect", label = "Search for gene:", value = "Actb"),
                                                         radioButtons("findMarkersViolinFeaturesSignature", label = "Select between gene or signature: ",
                                                                      choices = list("Gene" = "gene", 
                                                                                     "Gene signature" = "signature"
                                                                      ), 
                                                                      selected = "gene"),
                                                         selectizeInput(inputId = 'findMarkersGeneSelect2',
                                                                        label = 'Search for gene:',
                                                                        choices = NULL,
                                                                        selected = NULL,
                                                                        multiple = FALSE), # allow for multiple inputs
                                                         selectizeInput(inputId = 'findMarkersViolinSignatureSelect',
                                                                        label = 'Select signature/numeric variable:',
                                                                        choices = "-",
                                                                        selected = "-",
                                                                        multiple = FALSE),
                                                         actionButton(inputId = "findMarkersViolinConfirm", label = "Display plot")
                                                     ),
                                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                         div(class="ldBar", id="DEA6_loader", "data-preset"="circle"),
                                                         div(id="findMarkersViolinPlot_loader",
                                                             shinycssloaders::withSpinner(
                                                               plotlyOutput(outputId = "findMarkersViolinPlot", height = "800px")
                                                             )
                                                         )
                                                     )
                                                   )),
                                                   tabPanel("VolcanoPlot", fluidRow(
                                                     box(width = 3, status = "info", solidHeader = TRUE, title = "Cluster selection",
                                                         selectInput("findMarkersClusterSelect", "Cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                                                         actionButton(inputId = "findMarkersVolcanoConfirm", "Display volcano plot")
                                                         ),
                                                     
                                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Volcano plot",
                                                         div(class="ldBar", id="DEA7_loader", "data-preset"="circle"),
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
                                         sliderInput(inputId = "findMarkersLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.05),
                                         
                                         sliderInput(inputId = "findMarkersFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         actionButton(inputId = "findMarkersConfirmATAC", label = "OK"),
                                         
                                         tags$h3("Marker peaks"),
                                         tags$hr(),
                                         selectInput("findMarkersPeaksTestATAC", "Test used:",
                                                     c("Wilcoxon" = "wilcoxon",
                                                       "Binomial" = "binomial",
                                                       "T-test" = "ttest"
                                                     )),
                                         sliderInput(inputId = "findMarkersPeaksLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.05),
                                         
                                         sliderInput(inputId = "findMarkersPeaksFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         fileInput(inputId = "findMarkersPeaksCustomPeaks", label = "Please upload a .bed file", accept = ".bed"),
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
                                         
                                         actionButton(inputId = "findMarkersPeaksConfirmATAC", label = "OK"),
                                     ),
                                     
                                     box(
                                       tabsetPanel(type = "tabs", id = "ATAC_markers_tabs",
                                                   tabPanel("Marker genes (ATAC)", fluidRow(
                                                     tabsetPanel(type = "tabs", id = "marker_genes_tab_id",
                                                                 tabPanel("Marker genes table",
                                                                          div(class="ldBar", id="DEA8_loader", "data-preset"="circle"),
                                                                          
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
                                                                div(class="ldBar", id="DEA9_loader", "data-preset"="circle"),
                                                                
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
                                                   ),
                                                   tabPanel("Gene-score plot", fluidRow(
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
                                                         actionButton(inputId = "findMarkersFPConfirmATAC", label = "Display plot"),
                                                     ),
                                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Plot",
                                                         div(class="ldBar", id="DEA10_loader", "data-preset"="circle"),
                                                         div(id="findMarkersFeaturePlotATAC_loader",
                                                             shinycssloaders::withSpinner(
                                                               plotOutput(outputId = "findMarkersFeaturePlotATAC", height = "1100px")
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
      
      #Cell cycle phase analysis
      tabItem(tabName = "cellCycle",
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = T,
                  title = "Cell cycle phase analysis",
                  tabsetPanel(type = "tabs",
                              tabPanel("Dimensionality reduction plot", 
                                       selectInput("cellCycleReduction", "Plot type:", # TODO observe this selectbox instead, default value "-" -> do nothing
                                                   c("-" = "-")
                                       ),
                                       actionButton(inputId = "cellCycleRun", label = "Run cell cycle analysis"), # TODO remove this button
                                       div(class="ldBar", id="CC1_loader", "data-preset"="circle"),
                                       div(id="cellCyclePCA_loader",
                                           shinycssloaders::withSpinner(
                                             plotlyOutput(outputId = "cellCyclePCA", height = "700px")
                                           )
                                       )
                              ),
                              tabPanel("Barplot", # TODO this should appear together with cellCyclePCA output, by observed event cellCycleReduction
                                       div(class="ldBar", id="CC2_loader", "data-preset"="circle"),
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
                                         sliderInput("gProfilerSliderLogFC", "Log FC threshold:", min = 0, max = 3, value = 0.25, step = 0.05),
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
                                         actionButton(inputId = "gProfilerConfirm", label = "OK"),
                                         actionButton(inputId = "sendToFlame", label = "Send to Flame")
                                     ),
                                     box(
                                       width = 10, status = "info", solidHeader = TRUE, title = "Enrichment analysis results",
                                       tabsetPanel(type = "tabs",
                                                   tabPanel("Table of functional terms", 
                                                            div(class="ldBar", id="gprof1_loader", "data-preset"="circle"),
                                                            dataTableOutput(outputId = "gProfilerTable"),
                                                            downloadButton(outputId = "gProfilerRNAExport", label = "Save table")),
                                                   tabPanel("Manhatan plot", 
                                                            div(class="ldBar", id="gprof2_loader", "data-preset"="circle"),
                                                            div(id="gProfilerManhatan_loader",
                                                                shinycssloaders::withSpinner(
                                                                  plotlyOutput(outputId = "gProfilerManhatan")
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
                                         sliderInput(inputId = "findMotifsLogFCATAC", label = "Log2FC threshold:", min = 0, max = 3, value = 0.25, step = 0.05),
                                         sliderInput(inputId = "findMotifsFDRATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.01, step = 0.01),
                                         actionButton(inputId = "findMotifsConfirmATAC", label = "OK"),
                                     ),
                                     
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Motif enrichment analysis results",
                                         tabsetPanel(type = "tabs",
                                                     tabPanel("Table of enriched motifs", 
                                                                div(class="ldBar", id="motif_loader", "data-preset"="circle"),
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
                  actionButton(inputId = "annotateClustersConfirm", label = "OK"),
                ),
                box(
                  width = 9, status = "info", solidHeader = TRUE, title = "Cell type annotation",
                  tabsetPanel(type = "tabs",
                              tabPanel("Top-5 hits table", 
                                       div(class="ldBar", id="annot1_loader", "data-preset"="circle"),
                                       dataTableOutput(outputId="annotateClustersCIPRTable"),
                                       downloadButton(outputId = "annotationRNAExport", label = "Save table")
                                       ),
                              tabPanel("Top-5 hits dotplot", 
                                       div(class="ldBar", id="annot2_loader", "data-preset"="circle"),
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
      
      #Trajectory analysis
      tabItem(tabName = "trajectory",
              tabsetPanel(type="tabs", id = "trajectoryTabPanel",
                          tabPanel("scRNA-seq", 
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE,
                                       title = "Trajectory parameters",
                                       selectInput("trajectoryReduction", "Dimensionality reduction method:", choices=c("PCA"="pca","UMAP"="umap", "tSNE"="tsne", "Diffusion Map"="dfm"), selected = "PCA",
                                                   multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                                       sliderInput("trajectorySliderDimensions", "Number of dimensions to use :", min = 0, max = 100, value = 3, step = 1),
                                       selectInput("trajectoryStart", "Initial state:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                                       selectInput("trajectoryEnd", "Final state:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                                       actionButton(inputId = "trajectoryConfirm", label = "OK")
                                     ),

                                   box(
                                     width = 9, status = "info", solidHeader = TRUE, title = "Trajectory analysis results",
                                     tabsetPanel(type = "tabs",
                                                 tabPanel("Structure overview", 
                                                          div(class="ldBar", id="traj1_loader", "data-preset"="circle"),
                                                          div(id="trajectoryPlot_loader",
                                                              shinycssloaders::withSpinner(
                                                                plotOutput(outputId="trajectoryPlot", height = "1100px")
                                                                )
                                                              ),
                                                          verbatimTextOutput(outputId="trajectoryText")),
                                                 tabPanel("Lineage-Pseudotime view", fluidRow(
                                                   box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                                       selectInput(inputId = 'trajectoryLineageSelect', # TODO observe this instead of action button below, default = "-" and handle it
                                                                   label = 'Select lineage:',
                                                                   choices = c("Lineage1"),
                                                                   selected = "Lineage1",
                                                                   multiple = FALSE),
                                                       actionButton(inputId = "trajectoryConfirmLineage", label = "ok") # TODO remove this and observe selectbox above
                                                   ),
                                                   box(width = 9, status = "info", solidHeader = TRUE, title = "Pseudotime plot",
                                                       div(class="ldBar", id="traj2_loader", "data-preset"="circle"),
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
                             selectInput("trajectoryStartATAC", "Initial state:", choices=c("C1"="C1"), selected = "C1", multiple = F, selectize = F),
                             selectInput("trajectoryEndATAC", "Final state:", choices=c("C1"="C1"), selected = "C1", multiple = F, selectize = F),
                             actionButton(inputId = "trajectoryConfirmATAC", label = "OK")
                           ),
                           box(
                             width = 9, status = "info", solidHeader = TRUE, 
                             title = "Pseudotime plot",
                             selectInput(inputId = 'trajectoryLineageSelectATAC', # TODO observe this instead of action button below, default = "-" and handle it
                                         label = 'Select lineage:',
                                         choices = c("Lineage1"),
                                         selected = "Lineage1",
                                         multiple = FALSE),
                             actionButton(inputId = "trajectoryConfirmLineageATAC", label = "Display pseudotime ranking"),
                             div(class="ldBar", id="traj4_loader", "data-preset"="circle"),
                             div(id="trajectoryPseudotimePlotATAC_loader",
                                 shinycssloaders::withSpinner(
                                   plotOutput(outputId = "trajectoryPseudotimePlotATAC", height = "1100px")
                                   )
                                 ),
                             div(class="ldBar", id="traj3_loader", "data-preset"="circle"),
                             verbatimTextOutput(outputId="trajectoryTextATAC")
                           )
                         )
                )
              )
      ),
      
      #L-R analysis
      tabItem(tabName = "ligandReceptor",
               fluidRow(
                 box(
                   width = 3, status = "info", solidHeader = TRUE,
                   title = "L-R analysis parameters",
                   selectInput("ligandReceptorSender", "Ligand expressing cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                   selectInput("ligandReceptorReciever", "Receptor expressing cluster:", choices=c("-"="-"), multiple = F, selectize = F),
                   actionButton(inputId = "ligandReceptorConfirm", label = "OK")
                 ),
                 box(
                   width = 9, status = "info", solidHeader = TRUE, title = "L-R analysis results",
                   div(class="ldBar", id="lr_loader", "data-preset"="circle"),
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
              tabsetPanel(type="tabs", id = "grnTabPanel",
                          tabPanel("scRNA-seq", #TODO for pyscenic ctx minimun number of genes per module, AUC+NES thresholds [for the update]
                                   fluidRow(
                                     box(width = 3, status = "info", solidHeader = TRUE, title = "GRN input parameters",
                                         div(class="ldBar", id="loom_production_loader", "data-preset"="circle"),
                                         tags$h3("Prepare files for pyscenic"),
                                         tags$hr(),
                                         radioButtons("grnGenomeBuild", label = h3("Select genome build : "),
                                                      choices = list("Mus musculus (Mouse) - mm10" = "mm10", 
                                                                     "Homo sapiens (Human) - hg19" = "hg19",
                                                                     "Homo sapiens (Human) - hg38" = "hg38"
                                                      ), 
                                                      selected = "mm10"),
                                         
                                         actionButton(inputId = "grnProduceLoom", label = "Prepare .RDS and .loom files for download"),
                                         tags$br(),
                                         div(
                                           tags$br(),
                                           downloadButton(outputId = "grnDownloadLoom", label = "Export loom file"),
                                           downloadButton(outputId = "grnDownloadRDS", label = "Export RDS file"),
                                         ),
                                         
                                         div(class="ldBar", id="loom_analysis_loader", "data-preset"="circle"),
                                         tags$h3("Analyze pyscenic output"),
                                         tags$hr(),
                                         fileInput(inputId = "grnLoomInput", label = "Upload a .loom file", accept = ".loom"),
                                         fileInput(inputId = "grnRDSInput", label = "Upload an .RDS file", accept = ".RDS"),
                                         actionButton(inputId = "grnLoomAnalysis", label = "Analyze files"),
                                         
                                         tags$h3("Visualization options"),
                                         tags$hr(),
                                         selectInput(inputId = "grnMatrixSelectionRNA", label = "Regulons - display:", choices = c("Matrix of AUC values"="auc",
                                                                                                                                   "Matrix of RSS scores"="rss")),
                                         sliderInput(inputId = "grnTopRegulonsRNA", label = "Display top regulons:", min = 5, max = 100, value = 10, step = 1),
                                         actionButton(inputId = "grnConfirmVisualizationRNA", label = "Plot")
                                     ),
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "GRN output",
                                         div( dataTableOutput(outputId="grnMatrixRNA") ),
                                         tags$br(),
                                         tags$hr(),
                                         div(id="grnHeatmapRNA_loader",
                                             shinycssloaders::withSpinner(
                                               plotlyOutput(outputId = "grnHeatmapRNA", height = "800px")
                                             )
                                         )
                                     )
                                   )
                          ),
                          tabPanel("scATAC-seq",
                                   fluidRow(
                                     box(
                                       width = 3, status = "info", solidHeader = TRUE, title = "Analysis options",
                                       sliderInput(inputId = "grnFdrATAC", label = "FDR threshold:", min = 0, max = 1, value = 0.1, step = 0.01),
                                       sliderInput(inputId = "grnCorrlationATTAC", label = "Correllation threshold:", min = 0, max = 1, value = 0.7, step = 0.1),
                                       actionButton(inputId = "grnConfirmATAC", label = "OK")
                                     ),
                                     box(width = 9, status = "info", solidHeader = TRUE, title = "Gene regulatory networks results",
                                         tabsetPanel(type = "tabs",
                                                     tabPanel("Positive regulators table",
                                                              div(class="ldBar", id="grn2_loader", "data-preset"="circle"),
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
                    selectizeInput(inputId = 'visualizeTracksGene',
                                   label = 'Select a gene:',
                                   choices = NULL,
                                   selected = NULL,
                                   multiple = FALSE),
                    sliderInput("visualizeTracksBPupstream", "BP upstream :", min = 100, max = 100000, value = 50000, step = 1000),
                    sliderInput("visualizeTracksBPdownstream", "BP downstream :", min = 100, max = 100000, value = 50000, step = 1000),
                    actionButton(inputId = "visualizeTracksConfirm", label = "OK")
                    ),
                box(width = 9, status = "info", solidHeader = TRUE,
                    title = "scATAC-seq tracks",
                    div(class="ldBar", id="tracks_loader", "data-preset"="circle"),
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
                         tabPanel("Markers identification analysis"
                                  
                         ),
                         tabPanel("Cell Cycle phase"
                                  
                         ),
                         tabPanel("Functional/Motif Enrichment"
                                  
                         ),
                         tabPanel("Cluster Annotation"
                                  
                         ),
                         tabPanel("Trajectory Analysis"
                                  
                         ),
                         tabPanel("Ligand-Receptor Analysis"
                                  
                         ),
                         tabPanel("Gene regulatory networks analysis"
                                  
                         ),
                         tabPanel("Tracks")
                       )
                )
              ) #fluidRow end
      ),
      
      tabItem (tabName = "about",
               div(id = "about_div", class = "div_container",
                   h1(class = "container_title", "About scAnner"),
                   HTML("<p class=container_text> scAnner is actively developed and maintained by the Single-Cell Analysis Unit, 
                                                  and Bioinformatics and Integrative Biology Lab. </br></br> 
                        </p> 
                              <h2 class=sub_title> Developers </h2>
                              <ul>
                              <li> Christos Tzaferis, tzaferis[at]fleming[dot]com
                              <li> Evangelos Karatzas, karatzas[at]fleming[dot]gr
                              <li> Fotis Baltoumas, baltoumas[at]fleming[dot]gr
                              <li> Georgios A. Pavlopoulos, pavlopoulos[at]fleming[dot]gr 
                              <li> George Kollias, kollias[at]fleming[dot]gr
                              <li> Dimitris Konstantopoulos, konstantopoulos[at]fleming[dot]gr
                              </ul>
                              <footer>
                              &copy; 2021 <a href=\"https://fleming.gr/kollias-lab-single-cell-analysis-unit\" target=\"_blank\">Single Cell Analysis Unit</a> | 
                              <a href=\"https://sites.google.com/site/pavlopoulossite\" target=\"_blank\">Bioinformatics and Integrative Biology Lab</a> | 
                              <a href=\"https://www.fleming.gr\" target=\"_blank\">Biomedical Sciences Research Center \"Alexander Fleming\"</a>
                              </footer>"     
                   )
               )
      )
      
    )# tab item list
  )
)
