# TODO add loading bar divs and source js file
source("help.R", local=TRUE)

ui <- dashboardPage(
  title="scAnner",
  skin = "black",
  #------------------------------------------------------------Header
  dashboardHeader(
    titleWidth = "280px",
    title = tags$img(src='logo.png')
  ),
  #
  #------------------------------------------------------------Sidebar
  dashboardSidebar(
    width = "280px",
    
    sidebarMenu(id = "sidebarMenu",
      menuItem(text = "HOME", tabName = "home", icon = icon("home")),
      tags$hr(),
      menuItem(text = "DATA INPUT", tabName = "upload", icon = icon("upload")),
      menuItem(text = "UTILITY OPTIONS", tabName = "utilities", icon = icon("edit")),
      menuItem(text = "QUALITY CONTROL", tabName = "qc", icon = icon("check-circle")),
      menuItem(tags$div("DATA NORMALIZATION",
                        tags$br(),
                        "& SCALING", class = "menu_item_div"), tabName = "normalize", icon = icon("balance-scale")),
      tags$hr(),
      menuItem(tags$div("PRINCIPAL COMPONENT",
                        tags$br(),
                        "ANALYSIS", class = "menu_item_div"), tabName = "pca", icon = icon("chart-line")),
      menuItem(text = "CLUSTERING", tabName = "clustering", icon = icon("project-diagram")),
      menuItem(tags$div("ADDITIONAL DIMENSIONALITY",
                        tags$br(),
                        "REDUCTION METHODS", class = "menu_item_div"), tabName = "umap", icon = icon("draw-polygon")),
      menuItem(text = "MARKERS' IDENTIFICATION", tabName = "findMarkers", icon = icon("map-marker-alt")),
      menuItem(text = "CELL CYCLE PHASE ANALYSIS", tabName = "cellCycle", icon = icon("circle-notch")),
      tags$hr(),
      menuItem(tags$div("FUNCTIONAL ENRICHMENT",
                        tags$br(),
                        "ANALYSIS", class = "menu_item_div"), tabName = "gProfiler", icon=icon("chart-bar")),
      menuItem(text = "CLUSTERS' ANNOTATION", tabName = "annotateClusters", icon = icon("id-card")),
      menuItem(text = "TRAJECTORY ANALYSIS", tabName = "trajectory", icon = icon("route")),
      menuItem(tags$div("LIGAND - RECEPTOR",
                        tags$br(),
                        "ANALYSIS", class = "menu_item_div"), tabName = "ligandReceptor", icon = icon("satellite-dish")), #icon("satellite-dish"))
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
    tags$head(tags$script(src = "sliderfix.js")),
    useShinyjs(),
    extendShinyjs(text = js.enrich, functions = c("Enrich")),
    tabItems(
      #home tab
      tabItem(tabName = "home", 
              div(id = "home_div", class = "div_container",
                  h1(class = "container_title", "Welcome to scAnner"),
                  HTML("<p class=container_text> This is a web tool that handles the analysis of scRNAseq data, 
                  from quality control and normalization, to dimensionality reduction, differential expression analysis, clustering and visualization.
                  </br> Try out our sample data and visit the Help pages for guidance. </p>"
                  ),
              )
      ),
      
      #Upload tab
      tabItem(tabName = "upload", #TODO 3 tabs
              # two boxes inside upload tab
              fluidRow(
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "Select input files",
                  tabsetPanel(type = "tabs",
                              tabPanel("Gene-count matrix (scRNA-seq)",
                                       textInput(inputId = "uploadCountMatrixprojectID", label = "Project name : ", value = "Project1"),
                                       fileInput(inputId = "countMatrix", label = "1. Genes-Cells count matrix", accept = ".txt"),
                                       sliderInput(inputId = "uploadCountMatrixminCells", label = "Include features detected in at least this many cells :", min = 1, max = 20, value = 3, step = 1),
                                       sliderInput(inputId = "uploadCountMatrixminFeatures", label = "Include cells where at least this many features are detected :", min = 1, max = 200, value = 200, step = 1),
                                       radioButtons("uploadCountMatrixRadioSpecies", label = h3("Select organism : "),
                                                    choices = list("Mus musculus (Mouse)" = "mouse", 
                                                                   "Homo sapiens (Human)" = "human"
                                                    ), 
                                                    selected = "mouse"),
                                       actionButton(inputId = "uploadCountMatrixConfirm", label = "Submit")
                                       ),
                              tabPanel("10x input files (scRNA-seq)", 
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
                                       actionButton(inputId = "upload10xRNAConfirm", label = "Submit")
                                       ),
                              tabPanel("10x input files (scATAC-seq)", 
                                       actionButton(inputId = "upload10xATACConfirm", label = "Submit")
                              ),
                              tabPanel("Load example data", 
                                       actionButton(inputId = "upload10xExampleRNAConfirm", label = "Load PBMC 10x dataset (example scRNA-seq)"),
                                       actionButton(inputId = "upload10xExampleATACConfirm", label = "Load PBMC 10x dataset (example scATAC-seq)"))
                  )
                ),
                box(
                  width = 8, solidHeader = TRUE, status = "info",
                  title = "Metadata table",
                  div(class="ldBar", id="input_loader", "data-preset"="circle"),
                  dataTableOutput("metadataTable")
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
                  actionButton(inputId = "utilitiesConfirmRename", label = "Rename"),
                  tags$h3("Delete cluster"),
                  tags$hr(),
                  selectInput(inputId = "utilitiesDeleteCluster", label = "Cluster to be deleted:", choices = "-", multiple = F),
                  actionButton(inputId = "utilitiesConfirmDelete", label = "Delete"),
                  tags$h3("Export working object as .RDS file"),
                  tags$hr(),
                  downloadButton(outputId = "utilitiesConfirmExport", label = "Export .RDS"),
                  #tags$h3("Select palette for cluster colors"),
                  #tags$hr(),
                  #actionButton(inputId = "utilitiesColorPicker", label = "choose palette")
                  )
                )
              ),
      
      #QC tab
      tabItem(tabName = "qc",
              #two boxes inside QC tab
              fluidRow(
                box(
                  width = 3, status = "info", solidHeader = TRUE,
                  title = "Quality control",
                  tags$h3("1. Display unfiltered quality control plots"),
                  actionButton(inputId = "qcDisplay", label = "OK"),
                  tags$hr(),
                  tags$h3("2. Filter out low quality cells"),
                  tags$hr(),
                  #textInput(inputId = "minUniqueGenes", label = "Filter out cells that have unique feature counts less than :", value = "500"),
                  sliderInput(inputId = "minUniqueGenes", label = "Filter out cells that have unique feature counts less than :", min = 200, max = 2000, value = 500, step = 1),
                  sliderInput(inputId = "maxUniqueGenes", label = "Filter out cells that have unique feature counts over than :", min = 2001, max = 7000, value = 4500, step = 1),
                  #textInput(inputId = "maxUniqueGenes", label = "Filter out cells that have unique feature counts over than :", value = "4000"),
                  sliderInput(inputId = "maxMtReads", label = "Filter out cells with mitochondrial counts % over :", min = 1, max = 100, value = 10, step = 1),
                  #textInput(inputId = "maxMtReads", label = "Filter out cells with mitochondrial counts % over :", value = "10"),
                  selectInput("qcColorBy", "Color by:",
                              c("orig.ident" = "orig.ident")),
                  actionButton(inputId = "qcConfirm", label = "OK"),
                  ),
                box(
                  width = 9, status = "info", solidHeader = TRUE,
                  title = "Quality control plots",
                  div(class="ldBar", id="qc_loader", "data-preset"="circle"),
                  div(
                    column(tags$h3("Pre-filtering plots"), width=12),
                    column(tags$hr(), width = 12),
                    column(plotlyOutput(outputId = "nFeatureViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "totalCountsViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "mitoViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "genesCounts", height= "100%"), width = 6),
                    column(plotlyOutput(outputId = "mtCounts", height= "100%"), width = 6),
                    column(verbatimTextOutput(outputId = "cellStats"), width = 4),
                    column(tags$h3("Post-filtering plots"), width=12),
                    column(tags$hr(), width = 12),
                    column(plotlyOutput(outputId = "filteredNFeatureViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "filteredTotalCountsViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "filteredMitoViolin", height = "100%"), width = 4),
                    column(plotlyOutput(outputId = "filteredGenesCounts", height= "100%"), width = 6),
                    column(plotlyOutput(outputId = "filteredMtCounts", height= "100%"), width = 6),
                    column(verbatimTextOutput(outputId = "filteredCellStats"), width = 4)
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
                  sliderInput(inputId = "normScaleFactor", label = "Scale factor :", min = 1000, max = 1000000, value = 10000, step = 1000),
                  tags$h3("2. Identification of highly variable features"),
                  tags$hr(),
                  radioButtons("radioHVG", label = h3("Select one of the following methods : "),
                               choices = list("Variance Stabilizing Transformation method" = "vst", 
                                              
                                              "Mean-Variance method" = "mvp", 
                                              
                                              "Dispersion method" = "disp"), 
                               selected = "vst"),
                  #textInput(inputId = "nHVGs", label = "Number of genes to select as top variable genes (applicable only to the first and third option) :", value = "2000"),
                  sliderInput(inputId = "nHVGs", label = "Number of genes to select as top variable genes (applicable only to the first and third option) :", min = 200, max = 4000, value = 2000, step = 100),
                  tags$h3("3. Scaling the data"),
                  tags$hr(),
                  tags$p("Scaling transformation is implemented before dimensionality reduction of the dataset and performs the following steps : "),
                  tags$ul(),
                  tags$li("Shifts the expression of each gene, so that the mean expression across cells is 0"),
                  tags$li("Scales the expression of each gene, so that the variance across cells is 1(this step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate)"),
                  tags$br(),
                  selectInput("normalizeRegressColumns", "Select variables to regress out", list(), selected = NULL, multiple = TRUE, selectize = TRUE, width = NULL, size = NULL),
                  
                  actionButton(inputId = "normalizeConfirm", label = "OK"),
                ),
                box(
                  width = 8, status = "info", solidHeader = TRUE,
                  title = "Highly variable genes",
                  div(class="ldBar", id="normalize_loader", "data-preset"="circle"),
                  plotlyOutput(outputId = "hvgScatter", height = "800px")
                ),  
              )
      ),
      
      #PCA tab
      tabItem(tabName = "pca", 
              fluidRow(
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "PCA results", height = "990px",
                  column(actionButton(inputId = "PCrunPCA", label = "Run PCA"), width = 12),
                  column(selectInput("pcaColorBy", "Color by:", c("orig.ident" = "orig.ident")), width = 12),
                  div(class="ldBar", id="PCA1_loader", "data-preset"="circle"),
                  div(
                    column(plotlyOutput(outputId = "elbowPlotPCA", height = "790px"), width = 6),
                    column(plotlyOutput(outputId = "PCAscatter", height = "790px"), width = 6)
                  )
                ),
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Explore particular principal components", height = "990px",
                  #column(textInput(inputId = "PCin", label = "Select a principal component :", value = "1"), width = 6),
                  selectInput("PCin", "Select a principal component :", choices=1:50, selected = 1, multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                  column(actionButton(inputId = "PCconfirm", label = "OK"), width = 12),
                  div(class="ldBar", id="PCA2_loader", "data-preset"="circle"),
                  div(
                    column(plotlyOutput(outputId = "PCAloadings", height = "790px"), width = 6),
                    column(plotlyOutput(outputId = "PCAheatmap", height = "790px"), width = 6)
                  )
                )
              )
      ),
      
      #Clustering tab
      tabItem(tabName = "clustering", 
              fluidRow(
                box(
                  width = 4, status = "info", solidHeader = TRUE,
                  title = "k-NN & Clustering parameters",
                  tags$h3("1. Construction of the shared nearest neighbour"),
                  tags$hr(),
                  #textInput(inputId = "snnK", label = "k-nearest neighbours for each cell :", value = "20"),
                  sliderInput(inputId = "snnK", label = "k-nearest neighbours for each cell :", min = 1, max = 200, value = 20, step = 1),
                  #textInput(inputId = "snnPCs", label = "Number of principal components to use :", value = "15"),
                  sliderInput(inputId = "snnPCs", label = "Number of principal components to use :", min = 1, max = 50, value = 15, step = 1),
                  tags$h3("2. Clustering of the cells"),
                  tags$hr(),
                  #textInput(inputId = "clusterRes", label = "Clustering resolution :", value = "0.6"),
                  sliderInput(inputId = "clusterRes", label = "Clustering resolution :", min = 0.1, max = 60, value = 0.8, step = 0.1),
                  #textInput(inputId = "clusterPCs", label = "Number of principal components to use :", value = "15"),
                  sliderInput(inputId = "clusterPCs", label = "Number of principal components to use :", min = 1, max = 50, value = 15, step = 1),
                  actionButton(inputId = "snnConfirm", label = "OK"),
                ),
                box(
                  width = 8, status = "info", solidHeader = TRUE, title = "k-NN graph & clusters", height = "1500px",
                  tabsetPanel(type = "tabs",
                              tabPanel("Clustering results", 
                                       div(class="ldBar", id="clust1_loader", "data-preset"="circle"),
                                       dataTableOutput(outputId="clusterTable", height = "500px"),
                                       selectInput("clusterGroupBy", "Grouping variable:",
                                                   c("orig.ident" = "orig.ident")),
                                       actionButton(inputId = "clusterBarplotConfirm", label = "Display barchart"),
                                       div(class="ldBar", id="clust2_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "clusterBarplot", height = "700px")),
                              tabPanel("Shared Nearest Neighbour Graph", visNetworkOutput(outputId="snnSNN", height = "1300px"))
                  ),
                ),
              )
      ),
      
      #UMAP tab
      tabItem(tabName = "umap", 
              fluidRow(
                box(width = 3, status = "info", solidHeader = TRUE,
                    title = "Cell visualization options",
                    #textInput(inputId = "umapPCs", label = "Number of principal components to use:", value = "15"),
                    sliderInput(inputId = "umapPCs", label = "Number of principal components to use :", min = 1, max = 50, value = 15, step = 1),
                    #textInput(inputId = "umapOutComponents", label = "Number of principal components to fit output:", value = "3"),
                    sliderInput(inputId = "umapOutComponents", label = "Number of dimensions to fit output:", min = 1, max = 50, value = 15, step = 1),
                    actionButton(inputId = "umapRunUmap", label = "Run UMAP"),
                    actionButton(inputId = "umapRunTsne", label = "Run tSNE"),
                    actionButton(inputId = "umapRunDFM", label = "Run Diffusion Map"),
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
                    sliderInput("umapDotBorder", "Border width:", min = 0, max = 10, value = 1, step = 0.1),
                    actionButton(inputId = "umapConfirm", label = "Display plot")
                    ),

                box(width = 9, status = "info", solidHeader = TRUE, title = "Plot", height = "1200px",
                    div(class="ldBar", id="dim_red2_loader", "data-preset"="circle"),
                    plotlyOutput(outputId = "umapPlot", height = "1100px")
                )
              )
      ),
      
      #DEA tab
      tabItem(tabName = "findMarkers", 
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
                    #textInput(inputId = "findMarkersMinPct", label = "Only test genes that are detected in a minimum fraction of cells in either of the two populations :", value = "0.1"),
                    sliderInput(inputId = "findMarkersMinPct", label = "Only test genes that are detected in a minimum fraction of cells in either of the two populations :", min = 0, max = 1, value = 0.1, step = 0.05),
                    #textInput(inputId = "findMarkersLogFC", label = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells :", value = "0.25"),
                    sliderInput(inputId = "findMarkersLogFC", label = "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells :", min = 0, max = 3, value = 0.25, step = 0.05),
                    #textInput(inputId = "findMarkersPval", label = "Only return markers that have a p-value < slected threshold, or a power > selected threshold (if the test is ROC) :", value = "0.01"),
                    sliderInput(inputId = "findMarkersPval", label = "Only return markers that have a p-value < slected threshold, or a power > selected threshold (if the test is ROC) :", min = 0, max = 1, value = 0.01, step = 0.01),
                    actionButton(inputId = "findMarkersConfirm", label = "OK")
                ),
                
                box(
                  width = 9, status = "info", solidHeader = TRUE, title = "DEA results", height = "1500px",
                  tabsetPanel(type = "tabs",
                              tabPanel("Marker genes", 
                                       div(class="ldBar", id="DEA1_loader", "data-preset"="circle"),
                                       dataTableOutput(outputId="findMarkersTable", height = "1300px")),
                              tabPanel("Heatmap", 
                                       div(class="ldBar", id="DEA2_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "findMarkersHeatmap", height = "1300px")),
                              
                              tabPanel("Dotplot", 
                                       div(class="ldBar", id="DEA3_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "findMarkersDotplot", height = "1300px")),
                              
                              tabPanel("Feature plot", fluidRow(
                                box(width = 3, status = "info", solidHeader = TRUE, title = "Options",
                                    radioButtons("findMarkersFeatureSignature", label = "Select between gene or signature: ",
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
                                    #div(class="ldBar", id="DEA4_loader", "data-preset"="circle"),
                                    plotlyOutput(outputId = "findMarkersFeaturePlot", height = "1300px")
                                    )
                              )),
                              tabPanel("Gene-pair expression", fluidRow(
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
                                    div(
                                    column(plotlyOutput(outputId="findMarkersFPfeature1", height = "650px"), width = 6),
                                    column(plotlyOutput(outputId="findMarkersFPfeature2", height = "650px"), width = 6),
                                    column(plotlyOutput(outputId="findMarkersFPfeature1_2", height = "650px"), width = 6),
                                    column(plotlyOutput(outputId="findMarkersFPcolorbox", height = "650px"), width = 6),
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
                                    #div(class="ldBar", id="DEA5_loader", "data-preset"="circle"),
                                    plotlyOutput(outputId = "findMarkersViolinPlot", height = "1300px")
                                )
                              )),
                              tabPanel("VolcanoPlot", fluidRow(
                                box(width = 3, status = "info", solidHeader = TRUE, title = "Cluster selection",
                                    textInput(inputId = "findMarkersClusterSelect", label = "Cluster:", value = "0")),
                                box(width = 9, status = "info", solidHeader = TRUE, title = "Volcano plot",
                                    div(class="ldBar", id="DEA6_loader", "data-preset"="circle"),
                                    plotlyOutput(outputId = "findMarkersVolcanoPlot", height = "1300px"))
                              ))
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
                              tabPanel("PCA plot", 
                                       actionButton(inputId = "cellCycleRun", label = "Run cell cycle analysis"), # TODO remove this button
                                       selectInput("cellCycleReduction", "Plot type:", # TODO observe this selectbox instead, default value "-" -> do nothing
                                                   c("PCA" = "pca",
                                                     "UMAP" = "umap",
                                                     "tSNE" = "tsne"), selected ="pca"),
                                       div(class="ldBar", id="CC1_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "cellCyclePCA", height = "1100px")),
                              tabPanel("Barplot", # TODO this should appear together with cellCyclePCA output, by observed event cellCycleReduction
                                       div(class="ldBar", id="CC2_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "cellCycleBarplot", height = "1100px")))
                )
              )
      ),
      
      #Enrichment analysis -gProfiler
      tabItem(tabName = "gProfiler", 
              fluidRow(
                box(width = 2, status = "info", solidHeader = TRUE,
                    title = "gProfiler options",
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
                    sliderInput("gProfilerSliderLogFC", "Log FC threshold:", min = 0, max = 3, value = 0.35, step = 0.05),
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
                  width = 10, status = "info", solidHeader = TRUE, title = "gProfiler results",
                  tabsetPanel(type = "tabs",
                              tabPanel("Table of functional terms", 
                                       div(class="ldBar", id="gprof1_loader", "data-preset"="circle"),
                                       dataTableOutput(outputId = "gProfilerTable")),
                              tabPanel("Manhatan plot", 
                                       div(class="ldBar", id="gprof2_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId = "gProfilerManhatan"))
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
                                              "Hematopoietic diff (human)" = "hsrnaseq"
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
                               selected = "logfc_dot_product"
                  ),
                  tags$hr(),
                  actionButton(inputId = "annotateClustersConfirm", label = "OK"),
                ),
                box(
                  width = 9, status = "info", solidHeader = TRUE, title = "Cell type annotation",
                  tabsetPanel(type = "tabs",
                              tabPanel("Top-5 hits table", 
                                       div(class="ldBar", id="annot1_loader", "data-preset"="circle"),
                                       dataTableOutput(outputId="annotateClustersCIPRTable")),
                              tabPanel("Top-5 hits dotplot", 
                                       div(class="ldBar", id="annot2_loader", "data-preset"="circle"),
                                       plotlyOutput(outputId="annotateClustersCIPRDotplot", height = "1100px"))
                  )
                )
              )
      ),
      
      #Trajectory analysis
      tabItem(tabName = "trajectory",
              fluidRow(
                box(
                  width = 3, status = "info", solidHeader = TRUE,
                  title = "Trajectory parameters",
                  selectInput("trajectoryReduction", "Dimensionality reduction method:", choices=c("PCA"="pca","UMAP"="umap", "tSNE"="tsne", "Diffusion Map"="dfm"), selected = "PCA",
                              multiple = FALSE,selectize = TRUE, width = NULL, size = NULL),
                  sliderInput("trajectorySliderDimensions", "Number of dimensions to use :", min = 0, max = 50, value = 10, step = 1),
                  selectInput("trajectoryStart", "Initial state:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                  selectInput("trajectoryEnd", "Final state:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                  actionButton(inputId = "trajectoryConfirm", label = "OK")
                ),
                box(
                  width = 9, status = "info", solidHeader = TRUE, title = "Trajectory analysis results",
                  tabsetPanel(type = "tabs",
                              tabPanel("Structure overview", 
                                       div(class="ldBar", id="traj1_loader", "data-preset"="circle"),
                                       plotOutput(outputId="trajectoryPlot", height = "1100px"),
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
                                    plotOutput(outputId = "trajectoryPseudotimePlot", height = "1100px"))
                              ))
                  ),
                )
              )
      ),
      
      #L-R analysis
      tabItem(tabName = "ligandReceptor",
               fluidRow(
                 box(
                   width = 3, status = "info", solidHeader = TRUE,
                   title = "L-R analysis parameters",
                   selectInput("ligandReceptorSender", "Ligand expressing cluster:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                   selectInput("ligandReceptorReciever", "Receptor expressing cluster:", choices=c("0"="0"), selected = "0", multiple = F, selectize = F),
                   actionButton(inputId = "ligandReceptorConfirm", label = "OK")
                 ),
                 box(
                   width = 9, status = "info", solidHeader = TRUE, title = "L-R analysis results",
                   div(class="ldBar", id="lr_loader", "data-preset"="circle"),
                   div(
                     tabsetPanel(type = "tabs",
                                 tabPanel("All interactions",
                                          plotlyOutput(outputId="ligandReceptorFullHeatmap", height = "1100px")),
                                 tabPanel("Curated interactions (documented in literature and publicly available databases)",
                                          plotlyOutput(outputId="ligandReceptorCuratedHeatmap", height = "1100px"))
                                 )
                   )
                 )
               )
      ),
      
      tabItem(tabName = "help",
              fluidRow(
                column(12, 
                       tabsetPanel(
                         tabPanel("Examples",
                                  div(class = "div_container",
                                      examples_help
                                  )
                         ),
                         tabPanel("Data Input"
                                  
                         ),
                         tabPanel("Utility options"
                                  
                         ),
                         tabPanel("Quality Control"
                                  
                         ),
                         tabPanel("Normalization"
                                  
                         ),
                         tabPanel("PCA"
                                  
                         ),
                         tabPanel("Clustering"
                                  
                         ),
                         tabPanel("UMAP"
                                  
                         ),
                         tabPanel("Markers"
                                  
                         ),
                         tabPanel("Cell Cycle"
                                  
                         ),
                         tabPanel("Enrichment"
                                  
                         ),
                         tabPanel("Cluster Annotation"
                                  
                         ),
                         tabPanel("Trajectory Analysis"
                                  
                         ),
                         tabPanel("Ligand-Receptor Analysis"
                                  
                         )
                       )
                )
              ) #fluidRow end
      ),
      
      tabItem (tabName = "about",
               div(id = "about_div", class = "div_container",
                   h1(class = "container_title", "About scAnner"),
                   HTML("<p class=container_text> scAnner is actively developed and maintained by the 
                              Bioinformatics and Integrative Biology Lab. </br></br> </p> 
                              <h2 class=sub_title> Developers </h2>
                              <ul>
                              <li> Christos Tzaferis, tzaferis[at]gmail[dot]com
                              <li> Evangelos Karatzas, karatzas[at]fleming[dot]gr
                              <li> George Kollias, kollias[at]fleming[dot]gr
                              <li> Georgios A. Pavlopoulos, pavlopoulos[at]fleming[dot]gr 
                              </ul>
                              <footer>
                              &copy; 2021 <a href=\"https://sites.google.com/site/pavlopoulossite\" target=\"_blank\">Bioinformatics and Integrative Biology Lab</a> | 
                              <a href=\"https://www.fleming.gr\" target=\"_blank\">Biomedical Sciences Research Center \"Alexander Fleming\"</a>
                              </footer>"
                   )
               )
      )
      
    )# tab item list
  )
)
