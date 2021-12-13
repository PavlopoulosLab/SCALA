#SET THE NUMBER OF CPUs TO USE FOR PACKAGE INSTALLATION
cpus=8

# Step 0: setup environmental values. ####
#0.1 Set CRAN mirror: This is required if you want to run the script in the CMD/terminal with Rscript, instead of R-studio
local({r <- getOption("repos")
r["CRAN"] <- "http://cran.r-project.org" 
options(repos=r)
})
#0.2 Set number of cpus for compilation. This uses the 'cpus' variable defined above
local({options(Ncpus=cpus)})

# Step 1: check if RCurl remotes, devtools, BiocManager and reticulate are available. If not, install them, and load them with require ####
if ("RCurl" %in% rownames(installed.packages()))
{
  print(sprintf("RCurl is already installed."))
} else {
  install.packages("RCurl")
}
if("remotes" %in% rownames(installed.packages()))
{
  require(remotes)
} else {
  install.packages("remotes")
  require(remotes)
}
if("devtools" %in% rownames(installed.packages()))
{
  require(devtools)
} else {
  install.packages("devtools")
  require(devtools)
}
if("BiocManager" %in% rownames(installed.packages()))
{
  require(BiocManager)
} else {
  install.packages("BiocManager")
  require(BiocManager)
}
if("reticulate" %in% rownames(installed.packages()))
{
  require(BiocManager)
} else {
  install.packages("reticulate")
  require(reticulate)
}

# Step 2: Check and install (if required) standard libraries (i.e. libraries that can be installed through CRAN without any other dependencies) ####
libraries_CRAN <- c(
  "shiny",
  "shinyjs",
  "shinydashboard",
  "shinycssloaders",
  "shinyalert",
  "bsplus"
  "tidyverse",
  "DT",
  "ggplot2",
  "ggpubr",
  "plotly",
  "igraph",
  "rgl",
  "RColorBrewer",
  "colorspace",
  "dplyr",
  "visNetwork",
  "gprofiler2",
  "heatmaply",
  "Seurat",
  "SeuratObject",
  "readr",
  "stringr",
  "pheatmap",
  "hdf5r",
  "missMDA",
  "dismo"
)
for (i in 1:length(libraries_CRAN))
{
  if (libraries_CRAN[i] %in% rownames(installed.packages()))
  {
    print(sprintf("%s is already installed.", libraries_CRAN[i]))
  }
  else
  {
    install.packages(libraries_CRAN[i])
  }
}


# Step 3: Check and install (if required) libraries that need to be installed from GitHub with remotes or devtools ####
if("chromVARmotifs" %in% rownames(installed.packages()))
{
  print("chromVARmotifs is installed")
} else {  
  devtools::install_github("GreenleafLab/chromVARmotifs", upgrade = c("never"))
}
if("SCopeLoomR" %in% rownames(installed.packages()))
{
  print("SCopeLoomR is installed")
} else {  
  devtools::install_github("aertslab/SCopeLoomR", upgrade=c("never"))
}
if("CIPR" %in% rownames(installed.packages()))
{
  print("CIPR is installed")
} else {  
  devtools::install_github("atakanekiz/CIPR-Package", build_vignettes = F, upgrade=c("never"))
}
if("mrlMBO" %in% rownames(installed.packages()))
{
  print("mrlMBO is installed")
} else {  
  remotes::install_github("mlr-org/mlrMBO", upgrade=c("never"))
}
if("nichenetr" %in% rownames(installed.packages()))
{
  print("nichenetr is installed")
} else {  
  devtools::install_github("saeyslab/nichenetr", build_vignettes = F, upgrade=c("never"))
}
if("destiny" %in% rownames(installed.packages()))
{
  print("destiny is installed")
} else {  
  remotes::install_github("theislab/destiny", upgrade=c("never"))
}
if("UCell" %in% rownames(installed.packages()))
{
  print("UCell is installed")
} else {  
  remotes::install_github("carmonalab/UCell", upgrade=c("never"))
}

# 

# Step 4: Check and install (if required) libraries from BioConductor through BiocManager
libraries_bioconductor <- c(
  "DirichletMultinomial",
  "TFBSTools",
  "motifmatchr",
  "limma",
  "AUCell",
  "RcisTarget",
  "GENIE3",
  "GSEABase",
  "RcisTarget",
  "dittoSeq",
  "slingshot",
  "JASPAR2020",
  "JASPAR2018",
  "JASPAR2016",
  "MAST",
  "BSgenome.Mmusculus.UCSC.mm9",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38"
)
for (i in 1:length(libraries_bioconductor))
  {
  if (libraries_bioconductor[i] %in% rownames(installed.packages()))
  {
    print(sprintf("%s is already installed.", libraries_bioconductor[i]))
  }
  else
  {
    BiocManager::install(libraries_bioconductor[i], update=F)
  }
}


# Step 5: Install ArchR, SCENIC and other very specific ones ####

# 5.1 SCENIC ####
if("SCENIC" %in% rownames(installed.packages()))
{
  print("SCENIC is installed")
}  else {
  bmv = BiocManager::version() #first, get the bioconductor version
  #if bioconductor 4 or newer
  if(unlist(bmv)[1] >= 4)
  {
    devtools::install_github("aertslab/SCENIC", upgrade = c("never")) 
  }
  else
  {
    devtools::install_github("aertslab/SCENIC@v1.1.2", upgrade = c("never"))
  }
}
# 5.2 PhateR
if("PhateR" %in% rownames(installed.packages()))
{
  
  print("PhateR is installed")
}  else {
  #reticulate::py_install("phate", pip=TRUE) #uncomment this if you don't already have phate installed
  devtools::install_github("KrishnaswamyLab/phateR", upgrade=c("never"))
  devtools::install_github("scottgigante/seurat", ref="patch/add-PHATE-again", upgrade=c("never"))
}
# 5.3 ArchR ####
if("ArchR" %in% rownames(installed.packages()))
{
  print("ArchR is installed")
}  else { devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories(), upgrade = c("never"))
  library(ArchR)
  ArchR::installExtraPackages() 
}

