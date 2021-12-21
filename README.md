# scAnner

A web application for multimodal analysis of single cell next generation sequencing data


----

## Table of Contents

1 [Overview](#overview)

2 [Requirements](#requirements)
  * [System dependencies](#system-dependencies)
  * [List of required R libraries](#list-of-required-r-libraries)

3 [Installation Instructions](#installation-instructions)
  * [On Linux](#on-linux)
  * [On Windows](#on-windows)


----
## Overview
A web tool that handles the analysis of scRNAseq data, from quality control and normalization, to dimensionality reduction, differential expression analysis, clustering and visualization.

Online version: http://bib.fleming.gr:8084/SCANNER/, http://scanner.pavlopouloslab.info


----
## Requirements
### System dependencies

- Operating System: Linux (any Ubuntu/Debian-based distribution). Instructions for installing and using on Windows, MacOS, and other Linux distros will be available soon.
- R version 4 or newer (preferrably 4.1.2). Older R versions (i.e. 3.6.x) **WILL NOT WORK**.
- R-studio or, alternatively, shiny-server for deployment.
- A number of R libraries from CRAN, Bioconductor and some GitHub repositories. Lists are given below. An installation script is also offered.
- A Conda installation (Anaconda or Miniconda)
- Python 3.7 or newer (preferrably through the Conda environment)
- numpy
- Pyscenic
- MACS2
- phate
- The GNU C, C++ and Fortran compilers and libraries (gcc, g++, gfortran)
- gdal-bin, libgdal-dev
- The Arrow C++ library (libarrow-dev, libarrow-glib-dev)
- hdf5-tools, hdf5-helpers, libhdf5-dev
- gsl-bin, libgsl23, libgslcblas0, libgsl-dev
- curl, libcurl-dev

**Note:** Almost all of the above can be installed through your Linux distribution's  package manager (apt, zypper etc).  **An installation bash script ("install_dependencies.sh") is offered to automate setup in Debian and Debian-based (Ubuntu, Mint, etc) distributions**.

### List of required R libraries
#### CRAN packages (can be installed directly from R/CRAN)
- devtools
- remotes
- reticulate
- BiocManager
- RCurl
- shiny
- shinyjs
- shinydashboard
- shinycssloaders
- shinythemes
- shinyalert
- bsplus
- tidyverse
- DT
- ggplot2
- ggpubr
- plotly
- igraph
- rgl
- RColorBrewer
- colorspace
- dplyr
- visNetwork
- gprofiler2
- heatmaply
- Seurat
- SeuratObject
- readr
- stringr
- pheatmap
- hdf5r
- missMDA
- dismo
- Cairo
- mlrMBO
- rhandsontable


#### Bioconductor packages (can be installed using BiocManager)
- DirichletMultinomial
- TFBSTools
- motifmatchr
- limma
- AUCell
- RcisTarget
- GENIE3
- GSEABase
- dittoSeq
- DelayedMatrixStats
- JASPAR2020
- JASPAR2018
- JASPAR2016
- MAST
- BSgenome.Mmusculus.UCSC.mm9
- BSgenome.Mmusculus.UCSC.mm10
- BSgenome.Hsapiens.UCSC.hg19
- BSgenome.Hsapiens.UCSC.hg38
- slingshot **version 2.1.0 or newer** (**Note:** Slingshot needs to be installed **AFTER** phateR has been installed and Seurat has been patched - see below)

#### Other third-party packages (require installation from GitHub repositories with devtools/remotes)
- [SCopeLoomR](https://github.com/aertslab/SCopeLoomR)
- [CIPR](https://github.com/atakanekiz/CIPR-Package)
- [nichenetr](https://github.com/saeyslab/nichenetr)
- [destiny](https://github.com/theislab/destiny)
- [UCell](https://github.com/carmonalab/UCell)
- [SCENIC](https://github.com/aertslab/SCENIC)
- [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs)
- [ArchR](https://github.com/GreenleafLab/ArchR)
- [harmony](https://github.com/immunogenomics/harmony)
- [presto](https://github.com/immunogenomics/presto)
- [phateR] (**Note**: SCANNER requires the GitHub version of phateR to perform phate runs with Seurat.  The version of phateR that is available in CRAN **WILL NOT WORK**)
- The Seurat patch for support of PHATE (**Note**: This needs to be applied **after** the original Seurat, phate and phateR have been installed)

**Note 1:** All of the aforementioned packages can be installed through R, Rscript or R-studio. **A script ("install_libraries.R")** is included to automate the installation process.
**Note 2:** Some of the packages require the existence of external libraries in the system, prior to their installation. Make sure to **install the system dependencies prior to installing the R libraries**.
**Note 3:** The order of installing is important for some of the packages. If you perform the installation by hand, make sure you follow the precise order given in **install_libraries.R**.
**Note 4:** SCANNER requires modifications to the Seurat package, in order for it to be compatible with phate and phateR (see above). These modifications may lead to inconsistencies with other packages in your R installation.

----

## Installation Instructions

### On Linux
Linux is the native environment SCANNER is designed to operate in.  The main steps that need to be followed are:
1. Install all software dependencies and packages
2. Install R and R-studio
3. Install all required R libraries
4. Retrieve and install the required reference datasets
4. Open the tool's project file (SCANNER.rproj) in R-studio, select **ui.R**, **server.R** or **global.R** and click "Run App" **or** (alternatively), configure shiny-server and setup SCANNER as a service.

To aid you in installing and configuring SCANNER, we provide three installation scripts, "install_dependencies.sh",  "install_libraries.R" and "install_datasets.sh". 

1. The first script ("install_dependencies.sh") installs and configures all required software in your Linux distribution.  The file provided is written for Debian and Debian-based (Ubuntu, Debian, Linux Mint etc) and can be run as follows:

>     sudo bash install_dependencies.sh

or

>     chmod +x install_dependencies.sh
>     sudo ./install_dependencies.sh

**Note 1:** To run this script, you **will** need an account with administrative ("sudo") privileges. If you are working on a system without sudo access, please consult your system administrator.
**Note 2:** Users of other Linux distributions, like SUSE, Fedora, Arch etc, should edit this script and replace "apt-get" with the analogous package manager of their system (e.g. "zypper").

2. The second script ("install_libraries.R") will install all required libraries in your R environment.  It can be run as follows:

>     Rscript install_libraries.R

(local installation of R packages in the user's home environment)
or

>     sudo Rscript install_libraries.R

(system-wide installation)

or it can be loaded and executed in R-studio.

3. The third script ("install_datasets.sh") utilizes wget to retrieve a number of helper files that are required by several of the analyses in SCANNER. These need to be downloaded, unzipped and placed in specific directories. You can run the script as follows:

>     bash install_datasets.sh

or

>     chmod +x install_datasets.sh
>     ./install_ddatasets.sh


### On Windows
- To be added



----


