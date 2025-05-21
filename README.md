# SCALA

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

4 [Advanced Operations](#advanced-operations)
  * [Creating and using arrow files for scATAC analysis](#creating-and-using-arrow-files-for-scatac-analysis)
  * [Run pyscenic to perform analyses with SCALA](#run-pyscenic-to-perform-analyses-with-SCALA)

----
## Overview
A web tool that handles the analysis of scRNAseq data, from quality control and normalization, to dimensionality reduction, differential expression analysis, clustering and visualization.

Online version: http://scala.pavlopouloslab.info


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
- optparse


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
- [phateR](https://github.com/KrishnaswamyLab/phateR) (**Note**: SCALA requires the GitHub version of phateR to perform phate runs with Seurat.  The version of phateR that is available in CRAN **WILL NOT WORK**)
- [The Seurat patch for support of PHATE](https://github.com/scottgigante/seurat) (**Note**: This needs to be applied **after** the original Seurat, phate and phateR have been installed)

**Note 1:** All of the aforementioned packages can be installed through R, Rscript or R-studio. **A script ("install_libraries.R")** is included to automate the installation process.

**Note 2:** Some of the packages require the existence of external libraries in the system, prior to their installation. Make sure to **install the system dependencies prior to installing the R libraries**.

**Note 3:** The order of installing is important for some of the packages. If you perform the installation by hand, make sure you follow the precise order given in **install_libraries.R**.

**Note 4:** SCALA requires modifications to the Seurat package, in order for it to be compatible with phate and phateR (see above). These modifications may lead to inconsistencies with other packages in your R installation.

----

## Installation Instructions

### On Linux
Linux is the native environment SCALA is designed to operate in.  The main steps that need to be followed are:
1. Install all software dependencies and packages
2. Install R and R-studio
3. Install all required R libraries
4. Retrieve and install the required reference datasets
4. Open the tool's project file (SCALA.rproj) in R-studio, select **ui.R**, **server.R** or **global.R** and click "Run App" **or** (alternatively), configure shiny-server and setup SCALA as a service.

To aid you in installing and configuring SCALA, we provide three installation scripts, "install_dependencies.sh",  "install_libraries.R" and "install_datasets.sh". 

1. The first script ("install_dependencies.sh") installs and configures all required software in your Linux distribution.  The file provided is written for Debian and Debian-based (Ubuntu, Debian, Linux Mint etc) and can be run as follows:

     sudo bash install_dependencies.sh

or

     chmod +x install_dependencies.sh
     sudo ./install_dependencies.sh

**Note 1:** To run this script, you **will** need an account with administrative ("sudo") privileges. If you are working on a system without sudo access, please consult your system administrator.
**Note 2:** Users of other Linux distributions, like SUSE, Fedora, Arch etc, should edit this script and replace "apt-get" with the analogous package manager of their system (e.g. "zypper").

2. The second script ("install_libraries.R") will install all required libraries in your R environment.  It can be run as follows:

     Rscript install_libraries.R

(local installation of R packages in the user's home environment)
or

     sudo Rscript install_libraries.R

(system-wide installation)

or it can be loaded and executed in R-studio.

3. The third script ("install_datasets.sh") utilizes wget to retrieve a number of helper files that are required by several of the analyses in SCALA. These need to be downloaded, unzipped and placed in specific directories. You can run the script as follows:

     bash install_datasets.sh

or

     chmod +x install_datasets.sh
     ./install_ddatasets.sh


### On Windows
- To be added

## Advanced Operations

### Creating and using arrow files for scATAC analysis

We provide a script ("create_arrow_file.R"), in order to easily create scATAC-seq *.arrow files. The particular file type is compatible with ArchR package, which is the essential component of SCALA's scATAC-seq data analysis pipeline.

create_arrow_file.R requires R version 4, with "optparse" and "ArchR" packages installed. Since SCALA is compatible with hg19, hg38, mm9, and mm10 UCSC genome builds, the respective packages should be also installed accordingly. All of the aforementioned packages are part of SCALA's dependencies, but can also be installed as follows:

    install.package("optparse")
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
    BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

An example of how to run create_arrow_file.R using an mm10 mouse scATAC-seq fragments file (see the instruction of how to run cellranger-atac count) dataset is shown below:

    Rscript create_arrow_file.R -s mm10 -a wt4 -f wt4.fragments.ucsc.tsv.gz -t 4 -m 1000 -p 1

Parameter description:

    "-s" or  "--species": The genome buid used during the alignment and counting step. Valid options are mm9, mm10, hg19, hg38.
    "-a" or "--arrow_prefix": Arrow file prefix. The default name will be set as myProject.
    "-f" or "--fragmentfile": The scATAC-seq fragment file is located in the outs/ folder of the cellranger-atac count run. See https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments. for more details. The default file path will be set fragments.tsv.gz.
    "-t" or "--minTSS": The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering for use in downstream analyses. The default value will be set to 4.
    "-m" or "--minFrags": The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. The default value will be set to 1000.
    "-p" or "--threads": The number of threads to be used for parallel computing. The default value will be set to 8. For windows users, it will be automatically set to 1. 

create_arrow_file.R run will generate 1 arrow file per fragment file in the current directory. Each of these files can be uploaded to SCALA, under the "DATA INPUT" and "Arrow input files (scATAC-seq)" tabs, in order to initialize a SCALA scATAC-seq analysis.


### Run pyscenic to perform analyses with SCALA

We provide a convenient wrapper called "pyscenic_local.py" in order to easily run pyscenic with default parameters. If the user desires to tweak pyscenic parameters, he/she is encouraged to visit pyscenic site (https://pyscenic.readthedocs.io/) and run the pipeline accordingly.

The particular wrapper is applicable to scRNA-seq loom files generated by the R package SCopeLoomR. Ideally, the basic steps of a typical scRNA-seq analysis should have been applied.
Gene filtering is considered as a vital step. The particular wrapper currently supports only hg19, hg38, mm9 and mm10 genome builds.

We encourage the user to first run all the steps provided in our SCALA web server.


The particular script should be run by using python v.3.7 or newer.

We recommend to install pyscenic under a conda virtual environment, for example: 

    conda create -n pyscenic python=3.7
    conda activate pyscenic
    conda install numpy
    conda install -c anaconda cytoolz
    pip install pyscenic


An example of how to run pyscenic_local.py for an analyzed mm10 mouse scRNA-seq dataset is shown below:

    python pyscenic_local.py -l /media/hdd1/Single_Cell_Unit/Lyda_DARE/SCENIC/combined_myeloid_WT_only.gene.filtered2.loom -o test_pyscenic_local -g mm10 -a /media/hdd1/Single_Cell_Unit/Lyda_DARE/SCENIC/ -t 2 

Parameter description:

    '--loomFile' or '-l' : A valid loom file generated by SCopeLoomR R package.
    '--outputFolder' or '-o' : The output folder path where the pyscenic results will be stored.',
    --genomeBuild' or '-g' : The genome buid used during the alignment and counting step. Valid options are mm9, mm10, hg19, hg38'
    '--annFolder' or '-a' : 'A folder containing all the mandatory annotation files, including ranking databases *.feather files, TF names file and motif annotation database *.tbl file.
    '--threads' or '-t' :The number of threads to be used for each step of the analysis.

pyscenic_local.py run will generate 3 files, located in the predefined output folder ('--outputFolder' or '-o' ): adjacencies.*.tsv, auc.*.loom and nes.score.*.csv. 
auc.*.loom contains all the necessary information for visualizing the inferred Regulons in SCALA. This loom can be finally uploaded to http://SCALA.pavlopouloslab.info/ , and in particular under the scRNA-seq Gene Regulatory Network analysis tab, in order to proceed with regulon visualization.

The required files can be found in the following links.  They can also be retrieved automatically during the installation of SCALA, using the "install_datasets.sh" script.

Feather Files:

mm9:
 https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-10species.mc9nr.feather
 https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-10species.mc9nr.feather

mm10:
https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
 https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather

hg19:
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.feather

hg38:
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

Tbl Files:

Mouse:
https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl

Human:
https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

TF lists:

Mouse:
https://github.com/aertslab/pySCENIC/blob/master/resources/mm_mgi_tfs.txt

Human:
https://github.com/aertslab/pySCENIC/blob/master/resources/hs_hgnc_curated_tfs.txt



---

## 📚 Publication

**SCALA: A complete solution for multimodal analysis of single-cell Next Generation Sequencing data**
Tzaferis C., Karatzas E., Baltoumas F.A., Pavlopoulos G.A., Kollias G., Konstantopoulos D.
*Computational and Structural Biotechnology Journal*. 2023 Oct 20;21:5382–5393. doi: [10.1016/j.csbj.2023.10.032](https://doi.org/10.1016/j.csbj.2023.10.032)
PMID: [38022693](https://pubmed.ncbi.nlm.nih.gov/38022693/)  PMCID: [PMC10651449](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10651449/)

---

## 📄 License

This project is licensed under the **MIT License**. See [LICENSE](LICENSE) for details.







