###############################Example files####################################
examples_help <- HTML(
                  '<h4>Example matrix input file for scRNA-seq:</h4>
                  <li><a href="example_files/exampleMatrix.zip" download> PBMC 3K Cell-Gene count matrix (*.txt) </a></li>
                  <h4>Example 10x input files for scRNA-seq:</h4>
                  <li><a href="example_files/barcodes.tsv.gz" download> PBMC 3K barcodes file (*.tsv.gz) </a></li>
                  <li><a href="example_files/features.tsv.gz" download> PBMC 3K features file (*.tsv.gz) </a></li>
                  <li><a href="example_files/matrix.mtx.gz" download> PBMC 3K matrix file (*.mtx.gz) </a></li>
                  
                  <h4>Example files for GRN analysis:</h4>
                  <li><a href="example_files/processed_seurat_object-2021-12-19.zip" download> PBMC 3K seurat object (*.RDS) </a></li>
                  <li><a href="example_files/auc.hg19.zip" download> Pyscenic output file (*.loom) </a></li>
                  
                  <h4> Example files for scATAC-seq: </h4>
                  <li><a href="example_files/arrowFile.zip" download> PBMC arrow file (*.arrow) </a></li>
                  <li><a href="example_files/PBMCs_human_signac_peaks.zip" download> PBMC peakset file (*.bed) </a></li>
                 ')

###########################Upload###############################################
file_upload_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>Data upload tab</b> enables 
the user to upload scRNA-seq or scATAC-seq datasets to scAnner, 
in order to initialize single-cell analysis.<br> ScAnner provides the option to upload different file input formats. </h4>')

file_upload_txt <- HTML('<div class="col-md-4 scrollable">
                          <img src = "images/help_page/DATA_INPUT_scRNA_count_matrix.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 1: </b> The File Upload form for cell-gene count matrices </figcaption>
                          </div>
                          
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3> 1. Upload a cell-gene count matrix </h3>
                        <p>
                        <ol>
                        <li> Project name: User defined project name.
                        <li> File Upload: Click the <b>Browse</b> button of the upload form to select a single-cell RNA-seq gene/feature by cell/barcode matrix.<br>
                             <b>Notes: The file should not exceed 500 MB.</b>

                        <li> Gene and cell filtering parameters:<br>
                             <ul>
                                <li> The first sliding bar will subset the count matrix in order to include features/genes detected in at least this many cells.
                                <li> The second sliding bar will subset the count matrix in order to include cells where at least this many features are detected.
                             </ul>
                        <li> Select organism: <br>
                             <ul>
                                <li> If the counting procedure was accomplished by using UCSC mm9 or mm10 genome builds as a reference, the user should choose "Mouse".
                                <li> If the counting procedure was accomplished by using UCSC hg19 or hg38 genome builds as a reference, the user should choose "Human".
                             </ul>
                        </ol>
                        </p>
                        <hr>
                            <h3>2. Load the example count matrix </h3>
                            <p> 
                              By pressing the button <b>Load example</b> a scRNA-seq human example of Peripheral Blood Mononuclear Cells (PBMC),freely available from 10X Genomics, will be loaded. 
                            </p>
                        <hr>
                            <h3>3. Export working object as .RDS file </h3>
                            <p> 
                              At any point of the analysis the button <b>Export .RDS</b> exports a scRNA-seq Seurat Object as an RDS file compatible to R language. 
                            </p>    
                        </div>
                        ')

file_upload_10x <- HTML('<div class="col-md-4 scrollable">
                          <img src = "images/help_page/DATA_INPUT_scRNA_10x.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 2: </b> The File Upload form for 10x input files </figcaption>
                          </div>
                          
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3> 1. Upload 10x files </h3>
                        <p>
                        <ol>
                        <li> Project name: User defined project name.
                        <li> File Upload: Click the <b>Browse</b> buttons of the upload form to select "cellranger count" single-cell RNA-seq files
                          <ul>  
                            <li> A gzipped tsv file the contains only detected (filtered by cellranger count pipeline) cellular barcodes. The file should have the name "barcodes.tsv.gz".
                            <li> A gzipped tsv file of features (genes) that correspond to row indices. Feature ID,  feature name and feature type (Gene expression) are stored in the first, 
                                  second, and third column of the file, respectively. The file should have the name "features.tsv.gz".
                            <li> A gzipped Market Exchange Format (MEX) featutre-barcode count matrix. The file should have the name "matrix.mtx.gz". <br>
                                <b>Notes: The files should not exceed 500 MB in total. Additionally, the files should be named as described above in order to be uploaded successfully.</b>
                          </ul>
                        <li> Gene and cell filtering parameters:<br>
                             <ul>
                                <li> The first sliding bar will subset the count matrix in order to include features/genes detected in at least this many cells.
                                <li> The second sliding bar will subset the count matrix in order to include cells where at least this many features are detected.
                             </ul>
                        <li> Select organism: <br>
                             <ul>
                                <li> If the counting procedure was accomplished by using UCSC mm9 or mm10 genome builds as a reference, the user should choose "Mouse".
                                <li> If the counting procedure was accomplished by using UCSC hg19 or hg38 genome builds as a reference, the user should choose "Human".
                             </ul>
                        </ol>
                        </p>
                        <hr>
                            <h3>2. Load the example 10x files </h3>
                            <p> 
                              By pressing the button <b>Load example</b> a scRNA-seq human example of Peripheral Blood Mononuclear Cells (PBMC),freely available from 10X Genomics, will be loaded. 
                            </p>
                        <hr>
                            <h3>3. Export working object as .RDS file </h3>
                            <p> 
                              At any point of the analysis the button <b>Export .RDS</b> exports a scRNA-seq Seurat Object as an RDS file compatible to R language. 
                            </p>    
                        </div>
                        ')

file_upload_arrow <- HTML('<div class="col-md-4 scrollable">
                          <img src = "images/help_page/DATA_INPUT_scATAC-seq.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 3: </b> The File Upload form for arrow files </figcaption>
                          </div>
                          
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3> 1. Upload an arow file </h3>
                        <p>
                        <ol>
                        <li> Project name: User defined project name.
                        <li> File Upload: Click the <b>Browse</b> button of the upload form to select a single-cell ATAC-seq arrow file. Arrow files is a file type that stores all of the scATAC-seq data 
                        associated with an individual sample. Arrow files can be created by coupling scATAC-seq fragment files generated by cellranger-atac, and create_arrow_file.R Rscript (provided), as described in 
                        scAnner\'s github page. The particular operation should be applied locally.<br>
                        <b>Notes: The file should not exceed 2 GB.</b>
                        <li> Select organism: <br>
                             <ul>
                                <li> If the counting procedure was accomplished by using UCSC mm10 genome builds as a reference, the user should choose "mm10".
                                <li> If the counting procedure was accomplished by using UCSC hg19 or hg38 genome builds as a reference, the user should choose "hg19" or "hg38".
                             </ul>
                        <li> Threads to be used: The number of CPU threads to be used during the analysis.
                        </ol>
                        </p>
                        <hr>
                            <h3>2. Load an example arrow file </h3>
                            <p> 
                              By pressing the button <b>Load example</b> a scRNA-seq human example of Peripheral Blood Mononuclear Cells (PBMC),freely available from 10X Genomics, will be loaded. 
                            </p>
                        </div>
                        ')

file_upload_tab_new_project <- HTML('<div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px;">
<h4 style = "line-height: 1.5;">
                                      When you try to upload a new dataset at the same time that another one is already loaded, you will be prompted to discard the current working dataset. If you select <b>Yes</b> your 
                                      progress will be deleted, note that this is an irreversible action. After that you can start a new project by submitting the new input files.<br>
                                      If you wish to continue working on your current project slect <b>No</b>. 
                                    </h4>
                                    </div>')

file_upload_metadata_RNA <- HTML('
                        <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3><u> scRNA-seq metadata table </u></h3>
                        <p>
                          If the data upload was successfull, a scRNA-seq cell metadata table will appear. 
                          Project name, number of Unique Moleculare Identifiers (UMIs) per cell, number of detected features/genes per cell, percentage of mitochondrial RNA (percent.mt) 
                          and cell id are stored in the 1st, 2nd, 3rd, 4th and 5th column respectively. Additional columns will be created during the next steps of the analysis, and the cell 
                          metadata table will be updated automatically.The user is also able to download the metadata table using the <b>Save table</b> button.
                        </p>
                        
                        </div>
                        <div class="col-md-12 scrollable">
                          <img src = "images/help_page/DATA_INPUT_scRNA_Metedata_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 4: </b> RNA metadata table </figcaption>
                          </div>
                        ')

file_upload_metadata_ATAC <- HTML('
                        <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3><u> scATAC-seq metadata table </u></h3>
                        <p>
                          If the data upload was successfull, a scATAC-seq cell metadata table will appear. Project name, per-cell TSS enrichment, per-cell Unique Molecular Identifiers (UMIs) in Transcription Start Sites (TSSs), 
                          per-cell UMIs in promoters, per-cell UMIs in ENCODE black listed regions, per-cell fraction of accessible fragments that overlap promoters (PromoterRatio), 
                          information about high quality cell (PassQC), per-cell nucleosome ratio, per-cell number of fragments with size greater than 294 bp (nMultiFrags), per-cell 
                          number of fragments with size less than 147 bp (nMonoFrags), per-cell number of fragments (nFrags), per-cell number of fragments with size between 147 bp 
                          and 294 bp (diFrags),  per-cell ratio of fragments in ENCODE black listed regions, and cell id are stored in the 1st, 2nd, 3rd, 4th, 5th, 6th, 7th, 8th, 9th, 
                          10th, 11th, 12th, 13th, and 14th columns respectively. Additional column will be created during the next steps of the analysis, and the cell metadata table 
                          will be updated automatically. The user is also able to download the metadata table using the <b>Save table</b> button.
                        </p>
                        
                        </div>
                        <div class="col-md-12 scrollable">
                          <img src = "images/help_page/DATA_INPUT_scATAC-seq_Metadata_table_merged.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 5: </b> ATAC metadata table </figcaption>
                          </div>
                        ')

#######################################QC#######################################
qc_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>
                          Quality control tab</b> enables the generation of Quality Control (QC) plots for scRNA-seq and scATAC-seq datasets,<br> 
                          in order to examine the quality of the single-cell experiment. </h4>')

rna_qc <- HTML('
                          
                <div class="col-md-12 scrollable"> 
                  <img src = "images/help_page/QC_pre_params.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 6: </b> Quality control parameters and plots before filtering </figcaption>
                </div>
                 
                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                 <h3><u> Quality control parameters </u></h3>
                 <p>
                 <ol>
                  <li> Display quality control plots before filtering: The user is able to visualize QC plots, before the application of any cell-specific filter.<br>
                       The particular visualization allows the user to explore cell count metrics and filter-out low-quality cells and non-informative genes in the next steps of this analysis. 
                       The QC metrics included in this exploration are:
                       <ul>
                          <li> The number of detected features in each cell. Low-quality barcodes or empty droplets will exhibit low number of features, while multiplets may have a very high number of genes detected.
                          <li> The total number of UMIs detected in each cell. This metric correlates strongly with the feature-per-cell metric.
                          <li> The percentage of mitochondrial UMIs. Low-quality barcodes and dying cells contain high number of mitochondrial reads.
                          <li> The correlation between detected features and total counts per-cell. The particular metrics should exhibit high correlation.
                          <li> The correlation between mitochondrial counts and total counts per-cell. The particular metrics should exhibit negligible correlation.
                          <li> The total number of non-filtered cells.
                      </ul>
                  <li> Filter out low quality cells: After exploring the QC metrics of the unfiltered cells, the user is enabled to apply his/her filtering 
                      <ul>
                        <li> Minimum features detected: Filter out all cells with less than the user-defined number of detected features.
                        <li> Maximum features detected: Filter out all cells with more than the user-defined number of detected features.
                        <li> Mitochondrial %: Filter out cells that have greater than the user-defined percentage of mitochondrial counts.
                     </ul>
                 </ol>
                 </p>
                 </div>
               ')

rna_qc_pf <- HTML('
                          
                <div class="col-md-12 scrollable"> 
                  <img src = "images/help_page/QC_post.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 7: </b> Quality control plots post filtering </figcaption>
                </div>
                 
                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                 <h3><u> Inspection of filtering criteria </u></h3>
                 <p>
                    The particular visualization allows the user to explore cell count metrics after the filtering procedure. The user is able to readjust the defined criteria if the results 
                    are not satisfying. The QC metrics included in this figure are described in the previous tab.
                 </p>
                 </div>
               ')

atac_qc <- HTML('
                <div class="col-md-12 scrollable"> 
                  <img src = "images/help_page/QC_all_ATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 8: </b> Quality control plots after soft filtering </figcaption>
                </div>
                 
                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                 <h3><u> Inspection of filtering criteria for ATAC data </u></h3>
                 <p>
                    <h4>- Display soft filtered quality control plots: Visualization of soft filtered QC plots. Soft cell-filtering has been applied during the arrow file creation, 
                    locally, by using the create_arrow_file.R Rscript (provided in github). The filtering parameters and their default values are:</h4>
                    <ul>
                      <li> "-t" or "--minTSS": The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering for use in downstream analyses. The default value will be set to 4.
                      <li> "-m" or "--minFrags": The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses. The default value will be set to 1000.
                    </ul>
                    Note:	If the user wants to adjust additional parameters included in the cell metadata, the ArchR function createArrowFiles() should be used in an R environment.
                 </p>
                 
                 <p>
                    <h4>- The particular visualization allows the user to explore cell metrics after the soft-filtering procedure. The user is able to readjust the defined criteria if the 
                    results are not satisfying, by using the create_arrow_file.R Rscript and by tweaking the "--minTSS" and "--minFrags" parameters. The QC metrics included in this illustration are:</h4>
                    <ul>
                      <li> Violin plot of the per-cell TSS enrichment scores.
                      <li> Ridge plot of the per cell logged unique nuclear fragments number.
                      <li> Scatter plot of the logged unique nuclear fragments number versus the Transcription Start Site (TSS) enrichment score. Dashed lines indicate the thresholds used.
                    </ul>
                 </p>
                 </div>
               ')

############################Normalization#######################################

norm_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>DATA NORMALIZATION</b> tab enables the user to 
                          perform cell-specific normalization and gene-specific scaling of scRNA-seq count data. The aforementioned procedures are essential in order to perform the downstream steps of 
                          the analysis, including dimensionality reduction and differential expression. After these steps, the features that exhibit the higher variability in the dataset are detected.
                          </h4>')

rna_normalization_param <- HTML('<div class="col-md-4 scrollable">
                          <img src = "images/help_page/DATA_NORMALIZATION_AND_SCALING_parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 9: </b> Parameters for the normalization and scaling procedure in RNA seq data </figcaption>
                          </div>
                          
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        
                        <h3><u>Normalization</u></h3>
                          <li> A global-scaling normalization is applied. The particular methodology normalizes the per-cell gene expression counts by the total 
                               cell counts, multiplies this value by a scale factor (10,000 by default), and finally performs log-transformation. The user is able to 
                               alter the scaling factor.
                        
                        <h3><u> Identification of highly variable features </u></h3>
                        Detection of a set of genes that show high cell-to-cell variation in the scRNA-seq count matrix. One of the following methods can be selected:
                             <ul>
                                <li> Variance Stabilizing Transformation (vst) method. The particular method fits a line to the relationship of logged variance and logged 
                                     mean using local polynomial regression. Consequently, standardization of feature values using the observed mean and expected variance is performed. 
                                     Feature variance is finally calculated on standardized values, after clipping to a maximum. A fixed number of variable features is returned 
                                     (default: 2,000 features. Recommended values range: 1000 - 8000).
                                <li> Mean-Variance method. The particular method uses a function to calculate average gene counts and gene dispersions. All genes are separated into 20 
                                     bins according to their average counts. Finally, dispersion z-scores are calculated in each gene group.
                                <li> Dispersion method. Feature selection according to the highest dispersion values.
                             </ul>
                        
                        <h3><u> Scaling the data </u></h3>
                        <p>
                          A linear transformation that shifts the counts of each feature, so that the mean counts across cells is 0, and the variance across cells is 1. This step ensures that 
                          highly-expressed genes will not introduce biases during the downstream analysis.Additionally, in order to remove unwanted sources of variation, the user is able to 
                          select which metadata values would like to regress out during the scaling procedure. Typically, mitochondrial percentage is usually regressed out.
                        </p>
                          
                        </div>
                        ')

rna_normalization_output <- HTML('
                          
                <div class="col-md-12 scrollable"> 
                  <img src = "images/help_page/DATA_NORMALIZATION_AND_SCALING_mvgs.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 10: </b> Most variable genes </figcaption>
                </div>
                 
                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                 <h3><u> Exploration of MVGs </u></h3>
                 <p>
                    The particular visualization is a scatter plot that depicts the standardized variance versus the average expression of all features (vst method applied).
                    If one of the other two methods is selected dispersion values are shown in the Y axis. The number of highly variable genes is also reported (red dots).
                 </p>
                 </div>
               ')

##########################PCA/LSI###############################################

pca_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>PCA/LSI</b> tab enables the user 
                                to perform Principal Component Analysis (PCA) to scRNA-seq datasets and Latent Semantic Indexing analysis to scATAC-seq datasets.
                          </h4>')

pca_optimal_pcs <- HTML('
                        <div class="col-md-12 scrollable"> 
                          <img src = "images/help_page/PCA_selection.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <img src = "images/help_page/PCA_slow_results.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 11:</b> Determination of optimal number of principal components </figcaption>
                       </div>
                       
                       <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <h3>The user has two options during this analysis step:</h3>
                          <ol>
                            <li> Enable the automatic identification of the optimal number of Principal Components (PCs) using 1-fold SVA-CV. This option is significantly slower. 
                                 The optimal number of PCs is indicated by the red dotted line.   
                            <li> Perform PCA without automatic identification of the optimal number of Principal Components (PCs). The user will then decide about the dimensionality 
                                 of the dataset (the number of most informative PCS), based on the generated "elbow" plot.
                          <ol>
                          <p>
                            In the particular illustration,<br> (a) an elbow plot depicting the ranking of PCs based on the percentage of variance explained by each of them is illustrated, 
                            as also<br> (b) a scatter plot of cells in 2D PCA space, using the first two PCs.
                          </p>
                       </div>
                        ')

pca_explore_pcs <- HTML('
                        <div class="col-md-12 scrollable"> 
                          <img src = "images/help_page/PCA_exploration.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 12:</b> Exploration of principal components </figcaption>
                       </div>
                       
                       <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <p>
                            The particular visualization depicts<br> (a) the loading scores of the top genes of a PC of interest (top 30 features for the particular example), and<br> 
                            (b) a heatmap of scaled counts of the top loadings of the PC of interest, across cells.
                          </p>
                       </div>
                        ')

pca_lsi <- HTML('
                <div class="col-md-6"> 
                          <img src = "images/help_page/PCA_lsi.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 13:</b> LSI options </figcaption>
                </div>
                
                <div class="col-md-6" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <p>
                            In scAnner, LSI is run in an iterative manner (number of iterations). A first LSI transformation run is applied using the most accessible features. 
                            This procedure identifies lower resolution clusters that are not batch confounded. Consequently, average accessibility for each of these clusters calculated across all features. 
                            Finally, the most variable features are identified across the low resolution clusters, and are used as input for the next LSI iteration. 
                            The parameters for identifying low resolution clusters are:
                            <ul>
                              <li> Number of variable features. Defaults to 25,000 features.
                              <li> Number of LSI dimensions to use. Defaults to 30 dimensions.
                              <li> Cluster resolution. Defaults to 1.
                            </ul>
                          </p>
                       </div>
                ')

#########################Clustering#############################################
clustering_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                                The <b>CLUSTERING</b> tab enables the user to perform graph-based clustering in scRNA-seq and scATAC-seq datasets, in order to define cell types and cellular states. 
                             </h4>')

clustering_rna_input <- HTML('
                             <div class="col-md-4 scrollable"> 
                              <img src = "images/help_page/Clustering Parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 14:</b> Clustering options </figcaption>
                            </div>
                
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                           <p>
                              <h3><u>Clustering steps</u></h3>
                              <ol>
                                <li> Construction of the shared nearest neighbour (SNN) graph<br>
                                     Initially, cells are embedded in a K-nearest neighbor (KNN) graph structure based on Euclidean distances in the PCA space. 
                                     Cells that exhibit similar gene expression profiles are connected with edges. 
                                     The user can define (a) the maximum number of neighbors of each cell (defaults to 20), as also the number of principal components 
                                     to use (defaults to 10). 
                               <li> Communities\' detection (Louvain algorithm)<br>
                                    The formed graph of the previous step is partitioned into highly interconnected communities using the Louvain algorithm.
                                    The user can define (a) the desired clustering resolution (defaults to 0.5), as also the number of principal components to 
                                    use (defaults to 10). Higher values of these parameters will result to an increased formation of communities/cell clusters.
                              </ol>
                          </p>
                       </div>
                             ')

clustering_rna_output <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering Results_1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 15:</b> Table of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results can be explored by downloading the respective data table, which encapsulates information regarding the cluster name, 
                                    the number of cells in each cluster, and the percentage of total cells that is included in each cluster. The particular information is stored 
                                    in the 1st, 2nd and 3rd columns of the aforementioned data table respectively.
                                  </p>
                                  <br>
                              </div>
                              
                              <br>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering Results_2.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 16:</b> Barplot of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results are also depicted in a stacked barplot of percentages of cells in each of the identified clusters.
                                  </p>
                                  <br>
                              </div>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering SNN_graph.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 17:</b> The shared nearest neighbor graph (SNN) used in clustering </figcaption>
                                  <br>
                                  <p>
                                    The user is also able to explore the formation of the SNN graph that led to the identification of the final cluster set.
                                  </p>
                                  <br>
                              </div>
                              ')

clustering_atac_input <- HTML('
                             <div class="col-md-4 scrollable"> 
                              <img src = "images/help_page/ClusteringOptionsATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 18:</b> Clustering options </figcaption>
                            </div>
                
                        <div class="col-md-5" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                           <p>
                              <h3>scATAC-seq clustering, enables the same methodology as in described in scRNA-seq.</h3><br> 
                              The user is able to define<br> (a) the number of LSI dimensions to use 
                              (defaults to 30), as also<br> (b) the clustering resolution (defaults to 0.6).
                          </p>
                       </div>
                             ')

clustering_atac_output <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/ClusteringResultsATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 19:</b> Table of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results can be explored by downloading the respective data table, which encapsulates information regarding the cluster name, 
                                    the number of cells in each cluster, and the percentage of total cells that is included in each cluster. The particular information is stored 
                                    in the 1st, 2nd and 3rd columns of the aforementioned data table respectively.
                                  </p>
                                  <br>
                              </div>
                              
                              <br>
                              <br>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/ClusteringBarplotATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 20:</b> Table of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results are also depicted in a barplot of percentages of cells in each of the identified clusters.
                                  </p>
                                  <br>
                              </div>
                              ')

################################Umap############################################
umap_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                                The <b>ADDITIONAL DIMENSIONALITY REDUCTION METHODS</b> tab enables the user to perform nonlinear dimensionality 
                                reduction methodologies to scRNA-seq and scATAC-seq, in order to uncover patterns of cell similarity and differentiation. 
                             </h4>')

umap_rna_input <- HTML('
                      <div class="col-md-4 scrollable"> 
                              <img src = "images/help_page/UMAPparameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 21:</b> Options </figcaption>
                            </div>
                        
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                           <p>
                              <h3><u>Options for non linear dimensionality reductions</u></h3>
                              <ol>
                                <li> User-defined parameters for nonlinear dimensionality reduction application
                                  <ul>
                                    <li> Number of principal components to use. The user is able to define the number of the most informative principal components to be used for nonlinear dimensionality reduction.
                                    <li> Number of dimensions to fit output. The user is able to define the number of output dimensions for each of the available nonlinear dimensionality reduction methods.
                                  </ul>
                                <li> Available dimensionality reduction methodologies. The user can choose among: 
                                  <ul>
                                    <li> Uniform Manifold Approximation and Projection (UMAP)
                                    <li> t-distributed stochastic neighbor embedding (tSNE)
                                    <li> Diffusion Maps
                                    <li> PHATE
                                  </ul>   
                                <li> Display settings
                                  <ul>
                                    <li> Plot type. The user is able to choose which of the available nonlinear dimensionality reduction spaces to visualize. Defaults to pca.
                                    <li> Dimensions. The number of dimensions to be plotted. Defaults to 2D.
                                    <li> Color by. Apply different coloring based on the available cell metadata columns. Defaults to orig_ident.
                                    <li> Size. The size of the plotted cells/dots. Defaults to 5.
                                    <li> Opacity. Define the level of the transparency of the plotted cells/dots. Defaults to 1.
                                    <li> Border width. The dot/cell perimeter thickness. Defaults to 0.5.
                                  </ul>   
                              </ol>
                          </p>
                       </div>
                   ')

umap_rna_output <- HTML('
                        <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <br>
                          <img src = "images/help_page/UMAP.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 22:</b> Output example </figcaption>
                          <br>
                          <p>
                            Each nonlinear dimensionality reduction space is visualized using a scatter plot. Each dot represents a cell, and each axis a nonlinear dimension. 
                            The particular illustration hosts a UMAP visualization in 2D.
                          </p>
                          <br>
                      </div>
                   ')

umap_atac_input <- HTML('
                        <div class="col-md-4 scrollable"> 
                              <img src = "images/help_page/UMAPparameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 23:</b> Options </figcaption>
                            </div>
                        
                        <div class="col-md-8 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                           <p>
                              <h3><u>Options for non linear dimensionality reductions</u></h3>
                              <ol>
                                <li> User-defined parameters for nonlinear dimensionality reduction application
                                  <ul>
                                    <li> Number of LSI dimensions to use. The user is able to define the number of the most informative LSI dimensions to be used for nonlinear dimensionality reduction.
                                    <li> Number of dimensions to fit output. The user is able to define the number of output dimensions for each of the available nonlinear dimensionality reduction methods.
                                  </ul>
                                <li> Available dimensionality reduction methodologies. The user can choose among: 
                                  <ul>
                                    <li> Uniform Manifold Approximation and Projection (UMAP)
                                    <li> t-distributed stochastic neighbor embedding (tSNE)
                                  </ul>   
                                <li> Display settings
                                  <ul>
                                    <li> Plot type. The user is able to choose which of the available nonlinear dimensionality reduction spaces to visualize. Defaults to pca.
                                    <li> Dimensions. The number of dimensions to be plotted. Defaults to 2D.
                                    <li> Color by. Apply different coloring based on the available cell metadata columns. Defaults to orig_ident.
                                    <li> Size. The size of the plotted cells/dots. Defaults to 5.
                                    <li> Opacity. Define the level of the transparency of the plotted cells/dots. Defaults to 1.
                                    <li> Border width. The dot/cell perimeter thickness. Defaults to 0.5.
                                  </ul>   
                              </ol>
                          </p>
                       </div>
                   ')

umap_atac_output <- HTML('
                     <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <br>
                          <img src = "images/help_page/UMAPATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 24:</b> Output example </figcaption>
                          <br>
                          <p>
                            Each nonlinear dimensionality reduction space is visualized using a scatter plot. Each dot represents a cell, and each axis a nonlinear dimension. 
                            The particular illustration hosts a UMAP visualization in 2D.
                          </p>
                          <br>
                      </div>
                   ')

##############################D.E.A.############################################
dea_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                                The <b>MARKERS\' IDENTIFICATION</b> tab enables the user to identify marker genes (scRNA-seq and scATAC-seq) and marker peaks (scATAC-seq) by 
                                applying differential expression and differential accessibility analysis techniques respectively. This is a very crucial analysis step 
                                in single-cell analysis, since the respective results may guide procedures like cell type/state annotation, and identification of key transcriptional 
                                and regulatory programs that drive pathogenicity and/or development. 
                             </h4>')

dea_rna_input <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> Input parameters for marker genes identification </u><h3>
                                  <br>
                                  <img src = "images/help_page/DEA_RNA_parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 25:</b> DEA parameters for RNA-seq data </figcaption>
                                  <br>
                                  <p>
                                    <ol>
                                      <li> Test used: The user is able to choose amongst an extensive list of statistical tests and DEA methods. 
                                           The analysis is performed in a cluster-specific manner, where each cell-cluster\'s cells are tested against 
                                           all the other cells of the dataset. Defaults to Wilcoxon rank sum test.
                                        <ul>
                                          <li> Wilcoxon rank sum test. Identifies DEGs between two cell clusters by using a Wilcoxon Rank Sum test.
                                          <li> Likelihood-ratio test for single cell feature expression. Identifies DEGs between two cell clusters using a Likelihood-ratio test.
                                          <li> Standard AUC classifier. Identifies DEGs between two cell clusters by using Receiver operating characteristic (ROC) analysis. Creates 
                                          an AUC classifier for each gene, and tests the ability of the particular classifiew to separate two groups of cells. A value of 1 characterizes 
                                          the particular gene as a perfect classifier, implying upregulation. A value of 0 characterizes the particular gene as a perfect classifier, 
                                          implying downregulation. A value of 0.5 means that the gene has no predictive power.
                                          <li> Student\'s t-test. Identifies DEGs between two cell clusters using the Student\'s t-test.
                                          <li> MAST. Identifies DEGs between two cell clusters by using a hurdle model tailored to scRNA-seq data, by utilizing the MAST package.
                                          <li> DESeq2. Identifies DEGs between two cell clusters based a negative binomial distribution model, by utilizing the DESeq2 package. (slow operation)
                                        </ul>
                                      <li> DEA parameters
                                        <ul>
                                          <li> Base used for average logFC calculation. The user is able to choose what kind of log transformation will be applied to the average cell Fold Changes of 
                                          the DEA comparisons. Choices include (a) log with base e and (b) log with base 2. Defaults to log with base e.
                                          <li> Minimum % of expression. Only test genes that are detected to a minimum fraction of cells in either two groups. Defaults to 0.25.
                                          <li> Avg Log FC threshold. Limit testing to genes which show, on average, at least X-Fold change (logged) between the two groups of cells. Defaults to 0.25.
                                          <li> P-value threshold. Only report genes with a p-value lower than the user-defined limit. Defaults to 0.01.
                                        </ul>
                                    </ol>
                                  </p>
                                  <br>
                              </div>
                              
                              <br>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> DEA results </u><h3>
                                  <p>
                                      Marker genes data table. The particular, downloadable data table contains all the cluster-specific DEA results. 
                                      The information is stored as follow: <br>1st column: P value of the statistical test, 
                                      <br>2nd column: average log2 (or log(e)) Fold Change between the two groups of cells, 
                                      <br>3rd column: Percentage of cells of group1 that each gene is detected, 
                                      <br>4th column: Percentage of cells of group2 that each gene is detected, 
                                      <br>5th column: P adjusted value of the statistical test, 
                                      <br>6th column: cluster name, 7th column: gene name.
                                  </p>
                                  <img src = "images/help_page/DEA_RNA_results.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 26:</b> DEA results: table of marker genes </figcaption>
                                  <br>
                              </div>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <p>
                                    <h3><u>Heatmap visualization.<h3></u>
                                  <p>
                                  <p>
                                    A heatmap visualization of scaled expression values for the top 10 markers of each cell cluster. X-axis represents genes, 
                                    while y-axis represents cells. Cells are sorted by their cluster name. 
                                  <p>
                                  <br>
                                  <img src = "images/help_page/DEA_RNA_results_heat.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 26:</b> DEA results: heatmap of top-10 marker genes per cluster </figcaption>
                                  <br>
                              </div>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <p>
                                    <h3><u>Dotplot visualization.<h3></u>
                                  <p>
                                    Dot plot illustration of average scaled gene expression for the top 10 markers of each cell cluster. X-axis represents cell clusters, while y-axis represents gene markers. 
                                    The size of each dot indicates the gene detection percentage in each cluster.
                                  </p>
                                  <br>
                                  <img src = "images/help_page/DEA_RNA_results_dot1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 27:</b> DEA results: Dotplot of top-10 marker genes per cluster </figcaption>
                                  <br>
                              </div>
                              
                              <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <p>
                                    <h3><u>Volcano plot visualization.<h3></u>
                                  <p>
                                    Volcano plots is a very convenient way to visually summarize the DEA results between to compared cell groups. scanner implementation allows cluster specific Volcano plots.
                                    In the particular example, a volcano plot of differential expressed genes of cluster 2 is generated. X-axis depicts the average log2 Fold Change of each gene, while the y-axis 
                                    the respective -log10 p-values. Hovering over the visualization results allows the user to detect genes of interest, like for example CD14, which is a known marker of the 
                                    particular cell-type (CD14-Monocytes).
                                  </p>
                                  <br>
                                  <img src = "images/help_page/DEA_RNA_results_volcano.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 28:</b> DEA results: Volcano plot depicting the significant up and down regulated genes of the selected cluster </figcaption>
                                  <br>
                              </div>
                              ')

dea_rna_signature <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                            <br>
                            <h3><u> Feature plot visualization </u><h3>
                             A cell scatter plot in 2D reduced dimensional space that summarizes the average gene expression of particular markers or gene signatures
                              <ul>
                                <li> <b>Gene selection mode.</b> 
                                     <br>Visualization of particular genes. The user is able to adjust several parameters like <br>(a) the gene name to be plotted, <br>(b) which 2D reduced 
                                     space to use for plotting (available choices are describe in PCA/LSI and ADDITIONAL DIMENSIONALITY REDUCTION METHODS tabs), (c) showing or hiding the cluster names, 
                                     setting the <br>(d) maximum and the <br>(e) minimum expression value of the color scale. The plot will be generated after pressing the "Display Plot" button.
                                     
                                
                            <br>
                            <br>
                            <img src = "images/help_page/DEA_RNA_results_feat_genewise.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                            <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 29:</b> Parameters for feature/signature plot visualization </figcaption>
                            <br>
                            In the particular feature plot, a classical Monocyte marker (CD14) marks cluster 2, guiding the annotation of the particular cell group.
                            <br>
                                  <img src = "images/help_page/DEA_RNA_results_feat_gene1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 30:</b> Feature plot visualization of CD14 gene </figcaption>
                            <br>
                            
                               <li> <b>Gene signature mode</b>. <br>Visualization of sets of genes, as a combined gene expression signature. The user should initially define a gene signature name (Gene signature name), 
                               as also the members of this signature (Paste a list of genes). <b>Note: Gene names should be separated by entering a new line.</b> In the particular example, 
                               classical Bcell markers are defined. To calculate the signature score, the user should press the "Calculate signature score" button (Fig. 31). <br><br>

                               <img src = "images/help_page/DEA_RNA_results_feat_signature1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 31:</b> Creation of  B-cell markers specificsignature </figcaption>
                                  
                              To visualize the Gene signature, the user should activate the "Gene signature" button under the first option of the particular tab. Consequently, the user 
                              is able to choose his signature by using the "Select signature/numeric variable" sliding window. The rest of the options are described in the previous section. 
                              
                              <img src = "images/help_page/DEA_RNA_results_feat_signature_viz1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 31:</b> Vizualization of the signature  </figcaption>
                               
                              </ul>
                            </div>
                            
                            <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                              <h3><u> Violin plot vizualization </u></h3>
                              <p>
                              Violin plots of normalized gene expression for either individual genes ("Gene" button), or gene signatures ("Gene signature" button).
                              <br>
                              <img src = "images/help_page/DEA_RNA_results_violin1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 32:</b> Violin plots parameters </figcaption>
                              <br>
                              The user can choose either a gene ("Search for gene" sliding window), or a predefined gene signature of interest ("Select signature/numeric variable" sliding window). 
                              To generate the plot, the "Display plot" button should be pressed.
                              
                              <img src = "images/help_page/DEA_RNA_results_violin2.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 33:</b> Violin plots visualization </figcaption>

                            </div>
                            
                            <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                            <br>
                            <h3><u> Multi-feature visualization </u><h3>
                             Average expression visualization of couples of genes. This functionality allows the detection of marker co-occurrence and co-expression across cells and cell clusters.
                            <br>
                            <img src = "images/help_page/DEA_RNA_results_multi.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                            <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 34:</b> Parameters for multi-feature plot visualization </figcaption>
                            <br>
                            <p>
                              The user is able to select a pair of genes to be plotted together ("Select 1st feature" and "Select 2nd feature respectively"), 
                              while tweaking the color blending behavior is also applicable ("Select threshold for blending", defaults to 0.5). The rest of the 
                              user-defined parameters are the same as in "Feature plot" tab (previous section).
                            </p>
                            <img src = "images/help_page/DEA_RNA_results_multi_viz_combined.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                            <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 35:</b> Example of multi-feature plot visualization </figcaption>
                            <p>
                              In the particular example, CD3D, a classical T-cell marker, and GZMB, a classical Natural Killer (NK) cell marker are visualized using the multi-feature functionality. 
                              The first two Feature plots (first row) show the individual activity of each marker, the 3rd plot (lower-left) show the co-embedding/occurence of the two markers, which is 
                              typical for T-cells and NK cells, while  the forth plot (lower-right) depict the color scale values and mixed colors references.
                            </p>
                           </div>
                      ')

dea_atac_genes <- HTML('
                       <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> Differential Accessibility Analysis options (scATAC-seq) </u><h3>
                                  <br>
                                  <img src = "images/help_page/ATAC_markers_menu_genes.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 36:</b> DAA parameters for ATAC-seq data </figcaption>
                                  <br>
                                  <p>
                                   <u>Marker genes</u><br> 
                                   Marker gene detection using gene -score activity scores. The particular gene activity values are computed by aggregating the 
                                   accessibility signal along the regulatory space of each gene, and are considered very good predictors of gene expression. Based on this assumption, 
                                   differential expression analysis is performed as described above.
                                   <br>
                                   <br>
                                   The user is able to select among 3 statistical methods for differential expression testing:
                                        <ul>
                                          <li> Wilcoxon. Identifies DEGs between two cell clusters by using a Wilcoxon Rank Sum test.
                                          <li> Binomial. Identifies DEGs between two cell clusters using a binomial statistical test.
                                          <li> T-test. Identifies DEGs between two cell clusters using the Student\'s t-test.  
                                        </ul>
                                  </p>
                                  
                                  <br>
                                  
                                  <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> DAA - marker genes results </u><h3>
                                  <p>
                                      After the detection of cluster-specific marker genes, a downloadable data matrix depicting the DEA results is produced. 
                                      The 1st	column reports the chromosome location of each marker, the 2nd	and 3rd	columns report the start and end chromosomal 
                                      positions of each reported gene respectively, the 4th  column reports the gene strand, the 5th column reports the gene name, 
                                      the 6th column reports the gene index, the 7th	column reports the average log2  Fold Change, the 8th column reports the False 
                                      Discovery Rate of the applied statistical test, the 9th column reports the mean difference between the average normalized gene 
                                      activity scores of the two compared groups, and the 10th	column reports the marker cluster name.
                                  </p>
                                  <img src = "images/help_page/ATAC_marker_genes_results_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 37:</b> DAA results: table of marker genes </figcaption>
                                  <br>
                                 </div>
                                 
                                 <br>
                                 
                                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <h3><u> Heatmap visualization <h3></u>
                                  <p>
                                    The DEA results are visually summarized using a heatmap of average normalized gene activity scores of the top 10 gene markers of each cluster. 
                                    Both rows (clusters) and columns (genes) are clustered using a binary sorting procedure.
                                  </p>
                                  <img src = "images/help_page/ATAC_marker_genes_results_heatmap.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 37:</b> Heatmap of top-10 marker genes per cluster </figcaption>
                                 </div>
                        </div>
                       ')

dea_atac_peaks <- HTML('
                       <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> Differential Peak Accessibility Analysis options (scATAC-seq) </u><h3>
                                  <br>
                                  <img src = "images/help_page/ATAC_markers_menu_peaks.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 38:</b> DAA parameters for ATAC-seq data </figcaption>
                                  <br>
                                  <p>
                                   <u>Marker Peaks</u><br> 
                                   Marker peak detection using summarized chromatin accessibility counts in user-defined peak sets. 
                                   Marker peaks are calculated in a same fashion as described in marker gene detection in the section above. The main difference is that the 
                                   aggregation of the accessibility counts is applied in user defined peaks. Examples of such genomic regions are 
                                   <br>(a) the peaks generated by cellranger-atac count ("peaks.bed" file, in the outs/ directory of cellranger-atac count run), 
                                   <br>(b) collections of regulatory regions such as promoters and enhancers, derived by ATAC-seq, DNase-seq, FAIRE-seq, MNase-seq 
                                   and ChIP-seq datasets, or 
                                   <br>(c) fixed genomic segments called bins/tiles.
                                   <br>
                                   The particular input file should be stored in "bed" format, including the chromosome name, and the start and end genomic positions for each of the stored peak in the first 3 columns.
                                   
                                   <div class="col-md-12 scrollable">
                                   <br>
                                   <img src = "images/help_page/peak_file_head.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                   <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 39:</b> BED example file </figcaption>
                                   <br>
                                   </div>
                                   
                                  <br>
                                  
                                  <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <h3><u> DAA - marker peaks results </u><h3>
                                  <p>
                                      After the detection of cluster-specific marker peaks, a downloadable data matrix depicting the differential accessibility analysis results is produced. 
                                      The 1st	column reports the chromosome location of each marker, the 2nd column reports the peak index, the 3rd		and 4th	columns report the start and end 
                                      chromosomal positions of each reported peak respectively, the 5th	column reports the average log2  Fold Change, the 6th column reports the False Discovery 
                                      Rate of the applied statistical test, the 7th column reports the mean difference between the average normalized chromatin accessibility values of the two 
                                      compared groups, and the 8th	column reports the marker cluster name
                                  </p>
                                  <img src = "images/help_page/ATAC_marker_peaks_results_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 40:</b> DAA results: table of marker peaks </figcaption>
                                  <br>
                                 </div>
                                 
                                 <br>
                                 
                                 <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <h3><u> Heatmap visualization <h3></u>
                                  <p>
                                    The differential accessibility analysis results are visually summarized using a heatmap of average normalized peak accessibility of the top 10 markers peaks of each cluster. 
                                    Both rows (clusters) and columns (genes) are clustered using a binary sorting procedure.
                                  </p>
                                  <img src = "images/help_page/ATAC_marker_peaks_results_heatmap.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 41:</b> Heatmap of top-10 marker peaks per cluster </figcaption>
                                  <br>
                                 </div>
                        </div>
                       ')

dea_atac_activity <- HTML('
                          <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <br>
                          <h3><u>Feature plots of gene activity scores</u></h3>
                          
                          <img src = "images/help_page/CD14.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 42:</b> Feature plot depicting gene activity scores per cell for a selected gene. </figcaption>
                          <br>
                          <p>
                            The user is able to determine the gene of interest ("Select a gene" sliding window), and in which reduced dimensional space the visualization will take place 
                            ("Plot type" sliding window, defaults to UMAP). To generate the visualization, a respective "Display plot" button is available.
                            In the particular example, CD14 (a classical Monocyte marker) gene activity scores are plotted, revealing that C2 cluster could be annotated as Monocytes.
                          </p>
                          </div>
                          ')

##############################Cell cycle########################################
cellCycle_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                                The <b>CELL CYCLE PHASE ANALYSIS</b> tab enables the user identify effects of cell cycle heterogeneity in 
                                a scRNA-seq dataset, by calculating cell cycle phase scores based on canonical markers. If the analysis shows 
                                a clear cell grouping based on these markers, then the particular effect should be considered a putative unwanted bias.
                             </h4>')

cell_cycle_rna <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <br>
                        <p>
                          The user can choose in which reduced space the cell cycle phase scores will be visualized ("Plot type" sliding window, defaults to PCA).
                        </p>
                        
                        <img src = "images/help_page/CellCycleRNA.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                        <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 43:</b> Cells visualized PCA space and colored by cell phase identity </figcaption>
                        <br>
                        <p>
                          PCA scatter plot showing a reasonable mix of cell-cycle phases along the first two principal components.
                        </p>
                        <br>
                        <img src = "images/help_page/CellCycleRNA_barplot.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                        <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 44:</b> Per-cluster bar plots showing the percentages of the 3 basic cell-cycle phases. </figcaption>
                        <br>
                        </div>
                       ')

##############################Functional########################################

functional_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                                <b>The FUNCTIONAL/MOTIF ENRICHMENT ANALYSIS</b> tab enables the user to perform cluster-specific GO and pathway enrichment analysis 
                                based in scRNA-seq gene markers, as also cluster-specific Transcription Factor (TF) motif enrichment analysis, based in scATAC-seq peak markers. 
                                The particular analysis determines biological processes and regulatory programs that define cell-type and cell-state identity.
                             </h4>
                             ')

grpofiler_tab_rna <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Functional enrichment analysis (with gProfiler) </u><h3>
                                <p>
                                <ol>
                                  <li><b> Options for input list </b><br>
                                       The user is enabled to choose a cluster of interest to perform functional enrichment analysis ("Input list of genes" sliding window). 
                                       In order to subset the marker gene list to the most significant features, a significant method and threshold is applied. In particular, 
                                       the user chooses a statistical method (p-value, adjusted p-value, or power, defaults to "p-value"), a significant threshold (defaults to 0.01), 
                                       the "Direction of deregulation" (defaults to "up-regulated"), and the average log Fold Change threshold ("Log FC threshold", defaults to 0.25).
                                  <li><b> Options for enrichment analysis </b><br>
                                       The user is able to select among a list of Gene Ontologies (Molecular Function, Cellular Component, and Biological Process), 
                                       Biological Pathways (KEGG Pathways, REACTOME, WikiPathways), Regulatory Motifs (TRANSFAC, miRTarBase), Protein Databases (CORUM, Human Protein Atlas (HPA)), 
                                       and Human Phenotype Ontologies (Human Phenotype Ontology) by using the "Select datasources" sliding window.
                                       Consequently, the user should also choose the under study organism ("Select organism", defaults to Homo sapiens (Human)), 
                                       as also the statistical method of the enrichment analysis ("Correction method for multiple testing", defaults to Bonferoni correction.
                                       Finally, a statistical threshold should be also applied ("Significane for enriched terms", defaults to 0.05).
                                       <br>By pressing the "OK" button, a functional enrichment analysis is initialized. An additional feature of this module, is to perform the 
                                       analysis to Flame web application, a web tool for Functional and Literature Enrichment analysis of multiple sets.
                                </ol>
                                
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/gprofiler_parameters1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <br>
                               <img src = "images/help_page/gprofiler_parameters2.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 45:</b> Form containing the input parameters for functional enrichment analysis. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                                <h3><u> Functional enrichment analysis results </u><h3>
                                <img src = "images/help_page/gprofiler_results_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 46:</b> Table of enriched terms. </figcaption>
                                <br>
                                <p>
                                  The particular analysis generates a downloadable data table, that depicts cluster-specific significant enriched terms. The table 
                                  encapsulates information regarding the cluster name (1st column), the level of significance (2nd column), the term size (3rd column), 
                                  the query size (4th column), the number of marker genes present in the particular term (5th column), the term id (6th column), the data 
                                  resource for the enrichment analysis (7th column)  and a comma-separated list of marker genes present in the particular term (8th column).
                                </p>
                                
                                <img src = "images/help_page/gprofiler_manhattan.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 47:</b> Manhattan plot of enriched terms. </figcaption>
                                <br>
                                <p>
                                The particular analysis results are also visually summarized in an interactive Manhattan plot. The user is able to hover over the results, in order to identify the most significantly enriched terms.
                                </p>
                              </div>
                              
                              </div>') 


motif_tab_atac <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Motif enrichment analysis</u><h3>
                                <p>
                                The user is enabled to choose motif database, in order to identify enriched motifs in all the identified cell clusters ("Motif set" sliding window, defaults to Cisbp). 
                                Options include: "JASPAR2016", "JASPAR2018", "JASPAR2020" which correspond to the 2016, 2018 or 2020 version of JASPAR motifs, "cisbp", "encode", and "homer".
                                The user should also define the Log2 Fold Change and False Discovery Rate thresholds (default to 0.25 and 0.01 respectively) that will limit the analysis output 
                                to the most significant results.
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/Motif_Enrichment_Menu.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 48:</b> Form containing the input parameters for motif enrichment analysis. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                              <p>
                                <h3><u> Motif enrichment analysis results </u><h3>
                                <img src = "images/help_page/Motif_Enrichment_results_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 49:</b> Table of enriched motifs </figcaption>
                                <br>
                                
                                <img src = "images/help_page/Motif_Enrichment_results_heatmap.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 50:</b> Manhattan plot of enriched terms. </figcaption>
                                <br>
                                <p>
                                The particular analysis results are visually summarized in a heatmap illustration. For each cluster (rows), the top-10 most significant enriched motifs are depicted.
                                </p>
                              </div>
                              
                              </div>
                       ')

#############################CIPR###############################################
annot_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                            The <b>CLUSTERS\' ANNOTATION</b> tab enables the user to perform scRNA-seq cell-cluster annotation and to assign known cell-type identity to 
                            the analyzed dataset. This analysis pipeline utilizes the most variable features of the dataset, scores the gene expression profiles 
                            of the identified clusters against single-cell data resources, and finally assigns known cluster-specific cell-type legends.    
                         </h4>
                      ')

annot_cipr_rna <- HTML('<div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Cluster annotation (with CIPR)</u><h3>
                                <p>
                                Annotation parameters
                                <ol>
                                  <li> Reference dataset. A list of available cell-type annotation references. Available choices are: 
                                  (a) "Blueprint-Encode", (b) "Primary Cell Atlas", (c) "DICE", (d) "Hematopoietic diff" and (e) "Presorted RNA seq" for human, 
                                  and  (a) "ImmGen and (b) "Presorted RNAseq" for mouse.
                                  <li> Keep top Nth % of variable genes in reference. The percentage of most variable features to retain for the cluster annotation analysis. Defaults to 100%.
                                  <li> Select method for comparisons. The comparison method to be applied for matching the cell markers with the matching genes in the reference dataset. 
                                       The available options are (a) LogFC dot product, (b) LogFC Spearman\'s correlation and (c) LogFC Pearson\'s correlation. Alternatively, the nonlinear and 
                                       linear correlations between the input and reference datasets can be drawn by considering all genes regardless of differential expression status. This can 
                                       be applied by using (d) Pearson and (e) Spearman correlation in all the most variable genes respectively.
                                </ol>
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/cipr_parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 51:</b> Form containing parameters for the annotation. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                              <p>
                                <h3><u> Cluster\'s annotation analysis results </u><h3>
                                <img src = "images/help_page/cipr_results_table_top5.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 52:</b> Table of top preddicted annotations. </figcaption>
                                <br>
                                <p>
                                  The particular analysis generates a downloadable data table, that depicts thje top 5 cluster-specific cell-type annotations. The table encapsulates 
                                  information regarding the cluster name (1st column), the assigned cell-type annotation (2nd column), the assigned cell-type annotation id (3rd column), 
                                  the assigned cell-type annotation log name (4th column), assigned cell-type annotation description (5th column), the assigned cell-type annotation identity 
                                  score (6th column), the annotation index (7th column) and the assigned cell-type annotation z-score (8th column).
                                </p>
                                <br>
                                <img src = "images/help_page/cipr_results_dotplot_top5.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 53:</b> Dotplot of top preddicted annotations. </figcaption>
                                <br>
                                <p>
                                The particular analysis results are also visually summarized in an interactive Dot plot. The user is able to hover over the results, in order to 
                                identify the top 5 annotated cell-types for each of the analyzed clusters.
                                </p>
                              </div>
                              
                              </div>
                       ')

#######################Trajectory###############################################

traj_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                            The <b>TRAJECTORY ANALYSIS</b> tab enables the user to perform pseudotemporal ordering of scRNA-seq and scATAC-seq cells, in order to identify developmental traces and lineages of cell differentiation.    
                         </h4>
                      ')

traj_rna_slingshot <- HTML('
                          <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Trajectory analysis (with Slingshot)</u><h3>
                                <p>
                                Trajectory parameters
                                <ul>
                                  <li> Dimensionality reduction method. A sliding window of available dimensionality reduction spaces to use for the trajectory 
                                       inference analysis. Available options are derived by the previous analysis steps in the "ADDITIONAL DIMENSIONALITY 
                                       REDUCTION METHODS" tab.
                                  <li> Number of dimensions to use. The number of dimensions to incorporate during the trajectory inference procedure. Defaults to 3.
                                  <li> Initial and final state selection. The user should set the root and the terminal cell-cluster of the underlying differentiation process.
                                </ul>
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/traj_RNA_Parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 54:</b> Form containing parameters for the trajectory analysis. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                              <p>
                                <h3><u> Trajectory analysis results </u><h3>
                                Stracture overview. <br>
                                <img src = "images/help_page/traj_RNA_Structure_view.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 55:</b> A cell scatter plot in the UMAP space that depicts the identified 
                                lineages of the underlying differentiation process. </figcaption>
                                <br>
                              </p>
                                 <p>
                                  Lineage view. <br>
                                  The user chooses an inferred lineage ("Select lineage" sliding window and press "OK" button) in order to visualize the pseudotemporally-ordered cells. 
                                  The pseudotime values indicate the cell ranking along the underlying differentiation process.
                                </p>
                                <br>
                                <img src = "images/help_page/traj_RNA_Lineage_view.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 56:</b> Scatterplot with cells colored by pseudotime values. </figcaption>
                                <br>
                              </div>
                              
                              </div>
                       ')

traj_atac_slingshot <- HTML('
                          <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Trajectory analysis (with Slingshot)</u><h3>
                                <p>
                                Trajectory parameters
                                <ul>
                                  <li> Number of UMAP dimensions to use. The number of dimensions to incorporate during the trajectory inference procedure. Defaults to 3.
                                  <li> Initial and final state selection. The user should set the root and the terminal cell-cluster of the underlying differentiation process.
                                </ul>
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/atac_slingshot_Parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 57:</b> Form containing parameters for the trajectory analysis. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                                 <p>
                                  Lineage view. <br>
                                  For the particular visualization, the user chooses an inferred lineage ("Select lineage" sliding window and press "Display pseudotime ranking" button) 
                                  in order to visualize the pseudotemporally-ordered cells. The pseudotime values indicate the cell ranking along the underlying differentiation process.
                                </p>
                                <br>
                                <img src = "images/help_page/atac_slingshot_lineage1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 58:</b> Scatterplot with cells colored by pseudotime values. </figcaption>
                                <br>
                              </div>
                              
                              </div>
                       ')

##############################L-R###############################################

lr_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                            The <b>LIGAND - RECEPTOR ANALYSIS</b> tab enables the user to perform ligand activity analysis and target receptor predictions on cell-clusters of interest. 
                            Consequently, predictions of active ligands and their target features is performed in a pairwise manner (cluster-to-cluster).    
                         </h4>
                      ')

lr_rna_nichnet <- HTML('
                          <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                               <br>
                               
                               <div class="col-md-6">
                               <h3><u> Ligand - Receptor analysis (with nichenetR)</u><h3>
                                <p>
                                L-R analysis parameters<br>
                                The user defines a cell-cluster pair in order to perform ligand - receptor analysis ("Ligand expression cluster" and "Receptor expressing cluster" sliding windows).
                                </p>
                               </div>
                               
                               <div class="col-md-6">
                               <img src = "images/help_page/lr_parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                               <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 59:</b> Form containing parameters for L-R analysis. </figcaption>
                               <br>
                               </div>
                              <br>
                              
                              <div style="font-size:16px;">
                              <p>
                                <h3><u> L-R analysis results </u><h3>
                                All interactions. The results of the particular analysis are visually summarized by generating a heatmap of prior interaction potential between 
                                the inferred active ligands (y-axis) and paired receptors (x-axis).
                                
                                <img src = "images/help_page/lr_all_interactions_results.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 60:</b> Heatmap of all preddicted interactions. </figcaption>
                                <br>
                              </p>
                                 <p>
                                  Curated interactions. The results that include interactions documented in the literature and in publicly available databases are visually 
                                  summarized by generating a heatmap of "bona fide" prior interaction potential between the inferred active ligands (y-axis) and paired receptors (x-axis).
                                </p>
                                <br>
                                <img src = "images/help_page/lr_Curated_interactions.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 61:</b> Heatmap containing only "bona fide" interactions. </figcaption>
                                <br>
                              </div>
                              
                              </div>
                       ')

#########################GRN####################################################
grn_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                       The <b>GENE REGULATORY NETWORK ANALYSIS</b> tab enables the user to perform GRN, positive regulator inference analysis and visualization using scRNA-seq and scATAC-seq datasets, in single-cell resolution.      
                       </h4>
                      ')

grn_tab_rna <- HTML('
                    <div class="col-md-12 scrollable" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                      <h3><u> Gene Regulatory Network analysis</u><h3>
                      GRN input parameters<br>
                      <p>
                        Prepare files for pyscenic. The user should select the appropriate genome build used during the counting procedure of the scRNA-seq raw files. The available choices are: mm10, hg19, and hg38.
                        After scAnner prepares the respective RDS and loom files, the user is able to download them in order to perform a local run of pyscenic pipeline by using the provided pyscenic_local.py python 
                        script (see the respective help page).
                      </p>
                      <img src = "images/help_page/GRN_menu_prepare_input_pyscenic.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 62:</b> Preparation of input files for pyscenic. </figcaption>
                      <br>
                      <br>
                      <p>
                        Analyze pyscenic output. After the local pyscenic run is completed (this might take a while), the user is able to upload the respective pyscenic results 
                        (auc loom file), as also the RDS file used in the above step in order to proceed to regulon analysis and visualization.
                      </p>
                      <img src = "images/help_page/GRN_anlyze_pyscenic_files.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 63:</b> Analysis of pyscenic results. </figcaption>
                      <br>
                      <br>
                      <p>
                        Visualization options. The user is enabled to preview the generated regulon AUC and RSS values ("Regulons - display" sliding window), as also the 
                        number of the most specific regulons per analyzed cluster ("Display top regulons" options, defaults to 10). The results will be visualized 
                        after pressing the "Plot" button.
                      </p>
                      <img src = "images/help_page/GRN_menu_visualization.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 64:</b> Visualization of regulons. </figcaption>
                      <br>
                      <br>
                      <h3><u> Gene Regulatory Network output</u><h3>
                      <p>
                        The particular analysis generates a downloadable data table, that depicts the per cluster average AUC regulon values, that represent the activity of the inferred GRNs. 
                        Rows represent regulons, while columns represent cell clusters.Analysis results are also visually summarized in a heatmap illustration of z-scored AUC values. For each 
                        cluster (rows), the top-10 most specific regulons are depicted.
                        <br>
                        <img src = "images/help_page/GRN_results_full.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 65:</b> Table and heatmap of regulons depicting AUC values. </figcaption>
                      <br>
                      <br>
                      <p>
                    </div>
                    ')

grn_tab_atac <- HTML('<h3><u> Gene Regulatory Network analysis</u><h3>
                      Analysis options<br> 
                      The user is enabled to set the "FDR threshold" (defaults to 0.1) and the "correlation threshold" (defaults to 0.7) used during the procedure of identifying cluster-specific positive regulators. 
                      The particular values are used to determine significant correlation between motif accessibility and gene activity scores across cells.
                      <br>
                      <br>
                      <img src = "images/help_page/GRN_menu.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 66:</b> GRN analysis parameters. </figcaption>
                      <br>
                      <br>
                      <h3><u> Gene Regulatory Network results</u><h3>
                      <p>
                        The particular analysis generates a downloadable data table, that depicts the per cluster average motif deviation scores, that 
                        represent the activity of the inferred positive regulators. Rows represent positive regulators, while columns represent cell clusters. 
                        Positive regulator activity is also visually summarized in a heatmap illustration. The top-10 most active regulators are shown.
                      </p>
                      <br>
                      <img src = "images/help_page/GRN_Positive_Regulators_Table.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 67:</b> Positive regulators table. </figcaption>
                      <br>
                      <img src = "images/help_page/GRN_Positive_Regulators_Heatmap.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 68:</b> Heatmap of positive regulators table. </figcaption>
                                
                      <br>
                      <br>
                      <p>
                        Peak to gene links. The particular analysis also generates peak to gene linkages (P2G links), which are inferred by using correlation analysis between enhancer 
                        accessibility and gene activity scores. The results are included in a downloadable data table, that encapsulates information regarding the P2G link name (row names), 
                        the peak index (1st column), the gene index (2nd column), the correlation value (3rd column), correlation FDR (4th column), the chromatin accessibility variance (5th column), 
                        the gene activity variance (6th column), the peak name (7th column) and the gene name (8th column).
                      </p>
                      <br>
                      <img src = "images/help_page/GRN_Peak_To_Gene_Links_Table.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 69:</b> Table of peak to gene links. </figcaption>
                                
                      <br>
                      <br>
                      <p>
                        Peak x motif occurrence matrix. The particular analysis also generates a data table which informs the user for transcription factor motif (columns) presence in each peak (rows). 
                        This matrix can be combined with the P2G linkage, the marker peak and the marker gene matrix, in order to form GRNs.
                      </p>
                      <br>
                      <img src = "images/help_page/GRN_Peak_Motif_Table.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 70:</b> Peak x motif occurrence matrix. </figcaption>
                     ')

##############################Tracks############################################
tracks_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">
                        The <b>TRACKS</b> tab enables the user to generate cluster-specific genome browser track snapshots of scATAC-seq signal. Local chromatin accessibility signal is visualized in the regulatory space of particular genes.      
                       </h4>
                      ')

tracks_tab_atac <- HTML('
                        <h3><u> scATAC-seq tracks options</u><h3>
                        The user is able to select a gene ("Select a gene" ), as also the regulatory space upstream ("BP upstream") and downstream ("BP downstream") 
                        relative to the selected gene genomic position, for visualizing the scATAC-seq signal.
                        <br>
                        <img src = "images/help_page/Tracks_Menu.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 71:</b> scATAC-seq track options. </figcaption>
                        <br>
                        <h3><u> scATAC-seq tracks visualization</u><h3>
                        <p>
                        A genome browser tracks is generated, based on the user options. The chromatic accessibility signal is aggregated and organized per cluster. 
                        Additionally, Peaks, co-accessibility links, and gene annotations are also generated. In the particular example, the 50kb regulatory space of 
                        CCL2 is visualized. scATAC-seq is specifically enriched in C2 cluster, in CCL2 promoters, as also in distal regulatory elements (peaks) thar putatively 
                        regulate the particular gene\'s promoter (Co-accessibility).
                        </p>
                        <br>
                        <img src = "images/help_page/Tracks_Plot.PNG">
                                <figcaption style = "font-size:14px" class="figure-caption text-left"><b>Figure 72:</b> scATAC-seq tracks for CCL2 gene. </figcaption>
                        ')

################################################################################
################################################################################