###############################Example files####################################
examples_help <- HTML(
                  '<h4>Example matrix input file for scRNA-seq:</h4>
                  <li><a href="example_files.../" download> PBMC 3K Cell-Gene count matrix (*.txt) </a></li>
                  <h4>Example 10x input files for scRNA-seq:</h4>
                  <li><a href="example_files.../" download> PBMC 3K barcodes file (*.tsv.gz) </a></li>
                  <li><a href="example_files/..." download> PBMC 3K features file (*.tsv.gz) </a></li>
                  <li><a href="example_files/..." download> PBMC 3K matrix file (*.mtx.gz) </a></li>
                  
                  <h4>Example files for GRN analysis:</h4>
                  <li><a href="example_files/..." download> PBMC 3K seurat object (*.RDS) </a></li>
                  <li><a href="example_files/..." download> Pyscenic output file (*.loom) </a></li>
                  
                  <h4> Example files for scATAC-seq: </h4>
                  <li><a href="example_files/..." download> PBMC arrow file (*.arrow) </a></li>
                  <li><a href="example_files/..." download> PBMC peakset file (*.bed) </a></li>
                 ')

###########################Upload###############################################
file_upload_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>Data upload tab</b> enables 
the user to upload scRNA-seq or scATAC-seq datasets to scAnner, 
in order to initialize single-cell analysis.<br> ScAnner provides the option to upload different file input formats. </h4>')

file_upload_txt <- HTML('<div class="col-md-4">
                          <img src = "images/help_page/DATA_INPUT_scRNA_count_matrix.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 1: </b> The File Upload form for cell-gene count matrices </figcaption>
                          </div>
                          
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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

file_upload_10x <- HTML('<div class="col-md-4">
                          <img src = "images/help_page/DATA_INPUT_scRNA_10x.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 2: </b> The File Upload form for 10x input files </figcaption>
                          </div>
                          
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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

file_upload_arrow <- HTML('<div class="col-md-4">
                          <img src = "images/help_page/DATA_INPUT_scATAC-seq.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 3: </b> The File Upload form for arrow files </figcaption>
                          </div>
                          
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3> 1. Upload an arow file </h3>
                        <p>
                        <ol>
                        <li> Project name: User defined project name.
                        <li> File Upload: Click the <b>Browse</b> button of the upload form to select a single-cell ATAC-seq arrow file. Arrow files is a file type that stores all of the scATAC-seq data 
                        associated with an individual sample. Arrow files can be created by coupling scATAC-seq fragment files generated by cellranger-atac, and create_arrow_file.R Rscript (provided), as described in 
                        scAnner’s github page. The particular operation should be applied locally.<br>
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

file_upload_tab_new_project <- HTML('<div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px;">
<h4 style = "line-height: 1.5;">
                                      When you try to upload a new dataset at the same time that another one is already loaded, you will be prompted to discard the current working dataset. If you select <b>Yes</b> your 
                                      progress will be deleted, note that this is an irreversible action. After that you can start a new project by submitting the new input files.<br>
                                      If you wish to continue working on your current project slect <b>No</b>. 
                                    </h4>
                                    </div>')

file_upload_metadata_RNA <- HTML('
                        <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        <h3><u> scRNA-seq metadata table </u></h3>
                        <p>
                          If the data upload was successfull, a scRNA-seq cell metadata table will appear. 
                          Project name, number of Unique Moleculare Identifiers (UMIs) per cell, number of detected features/genes per cell, percentage of mitochondrial RNA (percent.mt) 
                          and cell id are stored in the 1st, 2nd, 3rd, 4th and 5th column respectively. Additional columns will be created during the next steps of the analysis, and the cell 
                          metadata table will be updated automatically.The user is also able to download the metadata table using the <b>Save table</b> button.
                        </p>
                        
                        </div>
                        <div class="col-md-12">
                          <img src = "images/help_page/DATA_INPUT_scRNA_Metedata_table.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 4: </b> RNA metadata table </figcaption>
                          </div>
                        ')

file_upload_metadata_ATAC <- HTML('
                        <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                        <div class="col-md-12">
                          <img src = "images/help_page/DATA_INPUT_scATAC-seq_Metadata_table_merged.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 5: </b> ATAC metadata table </figcaption>
                          </div>
                        ')

#######################################QC#######################################
qc_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>
                          Quality control tab</b> enables the generation of Quality Control (QC) plots for scRNA-seq and scATAC-seq datasets,<br> 
                          in order to examine the quality of the single-cell experiment. </h4>')

rna_qc <- HTML('
                          
                <div class="col-md-12"> 
                  <img src = "images/help_page/QC_pre_params.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 6: </b> Quality control parameters and plots before filtering </figcaption>
                </div>
                 
                 <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                          
                <div class="col-md-12"> 
                  <img src = "images/help_page/QC_post.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 7: </b> Quality control plots post filtering </figcaption>
                </div>
                 
                 <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                 <h3><u> Inspection of filtering criteria </u></h3>
                 <p>
                    The particular visualization allows the user to explore cell count metrics after the filtering procedure. The user is able to readjust the defined criteria if the results 
                    are not satisfying. The QC metrics included in this figure are described in the previous tab.
                 </p>
                 </div>
               ')

atac_qc <- HTML('
                <div class="col-md-12"> 
                  <img src = "images/help_page/QC_all_ATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 8: </b> Quality control plots after soft filtering </figcaption>
                </div>
                 
                 <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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

rna_normalization_param <- HTML('<div class="col-md-4">
                          <img src = "images/help_page/DATA_NORMALIZATION_AND_SCALING_parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 9: </b> Parameters for the normalization and scaling procedure in RNA seq data </figcaption>
                          </div>
                          
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                        
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
                          
                <div class="col-md-12"> 
                  <img src = "images/help_page/DATA_NORMALIZATION_AND_SCALING_mvgs.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <br>
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 10: </b> Most variable genes </figcaption>
                </div>
                 
                 <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                        <div class="col-md-12"> 
                          <img src = "images/help_page/PCA_selection.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <img src = "images/help_page/PCA_slow_results.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 11:</b> Determination of optimal number of principal components </figcaption>
                       </div>
                       
                       <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                        <div class="col-md-12"> 
                          <img src = "images/help_page/PCA_exploration.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 12:</b> Exploration of principal components </figcaption>
                       </div>
                       
                       <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <p>
                            The particular visualization depicts<br> (a) the loading scores of the top genes of a PC of interest (top 30 features for the particular example), and<br> 
                            (b) a heatmap of scaled counts of the top loadings of the PC of interest, across cells.
                          </p>
                       </div>
                        ')

pca_lsi <- HTML('
                <div class="col-md-6"> 
                          <img src = "images/help_page/PCA_lsi.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 13:</b> LSI options </figcaption>
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
                             <div class="col-md-4"> 
                              <img src = "images/help_page/Clustering Parameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 14:</b> Clustering options </figcaption>
                            </div>
                
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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

clustering_rna_output <- HTML('<div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering Results_1.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 15:</b> Table of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results can be explored by downloading the respective data table, which encapsulates information regarding the cluster name, 
                                    the number of cells in each cluster, and the percentage of total cells that is included in each cluster. The particular information is stored 
                                    in the 1st, 2nd and 3rd columns of the aforementioned data table respectively.
                                  </p>
                                  <br>
                              </div>
                              
                              <br>
                              
                              <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering Results_2.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 16:</b> Barplot of the identified clusters </figcaption>
                                  <br>
                                  <p>
                                    The clustering results are also depicted in a stacked barplot of percentages of cells in each of the identified clusters.
                                  </p>
                                  <br>
                              </div>
                              
                              <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/Clustering SNN_graph.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 17:</b> The shared nearest neighbor graph (SNN) used in clustering </figcaption>
                                  <br>
                                  <p>
                                    The user is also able to explore the formation of the SNN graph that led to the identification of the final cluster set.
                                  </p>
                                  <br>
                              </div>
                              ')

clustering_atac_input <- HTML('
                             <div class="col-md-4"> 
                              <img src = "images/help_page/ClusteringOptionsATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 18:</b> Clustering options </figcaption>
                            </div>
                
                        <div class="col-md-5" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                           <p>
                              <h3>scATAC-seq clustering, enables the same methodology as in described in scRNA-seq.</h3><br> 
                              The user is able to define<br> (a) the number of LSI dimensions to use 
                              (defaults to 30), as also<br> (b) the clustering resolution (defaults to 0.6).
                          </p>
                       </div>
                             ')

clustering_atac_output <- HTML('<div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/ClusteringResultsATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 19:</b> Table of the identified clusters </figcaption>
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
                              
                              <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                                  <br>
                                  <img src = "images/help_page/ClusteringBarplotATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                                  <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 20:</b> Table of the identified clusters </figcaption>
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
                      <div class="col-md-4"> 
                              <img src = "images/help_page/UMAPparameters.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 21:</b> Options </figcaption>
                            </div>
                        
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                                    <li> Plot type. The user is able to choose which of the available nonlinear dimensionality reduction spaces to visualize. Default to pca.
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
                        <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <br>
                          <img src = "images/help_page/UMAP.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 22:</b> Output example </figcaption>
                          <br>
                          <p>
                            Each nonlinear dimensionality reduction space is visualized using a scatter plot. Each dot represents a cell, and each axis a nonlinear dimension. 
                            The particular illustration hosts a UMAP visualization in 2D.
                          </p>
                          <br>
                      </div>
                   ')

umap_atac_input <- HTML('
                        <div class="col-md-4"> 
                              <img src = "images/help_page/UMAPparametersATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                              <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 23:</b> Options </figcaption>
                            </div>
                        
                        <div class="col-md-8" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
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
                                    <li> Plot type. The user is able to choose which of the available nonlinear dimensionality reduction spaces to visualize. Default to pca.
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
                     <div class="col-md-12" style="background-color: #ffffff; border: 1px solid #222d32; border-radius: 15px; font-size:16px;">
                          <br>
                          <img src = "images/help_page/UMAPATAC.PNG" style="border: 1px solid #222d32; border-radius: 15px;">
                          <figcaption style = "font-size:14px" class="figure-caption text-center"><b>Figure 24:</b> Output example </figcaption>
                          <br>
                          <p>
                            Each nonlinear dimensionality reduction space is visualized using a scatter plot. Each dot represents a cell, and each axis a nonlinear dimension. 
                            The particular illustration hosts a UMAP visualization in 2D.
                          </p>
                          <br>
                      </div>
                   ')
##############################Tracks############################################

