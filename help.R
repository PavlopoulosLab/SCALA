###############################Example files####################################
examples_help <- HTML(
                  '<h4>Example matrix input file for scRNA-seq:</h4>
                  <li><a href="example_files.../" download> PBMC 3K Cell-Gene count matrix (*.txt) </a></li>
                  <h4>Example 10x input files for scRNA-seq:</h4>
                  <li><a href="example_files.../" download> PBMC 3K barcodes file (*.tsv.gz) </a></li>
                  <li><a href="example_files/..." download> PBMC 3K features file (*.tsv.gz) </a></li>
                  <li><a href="example_files/..." download> PBMC 3K matrix file (*.mtx.gz) </a></li>
                  
                  <h4>Example files for GRN analysis in scRNA-seq PBMC 3K:</h4>
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
file_qc_tab_intro <- HTML('<h4 style = "line-height: 1.5; text-align:center; background-color: #ffffff; border: 1px solid #222d32; border-radius: 5px;">The <b>
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