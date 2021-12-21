#!/bin/bash

wget -P scenic_helper_files https://zenodo.org/record/3260758/files/ligand_target_matrix.rds
wget -P scenic_helper_files https://zenodo.org/record/3260758/files/lr_network.rds
wget -P scenic_helper_files https://zenodo.org/record/3260758/files/weighted_networks.rds
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-10species.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-10species.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl
wget -P scenic_helper_files https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
wget -P scenic_helper_files https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/mm_mgi_tfs.txt
wget -P scenic_helper_files https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt


wget -P . http://bib.fleming.gr:8084/SCANNER/example_files/example_form_files.tgz
tar xvf example_form_files.tgz


wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/matrix.mtx.gz
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/barcodes.tsv.gz
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/exampleMatrix.zip
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/processed_seurat_object-2021-12-19.zip
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/arrowFile.zip
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/PBMCs_human_signac_peaks.zip
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/features.tsv.gz
wget -P www/example_files/ http://bib.fleming.gr:8084/SCANNER/example_files/auc.hg19.zip

