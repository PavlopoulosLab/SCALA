#!/bin/bash

#NOTE 1: RUN BEFORE installing the R libraries
#NOTE 2: RUN AS ROOT/SUDO


#Check if conda exists in path and if not, download and install miniconda3
CONDA_CHK=$(which conda)
if [[ $CONDA_CHK  == ""  ]]
then
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda_scanner/
	conda_path="/opt/conda_scanner/"
	rm Miniconda3-latest-Linux-x86_64.sh
else
	conda_path=$(which conda)
	conda_path=$(echo $conda_path | sed 's/condabin\/conda//g' | sed 's/bin\/conda//g')
fi


# CREATE A CONDA ENVIRONMENT FOR WHERE PYSCENIC, MACS2 AND PHATE WILL BE INSTALLED
$conda_path/bin/conda create -y -n pyscenic python=3.7
$conda_path/bin/conda activate pyscenic
$conda_path/bin/conda install -y -c anaconda cytoolz
$conda_path/envs/pyscenic/bin/pip install pyscenic
$conda_path/envs/pyscenic/bin/pip install macs2
$conda_path/envs/pyscenic/bin/pip install phate

#install some dependencies that are needed by a few of the R packages
apt-get install -y gsl-bin libgsl23 libgslcblas0 libgsl-dev
apt-get install -y hdf5-tools hdf5-helpers libhdf5-dev
apt-get install -y gdal-bin libgdal-dev


#finally, check if user has R installed, if not, install it
R_CHK=$(which R)
if [[ $R_CHK == "" ]]
then
	apt-get install -y r-base r-base-core
fi
RSTUDIO_CHK=$(which rstudio)
if [[ $RSTUDIO_CHK == "" ]]
then
	apt-get install dpkg-sig
	wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2021.09.1-372-amd64.deb
	gpg --keyserver keyserver.ubuntu.com --recv-keys 3F32EE77E331692F
	dpkg-sig --verify rstudio-2021.09.1-372-amd64.deb
	dpkg -i rstudio-2021.09.1-372-amd64.deb
	rm rstudio-2021.09.1-372-amd64.deb
fi
