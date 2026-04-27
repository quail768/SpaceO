# Xenium_SpatialProteomics

This repository contains the code for . 
A start to end approach for the registration and quantification of paired Xenium and Spatial proteomics data.

## Installation
1) Ensure that your system has the appropriate Java Development Kit (JDK) containing the Java Runtime Environment (JRE). If necessary check with 
2) Download this repository and `cd` into the directory.
3) Download and install [Valis](https://valis.readthedocs.io/en/latest/index.html) prerequisites: 1)[Maven](https://maven.apache.org/index.html) and 2) [Libvips](https://www.libvips.org/). This tool was built and tested with [Maven v3.9.11](https://archive.apache.org/dist/maven/maven-3/3.9.11/binaries/apache-maven-3.9.11-bin.zip) and [Libvips v8.17.3](https://archive.apache.org/dist/maven/maven-3/3.9.11/binaries/apache-maven-3.9.11-bin.zip). Add Maven and Libvips to your path.
4) `conda env create -f environment.yml`.
5) Activate your conda environment
6) Install [Mcquant](https://github.com/quail768/quantification/tree/master) 
