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


## Usage
This tool can be used simply by running the SpaceO.py script 

## Run script
`python SpaceO.py --Xenium_Image ./outs/morphology.ome.tif --Phenocycler_Directory_Home ./PhenocylerImage/CODEX.ome.tif --Zarr ./outs/cells_zarr  --Channels ./my_channels.csv`

**SpaceO** options:

* `--Xenium_Image` Path to morphology.ome.tif from Xenium ranger output (Ex: ./outs/morphology.ome.tif) 

* `--Phenocycler_Directory_Home` Path to a directory containing a single Phenocycler image.  (Ex: ./PhenocylerImage/CODEX.ome.tif)

* `--Zarr` Path to a directory containing the unzipped cells.zarr file from the Xenium ranger output. (Ex: ./outs/cells_zarr)

* `--Channels` Path to a .csv file with each row having the name of a single marker that was imaged (Ex: ./my_channels.csv)
