# SpaceO

This repository contains the code for SpaceO. 
A start to end approach for the registration and quantification of paired Xenium and Spatial proteomics data.

## Installation
1) Ensure that your system has the appropriate Java Development Kit (JDK) containing the Java Runtime Environment (JRE). If necessary check with echo `$JAVA_HOME`
2) Download this repository and `cd` into the directory.
3) Download and install [Valis](https://valis.readthedocs.io/en/latest/index.html) prerequisites: 1)[Maven](https://maven.apache.org/index.html) and 2) [Libvips](https://www.libvips.org/). This tool was built and tested with [Maven v3.9.11](https://archive.apache.org/dist/maven/maven-3/3.9.11/binaries/apache-maven-3.9.11-bin.zip) and [Libvips v8.17.3](https://github.com/libvips/libvips/releases/download/v8.17.3/vips-8.17.3.tar.xz). Add Maven and Libvips to your path.
4) `conda env create -f environment.yml`.
5) Activate your conda environment
6) Install [Mcquant](https://github.com/quail768/quantification/tree/master) 


## Usage
This tool can be used simply by running the SpaceO.py script 

## Run script
`python SpaceO.py --Xenium_Image ./outs/morphology.ome.tif --Phenocycler_Directory_Home ./PhenocylerImage/CODEX.ome.tif --Zarr ./outs/cells_zarr  --Channels ./my_channels.csv   --Spatial_Proteomics_Magnification 20`

**SpaceO** options:

* `--Xenium_Image` Path to morphology.ome.tif from Xenium ranger output (Ex: ./outs/morphology.ome.tif) 

* `--Phenocycler_Directory_Home` Path to a directory containing a single Phenocycler image.  (Ex: ./PhenocylerImage/CODEX.ome.tif). The pipeline will crash if this directory contains any other files.

* `--Zarr` Path to a directory containing the unzipped cells.zarr file from the Xenium ranger output. (Ex: ./outs/cells_zarr). You can create this directory with `unzip  ./outs/cells.zarr.zip -d ./outs/cells_zarr`

* `--Channels` Path to a .csv file with each row having the name of a single marker that was imaged (Ex: ./my_channels.csv)

* `--Spatial_Proteomics_Magnification` The magnification used for imaging spatial proteomics data (Ex: 20)


## Outputs
SpaceO has a number of different outputs that help to understand your data. There will be 2 directories created within the directory where the script is run.

**1)  Results**

Here, both the reference `morphology.ome.tif` image and the registered `Phenocycler.ome.tif` image from the registration procedure with Valis are stored in Registerd_slides.
The other directory will be named after the '--Phenocycler_Directory_Home' (Ex: ./P134_BL/CODEX.ome.tif -> ./Results/P134_BL) and will contain the outputs from the Valis Micro_Rigid_Registration module. You can monitor how well the registration is proceeding with the help of this directory in ./Results/P134_BL/overlaps

<img width="549" height="365" alt="image" src="https://github.com/user-attachments/assets/098b7a8a-3d73-471e-bf6e-a81d954a153e" /> <img width="582" height="390" alt="image" src="https://github.com/user-attachments/assets/7cee7c7c-6ec3-4203-aa4b-0b995647cd76" />


**2) Quantification**
* `CountMatrixwithLabels.csv` - Count Matrix from McQuant with appended Xenium cell ids. Use the cell ids to append  this data to the Xenium data in R/Python.
* `masks.geojson` - Geojson file containing masks that can be imported into QuPath with the registered Phenocycler image to validate registration.
* `out_of_bounds_cells.txt` - Cells with masks in the Xenium data that have been rasterized outside the limits of the Phenocycler image. This should be not more than 20-30 if your registraton has worked well.
* `ImageMasks.tif` - Masks drawn from geojson object for perfroming quantification. Can be imported with geojson object into QuPath to validate the success of the drawing.
* `CODEX_ImageMasks.csv` - Output from Mcquant
* `label_lookup.csv`- Dictionary for masks and corresponding Xenium label. Not critical for performing analysis since the workflow assembles the count matrix and the corresponding cell id.




