# Equine Asthma Cell Atlas
Analysis scripts for the Equine Asthma Cell Atlas project 

## Citation
This is a public repository containing scripts described in the publication:

Raine et al. "-----" https://doi.org/-----
## Data
Sequencing data for this project is available at -----. All data has been analyzed on a high performance cluster (HPC) using R. The resulting analysis datasets are available at -----.

## Instructions
The analysis directories are numbered in the order the steps should be executed. Working directory for each script should be the same as the script's location, to align with relative paths.

To run the other scripts in this repository, you will need to do the following:

Install:

- R 3.6.0 and an integrated environment, e.g. RStudio
- R packages: Seurat v.4.3 (!), DoubletFinder, ShinyCell
- Additional R packages required by ShinyCell: data.table, Matrix, hdf5r, reticulate, ggplot2, gridExtra, glue, readr, RColorBrewer, R.utils, shiny, shinyhelper, DT, magrittr, ggdendro

Download the files from ------ and place them in the following directories:

- raw_data/Dropseq : all files from Dropseq folder (n=20)
- raw_data/HIVE : all files from HIVE folder (n=13)

Exceptions:

01_preprocessing/custom_code contains two functions from DoubletFinder package modified manually due to incompatibility with Seurat v.4. These do not require any actions. 
