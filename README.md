# Equine-Asthma-Cell-Atlas
Analysis scripts for the Equine Asthma Cell Atlas project 

## Citation
This is a public repository containing scripts described in the publication:

Raine et al. "add-publication-name-here" https://doi.org/add-doi-here

## Data
Sequencing data for this project is available at ---figshare-add-more-details. All data has been analyzed on a high performance cluster (HPC) using R. The resulting analysis datasets are available at --add-more-details.

## Instructions
The analysis steps are numbered in the order they should be executed. Working directory for each script should be the same as the script's location, to align with relative paths.

To run the other scripts in this repository, you will need to do the following.

Install:

- R 3.6.0 and an integrated environment, e.g. RStudio
- R packages: Seurat v.4.3, DoubletFinder, ShinyCell
- Additional R packages required by ShinyCell: data.table, Matrix, hdf5r, reticulate, ggplot2, gridExtra, glue, readr, RColorBrewer, R.utils, shiny, shinyhelper, DT, magrittr, ggdendro

Download the files from ---figshare?-- and place them in the following directories:

- raw_data/Dropseq : all files from Dropseq folder (n=20)
- raw_data/HIVE : all files from HIVE folder (n=13)
