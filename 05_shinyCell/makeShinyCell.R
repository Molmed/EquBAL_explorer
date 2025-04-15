library(Seurat)
library(ShinyCell)

#----------------------------------------------------------------
print("1. Preparing FULL dataset obj")
print("Reading in  FULL dataset obj")
fullobj <- readRDS("Equine_Asthma_fulldataset.rds")
print("Making initial config files")
scConf_fullobj = createConfig(fullobj)
print("Changing default plot colors")
scConf_fullobj = modColours(scConf_fullobj, meta.to.mod = "CellType_lvl1", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#a0c4ff", "#bdb2ff", "#ffc6ff"))
scConf_fullobj = modColours(scConf_fullobj, meta.to.mod = "Clusters res.0.3", new.colours= c("#fabed4", "#ffd8b1", "#fffac8", "#aaffc3", "#dcbeff", "#800000", "#9A6324", "#808000", "#469990", "#000075", "#e6194B", "#f58231", "#ffe119", "#bfef45", "#3cb44b"))
scConf_fullobj = modColours(scConf_fullobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(fullobj, scConf_fullobj, gex.assay = "SCT",  shiny.prefix = "fullobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "RESF1", default.gene2 = "LAMP3", default.multigene = c("FTH1", "APOE", "CD163", "CD3E", "TOP2A", "CD1E1", "MS4A2", "SNX10", "MARCO", "DHRS7", "CD2", "TPX2", "LAMP3", "FCER1A", "CSF3R", "STK17B", "ADAMDEC1", "CCR7", "GATA2", "TREM1", "RESF1"))
#----------------------------------------------------------------
print("2. Preparing TCELLS subset obj")
print("Reading in TCELLS subset obj")
tcellsobj <- readRDS("Equine_Asthma_subset_tcells.rds")
print("Making initial config files")
scConf_tcellsobj = createConfig(tcellsobj)
print("Changing default plot colors")
scConf_tcellsobj = modColours(scConf_tcellsobj, meta.to.mod = "CellType_lvl2", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#a0c4ff", "#bdb2ff", "#ffc6ff"))
scConf_tcellsobj = modColours(scConf_tcellsobj, meta.to.mod = "Clusters res.0.7", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#d0f4de", "#b5ead7", "#9bf6ff", "#a0c4ff", "#bdb2ff", "#ffc6ff"))
scConf_tcellsobj = modColours(scConf_tcellsobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff",  "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(tcellsobj, scConf_tcellsobj, gex.assay = "SCT",  shiny.prefix = "tcellsobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "CD2", default.gene2 = "KLRD1", default.multigene = c("CD3E", "FGL2", "KLRD1", "CD40LG", "CD4", "CD8A", "CTSW", "KLRK1", "GZMA", "GZMK", "CCL5",  "KLRB1", "TRDC", "TRAT1", "ETS1", "IGF2", "CD2", "RPS12", "RPS6", "MX1", "ISG15", "SAMD9L", "IL7R"))
#----------------------------------------------------------------
print("3. Preparing MACROPHAGES subset obj")
print("Reading in MACROPHAGES subset obj")
macrophagesobj <- readRDS("Equine_Asthma_subset_macrophages.rds")
print("Making initial config files")
scConf_macrophagesobj = createConfig(macrophagesobj)
print("Setting default metadata fields")
scConf_macrophagesobj = modDefault(scConf_macrophagesobj, default1 = "Clusters res.0.3", default2 = "BAL phenotype") 
print("Changing default plot colors")
scConf_macrophagesobj = modColours(scConf_macrophagesobj, meta.to.mod = "Clusters res.0.3", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#d0f4de", "#b5ead7", "#9bf6ff", "#a0c4ff", "#bdb2ff", "#c3b1e1", "#ffc6ff"))
scConf_macrophagesobj = modColours(scConf_macrophagesobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(macrophagesobj, scConf_macrophagesobj, gex.assay = "SCT",  shiny.prefix = "macrophagesobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "CD163", default.gene2 = "MARCO", default.multigene = c("MARCO", "CD163", "CCL8", "LGALS3", "MX1", "IFIT2", "GPNMB", "SOD2", "CCL24", "CTSW", "APOE", "TXNIP", "MRC1", "C4BPA", "FABP5", "FCER1G", "CD74", "PSAP", "ISG15", "C1QB", "ADD3", "ENSECAG00000003925", "EQMHCB2", "ANXA2", "ENSECAG00000032733", "TMSB4X", "ENSECAG00000027699", "ENSECAG00000027676", "C1QC", "ENSECAG00000003925", "ENSECAG00000029423", "CCDC88A", "IFIT3", "CTSD", "CTSB", "ENSECAG00000039383", "ENSECAG00000016787", "CXCL14"))
#----------------------------------------------------------------
print("4. Preparing DENDRITIC CELLS subset obj")
print("Reading in DENDRITIC CELLS subset obj")
dendriticobj <- readRDS("Equine_Asthma_subset_dendritic.rds")
print("Making initial config files")
scConf_dendriticobj = createConfig(dendriticobj)
print("Changing default plot colors")
scConf_dendriticobj = modColours(scConf_dendriticobj, meta.to.mod = "CellType_lvl2", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#a0c4ff"))
scConf_dendriticobj = modColours(scConf_dendriticobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff",  "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(dendriticobj, scConf_dendriticobj, gex.assay = "SCT",  shiny.prefix = "dendriticobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "CD1E1", default.gene2 = "CCR7", default.multigene = c("SPP1", "CD1E1", "LAMP1", "LAMP3", "CCR7", "CPVL", "CD14", "APOE", "FTH1", "CD1E2", "CD74", "DRA", "DRB", "CST3", "CLEC10A", "IDO1"))
#----------------------------------------------------------------
print("5. Preparing NEUTROPHILS subset obj")
print("Reading in NEUTROPHILS subset obj")
neutrophilsobj <- readRDS("Equine_Asthma_subset_neutrophils.rds")
print("Making initial config files")
scConf_neutrophilsobj = createConfig(neutrophilsobj)
print("Setting defailt metadata fields")
scConf_neutrophilsobj = modDefault(scConf_neutrophilsobj, default1 = "Clusters res.0.2", default2 = "BAL phenotype")
print("Changing default plot colors")
scConf_neutrophilsobj = modColours(scConf_neutrophilsobj, meta.to.mod = "Clusters res.0.2", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6"))
scConf_neutrophilsobj = modColours(scConf_neutrophilsobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(neutrophilsobj, scConf_neutrophilsobj, gex.assay = "SCT",  shiny.prefix = "neutrophilsobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "SNX10", default.gene2 = "CSF3R", default.multigene = c("CSF3R", "IFIT2", "IFIT3", "ISG15", "SNX10", "TREM1", "FTH1", "SRGN", "TREM1", "CTSZ", "VIM", "CD163"))
#----------------------------------------------------------------
print("6. Preparing MAST CELLS subset obj")
print("Reading in MAST CELLS subset obj") 
mastcellsobj <- readRDS("Equine_Asthma_subset_mastcells.rds")
print("Making initial config files")
scConf_mastcellsobj = createConfig(mastcellsobj)
print("Setting defailt metadata fields")
scConf_mastcellsobj = modDefault(scConf_mastcellsobj, default1 = "Clusters res.0.2", default2 = "BAL phenotype")
print("Changing default plot colors")
scConf_mastcellsobj = modColours(scConf_mastcellsobj, meta.to.mod = "Clusters res.0.2", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf"))
scConf_mastcellsobj = modColours(scConf_mastcellsobj, meta.to.mod = "BAL phenotype", new.colours= c("#ffadad", "#ffd6a5", "#fdffb6", "#caffbf", "#9bf6ff", "#bdb2ff", "#ffc6ff"))
print("Making shiny files, selecting default genes to show")
makeShinyFiles(mastcellsobj, scConf_mastcellsobj, gex.assay = "SCT",  shiny.prefix = "mastcellsobj", shiny.dir = "shinyAppMulti/", 
               default.gene1 = "MS4A2", default.gene2 = "GATA2", default.multigene = c("MS4A2", "FKBP5", "GATA2", "FCER1G", "KIT", "IGF2"))
#----------------------------------------------------------------
print("7. BUILDING shiny app")
citation = list(
  volume = "An Interactive Web Resource for Exploring Equine BAL Cell scRNA-seq Data", 
  page = "(2025) Sulyaeva J., Fegraeus K., Riihimäki M., Nordlund J., Raine A.", 
  year = "Funding: FORMAS")
makeShinyCodesMulti(
  shiny.title = "EquBAL explorer", 
  shiny.footnotes = citation, 
  shiny.prefix = c("fullobj", "tcellsobj", "macrophagesobj", "dendriticobj", "neutrophilsobj", "mastcellsobj"),
  shiny.headers = c("Full dataset", "T-cells subset", "Macrophages subset", "Dendritic cells subset", "Neutrophils subset", 
                    "Mast cells subset"),
  shiny.dir = "shinyAppMulti/")

