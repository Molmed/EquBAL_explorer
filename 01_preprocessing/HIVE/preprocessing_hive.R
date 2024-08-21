# EquineAsthmaCellAtlas code can ONLY run on Seuratv4 (preferably v.4.3)

library(Seurat)
library(DoubletFinder)

# Reading in count matrices
AM <- read.table("../../raw_data/HIVE/FU190-AM_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
BE <- read.table("../../raw_data/HIVE/FU190-BE_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
DE <- read.table("../../raw_data/HIVE/FU190-DE_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
FA <- read.table("../../raw_data/HIVE/FU190-FA_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
FL <- read.table("../../raw_data/HIVE/FU190-FL_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
GN <- read.table("../../raw_data/HIVE/FU190-GN_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
KA20 <- read.table("../../raw_data/HIVE/FU190-KA20_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
KA30 <- read.table("../../raw_data/HIVE/FU190-KA30_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
LE <- read.table("../../raw_data/HIVE/FU190-LE_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
NAT <- read.table("../../raw_data/HIVE/FU190-NAT_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
QU15 <- read.table("../../raw_data/HIVE/FU190-QU15_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
QU30 <- read.table("../../raw_data/HIVE/FU190-QU30_20230505120209_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)
ST <- read.table("../../raw_data/HIVE/FU190-ST_20230505120208_TCM.tsv.gz", header=1,row.names=1,check.names=FALSE , sep="\t", quote = "", stringsAsFactors=0)

# Creating seurat objects
Obj_AM<-CreateSeuratObject(counts=AM, project= "AM") 
Obj_BE <-CreateSeuratObject(counts=BE, project="BE") 
Obj_DE <-CreateSeuratObject(counts=DE, project = "DE") 
Obj_FA <-CreateSeuratObject(counts=FA, project="FA") 
Obj_FL <-CreateSeuratObject(counts=FL, project = "FL") 
Obj_GN <-CreateSeuratObject(counts=GN, project="GN") 
Obj_KA20 <-CreateSeuratObject(counts=KA20, project="KA20") 
Obj_KA30 <-CreateSeuratObject(counts=KA30, project ="KA30") 
Obj_LE <-CreateSeuratObject(counts=LE, project ="LE") 
Obj_NAT <-CreateSeuratObject(counts=NAT, project ="NAT")
Obj_QU15 <-CreateSeuratObject(counts=QU15, project ="QU15")
Obj_QU30 <-CreateSeuratObject(counts=QU30, project ="QU30")
Obj_ST <-CreateSeuratObject(counts=ST, project ="ST")


# Merging seurat objects
mergedobj <- merge(Obj_AM, c(Obj_BE, Obj_DE, Obj_FA, Obj_FL, Obj_GN, 
                             Obj_KA20, Obj_KA30, Obj_LE, Obj_NAT, Obj_QU15, 
                             Obj_QU30, Obj_ST), 
                 add.cell.ids = c("AM", "BE", "DE", "FA", "FL", 
                                  "GN", "KA20", "KA30", "LE", "NAT", 
                                  "QU15", "QU30", "ST"))

# Adding mitochondrial and ribosomal proportion data
total_counts_per_cell <- colSums(mergedobj@assays$RNA@counts)
mito_genes <- rownames(mergedobj)[grep("^MT-", rownames(mergedobj))]
mergedobj$percent_mito <- colSums(mergedobj@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
ribo_genes <- rownames(mergedobj)[grep("^RP[SL]", rownames(mergedobj))]
mergedobj$percent_ribo <- colSums(mergedobj@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell
obj <- mergedobj

# Filtering out cells with >8000 UMI, >10% mitochondrial genome, >15% ribosomal genome
obj_filtered <- subset(obj, nCount_RNA <= 8000) 
obj_filtered <- subset(obj_filtered, percent_mito <= 0.10) 
obj_filtered <- subset(obj_filtered, percent_ribo <= 0.15)

# Filtering out genes found in <0.1% of cells (for each horse separately, to preserve biological variability)
subsets <- SplitObject(obj_filtered, split.by = "orig.ident") 
subsets_filtered <- c()
for (i in 1:13) {
  horse_subset <- subsets[[i]]
  total_cells <- dim(horse_subset)[2]
  rarity_threshold <- as.integer(total_cells/1000)
  non_rare_genes <- rownames(horse_subset)[Matrix::rowSums(horse_subset)>rarity_threshold]
  subset_filtered <- subset(horse_subset, features = non_rare_genes)
  subsets_filtered <- c(subsets_filtered, subset_filtered)
}
obj_filtered <- merge(subsets_filtered[[1]], c(subsets_filtered[[2]], 
                                               subsets_filtered[[3]], 
                                               subsets_filtered[[4]],
                                               subsets_filtered[[5]],
                                               subsets_filtered[[6]],
                                               subsets_filtered[[7]], 
                                               subsets_filtered[[8]],
                                               subsets_filtered[[9]],
                                               subsets_filtered[[10]], 
                                               subsets_filtered[[11]],
                                               subsets_filtered[[12]],
                                               subsets_filtered[[13]]),
                      add.cell.ids = c("AM", "BE", "DE", "FA", "FL", 
                                       "GN", "KA20", "KA30", "LE", "NAT", 
                                       "QU15", "QU30", "ST"))

# Filtering out cells with <400 genes
obj_filtered <- subset(obj_filtered, nFeature_RNA >= 400)
obj <- obj_filtered

# Loading local files for two DoubletFinder functions due to incompatibility with Seuratv4
# (see DoubletFinder github #180)
source("../custom_code/paramSweep.R")
source("../custom_code/doubletFinder.R") 

# Running normalization and detecting doublets for each horse
subsets <- SplitObject(obj, split.by = "orig.ident") 
for (i in 1:13) {
  
  # Running SCTransform, PCA
  horse_obj <- SCTransform(subsets[[i]], vars.to.regress= "nFeature_RNA")  
  horse_obj <- RunPCA(horse_obj, nfeatures.print = 10)
  
  # Calculating significant PCs
  stdv <- horse_obj[["pca"]]@stdev
  sum.stdv <- sum(horse_obj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  
  # Running UMAP, FindNeighbors, FindClusters
  horse_obj <- RunUMAP(horse_obj, dims = 1:pcs)
  horse_obj <- FindNeighbors(object = horse_obj, dims = 1:pcs)              
  horse_obj <- FindClusters(object = horse_obj, resolution = 0.1)
  
  # Running PK identification (no ground-truth)
  sweep.list <- paramSweep(horse_obj, PCs = 1:min.pc, num.cores = 1, sct=T) 
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Calculating optimal pk
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  # Estimating homotypic doublet proportion
  annotations <- horse_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(horse_obj@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # Running DoubletFinder
  horse_obj <- doubletFinder(seu = horse_obj, PCs = 1:min.pc, pK = optimal.pk, nExp = nExp.poi.adj, sct=T)
  metadata <- horse_obj@meta.data
  colnames(metadata)[11] <- "doublet_finder2"
  horse_obj@meta.data <- metadata
  
  # Subsetting singlets
  horse_obj_singlets <- subset(horse_obj, doublet_finder2 == "Singlet")
  subsets[[i]] <- horse_obj_singlets
  remove(horse_obj_singlets)
}

# Merging singlet horse-sets back into one dataset
obj_singlets <- merge(x = subsets[[1]],
                      y = c(subsets[[2]], subsets[[3]], subsets[[4]],
                            subsets[[5]], subsets[[6]], subsets[[7]],
                            subsets[[8]], subsets[[9]], subsets[[10]],
                            subsets[[11]], subsets[[12]], subsets[[13]]), 
                      project = "HIVE Equine Asthma")

# Selecting variable features for the merged SCT object (seurat github #4145)
obj_features <- SelectIntegrationFeatures(object.list = subsets, nfeatures = 2000)
VariableFeatures(obj_singlets[["SCT"]]) <- obj_features
obj <- obj_singlets

# Cell-Cycle Scoring
obj <- SetIdent(obj, value = obj$orig.ident)
obj <- CellCycleScoring(obj, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
obj <- SetIdent(obj, value = obj$orig.ident)

# Adding metadata
HIVE <- obj
cells <- (dim(HIVE)[2])
HIVE$method <- rep("HIVE", times = cells)
HIVE$Phenotype1 <- rep("case", times = cells)

subsets <- SplitObject(HIVE, split.by = "orig.ident")
total_cells <- c()
for (i in 1:13) {
  total_cells <- c(total_cells, dim(subsets[[i]])[2])
}
subsets[[1]]$BAL_phenotype <- rep("mastocytic", times = total_cells[1]) #AM
subsets[[2]]$BAL_phenotype <- rep("mastocytic", times = total_cells[2]) #BE
subsets[[3]]$BAL_phenotype <- rep("mastocytic", times = total_cells[3]) #DE
subsets[[4]]$BAL_phenotype <- rep("mastocytic", times = total_cells[4]) #FA
subsets[[5]]$BAL_phenotype <- rep("mastocytic", times = total_cells[5]) #FL
subsets[[6]]$BAL_phenotype <- rep("neutrophilic", times = total_cells[6]) #GN
subsets[[7]]$BAL_phenotype <- rep("mastocytic_eosinophilic", times = total_cells[7]) #KA20
subsets[[8]]$BAL_phenotype <- rep("mastocytic_eosinophilic", times = total_cells[8]) #KA30
subsets[[9]]$BAL_phenotype <- rep("mastocytic_eosinophilic", times = total_cells[9]) #LE
subsets[[10]]$BAL_phenotype <- rep("normal_BAL", times = total_cells[10]) #NAT
subsets[[11]]$BAL_phenotype <- rep("mastocytic", times = total_cells[11]) #QU15
subsets[[12]]$BAL_phenotype <- rep("mastocytic", times = total_cells[12]) #QU30
subsets[[13]]$BAL_phenotype <- rep("mastocytic", times = total_cells[13]) #ST

HIVE <- merge(subsets[[1]], c(subsets[[2]], 
                              subsets[[3]], 
                              subsets[[4]],
                              subsets[[5]],
                              subsets[[6]],
                              subsets[[7]], 
                              subsets[[8]],
                              subsets[[9]],
                              subsets[[10]], 
                              subsets[[11]],
                              subsets[[12]],
                              subsets[[13]]))

saveRDS(HIVE, "preprocessed_hive.rds")
