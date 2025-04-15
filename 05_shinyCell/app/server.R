library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
library(ggplot2) 
library(ggrepel) 
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
fullobjconf = readRDS("fullobjconf.rds")
fullobjdef  = readRDS("fullobjdef.rds")
fullobjgene = readRDS("fullobjgene.rds")
fullobjmeta = readRDS("fullobjmeta.rds")



tcellsobjconf = readRDS("tcellsobjconf.rds")
tcellsobjdef  = readRDS("tcellsobjdef.rds")
tcellsobjgene = readRDS("tcellsobjgene.rds")
tcellsobjmeta = readRDS("tcellsobjmeta.rds")



macrophagesobjconf = readRDS("macrophagesobjconf.rds")
macrophagesobjdef  = readRDS("macrophagesobjdef.rds")
macrophagesobjgene = readRDS("macrophagesobjgene.rds")
macrophagesobjmeta = readRDS("macrophagesobjmeta.rds")



dendriticobjconf = readRDS("dendriticobjconf.rds")
dendriticobjdef  = readRDS("dendriticobjdef.rds")
dendriticobjgene = readRDS("dendriticobjgene.rds")
dendriticobjmeta = readRDS("dendriticobjmeta.rds")



neutrophilsobjconf = readRDS("neutrophilsobjconf.rds")
neutrophilsobjdef  = readRDS("neutrophilsobjdef.rds")
neutrophilsobjgene = readRDS("neutrophilsobjgene.rds")
neutrophilsobjmeta = readRDS("neutrophilsobjmeta.rds")



mastcellsobjconf = readRDS("mastcellsobjconf.rds")
mastcellsobjdef  = readRDS("mastcellsobjdef.rds")
mastcellsobjgene = readRDS("mastcellsobjgene.rds")
mastcellsobjmeta = readRDS("mastcellsobjmeta.rds")



### Useful stuff 
# Colour palette 
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 
 
# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 
 
# Function to extract legend 
g_legend <- function(a.gplot){  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
  legend <- tmp$grobs[[leg]]  
  legend 
}  
 
# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
  oupTheme = theme( 
    text =             element_text(size = base_size, family = "Helvetica"), 
    panel.background = element_rect(fill = "white", colour = NA), 
    axis.line =   element_line(colour = "black"), 
    axis.ticks =  element_line(colour = "black", size = base_size / 20), 
    axis.title =  element_text(face = "bold"), 
    axis.text =   element_text(size = base_size), 
    axis.text.x = element_text(angle = Xang, hjust = XjusH), 
    legend.position = "bottom", 
    legend.key =      element_rect(colour = NA, fill = NA) 
  ) 
  if(!XYval){ 
    oupTheme = oupTheme + theme( 
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
  } 
  return(oupTheme) 
} 
 
### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "val", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Do factoring if required 
  if(!is.na(inpConf[UI == inp1]$fCL)){ 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$val) 
    ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
    ggData$val = factor(ggData$val, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
  } 
 
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  } else { 
    sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
    sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
    ggOut = ggOut + scale_color_manual("", values = ggCol) + 
      guides(color = guide_legend(override.aes = list(size = 5),  
                                  nrow = inpConf[UI == inp1]$fRow)) + 
      theme(legend.text = element_text(size = sListX[inpfsz])) 
    if(inplab){ 
      ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
      lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
      lListX = lList - (0.25 * floor(lListX/50)) 
      ggOut = ggOut + 
        geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                        color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                        size = lListX[inpfsz], seed = 42) 
    } 
  } 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("group", "sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Split inp1 if necessary 
  if(is.na(inpConf[UI == inp1]$fCL)){ 
    if(inpsplt == "Quartile"){nBk = 4} 
    if(inpsplt == "Decile"){nBk = 10} 
    ggData$group = cut(ggData$group, breaks = nBk) 
  } 
  
  # Actual data.table 
  ggData$express = FALSE 
  ggData[val2 > 0]$express = TRUE 
  ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
  ggData = ggData[, .(nCells = .N), by = "group"] 
  ggData = ggData1[ggData, on = "group"] 
  ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
  ggData[is.na(nExpress)]$nExpress = 0 
  ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
  ggData = ggData[order(group)] 
  colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
  return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val < 0]$val = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(val)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-val)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
   
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y, color = val)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
    scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
      guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
  oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
  oup = oup / (xy*xy) 
  return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "Y", "sub") 
  rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
  
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  bgCells = FALSE 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    bgCells = TRUE 
    ggData2 = ggData[!sub %in% inpsub2] 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Map colours 
  ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
  ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
  ggData$v0 = ggData$v1 + ggData$v2 
  ggData = gg[ggData, on = c("v1", "v2")] 
  if(inpord == "Max-1st"){ 
    ggData = ggData[order(v0)] 
  } else if(inpord == "Min-1st"){ 
    ggData = ggData[order(-v0)] 
  } else if(inpord == "Random"){ 
    ggData = ggData[sample(nrow(ggData))] 
  } 
  
  # Actual ggplot 
  ggOut = ggplot(ggData, aes(X, Y)) 
  if(bgCells){ 
    ggOut = ggOut + 
      geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + 
    geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
    xlab(inpdrX) + ylab(inpdrY) + 
    sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
    scale_color_gradientn(inp1, colours = cList[[1]]) + 
    guides(color = guide_colorbar(barwidth = 15)) 
  if(inpasp == "Square") { 
    ggOut = ggOut + coord_fixed(ratio = rat) 
  } else if(inpasp == "Fixed") { 
    ggOut = ggOut + coord_fixed() 
  } 
  return(ggOut) 
} 
 
scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
  # Generate coex color palette 
  cInp = strsplit(inpcol, "; ")[[1]] 
  if(cInp[1] == "Red (Gene1)"){ 
    c10 = c(255,0,0) 
  } else if(cInp[1] == "Orange (Gene1)"){ 
    c10 = c(255,140,0) 
  } else { 
    c10 = c(0,255,0) 
  } 
  if(cInp[2] == "Green (Gene2)"){ 
    c01 = c(0,255,0) 
  } else { 
    c01 = c(0,0,255) 
  } 
  c00 = c(217,217,217) ; c11 = c10 + c01 
  nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
  gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
  gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
  gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
  gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
  gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
  gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
  gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
  gg = gg[, c("v1", "v2", "cMix")] 
  
  # Actual ggplot 
  ggOut = ggplot(gg, aes(v1, v2)) + 
    geom_tile(fill = gg$cMix) + 
    xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
    scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
    sctheme(base_size = sList[inpfsz], XYval = TRUE) 
  return(ggOut) 
} 
 
scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
  colnames(ggData) = c("sub") 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
  ggData[val1 < 0]$val1 = 0 
  ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
  ggData[val2 < 0]$val2 = 0 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Actual data.table 
  ggData$express = "none" 
  ggData[val1 > 0]$express = inp1 
  ggData[val2 > 0]$express = inp2 
  ggData[val1 > 0 & val2 > 0]$express = "both" 
  ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
  ggData = ggData[, .(nCells = .N), by = "express"] 
  ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
  ggData = ggData[order(express)] 
  colnames(ggData)[1] = "expression > 0" 
  return(ggData) 
} 
 
# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                   with = FALSE] 
  colnames(ggData) = c("X", "sub") 
  
  # Load in either cell meta or gene expr
  if(inp2 %in% inpConf$UI){ 
    ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
  } else { 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val < 0]$val = 0 
    set.seed(42) 
    tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
    ggData$val = ggData$val + tmpNoise 
    h5file$close_all() 
  } 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$X) 
  ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
  ggData$X = factor(ggData$X, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "violin"){ 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
  } else { 
    ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
  } 
  if(inppts){ 
    ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
  } 
  ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) +
    theme(legend.position = "none")
  return(ggOut) 
} 
 
# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Prepare ggData 
  ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                       inpConf[UI == inpsub1]$ID),  
                   with = FALSE] 
  colnames(ggData) = c("X", "grp", "sub") 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
  ggData = ggData[, {tot = sum(nCells) 
                      .SD[,.(pctCells = 100 * sum(nCells) / tot, 
                             nCells = nCells), by = "grp"]}, by = "X"] 
  
  # Do factoring 
  ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
  names(ggCol) = levels(ggData$grp) 
  ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
  ggData$grp = factor(ggData$grp, levels = ggLvl) 
  ggCol = ggCol[ggLvl] 
  
  # Actual ggplot 
  if(inptyp == "Proportion"){ 
    ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
      geom_col() + ylab("Cell Proportion (%)") 
  } else { 
    ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
      geom_col() + ylab("Number of Cells") 
  } 
  if(inpflp){ 
    ggOut = ggOut + coord_flip() 
  } 
  ggOut = ggOut + xlab(inp1) + 
    sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
    scale_fill_manual("", values = ggCol) + 
    theme(legend.position = "right") 
  return(ggOut) 
} 
 
# Get gene list 
scGeneList <- function(inp, inpGene){ 
  geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                        present = TRUE) 
  geneList[!gene %in% names(inpGene)]$present = FALSE 
  return(geneList) 
} 
 
# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
  if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
  # Identify genes that are in our dataset 
  geneList = scGeneList(inp, inpGene) 
  geneList = geneList[present == TRUE] 
  shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
  shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
   
  # Prepare ggData 
  h5file <- H5File$new(inpH5, mode = "r") 
  h5data <- h5file[["grp"]][["data"]] 
  ggData = data.table() 
  for(iGene in geneList$gene){ 
    tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(tmp) = c("sampleID", "sub") 
    tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
    tmp$geneName = iGene 
    tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
    ggData = rbindlist(list(ggData, tmp)) 
  } 
  h5file$close_all() 
  if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
    ggData = ggData[sub %in% inpsub2] 
  } 
  shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
   
  # Aggregate 
  ggData$val = expm1(ggData$val) 
  ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                  by = c("geneName", "grpBy")] 
  ggData$val = log1p(ggData$val) 
   
  # Scale if required 
  colRange = range(ggData$val) 
  if(inpScl){ 
    ggData[, val:= scale(val), keyby = "geneName"] 
    colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
  } 
   
  # hclust row/col if necessary 
  ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
  tmp = ggMat$geneName 
  ggMat = as.matrix(ggMat[, -1]) 
  rownames(ggMat) = tmp 
  if(inpRow){ 
    hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
    ggRow = ggplot() + coord_flip() + 
      geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                         labels = unique(ggData$grpBy), expand = c(0, 0)) + 
      scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                         labels = hcRow$labels$label, expand = c(0, 0.5)) + 
      sctheme(base_size = sList[inpfsz]) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.y = element_blank(), 
            axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
    ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
  } else { 
    ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
  } 
  if(inpCol){ 
    hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
    ggCol = ggplot() + 
      geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
      scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                         labels = hcCol$labels$label, expand = c(0.05, 0)) + 
      scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                         labels = unique(ggData$geneName), expand=c(0,0)) + 
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      theme(axis.title = element_blank(), axis.line = element_blank(), 
            axis.ticks = element_blank(), axis.text.x = element_blank(), 
            axis.text.y = element_text(color = "white")) 
    ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
  } 
   
  # Actual plot according to plottype 
  if(inpPlt == "Bubbleplot"){ 
    # Bubbleplot 
    ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
      geom_point() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_size_continuous("proportion", range = c(0, 8), 
                            limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
      scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(color = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank(), legend.box = "vertical") 
  } else { 
    # Heatmap 
    ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
      geom_tile() +  
      sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
      scale_x_discrete(expand = c(0.05, 0)) +  
      scale_y_discrete(expand = c(0, 0.5)) + 
      scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
      guides(fill = guide_colorbar(barwidth = 15)) + 
      theme(axis.title = element_blank()) 
  } 
     
  # Final tidy 
  ggLeg = g_legend(ggOut) 
  ggOut = ggOut + theme(legend.position = "none") 
  if(!save){ 
    if(inpRow & inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                   layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                   layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                   layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                   layout_matrix = rbind(c(1),c(2)))  
    }  
  } else { 
    if(inpRow & inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                  layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
    } else if(inpRow){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                  layout_matrix = rbind(c(1,3),c(2,NA)))  
    } else if(inpCol){ggOut =  
      arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                  layout_matrix = rbind(c(3),c(1),c(2)))  
    } else {ggOut =  
      arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                  layout_matrix = rbind(c(1),c(2)))  
    }  
  } 
  return(ggOut) 
} 
 
 
 
 
 
### Start server code 
shinyServer(function(input, output, session) { 
  ### For all tags and Server-side selectize 
  observe_helpers() 
 optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "fullobja1inp2", choices = names(fullobjgene), server = TRUE, 
                       selected = fullobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "fullobja3inp1", choices = names(fullobjgene), server = TRUE, 
                       selected = fullobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "fullobja3inp2", choices = names(fullobjgene), server = TRUE, 
                       selected = fullobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "fullobjb2inp1", choices = names(fullobjgene), server = TRUE, 
                       selected = fullobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "fullobjb2inp2", choices = names(fullobjgene), server = TRUE, 
                       selected = fullobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "fullobjc1inp2", server = TRUE, 
                       choices = c(fullobjconf[is.na(fID)]$UI,names(fullobjgene)), 
                       selected = fullobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(fullobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$fullobja1sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobja1sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobja1sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobja1oup1 <- renderPlot({ 
    scDRcell(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp1,  
             input$fullobja1sub1, input$fullobja1sub2, 
             input$fullobja1siz, input$fullobja1col1, input$fullobja1ord1, 
             input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt, input$fullobja1lab1) 
  }) 
  output$fullobja1oup1.ui <- renderUI({ 
    plotOutput("fullobja1oup1", height = pList[input$fullobja1psz]) 
  }) 
  output$fullobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja1drX,"_",input$fullobja1drY,"_",  
                                   input$fullobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja1oup1.h, width = input$fullobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp1,   
                      input$fullobja1sub1, input$fullobja1sub2, 
                      input$fullobja1siz, input$fullobja1col1, input$fullobja1ord1,  
                      input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt, input$fullobja1lab1) ) 
  }) 
  output$fullobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja1drX,"_",input$fullobja1drY,"_",  
                                   input$fullobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja1oup1.h, width = input$fullobja1oup1.w, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp1,   
                      input$fullobja1sub1, input$fullobja1sub2, 
                      input$fullobja1siz, input$fullobja1col1, input$fullobja1ord1,  
                      input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt, input$fullobja1lab1) ) 
  }) 
  output$fullobja1.dt <- renderDataTable({ 
    ggData = scDRnum(fullobjconf, fullobjmeta, input$fullobja1inp1, input$fullobja1inp2, 
                     input$fullobja1sub1, input$fullobja1sub2, 
                     "fullobjgexpr.h5", fullobjgene, input$fullobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$fullobja1oup2 <- renderPlot({ 
    scDRgene(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp2,  
             input$fullobja1sub1, input$fullobja1sub2, 
             "fullobjgexpr.h5", fullobjgene, 
             input$fullobja1siz, input$fullobja1col2, input$fullobja1ord2, 
             input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt) 
  }) 
  output$fullobja1oup2.ui <- renderUI({ 
    plotOutput("fullobja1oup2", height = pList[input$fullobja1psz]) 
  }) 
  output$fullobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja1drX,"_",input$fullobja1drY,"_",  
                                   input$fullobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja1oup2.h, width = input$fullobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp2,  
                      input$fullobja1sub1, input$fullobja1sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja1siz, input$fullobja1col2, input$fullobja1ord2, 
                      input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt) ) 
  }) 
  output$fullobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja1drX,"_",input$fullobja1drY,"_",  
                                   input$fullobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja1oup2.h, width = input$fullobja1oup2.w, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja1drX, input$fullobja1drY, input$fullobja1inp2,  
                      input$fullobja1sub1, input$fullobja1sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja1siz, input$fullobja1col2, input$fullobja1ord2, 
                      input$fullobja1fsz, input$fullobja1asp, input$fullobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$fullobja2sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobja2sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobja2sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobja2oup1 <- renderPlot({ 
    scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp1,  
             input$fullobja2sub1, input$fullobja2sub2, 
             input$fullobja2siz, input$fullobja2col1, input$fullobja2ord1, 
             input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab1) 
  }) 
  output$fullobja2oup1.ui <- renderUI({ 
    plotOutput("fullobja2oup1", height = pList[input$fullobja2psz]) 
  }) 
  output$fullobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja2drX,"_",input$fullobja2drY,"_",  
                                   input$fullobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja2oup1.h, width = input$fullobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp1,   
                      input$fullobja2sub1, input$fullobja2sub2, 
                      input$fullobja2siz, input$fullobja2col1, input$fullobja2ord1,  
                      input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab1) ) 
  }) 
  output$fullobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja2drX,"_",input$fullobja2drY,"_",  
                                   input$fullobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja2oup1.h, width = input$fullobja2oup1.w, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp1,   
                      input$fullobja2sub1, input$fullobja2sub2, 
                      input$fullobja2siz, input$fullobja2col1, input$fullobja2ord1,  
                      input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab1) ) 
  }) 
   
  output$fullobja2oup2 <- renderPlot({ 
    scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp2,  
             input$fullobja2sub1, input$fullobja2sub2, 
             input$fullobja2siz, input$fullobja2col2, input$fullobja2ord2, 
             input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab2) 
  }) 
  output$fullobja2oup2.ui <- renderUI({ 
    plotOutput("fullobja2oup2", height = pList[input$fullobja2psz]) 
  }) 
  output$fullobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja2drX,"_",input$fullobja2drY,"_",  
                                   input$fullobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja2oup2.h, width = input$fullobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp2,   
                      input$fullobja2sub1, input$fullobja2sub2, 
                      input$fullobja2siz, input$fullobja2col2, input$fullobja2ord2,  
                      input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab2) ) 
  }) 
  output$fullobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja2drX,"_",input$fullobja2drY,"_",  
                                   input$fullobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja2oup2.h, width = input$fullobja2oup2.w, 
      plot = scDRcell(fullobjconf, fullobjmeta, input$fullobja2drX, input$fullobja2drY, input$fullobja2inp2,   
                      input$fullobja2sub1, input$fullobja2sub2, 
                      input$fullobja2siz, input$fullobja2col2, input$fullobja2ord2,  
                      input$fullobja2fsz, input$fullobja2asp, input$fullobja2txt, input$fullobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$fullobja3sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobja3sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobja3sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobja3oup1 <- renderPlot({ 
    scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp1,  
             input$fullobja3sub1, input$fullobja3sub2, 
             "fullobjgexpr.h5", fullobjgene, 
             input$fullobja3siz, input$fullobja3col1, input$fullobja3ord1, 
             input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) 
  }) 
  output$fullobja3oup1.ui <- renderUI({ 
    plotOutput("fullobja3oup1", height = pList[input$fullobja3psz]) 
  }) 
  output$fullobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja3drX,"_",input$fullobja3drY,"_",  
                                   input$fullobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja3oup1.h, width = input$fullobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp1,  
                      input$fullobja3sub1, input$fullobja3sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja3siz, input$fullobja3col1, input$fullobja3ord1, 
                      input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) ) 
  }) 
  output$fullobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja3drX,"_",input$fullobja3drY,"_",  
                                   input$fullobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja3oup1.h, width = input$fullobja3oup1.w, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp1,  
                      input$fullobja3sub1, input$fullobja3sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja3siz, input$fullobja3col1, input$fullobja3ord1, 
                      input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) ) 
  }) 
   
  output$fullobja3oup2 <- renderPlot({ 
    scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp2,  
             input$fullobja3sub1, input$fullobja3sub2, 
             "fullobjgexpr.h5", fullobjgene, 
             input$fullobja3siz, input$fullobja3col2, input$fullobja3ord2, 
             input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) 
  }) 
  output$fullobja3oup2.ui <- renderUI({ 
    plotOutput("fullobja3oup2", height = pList[input$fullobja3psz]) 
  }) 
  output$fullobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja3drX,"_",input$fullobja3drY,"_",  
                                   input$fullobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobja3oup2.h, width = input$fullobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp2,  
                      input$fullobja3sub1, input$fullobja3sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja3siz, input$fullobja3col2, input$fullobja3ord2, 
                      input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) ) 
  }) 
  output$fullobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobja3drX,"_",input$fullobja3drY,"_",  
                                   input$fullobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobja3oup2.h, width = input$fullobja3oup2.w, 
      plot = scDRgene(fullobjconf, fullobjmeta, input$fullobja3drX, input$fullobja3drY, input$fullobja3inp2,  
                      input$fullobja3sub1, input$fullobja3sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobja3siz, input$fullobja3col2, input$fullobja3ord2, 
                      input$fullobja3fsz, input$fullobja3asp, input$fullobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$fullobjb2sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobjb2sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobjb2sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobjb2oup1 <- renderPlot({ 
    scDRcoex(fullobjconf, fullobjmeta, input$fullobjb2drX, input$fullobjb2drY,   
             input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2sub1, input$fullobjb2sub2, 
             "fullobjgexpr.h5", fullobjgene, 
             input$fullobjb2siz, input$fullobjb2col1, input$fullobjb2ord1, 
             input$fullobjb2fsz, input$fullobjb2asp, input$fullobjb2txt) 
  }) 
  output$fullobjb2oup1.ui <- renderUI({ 
    plotOutput("fullobjb2oup1", height = pList2[input$fullobjb2psz]) 
  }) 
  output$fullobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjb2drX,"_",input$fullobjb2drY,"_",  
                                    input$fullobjb2inp1,"_",input$fullobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobjb2oup1.h, width = input$fullobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(fullobjconf, fullobjmeta, input$fullobjb2drX, input$fullobjb2drY,  
                      input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2sub1, input$fullobjb2sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobjb2siz, input$fullobjb2col1, input$fullobjb2ord1, 
                      input$fullobjb2fsz, input$fullobjb2asp, input$fullobjb2txt) ) 
  }) 
  output$fullobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjb2drX,"_",input$fullobjb2drY,"_",  
                                    input$fullobjb2inp1,"_",input$fullobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobjb2oup1.h, width = input$fullobjb2oup1.w, 
      plot = scDRcoex(fullobjconf, fullobjmeta, input$fullobjb2drX, input$fullobjb2drY,  
                      input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2sub1, input$fullobjb2sub2, 
                      "fullobjgexpr.h5", fullobjgene, 
                      input$fullobjb2siz, input$fullobjb2col1, input$fullobjb2ord1, 
                      input$fullobjb2fsz, input$fullobjb2asp, input$fullobjb2txt) ) 
  }) 
  output$fullobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2col1, input$fullobjb2fsz) 
  }) 
  output$fullobjb2oup2.ui <- renderUI({ 
    plotOutput("fullobjb2oup2", height = "300px") 
  }) 
  output$fullobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjb2drX,"_",input$fullobjb2drY,"_",  
                                    input$fullobjb2inp1,"_",input$fullobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2col1, input$fullobjb2fsz) ) 
  }) 
  output$fullobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjb2drX,"_",input$fullobjb2drY,"_",  
                                    input$fullobjb2inp1,"_",input$fullobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$fullobjb2inp1, input$fullobjb2inp2, input$fullobjb2col1, input$fullobjb2fsz) ) 
  }) 
  output$fullobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(fullobjconf, fullobjmeta, input$fullobjb2inp1, input$fullobjb2inp2, 
                         input$fullobjb2sub1, input$fullobjb2sub2, "fullobjgexpr.h5", fullobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$fullobjc1sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobjc1sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobjc1sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobjc1oup <- renderPlot({ 
    scVioBox(fullobjconf, fullobjmeta, input$fullobjc1inp1, input$fullobjc1inp2, 
             input$fullobjc1sub1, input$fullobjc1sub2, 
             "fullobjgexpr.h5", fullobjgene, input$fullobjc1typ, input$fullobjc1pts, 
             input$fullobjc1siz, input$fullobjc1fsz) 
  }) 
  output$fullobjc1oup.ui <- renderUI({ 
    plotOutput("fullobjc1oup", height = pList2[input$fullobjc1psz]) 
  }) 
  output$fullobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjc1typ,"_",input$fullobjc1inp1,"_",  
                                   input$fullobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobjc1oup.h, width = input$fullobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(fullobjconf, fullobjmeta, input$fullobjc1inp1, input$fullobjc1inp2, 
                      input$fullobjc1sub1, input$fullobjc1sub2, 
                      "fullobjgexpr.h5", fullobjgene, input$fullobjc1typ, input$fullobjc1pts, 
                      input$fullobjc1siz, input$fullobjc1fsz) ) 
  }) 
  output$fullobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjc1typ,"_",input$fullobjc1inp1,"_",  
                                   input$fullobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobjc1oup.h, width = input$fullobjc1oup.w, 
      plot = scVioBox(fullobjconf, fullobjmeta, input$fullobjc1inp1, input$fullobjc1inp2, 
                      input$fullobjc1sub1, input$fullobjc1sub2, 
                      "fullobjgexpr.h5", fullobjgene, input$fullobjc1typ, input$fullobjc1pts, 
                      input$fullobjc1siz, input$fullobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$fullobjc2sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobjc2sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobjc2sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$fullobjc2oup <- renderPlot({ 
  scProp(fullobjconf, fullobjmeta, input$fullobjc2inp1, input$fullobjc2inp2,  
         input$fullobjc2sub1, input$fullobjc2sub2, 
         input$fullobjc2typ, input$fullobjc2flp, input$fullobjc2fsz) 
}) 
output$fullobjc2oup.ui <- renderUI({ 
  plotOutput("fullobjc2oup", height = pList2[input$fullobjc2psz]) 
}) 
output$fullobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("fullobj",input$fullobjc2typ,"_",input$fullobjc2inp1,"_",  
                                 input$fullobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$fullobjc2oup.h, width = input$fullobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(fullobjconf, fullobjmeta, input$fullobjc2inp1, input$fullobjc2inp2,  
                  input$fullobjc2sub1, input$fullobjc2sub2, 
                  input$fullobjc2typ, input$fullobjc2flp, input$fullobjc2fsz) ) 
  }) 
output$fullobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("fullobj",input$fullobjc2typ,"_",input$fullobjc2inp1,"_",  
                                 input$fullobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$fullobjc2oup.h, width = input$fullobjc2oup.w, 
    plot = scProp(fullobjconf, fullobjmeta, input$fullobjc2inp1, input$fullobjc2inp2,  
                  input$fullobjc2sub1, input$fullobjc2sub2, 
                  input$fullobjc2typ, input$fullobjc2flp, input$fullobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$fullobjd1sub1.ui <- renderUI({ 
    sub = strsplit(fullobjconf[UI == input$fullobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("fullobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$fullobjd1sub1non, { 
    sub = strsplit(fullobjconf[UI == input$fullobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$fullobjd1sub1all, { 
    sub = strsplit(fullobjconf[UI == input$fullobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "fullobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$fullobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$fullobjd1inp, fullobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$fullobjd1oup <- renderPlot({ 
    scBubbHeat(fullobjconf, fullobjmeta, input$fullobjd1inp, input$fullobjd1grp, input$fullobjd1plt, 
               input$fullobjd1sub1, input$fullobjd1sub2, "fullobjgexpr.h5", fullobjgene, 
               input$fullobjd1scl, input$fullobjd1row, input$fullobjd1col, 
               input$fullobjd1cols, input$fullobjd1fsz) 
  }) 
  output$fullobjd1oup.ui <- renderUI({ 
    plotOutput("fullobjd1oup", height = pList3[input$fullobjd1psz]) 
  }) 
  output$fullobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjd1plt,"_",input$fullobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$fullobjd1oup.h, width = input$fullobjd1oup.w, 
      plot = scBubbHeat(fullobjconf, fullobjmeta, input$fullobjd1inp, input$fullobjd1grp, input$fullobjd1plt, 
                        input$fullobjd1sub1, input$fullobjd1sub2, "fullobjgexpr.h5", fullobjgene, 
                        input$fullobjd1scl, input$fullobjd1row, input$fullobjd1col, 
                        input$fullobjd1cols, input$fullobjd1fsz, save = TRUE) ) 
  }) 
  output$fullobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("fullobj",input$fullobjd1plt,"_",input$fullobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$fullobjd1oup.h, width = input$fullobjd1oup.w, 
      plot = scBubbHeat(fullobjconf, fullobjmeta, input$fullobjd1inp, input$fullobjd1grp, input$fullobjd1plt, 
                        input$fullobjd1sub1, input$fullobjd1sub2, "fullobjgexpr.h5", fullobjgene, 
                        input$fullobjd1scl, input$fullobjd1row, input$fullobjd1col, 
                        input$fullobjd1cols, input$fullobjd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "tcellsobja1inp2", choices = names(tcellsobjgene), server = TRUE, 
                       selected = tcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "tcellsobja3inp1", choices = names(tcellsobjgene), server = TRUE, 
                       selected = tcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "tcellsobja3inp2", choices = names(tcellsobjgene), server = TRUE, 
                       selected = tcellsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "tcellsobjb2inp1", choices = names(tcellsobjgene), server = TRUE, 
                       selected = tcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "tcellsobjb2inp2", choices = names(tcellsobjgene), server = TRUE, 
                       selected = tcellsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "tcellsobjc1inp2", server = TRUE, 
                       choices = c(tcellsobjconf[is.na(fID)]$UI,names(tcellsobjgene)), 
                       selected = tcellsobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(tcellsobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$tcellsobja1sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobja1sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobja1sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobja1oup1 <- renderPlot({ 
    scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp1,  
             input$tcellsobja1sub1, input$tcellsobja1sub2, 
             input$tcellsobja1siz, input$tcellsobja1col1, input$tcellsobja1ord1, 
             input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt, input$tcellsobja1lab1) 
  }) 
  output$tcellsobja1oup1.ui <- renderUI({ 
    plotOutput("tcellsobja1oup1", height = pList[input$tcellsobja1psz]) 
  }) 
  output$tcellsobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja1drX,"_",input$tcellsobja1drY,"_",  
                                   input$tcellsobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja1oup1.h, width = input$tcellsobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp1,   
                      input$tcellsobja1sub1, input$tcellsobja1sub2, 
                      input$tcellsobja1siz, input$tcellsobja1col1, input$tcellsobja1ord1,  
                      input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt, input$tcellsobja1lab1) ) 
  }) 
  output$tcellsobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja1drX,"_",input$tcellsobja1drY,"_",  
                                   input$tcellsobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja1oup1.h, width = input$tcellsobja1oup1.w, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp1,   
                      input$tcellsobja1sub1, input$tcellsobja1sub2, 
                      input$tcellsobja1siz, input$tcellsobja1col1, input$tcellsobja1ord1,  
                      input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt, input$tcellsobja1lab1) ) 
  }) 
  output$tcellsobja1.dt <- renderDataTable({ 
    ggData = scDRnum(tcellsobjconf, tcellsobjmeta, input$tcellsobja1inp1, input$tcellsobja1inp2, 
                     input$tcellsobja1sub1, input$tcellsobja1sub2, 
                     "tcellsobjgexpr.h5", tcellsobjgene, input$tcellsobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$tcellsobja1oup2 <- renderPlot({ 
    scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp2,  
             input$tcellsobja1sub1, input$tcellsobja1sub2, 
             "tcellsobjgexpr.h5", tcellsobjgene, 
             input$tcellsobja1siz, input$tcellsobja1col2, input$tcellsobja1ord2, 
             input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt) 
  }) 
  output$tcellsobja1oup2.ui <- renderUI({ 
    plotOutput("tcellsobja1oup2", height = pList[input$tcellsobja1psz]) 
  }) 
  output$tcellsobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja1drX,"_",input$tcellsobja1drY,"_",  
                                   input$tcellsobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja1oup2.h, width = input$tcellsobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp2,  
                      input$tcellsobja1sub1, input$tcellsobja1sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja1siz, input$tcellsobja1col2, input$tcellsobja1ord2, 
                      input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt) ) 
  }) 
  output$tcellsobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja1drX,"_",input$tcellsobja1drY,"_",  
                                   input$tcellsobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja1oup2.h, width = input$tcellsobja1oup2.w, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja1drX, input$tcellsobja1drY, input$tcellsobja1inp2,  
                      input$tcellsobja1sub1, input$tcellsobja1sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja1siz, input$tcellsobja1col2, input$tcellsobja1ord2, 
                      input$tcellsobja1fsz, input$tcellsobja1asp, input$tcellsobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$tcellsobja2sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobja2sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobja2sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobja2oup1 <- renderPlot({ 
    scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp1,  
             input$tcellsobja2sub1, input$tcellsobja2sub2, 
             input$tcellsobja2siz, input$tcellsobja2col1, input$tcellsobja2ord1, 
             input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab1) 
  }) 
  output$tcellsobja2oup1.ui <- renderUI({ 
    plotOutput("tcellsobja2oup1", height = pList[input$tcellsobja2psz]) 
  }) 
  output$tcellsobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja2drX,"_",input$tcellsobja2drY,"_",  
                                   input$tcellsobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja2oup1.h, width = input$tcellsobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp1,   
                      input$tcellsobja2sub1, input$tcellsobja2sub2, 
                      input$tcellsobja2siz, input$tcellsobja2col1, input$tcellsobja2ord1,  
                      input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab1) ) 
  }) 
  output$tcellsobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja2drX,"_",input$tcellsobja2drY,"_",  
                                   input$tcellsobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja2oup1.h, width = input$tcellsobja2oup1.w, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp1,   
                      input$tcellsobja2sub1, input$tcellsobja2sub2, 
                      input$tcellsobja2siz, input$tcellsobja2col1, input$tcellsobja2ord1,  
                      input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab1) ) 
  }) 
   
  output$tcellsobja2oup2 <- renderPlot({ 
    scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp2,  
             input$tcellsobja2sub1, input$tcellsobja2sub2, 
             input$tcellsobja2siz, input$tcellsobja2col2, input$tcellsobja2ord2, 
             input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab2) 
  }) 
  output$tcellsobja2oup2.ui <- renderUI({ 
    plotOutput("tcellsobja2oup2", height = pList[input$tcellsobja2psz]) 
  }) 
  output$tcellsobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja2drX,"_",input$tcellsobja2drY,"_",  
                                   input$tcellsobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja2oup2.h, width = input$tcellsobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp2,   
                      input$tcellsobja2sub1, input$tcellsobja2sub2, 
                      input$tcellsobja2siz, input$tcellsobja2col2, input$tcellsobja2ord2,  
                      input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab2) ) 
  }) 
  output$tcellsobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja2drX,"_",input$tcellsobja2drY,"_",  
                                   input$tcellsobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja2oup2.h, width = input$tcellsobja2oup2.w, 
      plot = scDRcell(tcellsobjconf, tcellsobjmeta, input$tcellsobja2drX, input$tcellsobja2drY, input$tcellsobja2inp2,   
                      input$tcellsobja2sub1, input$tcellsobja2sub2, 
                      input$tcellsobja2siz, input$tcellsobja2col2, input$tcellsobja2ord2,  
                      input$tcellsobja2fsz, input$tcellsobja2asp, input$tcellsobja2txt, input$tcellsobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$tcellsobja3sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobja3sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobja3sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobja3oup1 <- renderPlot({ 
    scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp1,  
             input$tcellsobja3sub1, input$tcellsobja3sub2, 
             "tcellsobjgexpr.h5", tcellsobjgene, 
             input$tcellsobja3siz, input$tcellsobja3col1, input$tcellsobja3ord1, 
             input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) 
  }) 
  output$tcellsobja3oup1.ui <- renderUI({ 
    plotOutput("tcellsobja3oup1", height = pList[input$tcellsobja3psz]) 
  }) 
  output$tcellsobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja3drX,"_",input$tcellsobja3drY,"_",  
                                   input$tcellsobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja3oup1.h, width = input$tcellsobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp1,  
                      input$tcellsobja3sub1, input$tcellsobja3sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja3siz, input$tcellsobja3col1, input$tcellsobja3ord1, 
                      input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) ) 
  }) 
  output$tcellsobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja3drX,"_",input$tcellsobja3drY,"_",  
                                   input$tcellsobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja3oup1.h, width = input$tcellsobja3oup1.w, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp1,  
                      input$tcellsobja3sub1, input$tcellsobja3sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja3siz, input$tcellsobja3col1, input$tcellsobja3ord1, 
                      input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) ) 
  }) 
   
  output$tcellsobja3oup2 <- renderPlot({ 
    scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp2,  
             input$tcellsobja3sub1, input$tcellsobja3sub2, 
             "tcellsobjgexpr.h5", tcellsobjgene, 
             input$tcellsobja3siz, input$tcellsobja3col2, input$tcellsobja3ord2, 
             input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) 
  }) 
  output$tcellsobja3oup2.ui <- renderUI({ 
    plotOutput("tcellsobja3oup2", height = pList[input$tcellsobja3psz]) 
  }) 
  output$tcellsobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja3drX,"_",input$tcellsobja3drY,"_",  
                                   input$tcellsobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobja3oup2.h, width = input$tcellsobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp2,  
                      input$tcellsobja3sub1, input$tcellsobja3sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja3siz, input$tcellsobja3col2, input$tcellsobja3ord2, 
                      input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) ) 
  }) 
  output$tcellsobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobja3drX,"_",input$tcellsobja3drY,"_",  
                                   input$tcellsobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobja3oup2.h, width = input$tcellsobja3oup2.w, 
      plot = scDRgene(tcellsobjconf, tcellsobjmeta, input$tcellsobja3drX, input$tcellsobja3drY, input$tcellsobja3inp2,  
                      input$tcellsobja3sub1, input$tcellsobja3sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobja3siz, input$tcellsobja3col2, input$tcellsobja3ord2, 
                      input$tcellsobja3fsz, input$tcellsobja3asp, input$tcellsobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$tcellsobjb2sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobjb2sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobjb2sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobjb2oup1 <- renderPlot({ 
    scDRcoex(tcellsobjconf, tcellsobjmeta, input$tcellsobjb2drX, input$tcellsobjb2drY,   
             input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2sub1, input$tcellsobjb2sub2, 
             "tcellsobjgexpr.h5", tcellsobjgene, 
             input$tcellsobjb2siz, input$tcellsobjb2col1, input$tcellsobjb2ord1, 
             input$tcellsobjb2fsz, input$tcellsobjb2asp, input$tcellsobjb2txt) 
  }) 
  output$tcellsobjb2oup1.ui <- renderUI({ 
    plotOutput("tcellsobjb2oup1", height = pList2[input$tcellsobjb2psz]) 
  }) 
  output$tcellsobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjb2drX,"_",input$tcellsobjb2drY,"_",  
                                    input$tcellsobjb2inp1,"_",input$tcellsobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobjb2oup1.h, width = input$tcellsobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(tcellsobjconf, tcellsobjmeta, input$tcellsobjb2drX, input$tcellsobjb2drY,  
                      input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2sub1, input$tcellsobjb2sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobjb2siz, input$tcellsobjb2col1, input$tcellsobjb2ord1, 
                      input$tcellsobjb2fsz, input$tcellsobjb2asp, input$tcellsobjb2txt) ) 
  }) 
  output$tcellsobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjb2drX,"_",input$tcellsobjb2drY,"_",  
                                    input$tcellsobjb2inp1,"_",input$tcellsobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobjb2oup1.h, width = input$tcellsobjb2oup1.w, 
      plot = scDRcoex(tcellsobjconf, tcellsobjmeta, input$tcellsobjb2drX, input$tcellsobjb2drY,  
                      input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2sub1, input$tcellsobjb2sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, 
                      input$tcellsobjb2siz, input$tcellsobjb2col1, input$tcellsobjb2ord1, 
                      input$tcellsobjb2fsz, input$tcellsobjb2asp, input$tcellsobjb2txt) ) 
  }) 
  output$tcellsobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2col1, input$tcellsobjb2fsz) 
  }) 
  output$tcellsobjb2oup2.ui <- renderUI({ 
    plotOutput("tcellsobjb2oup2", height = "300px") 
  }) 
  output$tcellsobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjb2drX,"_",input$tcellsobjb2drY,"_",  
                                    input$tcellsobjb2inp1,"_",input$tcellsobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2col1, input$tcellsobjb2fsz) ) 
  }) 
  output$tcellsobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjb2drX,"_",input$tcellsobjb2drY,"_",  
                                    input$tcellsobjb2inp1,"_",input$tcellsobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$tcellsobjb2inp1, input$tcellsobjb2inp2, input$tcellsobjb2col1, input$tcellsobjb2fsz) ) 
  }) 
  output$tcellsobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(tcellsobjconf, tcellsobjmeta, input$tcellsobjb2inp1, input$tcellsobjb2inp2, 
                         input$tcellsobjb2sub1, input$tcellsobjb2sub2, "tcellsobjgexpr.h5", tcellsobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$tcellsobjc1sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobjc1sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobjc1sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobjc1oup <- renderPlot({ 
    scVioBox(tcellsobjconf, tcellsobjmeta, input$tcellsobjc1inp1, input$tcellsobjc1inp2, 
             input$tcellsobjc1sub1, input$tcellsobjc1sub2, 
             "tcellsobjgexpr.h5", tcellsobjgene, input$tcellsobjc1typ, input$tcellsobjc1pts, 
             input$tcellsobjc1siz, input$tcellsobjc1fsz) 
  }) 
  output$tcellsobjc1oup.ui <- renderUI({ 
    plotOutput("tcellsobjc1oup", height = pList2[input$tcellsobjc1psz]) 
  }) 
  output$tcellsobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjc1typ,"_",input$tcellsobjc1inp1,"_",  
                                   input$tcellsobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobjc1oup.h, width = input$tcellsobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(tcellsobjconf, tcellsobjmeta, input$tcellsobjc1inp1, input$tcellsobjc1inp2, 
                      input$tcellsobjc1sub1, input$tcellsobjc1sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, input$tcellsobjc1typ, input$tcellsobjc1pts, 
                      input$tcellsobjc1siz, input$tcellsobjc1fsz) ) 
  }) 
  output$tcellsobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjc1typ,"_",input$tcellsobjc1inp1,"_",  
                                   input$tcellsobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobjc1oup.h, width = input$tcellsobjc1oup.w, 
      plot = scVioBox(tcellsobjconf, tcellsobjmeta, input$tcellsobjc1inp1, input$tcellsobjc1inp2, 
                      input$tcellsobjc1sub1, input$tcellsobjc1sub2, 
                      "tcellsobjgexpr.h5", tcellsobjgene, input$tcellsobjc1typ, input$tcellsobjc1pts, 
                      input$tcellsobjc1siz, input$tcellsobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$tcellsobjc2sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobjc2sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobjc2sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$tcellsobjc2oup <- renderPlot({ 
  scProp(tcellsobjconf, tcellsobjmeta, input$tcellsobjc2inp1, input$tcellsobjc2inp2,  
         input$tcellsobjc2sub1, input$tcellsobjc2sub2, 
         input$tcellsobjc2typ, input$tcellsobjc2flp, input$tcellsobjc2fsz) 
}) 
output$tcellsobjc2oup.ui <- renderUI({ 
  plotOutput("tcellsobjc2oup", height = pList2[input$tcellsobjc2psz]) 
}) 
output$tcellsobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("tcellsobj",input$tcellsobjc2typ,"_",input$tcellsobjc2inp1,"_",  
                                 input$tcellsobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$tcellsobjc2oup.h, width = input$tcellsobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(tcellsobjconf, tcellsobjmeta, input$tcellsobjc2inp1, input$tcellsobjc2inp2,  
                  input$tcellsobjc2sub1, input$tcellsobjc2sub2, 
                  input$tcellsobjc2typ, input$tcellsobjc2flp, input$tcellsobjc2fsz) ) 
  }) 
output$tcellsobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("tcellsobj",input$tcellsobjc2typ,"_",input$tcellsobjc2inp1,"_",  
                                 input$tcellsobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$tcellsobjc2oup.h, width = input$tcellsobjc2oup.w, 
    plot = scProp(tcellsobjconf, tcellsobjmeta, input$tcellsobjc2inp1, input$tcellsobjc2inp2,  
                  input$tcellsobjc2sub1, input$tcellsobjc2sub2, 
                  input$tcellsobjc2typ, input$tcellsobjc2flp, input$tcellsobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$tcellsobjd1sub1.ui <- renderUI({ 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("tcellsobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$tcellsobjd1sub1non, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$tcellsobjd1sub1all, { 
    sub = strsplit(tcellsobjconf[UI == input$tcellsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "tcellsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$tcellsobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$tcellsobjd1inp, tcellsobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$tcellsobjd1oup <- renderPlot({ 
    scBubbHeat(tcellsobjconf, tcellsobjmeta, input$tcellsobjd1inp, input$tcellsobjd1grp, input$tcellsobjd1plt, 
               input$tcellsobjd1sub1, input$tcellsobjd1sub2, "tcellsobjgexpr.h5", tcellsobjgene, 
               input$tcellsobjd1scl, input$tcellsobjd1row, input$tcellsobjd1col, 
               input$tcellsobjd1cols, input$tcellsobjd1fsz) 
  }) 
  output$tcellsobjd1oup.ui <- renderUI({ 
    plotOutput("tcellsobjd1oup", height = pList3[input$tcellsobjd1psz]) 
  }) 
  output$tcellsobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjd1plt,"_",input$tcellsobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$tcellsobjd1oup.h, width = input$tcellsobjd1oup.w, 
      plot = scBubbHeat(tcellsobjconf, tcellsobjmeta, input$tcellsobjd1inp, input$tcellsobjd1grp, input$tcellsobjd1plt, 
                        input$tcellsobjd1sub1, input$tcellsobjd1sub2, "tcellsobjgexpr.h5", tcellsobjgene, 
                        input$tcellsobjd1scl, input$tcellsobjd1row, input$tcellsobjd1col, 
                        input$tcellsobjd1cols, input$tcellsobjd1fsz, save = TRUE) ) 
  }) 
  output$tcellsobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("tcellsobj",input$tcellsobjd1plt,"_",input$tcellsobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$tcellsobjd1oup.h, width = input$tcellsobjd1oup.w, 
      plot = scBubbHeat(tcellsobjconf, tcellsobjmeta, input$tcellsobjd1inp, input$tcellsobjd1grp, input$tcellsobjd1plt, 
                        input$tcellsobjd1sub1, input$tcellsobjd1sub2, "tcellsobjgexpr.h5", tcellsobjgene, 
                        input$tcellsobjd1scl, input$tcellsobjd1row, input$tcellsobjd1col, 
                        input$tcellsobjd1cols, input$tcellsobjd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "macrophagesobja1inp2", choices = names(macrophagesobjgene), server = TRUE, 
                       selected = macrophagesobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "macrophagesobja3inp1", choices = names(macrophagesobjgene), server = TRUE, 
                       selected = macrophagesobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "macrophagesobja3inp2", choices = names(macrophagesobjgene), server = TRUE, 
                       selected = macrophagesobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "macrophagesobjb2inp1", choices = names(macrophagesobjgene), server = TRUE, 
                       selected = macrophagesobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "macrophagesobjb2inp2", choices = names(macrophagesobjgene), server = TRUE, 
                       selected = macrophagesobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "macrophagesobjc1inp2", server = TRUE, 
                       choices = c(macrophagesobjconf[is.na(fID)]$UI,names(macrophagesobjgene)), 
                       selected = macrophagesobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(macrophagesobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$macrophagesobja1sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobja1sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobja1sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobja1oup1 <- renderPlot({ 
    scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp1,  
             input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
             input$macrophagesobja1siz, input$macrophagesobja1col1, input$macrophagesobja1ord1, 
             input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt, input$macrophagesobja1lab1) 
  }) 
  output$macrophagesobja1oup1.ui <- renderUI({ 
    plotOutput("macrophagesobja1oup1", height = pList[input$macrophagesobja1psz]) 
  }) 
  output$macrophagesobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja1drX,"_",input$macrophagesobja1drY,"_",  
                                   input$macrophagesobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja1oup1.h, width = input$macrophagesobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp1,   
                      input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
                      input$macrophagesobja1siz, input$macrophagesobja1col1, input$macrophagesobja1ord1,  
                      input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt, input$macrophagesobja1lab1) ) 
  }) 
  output$macrophagesobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja1drX,"_",input$macrophagesobja1drY,"_",  
                                   input$macrophagesobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja1oup1.h, width = input$macrophagesobja1oup1.w, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp1,   
                      input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
                      input$macrophagesobja1siz, input$macrophagesobja1col1, input$macrophagesobja1ord1,  
                      input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt, input$macrophagesobja1lab1) ) 
  }) 
  output$macrophagesobja1.dt <- renderDataTable({ 
    ggData = scDRnum(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1inp1, input$macrophagesobja1inp2, 
                     input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
                     "macrophagesobjgexpr.h5", macrophagesobjgene, input$macrophagesobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$macrophagesobja1oup2 <- renderPlot({ 
    scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp2,  
             input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
             "macrophagesobjgexpr.h5", macrophagesobjgene, 
             input$macrophagesobja1siz, input$macrophagesobja1col2, input$macrophagesobja1ord2, 
             input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt) 
  }) 
  output$macrophagesobja1oup2.ui <- renderUI({ 
    plotOutput("macrophagesobja1oup2", height = pList[input$macrophagesobja1psz]) 
  }) 
  output$macrophagesobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja1drX,"_",input$macrophagesobja1drY,"_",  
                                   input$macrophagesobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja1oup2.h, width = input$macrophagesobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp2,  
                      input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja1siz, input$macrophagesobja1col2, input$macrophagesobja1ord2, 
                      input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt) ) 
  }) 
  output$macrophagesobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja1drX,"_",input$macrophagesobja1drY,"_",  
                                   input$macrophagesobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja1oup2.h, width = input$macrophagesobja1oup2.w, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja1drX, input$macrophagesobja1drY, input$macrophagesobja1inp2,  
                      input$macrophagesobja1sub1, input$macrophagesobja1sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja1siz, input$macrophagesobja1col2, input$macrophagesobja1ord2, 
                      input$macrophagesobja1fsz, input$macrophagesobja1asp, input$macrophagesobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$macrophagesobja2sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobja2sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobja2sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobja2oup1 <- renderPlot({ 
    scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp1,  
             input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
             input$macrophagesobja2siz, input$macrophagesobja2col1, input$macrophagesobja2ord1, 
             input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab1) 
  }) 
  output$macrophagesobja2oup1.ui <- renderUI({ 
    plotOutput("macrophagesobja2oup1", height = pList[input$macrophagesobja2psz]) 
  }) 
  output$macrophagesobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja2drX,"_",input$macrophagesobja2drY,"_",  
                                   input$macrophagesobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja2oup1.h, width = input$macrophagesobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp1,   
                      input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
                      input$macrophagesobja2siz, input$macrophagesobja2col1, input$macrophagesobja2ord1,  
                      input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab1) ) 
  }) 
  output$macrophagesobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja2drX,"_",input$macrophagesobja2drY,"_",  
                                   input$macrophagesobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja2oup1.h, width = input$macrophagesobja2oup1.w, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp1,   
                      input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
                      input$macrophagesobja2siz, input$macrophagesobja2col1, input$macrophagesobja2ord1,  
                      input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab1) ) 
  }) 
   
  output$macrophagesobja2oup2 <- renderPlot({ 
    scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp2,  
             input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
             input$macrophagesobja2siz, input$macrophagesobja2col2, input$macrophagesobja2ord2, 
             input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab2) 
  }) 
  output$macrophagesobja2oup2.ui <- renderUI({ 
    plotOutput("macrophagesobja2oup2", height = pList[input$macrophagesobja2psz]) 
  }) 
  output$macrophagesobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja2drX,"_",input$macrophagesobja2drY,"_",  
                                   input$macrophagesobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja2oup2.h, width = input$macrophagesobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp2,   
                      input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
                      input$macrophagesobja2siz, input$macrophagesobja2col2, input$macrophagesobja2ord2,  
                      input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab2) ) 
  }) 
  output$macrophagesobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja2drX,"_",input$macrophagesobja2drY,"_",  
                                   input$macrophagesobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja2oup2.h, width = input$macrophagesobja2oup2.w, 
      plot = scDRcell(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja2drX, input$macrophagesobja2drY, input$macrophagesobja2inp2,   
                      input$macrophagesobja2sub1, input$macrophagesobja2sub2, 
                      input$macrophagesobja2siz, input$macrophagesobja2col2, input$macrophagesobja2ord2,  
                      input$macrophagesobja2fsz, input$macrophagesobja2asp, input$macrophagesobja2txt, input$macrophagesobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$macrophagesobja3sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobja3sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobja3sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobja3oup1 <- renderPlot({ 
    scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp1,  
             input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
             "macrophagesobjgexpr.h5", macrophagesobjgene, 
             input$macrophagesobja3siz, input$macrophagesobja3col1, input$macrophagesobja3ord1, 
             input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) 
  }) 
  output$macrophagesobja3oup1.ui <- renderUI({ 
    plotOutput("macrophagesobja3oup1", height = pList[input$macrophagesobja3psz]) 
  }) 
  output$macrophagesobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja3drX,"_",input$macrophagesobja3drY,"_",  
                                   input$macrophagesobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja3oup1.h, width = input$macrophagesobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp1,  
                      input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja3siz, input$macrophagesobja3col1, input$macrophagesobja3ord1, 
                      input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) ) 
  }) 
  output$macrophagesobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja3drX,"_",input$macrophagesobja3drY,"_",  
                                   input$macrophagesobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja3oup1.h, width = input$macrophagesobja3oup1.w, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp1,  
                      input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja3siz, input$macrophagesobja3col1, input$macrophagesobja3ord1, 
                      input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) ) 
  }) 
   
  output$macrophagesobja3oup2 <- renderPlot({ 
    scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp2,  
             input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
             "macrophagesobjgexpr.h5", macrophagesobjgene, 
             input$macrophagesobja3siz, input$macrophagesobja3col2, input$macrophagesobja3ord2, 
             input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) 
  }) 
  output$macrophagesobja3oup2.ui <- renderUI({ 
    plotOutput("macrophagesobja3oup2", height = pList[input$macrophagesobja3psz]) 
  }) 
  output$macrophagesobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja3drX,"_",input$macrophagesobja3drY,"_",  
                                   input$macrophagesobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobja3oup2.h, width = input$macrophagesobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp2,  
                      input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja3siz, input$macrophagesobja3col2, input$macrophagesobja3ord2, 
                      input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) ) 
  }) 
  output$macrophagesobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobja3drX,"_",input$macrophagesobja3drY,"_",  
                                   input$macrophagesobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobja3oup2.h, width = input$macrophagesobja3oup2.w, 
      plot = scDRgene(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobja3drX, input$macrophagesobja3drY, input$macrophagesobja3inp2,  
                      input$macrophagesobja3sub1, input$macrophagesobja3sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobja3siz, input$macrophagesobja3col2, input$macrophagesobja3ord2, 
                      input$macrophagesobja3fsz, input$macrophagesobja3asp, input$macrophagesobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$macrophagesobjb2sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobjb2sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobjb2sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobjb2oup1 <- renderPlot({ 
    scDRcoex(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjb2drX, input$macrophagesobjb2drY,   
             input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2sub1, input$macrophagesobjb2sub2, 
             "macrophagesobjgexpr.h5", macrophagesobjgene, 
             input$macrophagesobjb2siz, input$macrophagesobjb2col1, input$macrophagesobjb2ord1, 
             input$macrophagesobjb2fsz, input$macrophagesobjb2asp, input$macrophagesobjb2txt) 
  }) 
  output$macrophagesobjb2oup1.ui <- renderUI({ 
    plotOutput("macrophagesobjb2oup1", height = pList2[input$macrophagesobjb2psz]) 
  }) 
  output$macrophagesobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjb2drX,"_",input$macrophagesobjb2drY,"_",  
                                    input$macrophagesobjb2inp1,"_",input$macrophagesobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobjb2oup1.h, width = input$macrophagesobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjb2drX, input$macrophagesobjb2drY,  
                      input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2sub1, input$macrophagesobjb2sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobjb2siz, input$macrophagesobjb2col1, input$macrophagesobjb2ord1, 
                      input$macrophagesobjb2fsz, input$macrophagesobjb2asp, input$macrophagesobjb2txt) ) 
  }) 
  output$macrophagesobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjb2drX,"_",input$macrophagesobjb2drY,"_",  
                                    input$macrophagesobjb2inp1,"_",input$macrophagesobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobjb2oup1.h, width = input$macrophagesobjb2oup1.w, 
      plot = scDRcoex(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjb2drX, input$macrophagesobjb2drY,  
                      input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2sub1, input$macrophagesobjb2sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, 
                      input$macrophagesobjb2siz, input$macrophagesobjb2col1, input$macrophagesobjb2ord1, 
                      input$macrophagesobjb2fsz, input$macrophagesobjb2asp, input$macrophagesobjb2txt) ) 
  }) 
  output$macrophagesobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2col1, input$macrophagesobjb2fsz) 
  }) 
  output$macrophagesobjb2oup2.ui <- renderUI({ 
    plotOutput("macrophagesobjb2oup2", height = "300px") 
  }) 
  output$macrophagesobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjb2drX,"_",input$macrophagesobjb2drY,"_",  
                                    input$macrophagesobjb2inp1,"_",input$macrophagesobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2col1, input$macrophagesobjb2fsz) ) 
  }) 
  output$macrophagesobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjb2drX,"_",input$macrophagesobjb2drY,"_",  
                                    input$macrophagesobjb2inp1,"_",input$macrophagesobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, input$macrophagesobjb2col1, input$macrophagesobjb2fsz) ) 
  }) 
  output$macrophagesobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjb2inp1, input$macrophagesobjb2inp2, 
                         input$macrophagesobjb2sub1, input$macrophagesobjb2sub2, "macrophagesobjgexpr.h5", macrophagesobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$macrophagesobjc1sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobjc1sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobjc1sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobjc1oup <- renderPlot({ 
    scVioBox(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc1inp1, input$macrophagesobjc1inp2, 
             input$macrophagesobjc1sub1, input$macrophagesobjc1sub2, 
             "macrophagesobjgexpr.h5", macrophagesobjgene, input$macrophagesobjc1typ, input$macrophagesobjc1pts, 
             input$macrophagesobjc1siz, input$macrophagesobjc1fsz) 
  }) 
  output$macrophagesobjc1oup.ui <- renderUI({ 
    plotOutput("macrophagesobjc1oup", height = pList2[input$macrophagesobjc1psz]) 
  }) 
  output$macrophagesobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjc1typ,"_",input$macrophagesobjc1inp1,"_",  
                                   input$macrophagesobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobjc1oup.h, width = input$macrophagesobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc1inp1, input$macrophagesobjc1inp2, 
                      input$macrophagesobjc1sub1, input$macrophagesobjc1sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, input$macrophagesobjc1typ, input$macrophagesobjc1pts, 
                      input$macrophagesobjc1siz, input$macrophagesobjc1fsz) ) 
  }) 
  output$macrophagesobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjc1typ,"_",input$macrophagesobjc1inp1,"_",  
                                   input$macrophagesobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobjc1oup.h, width = input$macrophagesobjc1oup.w, 
      plot = scVioBox(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc1inp1, input$macrophagesobjc1inp2, 
                      input$macrophagesobjc1sub1, input$macrophagesobjc1sub2, 
                      "macrophagesobjgexpr.h5", macrophagesobjgene, input$macrophagesobjc1typ, input$macrophagesobjc1pts, 
                      input$macrophagesobjc1siz, input$macrophagesobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$macrophagesobjc2sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobjc2sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobjc2sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$macrophagesobjc2oup <- renderPlot({ 
  scProp(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc2inp1, input$macrophagesobjc2inp2,  
         input$macrophagesobjc2sub1, input$macrophagesobjc2sub2, 
         input$macrophagesobjc2typ, input$macrophagesobjc2flp, input$macrophagesobjc2fsz) 
}) 
output$macrophagesobjc2oup.ui <- renderUI({ 
  plotOutput("macrophagesobjc2oup", height = pList2[input$macrophagesobjc2psz]) 
}) 
output$macrophagesobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("macrophagesobj",input$macrophagesobjc2typ,"_",input$macrophagesobjc2inp1,"_",  
                                 input$macrophagesobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$macrophagesobjc2oup.h, width = input$macrophagesobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc2inp1, input$macrophagesobjc2inp2,  
                  input$macrophagesobjc2sub1, input$macrophagesobjc2sub2, 
                  input$macrophagesobjc2typ, input$macrophagesobjc2flp, input$macrophagesobjc2fsz) ) 
  }) 
output$macrophagesobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("macrophagesobj",input$macrophagesobjc2typ,"_",input$macrophagesobjc2inp1,"_",  
                                 input$macrophagesobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$macrophagesobjc2oup.h, width = input$macrophagesobjc2oup.w, 
    plot = scProp(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjc2inp1, input$macrophagesobjc2inp2,  
                  input$macrophagesobjc2sub1, input$macrophagesobjc2sub2, 
                  input$macrophagesobjc2typ, input$macrophagesobjc2flp, input$macrophagesobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$macrophagesobjd1sub1.ui <- renderUI({ 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("macrophagesobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$macrophagesobjd1sub1non, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$macrophagesobjd1sub1all, { 
    sub = strsplit(macrophagesobjconf[UI == input$macrophagesobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "macrophagesobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$macrophagesobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$macrophagesobjd1inp, macrophagesobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$macrophagesobjd1oup <- renderPlot({ 
    scBubbHeat(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjd1inp, input$macrophagesobjd1grp, input$macrophagesobjd1plt, 
               input$macrophagesobjd1sub1, input$macrophagesobjd1sub2, "macrophagesobjgexpr.h5", macrophagesobjgene, 
               input$macrophagesobjd1scl, input$macrophagesobjd1row, input$macrophagesobjd1col, 
               input$macrophagesobjd1cols, input$macrophagesobjd1fsz) 
  }) 
  output$macrophagesobjd1oup.ui <- renderUI({ 
    plotOutput("macrophagesobjd1oup", height = pList3[input$macrophagesobjd1psz]) 
  }) 
  output$macrophagesobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjd1plt,"_",input$macrophagesobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$macrophagesobjd1oup.h, width = input$macrophagesobjd1oup.w, 
      plot = scBubbHeat(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjd1inp, input$macrophagesobjd1grp, input$macrophagesobjd1plt, 
                        input$macrophagesobjd1sub1, input$macrophagesobjd1sub2, "macrophagesobjgexpr.h5", macrophagesobjgene, 
                        input$macrophagesobjd1scl, input$macrophagesobjd1row, input$macrophagesobjd1col, 
                        input$macrophagesobjd1cols, input$macrophagesobjd1fsz, save = TRUE) ) 
  }) 
  output$macrophagesobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("macrophagesobj",input$macrophagesobjd1plt,"_",input$macrophagesobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$macrophagesobjd1oup.h, width = input$macrophagesobjd1oup.w, 
      plot = scBubbHeat(macrophagesobjconf, macrophagesobjmeta, input$macrophagesobjd1inp, input$macrophagesobjd1grp, input$macrophagesobjd1plt, 
                        input$macrophagesobjd1sub1, input$macrophagesobjd1sub2, "macrophagesobjgexpr.h5", macrophagesobjgene, 
                        input$macrophagesobjd1scl, input$macrophagesobjd1row, input$macrophagesobjd1col, 
                        input$macrophagesobjd1cols, input$macrophagesobjd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "dendriticobja1inp2", choices = names(dendriticobjgene), server = TRUE, 
                       selected = dendriticobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "dendriticobja3inp1", choices = names(dendriticobjgene), server = TRUE, 
                       selected = dendriticobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "dendriticobja3inp2", choices = names(dendriticobjgene), server = TRUE, 
                       selected = dendriticobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "dendriticobjb2inp1", choices = names(dendriticobjgene), server = TRUE, 
                       selected = dendriticobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "dendriticobjb2inp2", choices = names(dendriticobjgene), server = TRUE, 
                       selected = dendriticobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "dendriticobjc1inp2", server = TRUE, 
                       choices = c(dendriticobjconf[is.na(fID)]$UI,names(dendriticobjgene)), 
                       selected = dendriticobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(dendriticobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$dendriticobja1sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobja1sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobja1sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobja1oup1 <- renderPlot({ 
    scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp1,  
             input$dendriticobja1sub1, input$dendriticobja1sub2, 
             input$dendriticobja1siz, input$dendriticobja1col1, input$dendriticobja1ord1, 
             input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt, input$dendriticobja1lab1) 
  }) 
  output$dendriticobja1oup1.ui <- renderUI({ 
    plotOutput("dendriticobja1oup1", height = pList[input$dendriticobja1psz]) 
  }) 
  output$dendriticobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja1drX,"_",input$dendriticobja1drY,"_",  
                                   input$dendriticobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja1oup1.h, width = input$dendriticobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp1,   
                      input$dendriticobja1sub1, input$dendriticobja1sub2, 
                      input$dendriticobja1siz, input$dendriticobja1col1, input$dendriticobja1ord1,  
                      input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt, input$dendriticobja1lab1) ) 
  }) 
  output$dendriticobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja1drX,"_",input$dendriticobja1drY,"_",  
                                   input$dendriticobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja1oup1.h, width = input$dendriticobja1oup1.w, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp1,   
                      input$dendriticobja1sub1, input$dendriticobja1sub2, 
                      input$dendriticobja1siz, input$dendriticobja1col1, input$dendriticobja1ord1,  
                      input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt, input$dendriticobja1lab1) ) 
  }) 
  output$dendriticobja1.dt <- renderDataTable({ 
    ggData = scDRnum(dendriticobjconf, dendriticobjmeta, input$dendriticobja1inp1, input$dendriticobja1inp2, 
                     input$dendriticobja1sub1, input$dendriticobja1sub2, 
                     "dendriticobjgexpr.h5", dendriticobjgene, input$dendriticobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$dendriticobja1oup2 <- renderPlot({ 
    scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp2,  
             input$dendriticobja1sub1, input$dendriticobja1sub2, 
             "dendriticobjgexpr.h5", dendriticobjgene, 
             input$dendriticobja1siz, input$dendriticobja1col2, input$dendriticobja1ord2, 
             input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt) 
  }) 
  output$dendriticobja1oup2.ui <- renderUI({ 
    plotOutput("dendriticobja1oup2", height = pList[input$dendriticobja1psz]) 
  }) 
  output$dendriticobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja1drX,"_",input$dendriticobja1drY,"_",  
                                   input$dendriticobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja1oup2.h, width = input$dendriticobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp2,  
                      input$dendriticobja1sub1, input$dendriticobja1sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja1siz, input$dendriticobja1col2, input$dendriticobja1ord2, 
                      input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt) ) 
  }) 
  output$dendriticobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja1drX,"_",input$dendriticobja1drY,"_",  
                                   input$dendriticobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja1oup2.h, width = input$dendriticobja1oup2.w, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja1drX, input$dendriticobja1drY, input$dendriticobja1inp2,  
                      input$dendriticobja1sub1, input$dendriticobja1sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja1siz, input$dendriticobja1col2, input$dendriticobja1ord2, 
                      input$dendriticobja1fsz, input$dendriticobja1asp, input$dendriticobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$dendriticobja2sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobja2sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobja2sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobja2oup1 <- renderPlot({ 
    scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp1,  
             input$dendriticobja2sub1, input$dendriticobja2sub2, 
             input$dendriticobja2siz, input$dendriticobja2col1, input$dendriticobja2ord1, 
             input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab1) 
  }) 
  output$dendriticobja2oup1.ui <- renderUI({ 
    plotOutput("dendriticobja2oup1", height = pList[input$dendriticobja2psz]) 
  }) 
  output$dendriticobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja2drX,"_",input$dendriticobja2drY,"_",  
                                   input$dendriticobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja2oup1.h, width = input$dendriticobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp1,   
                      input$dendriticobja2sub1, input$dendriticobja2sub2, 
                      input$dendriticobja2siz, input$dendriticobja2col1, input$dendriticobja2ord1,  
                      input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab1) ) 
  }) 
  output$dendriticobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja2drX,"_",input$dendriticobja2drY,"_",  
                                   input$dendriticobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja2oup1.h, width = input$dendriticobja2oup1.w, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp1,   
                      input$dendriticobja2sub1, input$dendriticobja2sub2, 
                      input$dendriticobja2siz, input$dendriticobja2col1, input$dendriticobja2ord1,  
                      input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab1) ) 
  }) 
   
  output$dendriticobja2oup2 <- renderPlot({ 
    scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp2,  
             input$dendriticobja2sub1, input$dendriticobja2sub2, 
             input$dendriticobja2siz, input$dendriticobja2col2, input$dendriticobja2ord2, 
             input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab2) 
  }) 
  output$dendriticobja2oup2.ui <- renderUI({ 
    plotOutput("dendriticobja2oup2", height = pList[input$dendriticobja2psz]) 
  }) 
  output$dendriticobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja2drX,"_",input$dendriticobja2drY,"_",  
                                   input$dendriticobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja2oup2.h, width = input$dendriticobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp2,   
                      input$dendriticobja2sub1, input$dendriticobja2sub2, 
                      input$dendriticobja2siz, input$dendriticobja2col2, input$dendriticobja2ord2,  
                      input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab2) ) 
  }) 
  output$dendriticobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja2drX,"_",input$dendriticobja2drY,"_",  
                                   input$dendriticobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja2oup2.h, width = input$dendriticobja2oup2.w, 
      plot = scDRcell(dendriticobjconf, dendriticobjmeta, input$dendriticobja2drX, input$dendriticobja2drY, input$dendriticobja2inp2,   
                      input$dendriticobja2sub1, input$dendriticobja2sub2, 
                      input$dendriticobja2siz, input$dendriticobja2col2, input$dendriticobja2ord2,  
                      input$dendriticobja2fsz, input$dendriticobja2asp, input$dendriticobja2txt, input$dendriticobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$dendriticobja3sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobja3sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobja3sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobja3oup1 <- renderPlot({ 
    scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp1,  
             input$dendriticobja3sub1, input$dendriticobja3sub2, 
             "dendriticobjgexpr.h5", dendriticobjgene, 
             input$dendriticobja3siz, input$dendriticobja3col1, input$dendriticobja3ord1, 
             input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) 
  }) 
  output$dendriticobja3oup1.ui <- renderUI({ 
    plotOutput("dendriticobja3oup1", height = pList[input$dendriticobja3psz]) 
  }) 
  output$dendriticobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja3drX,"_",input$dendriticobja3drY,"_",  
                                   input$dendriticobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja3oup1.h, width = input$dendriticobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp1,  
                      input$dendriticobja3sub1, input$dendriticobja3sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja3siz, input$dendriticobja3col1, input$dendriticobja3ord1, 
                      input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) ) 
  }) 
  output$dendriticobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja3drX,"_",input$dendriticobja3drY,"_",  
                                   input$dendriticobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja3oup1.h, width = input$dendriticobja3oup1.w, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp1,  
                      input$dendriticobja3sub1, input$dendriticobja3sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja3siz, input$dendriticobja3col1, input$dendriticobja3ord1, 
                      input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) ) 
  }) 
   
  output$dendriticobja3oup2 <- renderPlot({ 
    scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp2,  
             input$dendriticobja3sub1, input$dendriticobja3sub2, 
             "dendriticobjgexpr.h5", dendriticobjgene, 
             input$dendriticobja3siz, input$dendriticobja3col2, input$dendriticobja3ord2, 
             input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) 
  }) 
  output$dendriticobja3oup2.ui <- renderUI({ 
    plotOutput("dendriticobja3oup2", height = pList[input$dendriticobja3psz]) 
  }) 
  output$dendriticobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja3drX,"_",input$dendriticobja3drY,"_",  
                                   input$dendriticobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobja3oup2.h, width = input$dendriticobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp2,  
                      input$dendriticobja3sub1, input$dendriticobja3sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja3siz, input$dendriticobja3col2, input$dendriticobja3ord2, 
                      input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) ) 
  }) 
  output$dendriticobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobja3drX,"_",input$dendriticobja3drY,"_",  
                                   input$dendriticobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobja3oup2.h, width = input$dendriticobja3oup2.w, 
      plot = scDRgene(dendriticobjconf, dendriticobjmeta, input$dendriticobja3drX, input$dendriticobja3drY, input$dendriticobja3inp2,  
                      input$dendriticobja3sub1, input$dendriticobja3sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobja3siz, input$dendriticobja3col2, input$dendriticobja3ord2, 
                      input$dendriticobja3fsz, input$dendriticobja3asp, input$dendriticobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$dendriticobjb2sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobjb2sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobjb2sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobjb2oup1 <- renderPlot({ 
    scDRcoex(dendriticobjconf, dendriticobjmeta, input$dendriticobjb2drX, input$dendriticobjb2drY,   
             input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2sub1, input$dendriticobjb2sub2, 
             "dendriticobjgexpr.h5", dendriticobjgene, 
             input$dendriticobjb2siz, input$dendriticobjb2col1, input$dendriticobjb2ord1, 
             input$dendriticobjb2fsz, input$dendriticobjb2asp, input$dendriticobjb2txt) 
  }) 
  output$dendriticobjb2oup1.ui <- renderUI({ 
    plotOutput("dendriticobjb2oup1", height = pList2[input$dendriticobjb2psz]) 
  }) 
  output$dendriticobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjb2drX,"_",input$dendriticobjb2drY,"_",  
                                    input$dendriticobjb2inp1,"_",input$dendriticobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobjb2oup1.h, width = input$dendriticobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(dendriticobjconf, dendriticobjmeta, input$dendriticobjb2drX, input$dendriticobjb2drY,  
                      input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2sub1, input$dendriticobjb2sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobjb2siz, input$dendriticobjb2col1, input$dendriticobjb2ord1, 
                      input$dendriticobjb2fsz, input$dendriticobjb2asp, input$dendriticobjb2txt) ) 
  }) 
  output$dendriticobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjb2drX,"_",input$dendriticobjb2drY,"_",  
                                    input$dendriticobjb2inp1,"_",input$dendriticobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobjb2oup1.h, width = input$dendriticobjb2oup1.w, 
      plot = scDRcoex(dendriticobjconf, dendriticobjmeta, input$dendriticobjb2drX, input$dendriticobjb2drY,  
                      input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2sub1, input$dendriticobjb2sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, 
                      input$dendriticobjb2siz, input$dendriticobjb2col1, input$dendriticobjb2ord1, 
                      input$dendriticobjb2fsz, input$dendriticobjb2asp, input$dendriticobjb2txt) ) 
  }) 
  output$dendriticobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2col1, input$dendriticobjb2fsz) 
  }) 
  output$dendriticobjb2oup2.ui <- renderUI({ 
    plotOutput("dendriticobjb2oup2", height = "300px") 
  }) 
  output$dendriticobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjb2drX,"_",input$dendriticobjb2drY,"_",  
                                    input$dendriticobjb2inp1,"_",input$dendriticobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2col1, input$dendriticobjb2fsz) ) 
  }) 
  output$dendriticobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjb2drX,"_",input$dendriticobjb2drY,"_",  
                                    input$dendriticobjb2inp1,"_",input$dendriticobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$dendriticobjb2inp1, input$dendriticobjb2inp2, input$dendriticobjb2col1, input$dendriticobjb2fsz) ) 
  }) 
  output$dendriticobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(dendriticobjconf, dendriticobjmeta, input$dendriticobjb2inp1, input$dendriticobjb2inp2, 
                         input$dendriticobjb2sub1, input$dendriticobjb2sub2, "dendriticobjgexpr.h5", dendriticobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$dendriticobjc1sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobjc1sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobjc1sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobjc1oup <- renderPlot({ 
    scVioBox(dendriticobjconf, dendriticobjmeta, input$dendriticobjc1inp1, input$dendriticobjc1inp2, 
             input$dendriticobjc1sub1, input$dendriticobjc1sub2, 
             "dendriticobjgexpr.h5", dendriticobjgene, input$dendriticobjc1typ, input$dendriticobjc1pts, 
             input$dendriticobjc1siz, input$dendriticobjc1fsz) 
  }) 
  output$dendriticobjc1oup.ui <- renderUI({ 
    plotOutput("dendriticobjc1oup", height = pList2[input$dendriticobjc1psz]) 
  }) 
  output$dendriticobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjc1typ,"_",input$dendriticobjc1inp1,"_",  
                                   input$dendriticobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobjc1oup.h, width = input$dendriticobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(dendriticobjconf, dendriticobjmeta, input$dendriticobjc1inp1, input$dendriticobjc1inp2, 
                      input$dendriticobjc1sub1, input$dendriticobjc1sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, input$dendriticobjc1typ, input$dendriticobjc1pts, 
                      input$dendriticobjc1siz, input$dendriticobjc1fsz) ) 
  }) 
  output$dendriticobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjc1typ,"_",input$dendriticobjc1inp1,"_",  
                                   input$dendriticobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobjc1oup.h, width = input$dendriticobjc1oup.w, 
      plot = scVioBox(dendriticobjconf, dendriticobjmeta, input$dendriticobjc1inp1, input$dendriticobjc1inp2, 
                      input$dendriticobjc1sub1, input$dendriticobjc1sub2, 
                      "dendriticobjgexpr.h5", dendriticobjgene, input$dendriticobjc1typ, input$dendriticobjc1pts, 
                      input$dendriticobjc1siz, input$dendriticobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$dendriticobjc2sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobjc2sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobjc2sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$dendriticobjc2oup <- renderPlot({ 
  scProp(dendriticobjconf, dendriticobjmeta, input$dendriticobjc2inp1, input$dendriticobjc2inp2,  
         input$dendriticobjc2sub1, input$dendriticobjc2sub2, 
         input$dendriticobjc2typ, input$dendriticobjc2flp, input$dendriticobjc2fsz) 
}) 
output$dendriticobjc2oup.ui <- renderUI({ 
  plotOutput("dendriticobjc2oup", height = pList2[input$dendriticobjc2psz]) 
}) 
output$dendriticobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("dendriticobj",input$dendriticobjc2typ,"_",input$dendriticobjc2inp1,"_",  
                                 input$dendriticobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$dendriticobjc2oup.h, width = input$dendriticobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(dendriticobjconf, dendriticobjmeta, input$dendriticobjc2inp1, input$dendriticobjc2inp2,  
                  input$dendriticobjc2sub1, input$dendriticobjc2sub2, 
                  input$dendriticobjc2typ, input$dendriticobjc2flp, input$dendriticobjc2fsz) ) 
  }) 
output$dendriticobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("dendriticobj",input$dendriticobjc2typ,"_",input$dendriticobjc2inp1,"_",  
                                 input$dendriticobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$dendriticobjc2oup.h, width = input$dendriticobjc2oup.w, 
    plot = scProp(dendriticobjconf, dendriticobjmeta, input$dendriticobjc2inp1, input$dendriticobjc2inp2,  
                  input$dendriticobjc2sub1, input$dendriticobjc2sub2, 
                  input$dendriticobjc2typ, input$dendriticobjc2flp, input$dendriticobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$dendriticobjd1sub1.ui <- renderUI({ 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("dendriticobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$dendriticobjd1sub1non, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$dendriticobjd1sub1all, { 
    sub = strsplit(dendriticobjconf[UI == input$dendriticobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "dendriticobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$dendriticobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$dendriticobjd1inp, dendriticobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$dendriticobjd1oup <- renderPlot({ 
    scBubbHeat(dendriticobjconf, dendriticobjmeta, input$dendriticobjd1inp, input$dendriticobjd1grp, input$dendriticobjd1plt, 
               input$dendriticobjd1sub1, input$dendriticobjd1sub2, "dendriticobjgexpr.h5", dendriticobjgene, 
               input$dendriticobjd1scl, input$dendriticobjd1row, input$dendriticobjd1col, 
               input$dendriticobjd1cols, input$dendriticobjd1fsz) 
  }) 
  output$dendriticobjd1oup.ui <- renderUI({ 
    plotOutput("dendriticobjd1oup", height = pList3[input$dendriticobjd1psz]) 
  }) 
  output$dendriticobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjd1plt,"_",input$dendriticobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$dendriticobjd1oup.h, width = input$dendriticobjd1oup.w, 
      plot = scBubbHeat(dendriticobjconf, dendriticobjmeta, input$dendriticobjd1inp, input$dendriticobjd1grp, input$dendriticobjd1plt, 
                        input$dendriticobjd1sub1, input$dendriticobjd1sub2, "dendriticobjgexpr.h5", dendriticobjgene, 
                        input$dendriticobjd1scl, input$dendriticobjd1row, input$dendriticobjd1col, 
                        input$dendriticobjd1cols, input$dendriticobjd1fsz, save = TRUE) ) 
  }) 
  output$dendriticobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("dendriticobj",input$dendriticobjd1plt,"_",input$dendriticobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$dendriticobjd1oup.h, width = input$dendriticobjd1oup.w, 
      plot = scBubbHeat(dendriticobjconf, dendriticobjmeta, input$dendriticobjd1inp, input$dendriticobjd1grp, input$dendriticobjd1plt, 
                        input$dendriticobjd1sub1, input$dendriticobjd1sub2, "dendriticobjgexpr.h5", dendriticobjgene, 
                        input$dendriticobjd1scl, input$dendriticobjd1row, input$dendriticobjd1col, 
                        input$dendriticobjd1cols, input$dendriticobjd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "neutrophilsobja1inp2", choices = names(neutrophilsobjgene), server = TRUE, 
                       selected = neutrophilsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "neutrophilsobja3inp1", choices = names(neutrophilsobjgene), server = TRUE, 
                       selected = neutrophilsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "neutrophilsobja3inp2", choices = names(neutrophilsobjgene), server = TRUE, 
                       selected = neutrophilsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "neutrophilsobjb2inp1", choices = names(neutrophilsobjgene), server = TRUE, 
                       selected = neutrophilsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "neutrophilsobjb2inp2", choices = names(neutrophilsobjgene), server = TRUE, 
                       selected = neutrophilsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "neutrophilsobjc1inp2", server = TRUE, 
                       choices = c(neutrophilsobjconf[is.na(fID)]$UI,names(neutrophilsobjgene)), 
                       selected = neutrophilsobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(neutrophilsobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$neutrophilsobja1sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobja1sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobja1sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobja1oup1 <- renderPlot({ 
    scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp1,  
             input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
             input$neutrophilsobja1siz, input$neutrophilsobja1col1, input$neutrophilsobja1ord1, 
             input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt, input$neutrophilsobja1lab1) 
  }) 
  output$neutrophilsobja1oup1.ui <- renderUI({ 
    plotOutput("neutrophilsobja1oup1", height = pList[input$neutrophilsobja1psz]) 
  }) 
  output$neutrophilsobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja1drX,"_",input$neutrophilsobja1drY,"_",  
                                   input$neutrophilsobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja1oup1.h, width = input$neutrophilsobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp1,   
                      input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
                      input$neutrophilsobja1siz, input$neutrophilsobja1col1, input$neutrophilsobja1ord1,  
                      input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt, input$neutrophilsobja1lab1) ) 
  }) 
  output$neutrophilsobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja1drX,"_",input$neutrophilsobja1drY,"_",  
                                   input$neutrophilsobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja1oup1.h, width = input$neutrophilsobja1oup1.w, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp1,   
                      input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
                      input$neutrophilsobja1siz, input$neutrophilsobja1col1, input$neutrophilsobja1ord1,  
                      input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt, input$neutrophilsobja1lab1) ) 
  }) 
  output$neutrophilsobja1.dt <- renderDataTable({ 
    ggData = scDRnum(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1inp1, input$neutrophilsobja1inp2, 
                     input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
                     "neutrophilsobjgexpr.h5", neutrophilsobjgene, input$neutrophilsobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$neutrophilsobja1oup2 <- renderPlot({ 
    scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp2,  
             input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
             "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
             input$neutrophilsobja1siz, input$neutrophilsobja1col2, input$neutrophilsobja1ord2, 
             input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt) 
  }) 
  output$neutrophilsobja1oup2.ui <- renderUI({ 
    plotOutput("neutrophilsobja1oup2", height = pList[input$neutrophilsobja1psz]) 
  }) 
  output$neutrophilsobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja1drX,"_",input$neutrophilsobja1drY,"_",  
                                   input$neutrophilsobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja1oup2.h, width = input$neutrophilsobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp2,  
                      input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja1siz, input$neutrophilsobja1col2, input$neutrophilsobja1ord2, 
                      input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt) ) 
  }) 
  output$neutrophilsobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja1drX,"_",input$neutrophilsobja1drY,"_",  
                                   input$neutrophilsobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja1oup2.h, width = input$neutrophilsobja1oup2.w, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja1drX, input$neutrophilsobja1drY, input$neutrophilsobja1inp2,  
                      input$neutrophilsobja1sub1, input$neutrophilsobja1sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja1siz, input$neutrophilsobja1col2, input$neutrophilsobja1ord2, 
                      input$neutrophilsobja1fsz, input$neutrophilsobja1asp, input$neutrophilsobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$neutrophilsobja2sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobja2sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobja2sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobja2oup1 <- renderPlot({ 
    scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp1,  
             input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
             input$neutrophilsobja2siz, input$neutrophilsobja2col1, input$neutrophilsobja2ord1, 
             input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab1) 
  }) 
  output$neutrophilsobja2oup1.ui <- renderUI({ 
    plotOutput("neutrophilsobja2oup1", height = pList[input$neutrophilsobja2psz]) 
  }) 
  output$neutrophilsobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja2drX,"_",input$neutrophilsobja2drY,"_",  
                                   input$neutrophilsobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja2oup1.h, width = input$neutrophilsobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp1,   
                      input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
                      input$neutrophilsobja2siz, input$neutrophilsobja2col1, input$neutrophilsobja2ord1,  
                      input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab1) ) 
  }) 
  output$neutrophilsobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja2drX,"_",input$neutrophilsobja2drY,"_",  
                                   input$neutrophilsobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja2oup1.h, width = input$neutrophilsobja2oup1.w, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp1,   
                      input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
                      input$neutrophilsobja2siz, input$neutrophilsobja2col1, input$neutrophilsobja2ord1,  
                      input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab1) ) 
  }) 
   
  output$neutrophilsobja2oup2 <- renderPlot({ 
    scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp2,  
             input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
             input$neutrophilsobja2siz, input$neutrophilsobja2col2, input$neutrophilsobja2ord2, 
             input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab2) 
  }) 
  output$neutrophilsobja2oup2.ui <- renderUI({ 
    plotOutput("neutrophilsobja2oup2", height = pList[input$neutrophilsobja2psz]) 
  }) 
  output$neutrophilsobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja2drX,"_",input$neutrophilsobja2drY,"_",  
                                   input$neutrophilsobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja2oup2.h, width = input$neutrophilsobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp2,   
                      input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
                      input$neutrophilsobja2siz, input$neutrophilsobja2col2, input$neutrophilsobja2ord2,  
                      input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab2) ) 
  }) 
  output$neutrophilsobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja2drX,"_",input$neutrophilsobja2drY,"_",  
                                   input$neutrophilsobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja2oup2.h, width = input$neutrophilsobja2oup2.w, 
      plot = scDRcell(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja2drX, input$neutrophilsobja2drY, input$neutrophilsobja2inp2,   
                      input$neutrophilsobja2sub1, input$neutrophilsobja2sub2, 
                      input$neutrophilsobja2siz, input$neutrophilsobja2col2, input$neutrophilsobja2ord2,  
                      input$neutrophilsobja2fsz, input$neutrophilsobja2asp, input$neutrophilsobja2txt, input$neutrophilsobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$neutrophilsobja3sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobja3sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobja3sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobja3oup1 <- renderPlot({ 
    scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp1,  
             input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
             "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
             input$neutrophilsobja3siz, input$neutrophilsobja3col1, input$neutrophilsobja3ord1, 
             input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) 
  }) 
  output$neutrophilsobja3oup1.ui <- renderUI({ 
    plotOutput("neutrophilsobja3oup1", height = pList[input$neutrophilsobja3psz]) 
  }) 
  output$neutrophilsobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja3drX,"_",input$neutrophilsobja3drY,"_",  
                                   input$neutrophilsobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja3oup1.h, width = input$neutrophilsobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp1,  
                      input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja3siz, input$neutrophilsobja3col1, input$neutrophilsobja3ord1, 
                      input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) ) 
  }) 
  output$neutrophilsobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja3drX,"_",input$neutrophilsobja3drY,"_",  
                                   input$neutrophilsobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja3oup1.h, width = input$neutrophilsobja3oup1.w, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp1,  
                      input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja3siz, input$neutrophilsobja3col1, input$neutrophilsobja3ord1, 
                      input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) ) 
  }) 
   
  output$neutrophilsobja3oup2 <- renderPlot({ 
    scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp2,  
             input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
             "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
             input$neutrophilsobja3siz, input$neutrophilsobja3col2, input$neutrophilsobja3ord2, 
             input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) 
  }) 
  output$neutrophilsobja3oup2.ui <- renderUI({ 
    plotOutput("neutrophilsobja3oup2", height = pList[input$neutrophilsobja3psz]) 
  }) 
  output$neutrophilsobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja3drX,"_",input$neutrophilsobja3drY,"_",  
                                   input$neutrophilsobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobja3oup2.h, width = input$neutrophilsobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp2,  
                      input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja3siz, input$neutrophilsobja3col2, input$neutrophilsobja3ord2, 
                      input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) ) 
  }) 
  output$neutrophilsobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobja3drX,"_",input$neutrophilsobja3drY,"_",  
                                   input$neutrophilsobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobja3oup2.h, width = input$neutrophilsobja3oup2.w, 
      plot = scDRgene(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobja3drX, input$neutrophilsobja3drY, input$neutrophilsobja3inp2,  
                      input$neutrophilsobja3sub1, input$neutrophilsobja3sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobja3siz, input$neutrophilsobja3col2, input$neutrophilsobja3ord2, 
                      input$neutrophilsobja3fsz, input$neutrophilsobja3asp, input$neutrophilsobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$neutrophilsobjb2sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobjb2sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobjb2sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobjb2oup1 <- renderPlot({ 
    scDRcoex(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjb2drX, input$neutrophilsobjb2drY,   
             input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2sub1, input$neutrophilsobjb2sub2, 
             "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
             input$neutrophilsobjb2siz, input$neutrophilsobjb2col1, input$neutrophilsobjb2ord1, 
             input$neutrophilsobjb2fsz, input$neutrophilsobjb2asp, input$neutrophilsobjb2txt) 
  }) 
  output$neutrophilsobjb2oup1.ui <- renderUI({ 
    plotOutput("neutrophilsobjb2oup1", height = pList2[input$neutrophilsobjb2psz]) 
  }) 
  output$neutrophilsobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjb2drX,"_",input$neutrophilsobjb2drY,"_",  
                                    input$neutrophilsobjb2inp1,"_",input$neutrophilsobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobjb2oup1.h, width = input$neutrophilsobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjb2drX, input$neutrophilsobjb2drY,  
                      input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2sub1, input$neutrophilsobjb2sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobjb2siz, input$neutrophilsobjb2col1, input$neutrophilsobjb2ord1, 
                      input$neutrophilsobjb2fsz, input$neutrophilsobjb2asp, input$neutrophilsobjb2txt) ) 
  }) 
  output$neutrophilsobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjb2drX,"_",input$neutrophilsobjb2drY,"_",  
                                    input$neutrophilsobjb2inp1,"_",input$neutrophilsobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobjb2oup1.h, width = input$neutrophilsobjb2oup1.w, 
      plot = scDRcoex(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjb2drX, input$neutrophilsobjb2drY,  
                      input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2sub1, input$neutrophilsobjb2sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                      input$neutrophilsobjb2siz, input$neutrophilsobjb2col1, input$neutrophilsobjb2ord1, 
                      input$neutrophilsobjb2fsz, input$neutrophilsobjb2asp, input$neutrophilsobjb2txt) ) 
  }) 
  output$neutrophilsobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2col1, input$neutrophilsobjb2fsz) 
  }) 
  output$neutrophilsobjb2oup2.ui <- renderUI({ 
    plotOutput("neutrophilsobjb2oup2", height = "300px") 
  }) 
  output$neutrophilsobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjb2drX,"_",input$neutrophilsobjb2drY,"_",  
                                    input$neutrophilsobjb2inp1,"_",input$neutrophilsobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2col1, input$neutrophilsobjb2fsz) ) 
  }) 
  output$neutrophilsobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjb2drX,"_",input$neutrophilsobjb2drY,"_",  
                                    input$neutrophilsobjb2inp1,"_",input$neutrophilsobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, input$neutrophilsobjb2col1, input$neutrophilsobjb2fsz) ) 
  }) 
  output$neutrophilsobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjb2inp1, input$neutrophilsobjb2inp2, 
                         input$neutrophilsobjb2sub1, input$neutrophilsobjb2sub2, "neutrophilsobjgexpr.h5", neutrophilsobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$neutrophilsobjc1sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobjc1sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobjc1sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobjc1oup <- renderPlot({ 
    scVioBox(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc1inp1, input$neutrophilsobjc1inp2, 
             input$neutrophilsobjc1sub1, input$neutrophilsobjc1sub2, 
             "neutrophilsobjgexpr.h5", neutrophilsobjgene, input$neutrophilsobjc1typ, input$neutrophilsobjc1pts, 
             input$neutrophilsobjc1siz, input$neutrophilsobjc1fsz) 
  }) 
  output$neutrophilsobjc1oup.ui <- renderUI({ 
    plotOutput("neutrophilsobjc1oup", height = pList2[input$neutrophilsobjc1psz]) 
  }) 
  output$neutrophilsobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjc1typ,"_",input$neutrophilsobjc1inp1,"_",  
                                   input$neutrophilsobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobjc1oup.h, width = input$neutrophilsobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc1inp1, input$neutrophilsobjc1inp2, 
                      input$neutrophilsobjc1sub1, input$neutrophilsobjc1sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, input$neutrophilsobjc1typ, input$neutrophilsobjc1pts, 
                      input$neutrophilsobjc1siz, input$neutrophilsobjc1fsz) ) 
  }) 
  output$neutrophilsobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjc1typ,"_",input$neutrophilsobjc1inp1,"_",  
                                   input$neutrophilsobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobjc1oup.h, width = input$neutrophilsobjc1oup.w, 
      plot = scVioBox(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc1inp1, input$neutrophilsobjc1inp2, 
                      input$neutrophilsobjc1sub1, input$neutrophilsobjc1sub2, 
                      "neutrophilsobjgexpr.h5", neutrophilsobjgene, input$neutrophilsobjc1typ, input$neutrophilsobjc1pts, 
                      input$neutrophilsobjc1siz, input$neutrophilsobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$neutrophilsobjc2sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobjc2sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobjc2sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$neutrophilsobjc2oup <- renderPlot({ 
  scProp(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc2inp1, input$neutrophilsobjc2inp2,  
         input$neutrophilsobjc2sub1, input$neutrophilsobjc2sub2, 
         input$neutrophilsobjc2typ, input$neutrophilsobjc2flp, input$neutrophilsobjc2fsz) 
}) 
output$neutrophilsobjc2oup.ui <- renderUI({ 
  plotOutput("neutrophilsobjc2oup", height = pList2[input$neutrophilsobjc2psz]) 
}) 
output$neutrophilsobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("neutrophilsobj",input$neutrophilsobjc2typ,"_",input$neutrophilsobjc2inp1,"_",  
                                 input$neutrophilsobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$neutrophilsobjc2oup.h, width = input$neutrophilsobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc2inp1, input$neutrophilsobjc2inp2,  
                  input$neutrophilsobjc2sub1, input$neutrophilsobjc2sub2, 
                  input$neutrophilsobjc2typ, input$neutrophilsobjc2flp, input$neutrophilsobjc2fsz) ) 
  }) 
output$neutrophilsobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("neutrophilsobj",input$neutrophilsobjc2typ,"_",input$neutrophilsobjc2inp1,"_",  
                                 input$neutrophilsobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$neutrophilsobjc2oup.h, width = input$neutrophilsobjc2oup.w, 
    plot = scProp(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjc2inp1, input$neutrophilsobjc2inp2,  
                  input$neutrophilsobjc2sub1, input$neutrophilsobjc2sub2, 
                  input$neutrophilsobjc2typ, input$neutrophilsobjc2flp, input$neutrophilsobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$neutrophilsobjd1sub1.ui <- renderUI({ 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("neutrophilsobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$neutrophilsobjd1sub1non, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$neutrophilsobjd1sub1all, { 
    sub = strsplit(neutrophilsobjconf[UI == input$neutrophilsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "neutrophilsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$neutrophilsobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$neutrophilsobjd1inp, neutrophilsobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$neutrophilsobjd1oup <- renderPlot({ 
    scBubbHeat(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjd1inp, input$neutrophilsobjd1grp, input$neutrophilsobjd1plt, 
               input$neutrophilsobjd1sub1, input$neutrophilsobjd1sub2, "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
               input$neutrophilsobjd1scl, input$neutrophilsobjd1row, input$neutrophilsobjd1col, 
               input$neutrophilsobjd1cols, input$neutrophilsobjd1fsz) 
  }) 
  output$neutrophilsobjd1oup.ui <- renderUI({ 
    plotOutput("neutrophilsobjd1oup", height = pList3[input$neutrophilsobjd1psz]) 
  }) 
  output$neutrophilsobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjd1plt,"_",input$neutrophilsobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$neutrophilsobjd1oup.h, width = input$neutrophilsobjd1oup.w, 
      plot = scBubbHeat(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjd1inp, input$neutrophilsobjd1grp, input$neutrophilsobjd1plt, 
                        input$neutrophilsobjd1sub1, input$neutrophilsobjd1sub2, "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                        input$neutrophilsobjd1scl, input$neutrophilsobjd1row, input$neutrophilsobjd1col, 
                        input$neutrophilsobjd1cols, input$neutrophilsobjd1fsz, save = TRUE) ) 
  }) 
  output$neutrophilsobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("neutrophilsobj",input$neutrophilsobjd1plt,"_",input$neutrophilsobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$neutrophilsobjd1oup.h, width = input$neutrophilsobjd1oup.w, 
      plot = scBubbHeat(neutrophilsobjconf, neutrophilsobjmeta, input$neutrophilsobjd1inp, input$neutrophilsobjd1grp, input$neutrophilsobjd1plt, 
                        input$neutrophilsobjd1sub1, input$neutrophilsobjd1sub2, "neutrophilsobjgexpr.h5", neutrophilsobjgene, 
                        input$neutrophilsobjd1scl, input$neutrophilsobjd1row, input$neutrophilsobjd1col, 
                        input$neutrophilsobjd1cols, input$neutrophilsobjd1fsz, save = TRUE) ) 
  }) 
   
   
   optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
  updateSelectizeInput(session, "mastcellsobja1inp2", choices = names(mastcellsobjgene), server = TRUE, 
                       selected = mastcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "mastcellsobja3inp1", choices = names(mastcellsobjgene), server = TRUE, 
                       selected = mastcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "mastcellsobja3inp2", choices = names(mastcellsobjgene), server = TRUE, 
                       selected = mastcellsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "mastcellsobjb2inp1", choices = names(mastcellsobjgene), server = TRUE, 
                       selected = mastcellsobjdef$gene1, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "mastcellsobjb2inp2", choices = names(mastcellsobjgene), server = TRUE, 
                       selected = mastcellsobjdef$gene2, options = list( 
                         maxOptions = 7, create = TRUE, persist = TRUE, render = I(optCrt))) 
  updateSelectizeInput(session, "mastcellsobjc1inp2", server = TRUE, 
                       choices = c(mastcellsobjconf[is.na(fID)]$UI,names(mastcellsobjgene)), 
                       selected = mastcellsobjconf[is.na(fID)]$UI[1], options = list( 
                         maxOptions = length(mastcellsobjconf[is.na(fID)]$UI) + 3, 
                         create = TRUE, persist = TRUE, render = I(optCrt))) 
 
 
  ### Plots for tab a1 
  output$mastcellsobja1sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobja1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobja1sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobja1sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobja1oup1 <- renderPlot({ 
    scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp1,  
             input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
             input$mastcellsobja1siz, input$mastcellsobja1col1, input$mastcellsobja1ord1, 
             input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt, input$mastcellsobja1lab1) 
  }) 
  output$mastcellsobja1oup1.ui <- renderUI({ 
    plotOutput("mastcellsobja1oup1", height = pList[input$mastcellsobja1psz]) 
  }) 
  output$mastcellsobja1oup1.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja1drX,"_",input$mastcellsobja1drY,"_",  
                                   input$mastcellsobja1inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja1oup1.h, width = input$mastcellsobja1oup1.w, useDingbats = FALSE, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp1,   
                      input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
                      input$mastcellsobja1siz, input$mastcellsobja1col1, input$mastcellsobja1ord1,  
                      input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt, input$mastcellsobja1lab1) ) 
  }) 
  output$mastcellsobja1oup1.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja1drX,"_",input$mastcellsobja1drY,"_",  
                                   input$mastcellsobja1inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja1oup1.h, width = input$mastcellsobja1oup1.w, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp1,   
                      input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
                      input$mastcellsobja1siz, input$mastcellsobja1col1, input$mastcellsobja1ord1,  
                      input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt, input$mastcellsobja1lab1) ) 
  }) 
  output$mastcellsobja1.dt <- renderDataTable({ 
    ggData = scDRnum(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1inp1, input$mastcellsobja1inp2, 
                     input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
                     "mastcellsobjgexpr.h5", mastcellsobjgene, input$mastcellsobja1splt) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("pctExpress"), digits = 2) 
  }) 
   
  output$mastcellsobja1oup2 <- renderPlot({ 
    scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp2,  
             input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
             "mastcellsobjgexpr.h5", mastcellsobjgene, 
             input$mastcellsobja1siz, input$mastcellsobja1col2, input$mastcellsobja1ord2, 
             input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt) 
  }) 
  output$mastcellsobja1oup2.ui <- renderUI({ 
    plotOutput("mastcellsobja1oup2", height = pList[input$mastcellsobja1psz]) 
  }) 
  output$mastcellsobja1oup2.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja1drX,"_",input$mastcellsobja1drY,"_",  
                                   input$mastcellsobja1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja1oup2.h, width = input$mastcellsobja1oup2.w, useDingbats = FALSE, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp2,  
                      input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja1siz, input$mastcellsobja1col2, input$mastcellsobja1ord2, 
                      input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt) ) 
  }) 
  output$mastcellsobja1oup2.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja1drX,"_",input$mastcellsobja1drY,"_",  
                                   input$mastcellsobja1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja1oup2.h, width = input$mastcellsobja1oup2.w, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja1drX, input$mastcellsobja1drY, input$mastcellsobja1inp2,  
                      input$mastcellsobja1sub1, input$mastcellsobja1sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja1siz, input$mastcellsobja1col2, input$mastcellsobja1ord2, 
                      input$mastcellsobja1fsz, input$mastcellsobja1asp, input$mastcellsobja1txt) ) 
  }) 
   
   
  ### Plots for tab a2 
  output$mastcellsobja2sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobja2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobja2sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobja2sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobja2oup1 <- renderPlot({ 
    scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp1,  
             input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
             input$mastcellsobja2siz, input$mastcellsobja2col1, input$mastcellsobja2ord1, 
             input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab1) 
  }) 
  output$mastcellsobja2oup1.ui <- renderUI({ 
    plotOutput("mastcellsobja2oup1", height = pList[input$mastcellsobja2psz]) 
  }) 
  output$mastcellsobja2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja2drX,"_",input$mastcellsobja2drY,"_",  
                                   input$mastcellsobja2inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja2oup1.h, width = input$mastcellsobja2oup1.w, useDingbats = FALSE, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp1,   
                      input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
                      input$mastcellsobja2siz, input$mastcellsobja2col1, input$mastcellsobja2ord1,  
                      input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab1) ) 
  }) 
  output$mastcellsobja2oup1.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja2drX,"_",input$mastcellsobja2drY,"_",  
                                   input$mastcellsobja2inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja2oup1.h, width = input$mastcellsobja2oup1.w, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp1,   
                      input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
                      input$mastcellsobja2siz, input$mastcellsobja2col1, input$mastcellsobja2ord1,  
                      input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab1) ) 
  }) 
   
  output$mastcellsobja2oup2 <- renderPlot({ 
    scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp2,  
             input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
             input$mastcellsobja2siz, input$mastcellsobja2col2, input$mastcellsobja2ord2, 
             input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab2) 
  }) 
  output$mastcellsobja2oup2.ui <- renderUI({ 
    plotOutput("mastcellsobja2oup2", height = pList[input$mastcellsobja2psz]) 
  }) 
  output$mastcellsobja2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja2drX,"_",input$mastcellsobja2drY,"_",  
                                   input$mastcellsobja2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja2oup2.h, width = input$mastcellsobja2oup2.w, useDingbats = FALSE, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp2,   
                      input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
                      input$mastcellsobja2siz, input$mastcellsobja2col2, input$mastcellsobja2ord2,  
                      input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab2) ) 
  }) 
  output$mastcellsobja2oup2.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja2drX,"_",input$mastcellsobja2drY,"_",  
                                   input$mastcellsobja2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja2oup2.h, width = input$mastcellsobja2oup2.w, 
      plot = scDRcell(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja2drX, input$mastcellsobja2drY, input$mastcellsobja2inp2,   
                      input$mastcellsobja2sub1, input$mastcellsobja2sub2, 
                      input$mastcellsobja2siz, input$mastcellsobja2col2, input$mastcellsobja2ord2,  
                      input$mastcellsobja2fsz, input$mastcellsobja2asp, input$mastcellsobja2txt, input$mastcellsobja2lab2) ) 
  }) 
   
   
  ### Plots for tab a3 
  output$mastcellsobja3sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja3sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobja3sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobja3sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobja3sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobja3sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobja3sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobja3oup1 <- renderPlot({ 
    scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp1,  
             input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
             "mastcellsobjgexpr.h5", mastcellsobjgene, 
             input$mastcellsobja3siz, input$mastcellsobja3col1, input$mastcellsobja3ord1, 
             input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) 
  }) 
  output$mastcellsobja3oup1.ui <- renderUI({ 
    plotOutput("mastcellsobja3oup1", height = pList[input$mastcellsobja3psz]) 
  }) 
  output$mastcellsobja3oup1.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja3drX,"_",input$mastcellsobja3drY,"_",  
                                   input$mastcellsobja3inp1,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja3oup1.h, width = input$mastcellsobja3oup1.w, useDingbats = FALSE, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp1,  
                      input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja3siz, input$mastcellsobja3col1, input$mastcellsobja3ord1, 
                      input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) ) 
  }) 
  output$mastcellsobja3oup1.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja3drX,"_",input$mastcellsobja3drY,"_",  
                                   input$mastcellsobja3inp1,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja3oup1.h, width = input$mastcellsobja3oup1.w, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp1,  
                      input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja3siz, input$mastcellsobja3col1, input$mastcellsobja3ord1, 
                      input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) ) 
  }) 
   
  output$mastcellsobja3oup2 <- renderPlot({ 
    scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp2,  
             input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
             "mastcellsobjgexpr.h5", mastcellsobjgene, 
             input$mastcellsobja3siz, input$mastcellsobja3col2, input$mastcellsobja3ord2, 
             input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) 
  }) 
  output$mastcellsobja3oup2.ui <- renderUI({ 
    plotOutput("mastcellsobja3oup2", height = pList[input$mastcellsobja3psz]) 
  }) 
  output$mastcellsobja3oup2.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja3drX,"_",input$mastcellsobja3drY,"_",  
                                   input$mastcellsobja3inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobja3oup2.h, width = input$mastcellsobja3oup2.w, useDingbats = FALSE, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp2,  
                      input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja3siz, input$mastcellsobja3col2, input$mastcellsobja3ord2, 
                      input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) ) 
  }) 
  output$mastcellsobja3oup2.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobja3drX,"_",input$mastcellsobja3drY,"_",  
                                   input$mastcellsobja3inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobja3oup2.h, width = input$mastcellsobja3oup2.w, 
      plot = scDRgene(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobja3drX, input$mastcellsobja3drY, input$mastcellsobja3inp2,  
                      input$mastcellsobja3sub1, input$mastcellsobja3sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobja3siz, input$mastcellsobja3col2, input$mastcellsobja3ord2, 
                      input$mastcellsobja3fsz, input$mastcellsobja3asp, input$mastcellsobja3txt) ) 
  }) 
     
   
  ### Plots for tab b2 
  output$mastcellsobjb2sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjb2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobjb2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobjb2sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobjb2sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjb2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjb2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobjb2oup1 <- renderPlot({ 
    scDRcoex(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjb2drX, input$mastcellsobjb2drY,   
             input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2sub1, input$mastcellsobjb2sub2, 
             "mastcellsobjgexpr.h5", mastcellsobjgene, 
             input$mastcellsobjb2siz, input$mastcellsobjb2col1, input$mastcellsobjb2ord1, 
             input$mastcellsobjb2fsz, input$mastcellsobjb2asp, input$mastcellsobjb2txt) 
  }) 
  output$mastcellsobjb2oup1.ui <- renderUI({ 
    plotOutput("mastcellsobjb2oup1", height = pList2[input$mastcellsobjb2psz]) 
  }) 
  output$mastcellsobjb2oup1.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjb2drX,"_",input$mastcellsobjb2drY,"_",  
                                    input$mastcellsobjb2inp1,"_",input$mastcellsobjb2inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobjb2oup1.h, width = input$mastcellsobjb2oup1.w, useDingbats = FALSE, 
      plot = scDRcoex(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjb2drX, input$mastcellsobjb2drY,  
                      input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2sub1, input$mastcellsobjb2sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobjb2siz, input$mastcellsobjb2col1, input$mastcellsobjb2ord1, 
                      input$mastcellsobjb2fsz, input$mastcellsobjb2asp, input$mastcellsobjb2txt) ) 
  }) 
  output$mastcellsobjb2oup1.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjb2drX,"_",input$mastcellsobjb2drY,"_",  
                                    input$mastcellsobjb2inp1,"_",input$mastcellsobjb2inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobjb2oup1.h, width = input$mastcellsobjb2oup1.w, 
      plot = scDRcoex(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjb2drX, input$mastcellsobjb2drY,  
                      input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2sub1, input$mastcellsobjb2sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, 
                      input$mastcellsobjb2siz, input$mastcellsobjb2col1, input$mastcellsobjb2ord1, 
                      input$mastcellsobjb2fsz, input$mastcellsobjb2asp, input$mastcellsobjb2txt) ) 
  }) 
  output$mastcellsobjb2oup2 <- renderPlot({ 
    scDRcoexLeg(input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2col1, input$mastcellsobjb2fsz) 
  }) 
  output$mastcellsobjb2oup2.ui <- renderUI({ 
    plotOutput("mastcellsobjb2oup2", height = "300px") 
  }) 
  output$mastcellsobjb2oup2.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjb2drX,"_",input$mastcellsobjb2drY,"_",  
                                    input$mastcellsobjb2inp1,"_",input$mastcellsobjb2inp2,"_leg.pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
      plot = scDRcoexLeg(input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2col1, input$mastcellsobjb2fsz) ) 
  }) 
  output$mastcellsobjb2oup2.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjb2drX,"_",input$mastcellsobjb2drY,"_",  
                                    input$mastcellsobjb2inp1,"_",input$mastcellsobjb2inp2,"_leg.png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = 3, width = 4, 
      plot = scDRcoexLeg(input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, input$mastcellsobjb2col1, input$mastcellsobjb2fsz) ) 
  }) 
  output$mastcellsobjb2.dt <- renderDataTable({ 
    ggData = scDRcoexNum(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjb2inp1, input$mastcellsobjb2inp2, 
                         input$mastcellsobjb2sub1, input$mastcellsobjb2sub2, "mastcellsobjgexpr.h5", mastcellsobjgene) 
    datatable(ggData, rownames = FALSE, extensions = "Buttons", 
              options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
      formatRound(columns = c("percent"), digits = 2) 
  }) 
     
   
  ### Plots for tab c1 
  output$mastcellsobjc1sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobjc1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobjc1sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobjc1sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjc1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobjc1oup <- renderPlot({ 
    scVioBox(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc1inp1, input$mastcellsobjc1inp2, 
             input$mastcellsobjc1sub1, input$mastcellsobjc1sub2, 
             "mastcellsobjgexpr.h5", mastcellsobjgene, input$mastcellsobjc1typ, input$mastcellsobjc1pts, 
             input$mastcellsobjc1siz, input$mastcellsobjc1fsz) 
  }) 
  output$mastcellsobjc1oup.ui <- renderUI({ 
    plotOutput("mastcellsobjc1oup", height = pList2[input$mastcellsobjc1psz]) 
  }) 
  output$mastcellsobjc1oup.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjc1typ,"_",input$mastcellsobjc1inp1,"_",  
                                   input$mastcellsobjc1inp2,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobjc1oup.h, width = input$mastcellsobjc1oup.w, useDingbats = FALSE, 
      plot = scVioBox(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc1inp1, input$mastcellsobjc1inp2, 
                      input$mastcellsobjc1sub1, input$mastcellsobjc1sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, input$mastcellsobjc1typ, input$mastcellsobjc1pts, 
                      input$mastcellsobjc1siz, input$mastcellsobjc1fsz) ) 
  }) 
  output$mastcellsobjc1oup.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjc1typ,"_",input$mastcellsobjc1inp1,"_",  
                                   input$mastcellsobjc1inp2,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobjc1oup.h, width = input$mastcellsobjc1oup.w, 
      plot = scVioBox(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc1inp1, input$mastcellsobjc1inp2, 
                      input$mastcellsobjc1sub1, input$mastcellsobjc1sub2, 
                      "mastcellsobjgexpr.h5", mastcellsobjgene, input$mastcellsobjc1typ, input$mastcellsobjc1pts, 
                      input$mastcellsobjc1siz, input$mastcellsobjc1fsz) ) 
  }) 
     
   
### Plots for tab c2 
  output$mastcellsobjc2sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc2sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobjc2sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobjc2sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobjc2sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjc2sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjc2sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
output$mastcellsobjc2oup <- renderPlot({ 
  scProp(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc2inp1, input$mastcellsobjc2inp2,  
         input$mastcellsobjc2sub1, input$mastcellsobjc2sub2, 
         input$mastcellsobjc2typ, input$mastcellsobjc2flp, input$mastcellsobjc2fsz) 
}) 
output$mastcellsobjc2oup.ui <- renderUI({ 
  plotOutput("mastcellsobjc2oup", height = pList2[input$mastcellsobjc2psz]) 
}) 
output$mastcellsobjc2oup.pdf <- downloadHandler( 
  filename = function() { paste0("mastcellsobj",input$mastcellsobjc2typ,"_",input$mastcellsobjc2inp1,"_",  
                                 input$mastcellsobjc2inp2,".pdf") }, 
  content = function(file) { ggsave( 
    file, device = "pdf", height = input$mastcellsobjc2oup.h, width = input$mastcellsobjc2oup.w, useDingbats = FALSE, 
    plot = scProp(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc2inp1, input$mastcellsobjc2inp2,  
                  input$mastcellsobjc2sub1, input$mastcellsobjc2sub2, 
                  input$mastcellsobjc2typ, input$mastcellsobjc2flp, input$mastcellsobjc2fsz) ) 
  }) 
output$mastcellsobjc2oup.png <- downloadHandler( 
  filename = function() { paste0("mastcellsobj",input$mastcellsobjc2typ,"_",input$mastcellsobjc2inp1,"_",  
                                 input$mastcellsobjc2inp2,".png") }, 
  content = function(file) { ggsave( 
    file, device = "png", height = input$mastcellsobjc2oup.h, width = input$mastcellsobjc2oup.w, 
    plot = scProp(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjc2inp1, input$mastcellsobjc2inp2,  
                  input$mastcellsobjc2sub1, input$mastcellsobjc2sub2, 
                  input$mastcellsobjc2typ, input$mastcellsobjc2flp, input$mastcellsobjc2fsz) ) 
  }) 
     
   
  ### Plots for tab d1 
  output$mastcellsobjd1sub1.ui <- renderUI({ 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjd1sub1]$fID, "\\|")[[1]] 
    checkboxGroupInput("mastcellsobjd1sub2", "Select which cells to show", inline = TRUE, 
                       choices = sub, selected = sub) 
  }) 
  observeEvent(input$mastcellsobjd1sub1non, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = NULL, inline = TRUE) 
  }) 
  observeEvent(input$mastcellsobjd1sub1all, { 
    sub = strsplit(mastcellsobjconf[UI == input$mastcellsobjd1sub1]$fID, "\\|")[[1]] 
    updateCheckboxGroupInput(session, inputId = "mastcellsobjd1sub2", label = "Select which cells to show", 
                             choices = sub, selected = sub, inline = TRUE) 
  }) 
  output$mastcellsobjd1oupTxt <- renderUI({ 
    geneList = scGeneList(input$mastcellsobjd1inp, mastcellsobjgene) 
    if(nrow(geneList) > 50){ 
      HTML("More than 50 input genes! Please reduce the gene list!") 
    } else { 
      oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
      if(nrow(geneList[present == FALSE]) > 0){ 
        oup = paste0(oup, "<br/>", 
                     nrow(geneList[present == FALSE]), " genes not found (", 
                     paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
      } 
      HTML(oup) 
    } 
  }) 
  output$mastcellsobjd1oup <- renderPlot({ 
    scBubbHeat(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjd1inp, input$mastcellsobjd1grp, input$mastcellsobjd1plt, 
               input$mastcellsobjd1sub1, input$mastcellsobjd1sub2, "mastcellsobjgexpr.h5", mastcellsobjgene, 
               input$mastcellsobjd1scl, input$mastcellsobjd1row, input$mastcellsobjd1col, 
               input$mastcellsobjd1cols, input$mastcellsobjd1fsz) 
  }) 
  output$mastcellsobjd1oup.ui <- renderUI({ 
    plotOutput("mastcellsobjd1oup", height = pList3[input$mastcellsobjd1psz]) 
  }) 
  output$mastcellsobjd1oup.pdf <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjd1plt,"_",input$mastcellsobjd1grp,".pdf") }, 
    content = function(file) { ggsave( 
      file, device = "pdf", height = input$mastcellsobjd1oup.h, width = input$mastcellsobjd1oup.w, 
      plot = scBubbHeat(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjd1inp, input$mastcellsobjd1grp, input$mastcellsobjd1plt, 
                        input$mastcellsobjd1sub1, input$mastcellsobjd1sub2, "mastcellsobjgexpr.h5", mastcellsobjgene, 
                        input$mastcellsobjd1scl, input$mastcellsobjd1row, input$mastcellsobjd1col, 
                        input$mastcellsobjd1cols, input$mastcellsobjd1fsz, save = TRUE) ) 
  }) 
  output$mastcellsobjd1oup.png <- downloadHandler( 
    filename = function() { paste0("mastcellsobj",input$mastcellsobjd1plt,"_",input$mastcellsobjd1grp,".png") }, 
    content = function(file) { ggsave( 
      file, device = "png", height = input$mastcellsobjd1oup.h, width = input$mastcellsobjd1oup.w, 
      plot = scBubbHeat(mastcellsobjconf, mastcellsobjmeta, input$mastcellsobjd1inp, input$mastcellsobjd1grp, input$mastcellsobjd1plt, 
                        input$mastcellsobjd1sub1, input$mastcellsobjd1sub2, "mastcellsobjgexpr.h5", mastcellsobjgene, 
                        input$mastcellsobjd1scl, input$mastcellsobjd1row, input$mastcellsobjd1col, 
                        input$mastcellsobjd1cols, input$mastcellsobjd1fsz, save = TRUE) ) 
  }) 
   
   
      
}) 
 
 
 
 