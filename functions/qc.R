############################################################
## qc.R
############################################################

# AddQC metric adds percent.mt and percent.ribo to seurat based on genes found 
# with mitoPattern and riboPattern
# zebrafish mitoPattern = "^mt-", riboPattern = "^rp[ls]"
# mouse mitoPattern = , riboPattern = 

AddQCmetrics <- function(
    seurat,
    mitoPattern,
    riboPattern){
  
  seurat[["percent.mt"]] <-
    PercentageFeatureSet(
      seurat,
      pattern=mitoPattern
    )
  
  seurat[["percent.ribo"]] <-
    PercentageFeatureSet(
      seurat,
      pattern=riboPattern
    )
  
  return(seurat)
}

# PlotQC returns my standard initial QC metric plots which are used for visualizing quality
# of the data and determining cell QC filtering

PlotQC <- function(
    seurat){
  v <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)
  f1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
  f2 <- FeatureScatter(seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt")
  f3 <- FeatureScatter(seurat, feature1 = "nFeature_RNA", feature2 = "percent.ribo")
  f4 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  qc1 <- ggplot(seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + geom_point(size=0.1) +
    scale_color_gradient(low="blue",high="red") + theme_classic()
  
  qc2 <- ggplot(seurat@meta.data, aes(x=nCount_RNA, y=percent.mt)) + geom_point(size=0.1) + scale_x_log10() +
    geom_density2d()
  
  qc3 <- ggplot(seurat@meta.data, aes(x=nCount_RNA, y=nFeature_RNA ,color=percent.mt)) + geom_point(size=0.1) +
    scale_x_log10() + scale_y_log10(breaks = c(100, 500, 750, 1000, 1500, 2000, 3000, 10000)) +  geom_density2d() +
    scale_color_gradient(low="gray",high="darkblue",limits = c(0,20)) + theme_classic()
  
  print(v)
  print(f1)
  print(f2)
  print(f3)
  print(f4)
  print(qc1)
  print(qc2)
  print(qc3)
}

# FilterCells filters the SeuratObject based on maximum percent mitochondrial and 
# minimum nFeatures (number of genes)

FilterCells <- function(seurat, percentMT, nFeature){
  cells.orig <- table(seurat@meta.data$orig.ident)[[1]]
  seurat <- subset(x = seurat,
                   subset = percent.mt < percentMT)
  cells.mt <- table(seurat@meta.data$orig.ident)[[1]]
  seurat <- subset(x = seurat,
                   subset = nFeature_RNA > nFeature)
  cells.filterlowGenes <- table(seurat@meta.data$orig.ident)[[1]]
  
  print(unique(seurat$orig.ident))
  print(paste("cells in original sample: ", cells.orig))
  print(paste("cells after percent.mt <", percentMT,":", cells.mt))
  print(paste("cells after genes > ", nFeature, cells.filterlowGenes))
  
  return(seurat)
}

