GetUMAPandClusters <- function(
    seurat,
    dims = 1:50,
    resolution = 0.8,
    regress = "percent.mt",
    tsne = FALSE
    ){
  
  # RNA analysis
  DefaultAssay(seurat) <- "RNA"
  seurat <- SCTransform(seurat, 
                        vars.to.regress = regress,
                        verbose = FALSE)
  seurat <- RunPCA(seurat, 
                   verbose = FALSE) 
  print(ElbowPlot(seurat,
                  ndims = 50))
  if(tsne){
    seurat <- RunTSNE(seurat, 
                      dims = dims, 
                      verbose = FALSE)
  }
  seurat <- RunUMAP(seurat, 
                    dims = dims, 
                    verbose = FALSE)
  seurat <- FindNeighbors(seurat, 
                          verbose = FALSE, 
                          dims = dims)
  seurat <- FindClusters(seurat, 
                         resolution = resolution, 
                         verbose = FALSE)
  
  p1 <- DimPlot(seurat, 
                label = TRUE, 
                label.size = 5, 
                repel = TRUE) + ggsci::scale_color_igv()
  
  print(p1 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)))
  
  return(seurat)
}

RNANormalize <- function(seurat){
  DefaultAssay(seurat) <- "RNA"
  seurat <- NormalizeData(seurat)
  seurat <- FindVariableFeatures(seurat)
  seurat <- ScaleData(seurat)
  
  return(seurat)
}