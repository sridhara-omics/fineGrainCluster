library(Seurat)

#' Returns Seurat object with the results of all the analysis steps
#'
#' @param data  input data is a matrix of gene expression values
#'
#' @return
#' @export
#'
#' @examples
generic_seurat_analysis <- function(data, scale.factor=1000, nfeatures=20, num.replicate=100, dims = 1:10, resolution=0.6, min.pct = 0.25, logfc.threshold = 0.25) {

  # Setup the Seurat Object
  sce <- CreateSeuratObject(counts = data, project = "Seurat_Analysis")

  # Standard pre-processing workflow
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = scale.factor)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = nfeatures)

  # Normalizing the data
  sce <- ScaleData(sce)

  # Identification of highly variable features (feature selection)
  sce <- RunPCA(sce, features = VariableFeatures(sce), pcs.print = 1:10)

  # Scaling the data
  sce <- JackStraw(sce, num.replicate = num.replicate)
  # sce <- ScoreJackStraw(sce)

  # Perform linear dimensional reduction
  sce <- FindNeighbors(sce, reduction = "pca", dims = dims)
  sce <- FindClusters(sce, resolution = resolution)

  # Determine the ‘dimensionality’ of the dataset
  sce@meta.data$dimension <- as.factor(unique(sce@meta.data$cluster))

  # Cluster the cells
  sce <- RunUMAP(sce, dims = dims)

  # Run non-linear dimensional reduction (UMAP/tSNE)
  sce <- FindAllMarkers(sce, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)

  # Finding differentially expressed features (cluster biomarkers)
  sce <- CellTypeScore(sce, dims = c(1,2), reduction = "umap")

  # Assigning cell type identity to clusters
  return(sce)
}
