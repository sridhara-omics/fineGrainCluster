library(Seurat)

#' Returns Seurat object with the results of all the analysis steps
#'
#' @param data  input data is a matrix of gene expression values
#'
#' @return
#' @export
#'
#' @examples
generic_seurat_analysis <- function(data) {

  # Setup the Seurat Object
  sce <- CreateSeuratObject(counts = data, project = "Seurat_Analysis")

  # Standard pre-processing workflow
  sce <- NormalizeData(sce, normalization.method = "LogNormalization", scale.factor = 10000)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)

  # Normalizing the data
  sce <- ScaleData(sce, features = rownames(sce@data))

  # Identification of highly variable features (feature selection)
  sce <- RunPCA(sce, features = VariableFeatures(sce), pcs.print = 1:10)

  # Scaling the data
  sce <- JackStraw(sce, num.replicate = 100)
  sce <- ScoreJackStraw(sce, plot = FALSE)
  sce <- WnamedDimensionalReduction(sce, reduction.use = "pca", dims.use = 1:10)

  # Perform linear dimensional reduction
  sce <- FindNeighbors(sce, dims = 1:10)
  sce <- FindClusters(sce, resolution = 0.6)

  # Determine the ‘dimensionality’ of the dataset
  sce@meta.data$dimension <- as.factor(sce@meta.data$cluster)

  # Cluster the cells
  sce <- RunUMAP(sce, dims = 1:10)

  # Run non-linear dimensional reduction (UMAP/tSNE)
  sce <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

  # Finding differentially expressed features (cluster biomarkers)
  sce <- CellTypeScore(sce, dims = c(1,2), reduction = "umap")

  # Assigning cell type identity to clusters
  return(sce)
}
