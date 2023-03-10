library(Seurat)

#' Returns Seurat object with the results of all the analysis steps
#'
#' @param data  input data is a matrix of gene expression values
#'
#' @return
#' @export
#'
#' @examples
generic_seurat_analysis <- function(data.dir, scale.factor=1000, nfeatures=2000, num.replicate=100, dims = 1:10, resolution=0.6, min.pct = 0.25, logfc.threshold = 0.25) {

  # Read 10X data
  data <- Read10X(data.dir = data.dir)

  # Create Seurat object
  sce <- CreateSeuratObject(counts = data, project = "Seurat_Analysis")

  # Pre-processing
  sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = scale.factor)
  sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = nfeatures)
  sce <- ScaleData(sce)

  # PCA and feature selection
  sce <- RunPCA(sce, features = VariableFeatures(sce), pcs.print = 1:10)
  sce <- JackStraw(sce, num.replicate = num.replicate)

  # Dimensionality reduction
  sce <- FindNeighbors(sce, reduction = "pca", dims = dims)
  sce <- FindClusters(sce, resolution = resolution)
  # sce@meta.data$dimension <- as.factor(unique(sce@meta.data$cluster))
  sce <- RunUMAP(sce, dims = dims)

  # Finding markers
  # sce@meta.data$all_markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = min.pct, logfc.threshold = logfc.threshold)

  return(sce)
}



