#' @title run_seurat_pipeline
#' @details
#'   this function can simply run base seurat pipeline from Seurat object to
#'   TSNE or UMAP. However, we do not provide data quality control services.
#'   So please do some common quality control on the data first,such as :
#'   nFeature RNA; percent mt; cellcycle; double cells and more.
#' @param object seurat object
#' @param runSCTransform TRUE or FALSE,Deciding whether to standardize the data with sctransform
#' @param runHarmony  TRUE or FALSE,Deciding whether to remove batch use harmony
#' @return a seurat object
#' @export
#' @import Seurat harmony dplyr
#' @examples
#'  \donttest{
#'    pbmc <- readRDS(system.file("data","pbmc.rda",package="scrnaVis"))
#'    pbmc_new <- run_seurat_pipeline(object=pbmc,runSCTransform=TRUE,runHarmony=TRUE)
#'  }
run_seurat_pipeline <- function(object = NULL,runSCTransform = TRUE,runHarmony = TRUE) {
  ## do SCTranfrom
  if (runSCTransform == TRUE) {
    scRNA <- object %>% SCTransform(return.only.var.genes = T) %>%
      RunPCA(assay.use = "SCT")
    # do harmony to remove batch
    if (runHarmony == TRUE) {
      scRNA %<>% RunHarmony(group.by.vars = "orig.ident", assay.use = "SCT") %>%
        RunTSNE(reduction = "harmony", dims = 1:20) %>%
        RunUMAP(reduction = "harmony", dims = 1:20) %>%
        FindNeighbors(reduction = "harmony") %>%
        FindClusters(resolution = seq(0, 1.2, .2))
      # just a sample, we do not need to remove batch
    } else {
        scRNA %<>% RunTSNE(reduction = "pca", dims = 1:20) %>%
          RunUMAP(reduction = "pca", dims = 1:20) %>%
          FindNeighbors(reduction = "pca") %>%
          FindClusters(resolution = seq(0, 1.2, .2))
      }

    ## do normalise
  } else {
      scRNA <- NormalizeData(object,normalization.method = "LogNormalize",scale.factor = 10000,assay = "RNA") %>%
        FindVariableFeatures(selection.method = "vst",nfeatures = 3000,assay = "RNA") 
      scRNA %<>% ScaleData(features = VariableFeatures(scRNA), assay = "RNA") %>%
        RunPCA(assay = "RNA")
    # do harmony to remove batch
      if (runHarmony == TRUE) {
        scRNA %<>% RunHarmony(group.by.vars = "orig.ident") %>%
          RunTSNE(reduction = "harmony", dims = 1:20) %>%
          RunUMAP(reduction = "harmony", dims = 1:20) %>%
          FindNeighbors(reduction = "harmony") %>%
          FindClusters(resolution = seq(0, 1.2, .2))
      # just a sample, we do not need to remove batch
      } else{
        scRNA %<>% RunTSNE(reduction = "pca", dims = 1:20) %>%
          RunUMAP(reduction = "pca", dims = 1:20) %>%
          FindNeighbors(reduction = "pca") %>%
          FindClusters(resolution = seq(0, 1.2, .2))
        }
    }
    return(scRNA)
  }


