#' @import reticulate
.setup_deconv_conda <- function() {
    envs <- conda_list()
    if(!"deconverse" %in% envs$name) {
        message("Creating `deconverse` conda environment...")
        conda_create("deconverse",
                     python_version = "3.8",
                     packages = c("numpy", "numexpr>2.7.3", "pandas",
                                  "scipy", "scikit-learn", "matplotlib",
                                  "anndata>0.8", "scanpy"), pip = TRUE, forge = FALSE)
        conda_install("deconverse", "anndata>0.8", pip = TRUE, forge = FALSE)
    } else {
        message("Using `deconverse` conda environment...")
    }
    use_miniconda("deconverse", required = TRUE)
}

#' Read a h5ad file and convert it to a Seurat object
#'
#' This function reads a h5ad file, extracts metadata, gene and cell names, and count data.
#' Then, it creates a Seurat object with the extracted information.
#'
#' Note that this version doesn't import unstructured data saved in `uns`.
#'
#' @param file A character string representing the path to the h5ad file.
#'
#' @return A Seurat object containing the data from the h5ad file.
#' @import Seurat SeuratObject
#' @importFrom dplyr relocate
#' @importFrom reticulate import
#' @importFrom stringr str_detect
#'
#' @export
read_h5ad <- function(file) {
    warning("This version doesn't import unstructured data saved in `uns`.")
    #-- Use conda h5 environment
    .setup_deconv_conda()
    #-- Get elements of object
    ad <- import("anndata")
    if(str_detect(file, "h5ad$")) {
        adata <- ad$read_h5ad(file)
        adata$obs_names_make_unique()
    } else {
        stop("`file` is not of format `h5ad`.")
    }

    #-- Get meta
    message("Getting metadata...")
    cell_meta <- adata$obs
    cell_meta$barcode <- rownames(cell_meta)
    cell_meta <- relocate(cell_meta, "barcode")
    gene_meta <- adata$var

    cell_names <- rownames(cell_meta)
    gene_names <- rownames(gene_meta)

    #-- Get data
    message("Getting counts...")
    layers <- adata$layers$as_dict()
    layers <- lapply(layers, function(mat) {
        # mat <- py_to_r(mat)
        mat <- t(mat)
        rownames(mat) <- gene_names
        colnames(mat) <- cell_names
        return(mat)
    })

    message("Creating Seurat object...")
    seu <- CreateSeuratObject(counts = layers$counts)
    seu@meta.data <- left_df_join(seu@meta.data, cell_meta)
    seu@assays$RNA@meta.features <- gene_meta
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    gc()
    return(seu)
}
