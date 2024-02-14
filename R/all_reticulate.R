.install_reticulate <- function() {
    if(!"reticulate" %in% installed.packages()) {
        message("R package reticulate not detected. Installing...")
        install.packages("reticulate")
    }
}

#' @import reticulate
.setup_deconv_conda <- function() {
    envs <- conda_list()
    if(!"deconverse" %in% envs$name) {
        message("Creating `deconverse` conda environment...")
        conda_create("deconverse",
                     python_version = "3.9",
                     packages = c("numpy", "numexpr>2.7.3", "pandas",
                                  "scipy", "scikit-learn", "matplotlib",
                                  "anndata>0.8", "scanpy"), pip = TRUE, forge = FALSE)
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
#'
#' @importFrom reticulate import
#'
#' @export
read_h5ad <- function(file, counts_layer = "counts") {
    .install_reticulate()
    #-- Use conda deconverse environment
    .setup_deconv_conda()
    #-- Get elements of object
    ad <- import("anndata")
    if(str_detect(file, "h5ad$")) {
        adata <- ad$read_h5ad(file)
    } else if(str_detect(file, "h5ad.zip$")) {
        message("Unzipping, reading, then removing unzipped file...")
        unzip_file <- unzip(file, exdir = "/tmp")
        adata <- ad$read_h5ad(unzip_file)
        file.remove(unzip_file)
    } else {
        stop("`file` is not of format `h5ad`.")
    }
    adata$obs_names_make_unique()

    seu <- adata_to_seurat(adata, counts_layer)
    return(seu)

}


#' Convert adata to a Seurat object
#'
#' This function reads a h5ad file, extracts metadata, gene and cell names, and count data.
#' Then, it creates a Seurat object with the extracted information.
#'
#' Note that this version doesn't import unstructured data saved in `uns`.
#'
#' @param adata An `anndata` python object
#' @param counts_layer layer of adata containing counts information
#'
#' @import Matrix SeuratObject
#' @importFrom Seurat CreateSeuratObject PercentageFeatureSet
#' @importFrom dplyr relocate
#' @importFrom stringr str_detect
#'
#' @return A converted Seurat object
#'
#' @export
adata_to_seurat <- function(adata, counts_layer) {
    .install_reticulate()
    warning("This version doesn't import unstructured data saved in `uns`.")
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
    layers <- dict(adata$layers)
    if(length(layers) == 0) {
        mat <- as.matrix(adata$X)
        mat <- base::t(mat)
        rownames(mat) <- gene_names
        colnames(mat) <- cell_names
        layers <- list(counts = mat)
    } else {
        layers <- lapply(layers, function(mat) {
            # mat <- py_to_r(mat)
            mat <- t(mat)
            rownames(mat) <- gene_names
            colnames(mat) <- cell_names
            return(mat)
        })
    }
    other_layers <- setdiff(names(layers), counts_layer)
    message("Creating Seurat object...")
    seu <- CreateSeuratObject(counts = layers[[counts_layer]])
    added_meta <- setdiff(colnames(cell_meta), colnames(seu@meta.data))
    seu@meta.data <- left_df_join(seu@meta.data, cell_meta[,added_meta])
    if(str_detect(as.character(packageVersion("Seurat")), "^5")) {
        seu@assays$RNA@meta.data <- gene_meta
    } else {
        seu@assays$RNA@meta.features <- gene_meta
    }

    for(ly in other_layers) {
        seu@assays[[ly]] <- layers[[ly]]
    }

    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    gc()
    return(seu)
}
#' Convert Seurat to a adata
#'
#' This function reads a h5ad file, extracts metadata, gene and cell names, and count data.
#' Then, it creates a Seurat object with the extracted information.
#'
#' Note that this version doesn't import unstructured data saved in `uns`.
#'
#' @param seurat_obj A `Seurat` object
#' @param counts_assay name of assay containing counts information. Defaults to
#' first assay found.
#'
#' @importFrom reticulate import
#' @importFrom stringr str_detect
#'
#' @return A converted adata object
#'
#' @export
seurat_to_adata <- function(seurat_obj, counts_assay = 1) {
    .install_reticulate()
    warning("This version doesn't import analysis results not saved as meta data")
    #-- Use conda deconverse environment
    .setup_deconv_conda()

    ad <- import("anndata")
    sp <- import("scipy.sparse")
    pd <- import("pandas")
    message("Converting Seurat to anndata...")
    if(str_detect(as.character(packageVersion("Seurat")), "^5")) {
        adata <- ad$AnnData(X = sp$csr_matrix(t(as.matrix(seurat_obj@assays[[counts_assay]]@counts))),
                            obs = r_to_py(seurat_obj@meta.data),
                            var = r_to_py(seurat_obj@assays[[1]]@meta.data))
    } else {
        adata <- ad$AnnData(X = sp$csr_matrix(t(as.matrix(seurat_obj@assays[[counts_assay]]@counts))),
                            obs = r_to_py(seurat_obj@meta.data),
                            var = r_to_py(seurat_obj@assays[[1]]@meta.features))
    }

    return(adata)

}
