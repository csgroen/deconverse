#' RCTD deconvolution using an `screference`
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay
#' @param scref an object of `screference`
#' @param gene_cutoff minimum normalized gene expression for genes to be included in the platform effect normalization step.
#' @param fc_cutoff minimum log-fold-change (across cell types) for genes to be included in the platform effect normalization step.
#' @param gene_cutoff_reg minimum normalized gene expression for genes to be included in the RCTD step.
#' @param fc_cutoff_reg minimum log-fold-change (across cell types) for genes to be included in the RCTD step.
#' @param UMI_min minimum UMI per pixel included in the analysis
#' @param UMI_max maximum UMI per pixel included in the analysis
#' @param counts_MIN minimum total counts per pixel of genes used in the analysis.
#' @param UMI_min_sigma minimum UMI per pixel for the \link{choose_sigma_c} function
#' @param max_cores for parallel processing, the number of cores used. If set to 1, parallel processing is not used. The system will additionally be checked for
#' number of available cores.
#' @param CELL_MIN_INSTANCE minimum number of cells required per cell type. Default 25, can be lowered if desired.
#' @param MAX_MULTI_TYPES max number of cell types per pixel
#' @param CONFIDENCE_THRESHOLD the minimum change in likelihood (compared to other cell types) necessary to determine a cell type identity with confidence
#' @param DOUBLET_THRESHOLD the penalty weight of predicting a doublet instead of a singlet for a pixel
#' @param cache_path a path to cache the results
#'
#' @return a tibble with deconvolution fractions
#'
#' @import tidyverse
#' @importFrom spacexr SpatialRNA Reference create.RCTD run.RCTD
#' @importFrom R.utils filePath
#'
#' @export
rctd_deconvolute <- function(spatial_obj, scref, ncores = 4, gene_cutoff = 0.000125, fc_cutoff = 0.5,
                             gene_cutoff_reg = 2e-4, fc_cutoff_reg = 0.75, UMI_min = 100,
                             UMI_max = 2e+7, counts_MIN = 10, UMI_min_sigma = 300, CELL_MIN_INSTANCE = 25,
                             MAX_MULTI_TYPES = 4, CONFIDENCE_THRESHOLD = 5, DOUBLET_THRESHOLD = 20,
                             cache_path = "rctd") {
    .install_rctd()
    assert(class(spatial_obj) == "Seurat")
    assert(class(scref) == "screference")

    # Cache --
    if(!is.null(cache_path)) {
        dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
        res_file <- filePath(cache_path, "deconverse_results.RDS")
    }

    if (file.exists(res_file)) {
        message("Results found in cache. Returning...")
        deconv_res <- readRDS(res_file)
        return(deconv_res)
    }

    query <- SpatialRNA(coords = spatial_obj@images[[1]]@coordinates[,2:3],
                        counts = spatial_obj@assays$Spatial@counts)
    reference <- Reference(counts = scref$seurat_obj@assays$RNA@counts,
                           cell_types = structure(factor(scref$seurat_obj@meta.data$annot_id),
                                                  names = colnames(scref$seurat_obj)))
    message("--- RCTD: creating object...")
    RCTD <- create.RCTD(query, reference, max_cores = ncores,
                        gene_cutoff = gene_cutoff, fc_cutoff = fc_cutoff,
                        gene_cutoff_reg = gene_cutoff_reg,
                        fc_cutoff_reg = fc_cutoff_reg, UMI_min = UMI_min,
                        UMI_max = UMI_max, counts_MIN = counts_MIN,
                        UMI_min_sigma = UMI_min_sigma,
                        CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
                        MAX_MULTI_TYPES = MAX_MULTI_TYPES,
                        CONFIDENCE_THRESHOLD = CONFIDENCE_THRESHOLD,
                        DOUBLET_THRESHOLD = DOUBLET_THRESHOLD)
    message("--- RCTD: running deconvolution...")
    RCTD <- run.RCTD(RCTD)

    deconv_res <- RCTD@results$weights / rowSums(RCTD@results$weights)
    deconv_res <- deconv_res %>%
        as.data.frame() %>%
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        mutate(method = "RCTD")

    if(!is.null(cache_path)) saveRDS(deconv_res, file = res_file)

    return(deconv_res)
}

#' @importFrom remotes install_github
.install_rctd <- function() {
    if(!"spacexr" %in% installed.packages()) {
        message("R package spacexr not detected. Installing...")
        install_github("dmcable/spacexr")
    }
}
