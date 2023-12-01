#' Select markers for SPOTlight using an `screference`
#'
#' @param scref an `screference` object
#' @param n_hvg number of highly variable genes to be evaluated
#' @param marker_selection_method a string with the name of a metric, indicating
#' how to select markers. See `scran::scoreMarkers` for available metrics.
#' @param marker_selection_params a list of parameters to select marker genes.
#' `min` and `max` are currently supported. For example, for
#' `marker_selection_method = "mean.AUC"`, `min = 0.8` indicates markers have
#' minimum mean AUC of 0.8, while `max = Inf` indicates no maximum value of
#' mean AUC.
#' @param downsample_n_cells an integer, how many cells of each type to keep
#' in the reference
#'
#' @return reference results to be used for SPOTlight deconvolution.
#'
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @importFrom scran modelGeneVar getTopHVGs scoreMarkers
#' @export
spotlight_scref <- function(scref,
                            n_hvg = 3000,
                            marker_selection_method = "mean.AUC",
                            marker_selection_params = list(
                                min = 0.8,
                                max = Inf),
                            downsample_n_cells = 100) {
    .install_spotlight()
    # Preprocess --
    sce <- as.SingleCellExperiment(scref$seurat_obj)
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce)
    hvg <- getTopHVGs(dec, n = 3000)

    # Select markers --
    colLabels(sce) <- colData(sce)$annot_id
    genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))
    marker_genes <- scoreMarkers(sce, subset.row = genes)
    marker_genes_df <- lapply(names(marker_genes), function(pop) {
        marker_genes[[pop]] %>%
            as.data.frame() %>%
            rownames_to_column("gene") %>%
            select(gene, method = {{ marker_selection_method }}) %>%
            filter(method >= marker_selection_params$min,
                   method <= marker_selection_params$max) %>%
            arrange(desc(method)) %>%
            mutate(population = pop)
    }) %>% bind_rows()
    colnames(marker_genes_df)[2] <- marker_selection_method

    # Downsample ---
    idx <- split(seq(ncol(sce)), sce$annot_id)
    cs_keep <- lapply(idx, function(i) {
        n <- length(i)
        n_cells <- if_else(n < downsample_n_cells, n, downsample_n_cells)
        sample(i, n_cells)
    })
    sce <- sce[, unlist(cs_keep)]

    spotlight_ref <- list(sce = sce, hvg = hvg, mgs = marker_genes_df,
                          weight_id = marker_selection_method,
                          group_id = "population", gene_id = "gene")

    return(spotlight_ref)
}


#' SPOTlight deconvolution of spatial data using an `screference`
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay
#' @param scref an object of `screference`
#' @param model character string indicating which model to use when running NMF.
#' Either "ns" (default) or "std".
#' @param n_top integer scalar specifying the number of markers to select per
#' group. By default NULL uses all the marker genes to initialize the model.
#' @param min_prop scalar in [0,1] setting the minimum contribution expected from
#' a cell type in x to observations in y. By default 0.
#' @param cache_path a path to cache the results
#'
#' @return a tibble with deconvolution fractions
#'
#' @importFrom SPOTlight SPOTlight
#' @import tidyverse
#' @export
spotlight_deconvolute <- function(spatial_obj, scref, model = "ns",
                                  n_top = NULL, min_prop = 0, cache_path = NULL) {
    .install_spotlight()
    # Checks --
    assert(class(spatial_obj) == "Seurat")
    assert(class(scref) == "screference")
    assert("spotlight" %in% names(scref$cached_results))

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

    # Run SPOTlight --
    spotlight_ref <- scref$cached_results$spotlight

    res <- SPOTlight::SPOTlight(
        x = scref$seurat_obj,
        y = as.SingleCellExperiment(spatial_obj),
        groups = scref$seurat_obj$annot_id,
        mgs = spotlight_ref$mgs,
        gene_id = spotlight_ref$gene_id,
        weight_id = spotlight_ref$weight_id,
        group_id = spotlight_ref$group_id,
        n_top = n_top,
        model = model,
        min_prop = min_prop)

    deconv_res <- res$mat %>%
        as.data.frame() %>%
        dplyr::rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        mutate(method = "SPOTlight") %>%
        tibble() %>%
        relocate(sample)

    # Return --
    if(!is.null(cache_path)) saveRDS(deconv_res, file = res_file)
    return(deconv_res)

}

#' @importFrom remotes install_bioc
.install_spotlight <- function() {
    if(! all(c("scuttle", "scran", "SPOTlight") %in% installed.packages())) {
        message("R package SPOTlight or one of its dependencies not detected. Installing...")
        install_bioc(c("SingleCellExperiment", "scuttle", "scran", "SPOTlight"))
    }
}
