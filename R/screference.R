# screference ------------------
#' Generate `screference` object
#'
#' @param seurat_obj a reference Seurat object
#' @param annot_id the column name of the Seurat object metadata containing
#' the reference annotation
#' @param batch_id NULL or a column name of the Seurat object metadata containing
#' information about batch to be used as batch information downstream. Examples
#' of batch are patient ID, study ID, etc.
#' @param project_name a string indicating the project name, used for caching
#' @param cache_path path to the directory where results will be cached
#' @param sample_cells NULL or a list with how many cells to randomly sub-sample of each
#' population. Sub-sampling can make the object significantly smaller and accelerate
#' downstream computation, but make the results less accurate.
#' @param seed an integer representing a "seed", for reproducibility of sub-sampling.
#'
#' @return an object of class `screference`
#'
#' @export
new_screference <- function(
        seurat_obj,
        annot_id,
        batch_id = NULL,
        project_name = "project",
        cache_path = "scref_cache",
        sample_cells = NULL,
        seed = NULL) {
    message("Generating new `screference` object...")
    #-- Check populations
    meta <- seurat_obj@meta.data
    colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == annot_id] <- "annot_id"

    #-- Get random cells
    if(is.null(sample_cells)) {
        scref <- seurat_obj
    } else {
        message("-- Subsampling cells...")
        assert(all(names(sample_cells) %in% unique(meta[,annot_id])))
        if(!is.null(seed)) set.seed(seed)
        sampled_cells <- lapply(names(sample_cells), function(pop) {
            all_pop_cells <- seurat_obj@meta.data %>%
                filter(annot_id == pop) %>%
                rownames()
            assert(length(all_pop_cells) > sample_cells[pop])
            sampled_cells <- sample(all_pop_cells, size = sample_cells[pop], replace = FALSE)
        }) %>% reduce(c)
        scref <- seurat_obj[,sampled_cells]
    }
    if(!is.null(batch_id)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == batch_id] <- "batch_id"
    }

    populations <- unique(scref@meta.data$annot_id)
    #-- Make new class
    scref <- list(seurat_obj = scref)
    class(scref) <- "screference"
    scref$project_name <- project_name
    scref$populations <- populations
    scref$cache <- cache_path
    return(scref)
}

#' Compute reference given a deconvolution method from an `screference` object
#'
#' @param scref an `screference` object
#' @param method a string, one of `deconvolution_methods()`. Only CIBERSORTx
#' and methods implemented in the DWLS package (`dwls`, `svr`, `ols`) need
#' to pre-compute reference.
#' @param ... passed to deconvolution method's reference computation function.
#' See the wrapper (`{method}_scref`) where `{method}` is the method name for
#' parameters.
#'
#' @return an object of class `screference`
#'
#' @export
compute_reference <- function(scref,
                              method = deconvolution_methods()[1],
                              ...) {
    #-- Checks
    assert(class(scref) == "screference")
    assert(method %in% deconvolution_methods())

    #-- Get cached results
    cache_path <- getAbsolutePath(filePath(scref$cache, scref$project_name))
    reference_res <- .scref_cache_check(cache_path, method)

    if(!is.null(reference_res)) {
        scref[["cached_results"]][[method]] <- reference_res
        return(scref)
    }
    #-- Compute reference
    if(method == "cibersortx") {
        message("CIBERSORTx: Building reference matrix...")
        cx_cache <- filePath(cache_path, "cibersortx")
        cx_ref <- filePath(cx_cache, "out_dir/CIBERSORTx_ref_data_inferred_phenoclasses.CIBERSORTx_ref_data_inferred_refsample.bm.K999.txt")
        out_path <- cibersortx_scref(scref, cache_path = cx_cache, ...)
        reference_res <- cx_ref
        saveRDS(reference_res, file = filePath(cx_cache, "reference_res.RDS"))
    } else if (method == "dwls" | method == "ols" | method == "svr") {
        message("DWLS/OLS/SVR: Building reference matrix using Seurat...")
        dwls_cache <- filePath(cache_path, "dwls")
        reference_res <- dwls_scref(scref, cache_path = dwls_cache, ...)
        saveRDS(reference_res, file = filePath(dwls_cache, "reference_res.RDS"))
    } else if(method == "bayesprism") {
        message("BayesPrism: Building reference matrix...")
        bp_cache <- filePath(cache_path, "bayesprism")
        reference_res <- bayesprism_scref(scref, cache_path = bayesprism_cache, ...)
        saveRDS(reference_res, file = filePath(bp_cache, "reference_res.RDS"))
    } else {
        message("No need to compute reference beforehand. Run `deconvolute` with the `scref` object directly.")
    }
    scref[["cached_results"]][[method]] <- reference_res
    return(scref)
}

# Helpers ----------------------------
.scref_cache_check <- function(cache_path, method) {
    cache_method <- filePath(cache_path, method)
    method_fname <- filePath(cache_method, "reference_res.RDS")
    if(file.exists(method_fname)) {
        message("Results found in cache, returning...")
        reference_res <- readRDS(method_fname)
    } else {
        dir.create(cache_method, showWarnings = FALSE, recursive = TRUE)
        reference_res <- NULL
    }
    return(reference_res)
}

#' @export
print.screference <- function(scref) {
    cat(str_glue("screference object named `{scref$project_name}` with {ncol(scref$seurat_obj)} cells and {length(scref$populations)} populations"))
    cat("\n")
    cat(paste0("populations: ", paste(scref$populations, collapse = ", ")))
    cat("\n")
    cat(paste0("cached results: "))
    cat(paste0(names(scref$cached_results), collapse = ", "))
}
