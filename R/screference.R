# screference ------------------

#' Generate single-cell reference (`screference`) object
#'
#' @param seurat_obj a reference Seurat object
#' @param annot_id the column name of the Seurat object metadata containing
#' the reference annotation
#' @param batch_id NULL or a column name of the Seurat object metadata containing
#' information about batch to be used as batch information downstream. Examples
#' of batch are patient ID, study ID, etc.
#' @param project_name a string indicating the project name, used for caching
#' @param cache_path path to the directory where results will be cached.
#' @param sample_cells NULL or a list with how many cells to randomly sub-sample of each
#' population. Sub-sampling can make the object significantly smaller and accelerate
#' downstream computation, but make the results less accurate.
#' @param seed an integer representing a "seed", for reproducibility of sub-sampling.
#'
#' @import tidyverse
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
    colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == annot_id] <- "annot_id"
    meta <- seurat_obj@meta.data
    populations <- unique(meta$annot_id)

    #-- Get random cells
    if(!is.null(sample_cells)) {
        message("-- Subsampling cells...")
        assert(all(names(sample_cells) %in% unique(meta[,annot_id])))
        if(!is.null(seed)) set.seed(seed)
        sampled_cells <- lapply(names(sample_cells), function(pop) {
            all_pop_cells <- meta %>%
                filter(annot_id == pop) %>%
                rownames()
            assert(length(all_pop_cells) > sample_cells[pop])
            sampled_cells <- sample(all_pop_cells, size = sample_cells[pop], replace = FALSE)
        }) %>% reduce(c)
        other_pops <- setdiff(populations, names(sample_cells))
        if(length(other_pops) > 0) {
            sampled_cells <- c(sampled_cells, rownames(meta)[meta[,"annot_id"] %in% other_pops])
        }
        seurat_obj <- seurat_obj[,sampled_cells]
    }
    if(!is.null(batch_id)) {
        colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) == batch_id] <- "batch_id"
    }

    #-- Make new class
    scref <- list(seurat_obj = seurat_obj)
    class(scref) <- "screference"
    scref$project_name <- project_name
    scref$populations <- populations
    scref$cache <- cache_path
    scref$metrics <- tibble(method = character(), level = character(),
                            type = character(),
                            time_elapsed_s = numeric(),
                            peak_ram_mib = numeric())
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
#' @importFrom peakRAM peakRAM
#' @return an object of class `screference`
#' @rdname compute_reference
#'
#' @export
compute_reference.screference <- function(scref,
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
    method_cache <- filePath(cache_path, method)

    #-- Compute reference
    if(method == "cibersortx") {
        message("CIBERSORTx: Building reference matrix...")
        ram_use <- peakRAM(
            out_path <- cibersortx_scref(scref, cache_path = method_cache, ...))
        reference_res <- cx_ref
    } else if (method == "dwls" | method == "ols" | method == "svr") {
        message("DWLS/OLS/SVR: Building reference matrix using Seurat...")
        ram_use <- peakRAM(
            reference_res <- dwls_scref(scref, cache_path = method_cache, ...))
    } else if(method == "bayesprism") {
        message("BayesPrism: Building reference matrix...")
        reference_res <- bayesprism_scref(scref, cache_path = method_cache, ...)
    } else if(method == "autogenes") {
        message("AutoGeneS: Building reference centroids...")
        reference_res <- autogenes_scref(scref, ...)
    } else if (method == "scaden") {
        message("scaden: Building reference model...")
        message("--- Note: this method might need to be retrained on features shared with the target bulk matrix")
        ram_use <- peakRAM(
            reference_res <- scaden_scref(scref, cache_path = method_cache, ...))
    } else {
        message("No need to compute reference beforehand. Run `deconvolute` with the `scref` object directly.")
    }
    saveRDS(reference_res, file = filePath(method_cache, "reference_res.RDS"))
    scref[["cached_results"]][[method]] <- reference_res
    return(scref)
}

# hscreference ----------------------------------
#' Generate multi-level hierarchical single-cell reference (`hscreference`) object
#'
#' @param seurat_obj a reference Seurat object
#' @param annot_ids in order, column names of the Seurat object metadata from
#' coarsest to finest grained annotation
#' @param batch_id NULL or a column name of the Seurat object metadata containing
#' information about batch to be used as batch information downstream. Examples
#' of batch are patient ID, study ID, etc.
#' @param project_name a string indicating the project name, used for caching
#' @param cache_path path to the directory where results will be cached
#' @param sample_cells NULL or a list with how many cells to randomly sub-sample of each
#' population in the coarsest grained annotation.
#' Sub-sampling can make the object significantly smaller and accelerate
#' downstream computation, but make the results less accurate.
#' @param seed an integer representing a "seed", for reproducibility of sub-sampling.
#'
#' @return an object of class `hscreference`
#' @importFrom data.tree as.Node
#' @import tidyverse
#' @export
new_hscreference <- function(
        seurat_obj,
        annot_ids,
        batch_id = NULL,
        project_name = "project",
        cache_path = "scref_cache",
        sample_cells = NULL,
        seed = NULL) {
    assert(length(annot_ids) > 1)
    message("Generating new `hscreference` object...")
    n_lvl <- length(annot_ids)

    #-- Subsample cells
    if(!is.null(sample_cells)) {
        meta <- seurat_obj@meta.data
        coarse_annot <- annot_ids[1]
        colnames(meta)[colnames(meta) == coarse_annot] <- "annot_id"
        pops <- unique(meta[,"annot_id"])
        if(any(str_detect(pops, c("\\/|\\-"), negate = TRUE))) {
            stop("Population names should not contain `/` or `-`")
        }

        message("-- Subsampling cells...")
        assert(all(names(sample_cells) %in% pops))

        if(!is.null(seed)) set.seed(seed)
        sampled_cells <- lapply(names(sample_cells), function(pop) {
            all_pop_cells <- meta %>%
                filter(annot_id == pop) %>%
                rownames()
            assert(length(all_pop_cells) > sample_cells[pop])
            sampled_cells <- sample(all_pop_cells, size = sample_cells[pop], replace = FALSE)
        }) %>% reduce(c)
        other_pops <- setdiff(pops, names(sample_cells))
        if(length(other_pops) > 0) {
            sampled_cells <- c(sampled_cells, rownames(meta)[meta[,"annot_id"] %in% other_pops])
        }
        seurat_obj <- seurat_obj[,sampled_cells]
    }

    screfs <- lapply(1:n_lvl, function(i) {
        l_scref <- new_screference(seurat_obj, annot_ids[i], batch_id = batch_id,
                                   project_name = paste0(project_name, "_l", i),
                                   cache_path = cache_path, seed = seed)
    })
    names(screfs) <- paste0("l", 1:n_lvl)

    #-- Get population hierarchy
    annots <- seurat_obj@meta.data[,annot_ids]
    colnames(annots) <- paste0("l", 1:n_lvl)
    hpops <- .get_pop_hierarchy(annots)

    hscref <- list(
        screfs = screfs,
        nlevels = n_lvl,
        project_name = project_name,
        hpop_table = hpops$hpop_table,
        hpop_tree = hpops$hpop_tree
    )
    class(hscref) <- "hscreference"
    return(hscref)
}

#' Compute reference given a deconvolution method from an `h-screference` object
#'
#' @param hscref an `h-screference` object
#' @param method a string, one of `deconvolution_methods()`
#' @param ... passed to deconvolution method's reference computation function.
#' See the wrapper (`(method)_scref`) where `(method)` is the method name for
#' parameters.
#'
#' @return an object of class `hscreference`
#'
#' @import tidyverse
#'
#' @export
#' @rdname compute_reference
compute_reference.hscreference <- function(hscref,
                                            method = deconvolution_methods()[1],
                                            ...) {
    assert(class(hscref) == "hscreference")
    assert(method %in% deconvolution_methods())
    hscref$screfs <- lapply(hscref$screfs, compute_reference, method = method, ...)
    return(hscref)
}

# Generics ------------------
#' @export
compute_reference <- function(x, ...) {
    UseMethod("compute_reference")
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

#' @importFrom data.tree as.Node
.get_pop_hierarchy <- function(annots) {
    for(i in 1:(ncol(annots)-1)) {
        l_coarser <- paste0("l", i)
        l_finer <- paste0("l", i+1)

        table <- table(annots[,l_coarser], annots[,l_finer]) %>%
            as.data.frame() %>%
            filter(Freq > 0) %>%
            select(-Freq) %>%
            arrange(Var1, Var2)

        colnames(table) <- c(l_coarser, l_finer)
        if (i == 1) {
            hpop_table <- table
        } else {
            hpop_table <- full_join(hpop_table, table)
        }
    }
    hpop_tree <- hpop_table %>%
        mutate(pathString = apply(hpop_table, 1, function(row) { paste(c("populations", row), collapse = "/") })) %>%
        as.Node()
    hpop_table
    return(list(hpop_tree = hpop_tree, hpop_table = hpop_table))
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

#' @export
print.hscreference <- function(hscref) {
    cat(str_glue("h-screference object named `{hscref$project_name}` with {hscref$nlevels} levels of annotation from {ncol(hscref$screfs[[1]]$seurat_obj)} cells"))
    cat("\n")
    cat("population annotation tree: ")
    cat("\n")
    print(hscref$hpop_tree, limit = 15)
    cat(paste0("cached results: "))
    cat("\n")
    methods <- lapply(hscref$screfs, function(scref) names(scref$cached_results))
    for(i in 1:length(methods)) {
        cat(paste0("l", i, ": "), paste0(methods[[1]], collapse = ", "))
        cat("\n")
    }
}


#' Delete cached data on an `screference` object
#' @param scref an `screference` object
#' @param which a string, which cached result to delete
#'
#' @return screference object
#' @export
delete_cached <- function(scref, which) {
    assert(class(scref) == "screference")
    scref$cached_results[[which]] <- NULL
    return(scref)
}
