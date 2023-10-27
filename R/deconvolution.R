# Deconvolution ----------------------------------
#' @export
deconvolution_methods <- function() {
    c("OLS" = "ols", "DWLS" = "dwls", "SVR" = "svr", "CIBERSORTx" =  "cibersortx",
      "MuSiC" = "music", "BayesPrism" = "bayesprism", "Bisque" = "bisque",
      "AutoGeneS" = "autogenes", "scaden" = "scaden", "CARD" = "card", "RCTD" = "rctd",
      "SPOTlight" = "spotlight")
}

#' @export
bulk_deconvolution_methods <- function() {
    setdiff(deconvolution_methods(), spatial_only_methods())
}

#' @export
spatial_only_methods <- function() {
    c("CARD" = "card", "RCTD" = "rctd", "SPOTlight" = "spotlight")
}


#' Deconvolute `scbench` object using a chosen method
#'
#' @param scbench an `scbench` object already processed by `pseudobulks`
#' @param scref an `screference` object containing single-cell reference data
#' and pre-computed references for methods that require them.
#' @param method a string with the name of the deconvolution method. For available
#' methods, consult `deconvolution_methods()`
#' @param type type of mixtures to deconvolute. One of: `"population"`,
#' `"spillover"` or `"lod"`
#' @param pseudobulk_norm the normalization method for pseudobulk counts. One of:
#' `"rpm"` for reads per million, `"none"` for raw counts, or `"proportional_fitting"`
#' for using the mean library size of the pseudobulk data to normalize.
#' @param correct_finer normalize the resulting fractions from finer-grained
#' annotation using coarser-grained results
#' @param ... other parameters, passed to the method wrapper, enables the user
#' to change the method parameters. See: `deconvolute_{method}` where `{method}`
#' is the method name in lowercase for method-specific parameters
#'
#' @import tidyverse
#' @importFrom peakRAM peakRAM
#'
#' @return an object of class `scbench`
#'
#' @rdname deconvolute
#' @export
deconvolute.scbench <- function(scbench,
                                scref,
                                method = deconvolution_methods()[1],
                                type = c("population", "spillover", "lod")[1],
                                pseudobulk_norm = c("rpm", "none", "proportional_fitting")[3],
                                correct_finer = TRUE,
                                ...) {
    #-- Error handling
    # assert(class(scbench) == "scbench")
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(!is.null(scbench$pseudobulk_counts[[type]]))
    assert(method %in% bulk_deconvolution_methods())

    #-- Reference picking
    if(class(scref) == "hscreference") {
        #-- Error handling
        assert(scref$nlevels == scbench$nlevels)
        levels <- paste0("l", 1:scbench$nlevels)
        missing_pops <- list()
        same_populations <- rep(FALSE, length(levels))
        for(lv in 1:scref$nlevels) {
            bench_pops <- as.character(unique(scbench$pop_hierarchy[[lv]]))
            ref_pops <- as.character(unique(scref$hpop_table[[lv]]))
            missing_pops[[lv]] <- setdiff(ref_pops, bench_pops)
            same_populations[lv] <- all(bench_pops %in% intersect(bench_pops, ref_pops))
        }
        #-- Checking for missing populations
        if(!all(same_populations)) {
            stop("Not all populations in `scbench` are contained in the reference `scref`.")
        }
        for(i in 1:length(missing_pops)) {
            lv_mp <- missing_pops[[i]]
            lv <- paste0("l", i)
            if(length(lv_mp) > 0) {
                warning("Missing populations in ", lv, ": ", paste(lv_mp, collapse = ", "), "\n")
            }
        }
    } else {
        if(scbench$nlevels > 1) {
            match <- sapply(1:scbench$nlevels, function(i) {
                all(scbench$pop_hierarchy[[i]] %in% scref$populations)
            }) %>% which()
            level_match <- paste0("l", match[length(match)])

            message(str_glue("`scref` contains reference for {level_match}."))
            levels <- level_match
        } else {
            levels <- "l1"
        }
    }
    for(level in levels) {
        if(class(scref) == "hscreference") {
            ref <- scref$screfs[[level]]
        } else {
            ref <- scref
        }
        message("Deconvoluting ", level, " populations in ", type, " analysis...")
        #-- Cache handling
        method_cache <- filePath(scbench$cache, scbench$project_name, level, type, method)
        deconv_res <- .scbench_cache_check(method_cache)
        if(!is.null(deconv_res)) {
            message("Found cached results. Returning...")
            scbench[["deconvolution"]][[level]][[type]][[method]] <- deconv_res
            metrics <- readRDS(filePath(method_cache, "metrics.RDS"))
            scbench$metrics <- bind_rows(scbench$metrics, metrics)
        } else {
            #-- Get pseudocounts and make linear transformation
            if(type == "population") {
                pb <- scbench$pseudobulk_counts[[type]]
            } else {
                pb <- scbench$pseudobulk_counts[[type]][[level]]
            }
            if(pseudobulk_norm == "rpm") {
                data <- librarySizeNormalization(pb, 10^6)
            } else if (pseudobulk_norm == "proportional_fitting") {
                data <- librarySizeNormalization(pb, mean(colSums(pb)))
            }
            #-- Compute
            if(method == "cibersortx") {
                message("Running CIBERSORTx...")
                assert(!is.null(ref$cached_results[["cibersortx"]]))
                ram_change <- peakRAM(
                    deconv_res <- cibersortx_deconvolute(data, ref, cache_path = method_cache, ...))
            } else if (method == "music") {
                message("Running MuSiC...")
                ram_change <- peakRAM(
                    deconv_res <- music_deconvolute(data, ref, ...))
            } else if (method == "dwls") {
                message("Running DWLS...")
                assert(!is.null(ref$cached_results[["dwls"]]))
                ram_change <- peakRAM(
                    deconv_res <- dwls_deconvolute(data, ref, ...))
            } else if(method == "ols") {
                assert(!is.null(ref$cached_results[["dwls"]]))
                message("Running OLS using the DWLS signature matrix...")
                ram_change <- peakRAM(
                    deconv_res <- ols_deconvolute(data, ref, ...))
            } else if(method == "svr") {
                assert(!is.null(ref$cached_results[["dwls"]]))
                message("Running SVR using the DWLS signature matrix...")
                ram_change <- peakRAM(
                    deconv_res <- svr_deconvolute(data, ref, ...))
            } else if(method == "bayesprism") {
                message("Running BayesPrism...")
                ram_change <- peakRAM(
                    deconv_res <- bayesprism_deconvolute(data, ref, cache_path = method_cache, ...))
            } else if(method == "bisque") {
                message("Running Bisque...")
                ram_change <- peakRAM(
                    deconv_res <- bisque_deconvolute(data, ref))
            } else if(method == "autogenes") {
                message("Running AutoGeneS...")
                ram_change <- peakRAM(
                    deconv_res <- autogenes_deconvolute(data, ref, ...))
            } else if(method == "scaden") {
                message("Running scaden...")
                ram_change <- peakRAM(
                    deconv_res <- scaden_deconvolute(data, ref, cache_path = method_cache, ...))
            } else {
                stop(paste0("`", method, "` is not supported for this application."))
            }
            #-- Get metrics
            metrics <- tibble(method = method,
                              level = level,
                              type = type,
                              time_elapsed_s = ram_change$Elapsed_Time_sec,
                              peak_ram_mib = ram_change$Peak_RAM_Used_MiB)
            #-- Save to cache
            saveRDS(deconv_res, file = filePath(method_cache, "deconv_res.RDS"))
            saveRDS(metrics, file = filePath(method_cache, "metrics.RDS"))
            #-- Return
            scbench$metrics <- bind_rows(scbench$metrics, metrics)
            scbench[["deconvolution"]][[level]][[type]][[method]] <- deconv_res
        }
    }
    lv_missing <- levels[sapply(missing_pops, length) > 0]
    names(missing_pops) <- levels
    if(length(lv_missing) > 0) {
        warning("Handling missing populations in benchmark...\n")
        for(lv in lv_missing) {
            to_remove <- paste0("frac_", missing_pops[[lv]])
            if(type == "population") {
                scbench[["deconvolution"]][[lv]][["complete_deconv"]][[method]] <-
                    scbench[["deconvolution"]][[lv]][[type]][[method]]
            }
            scbench[["deconvolution"]][[lv]][[type]][[method]] <-
                scbench[["deconvolution"]][[lv]][[type]][[method]] %>%
                dplyr::select(- {{to_remove}})
        }
    }

    if(length(levels) > 1 & type == "population" & correct_finer) {
        for(i in 2:length(levels)) {
            finer <- paste0("l", i); coarser <- paste0("l", i-1)
            finer_deconv <- scbench[["deconvolution"]][[finer]][["population"]][[method]]
            coarser_deconv <- scbench[["deconvolution"]][[coarser]][["population"]][[method]]
            pop_hierarchy <- scbench$pop_hierarchy
            hierarchy_list <- pop_hierarchy %>% pull(!! finer, !! coarser) %>% split(., names(.)) %>%
                sapply(unique)

            finer_deconv_norm <- .normalize_deconvolution_by_hierarchy(
                hierarchy_list, finer_deconv, coarser_deconv)

            scbench[["deconvolution"]][[finer]][[type]][[method]] <- finer_deconv_norm
        }

    }

    return(scbench)
}

#' Deconvolute bulk RNA-seq matrix using a chosen method
#'
#' @param bulk_data a count or linear-normalized (cpm, etc) matrix of genes-by-samples
#' @param scref an `screference` object containing single-cell reference data
#' and pre-computed references for methods that require them.
#' @param method a string with the name of the deconvolution method. For available
#' methods, consult `deconvolution_methods()`
#' @param bulk_norm the normalization to be applied in the bulk data. One of:
#' `"rpm"` for reads per million, `"none"` for raw counts, or `"proportional_fitting"`
#' for using the mean library size of the pseudobulk data to normalize. If already
#' normalized, make sure to choose `"none"` (default)
#' @param correct_finer normalize the resulting fractions from finer-grained
#' annotation using coarser-grained results
#' @param ... other parameters, passed to the method wrapper, enables the user
#' to change the method parameters. See: `deconvolute_{method}` where `{method}`
#' is the method name in lowercase for method-specific parameters
#'
#' @import tidyverse
#' @return a tibble with deconvolution results
#'
#' @rdname deconvolute
#' @export
deconvolute.matrix <- function(bulk_data, scref,
                               method = deconvolution_methods()[1],
                               bulk_norm = c("none", "rpm", "proportional_fitting")[1],
                               correct_finer = TRUE,
                               cache_path = ".",
                               ...) {
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(method %in% bulk_deconvolution_methods())

    if(class(scref) == "hscreference") {
        levels <- paste0("l", 1:scref$nlevels)
    } else {
        levels <- "l1"
    }
    if(bulk_norm == "rpm") {
        bulk_data <- librarySizeNormalization(pb, 10^6)
    } else if (bulk_norm == "proportional_fitting") {
        bulk_data <- librarySizeNormalization(pb, mean(colSums(pb)))
    }
    all_results <- list()
    for (level in levels) {
        if(class(scref) == "hscreference") {
            ref <- scref$screfs[[level]]
        } else {
            ref <- scref
        }
        #-- Compute
        if(method == "cibersortx") {
            message("Running CIBERSORTx...")
            assert(!is.null(ref$cached_results[["cibersortx"]]))
            if(cache_path == ".") cache_path <- "cibersortx"
            cache_path <- filePath(cache_path, level)
            deconv_res <- cibersortx_deconvolute(bulk_data, ref, cache_path = cache_path, ...)
        } else if (method == "music") {
            message("Running MuSiC...")
            deconv_res <- music_deconvolute(bulk_data, ref, ...)
        } else if (method == "dwls") {
            message("Running DWLS...")
            assert(!is.null(ref$cached_results[["dwls"]]))
            deconv_res <- dwls_deconvolute(bulk_data, ref, ...)
        } else if(method == "ols") {
            assert(!is.null(ref$cached_results[["dwls"]]))
            message("Running OLS using the DWLS signature matrix...")
            deconv_res <- ols_deconvolute(bulk_data, ref, ...)
        } else if(method == "svr") {
            assert(!is.null(ref$cached_results[["dwls"]]))
            message("Running SVR using the DWLS signature matrix...")
            deconv_res <- svr_deconvolute(bulk_data, ref, ...)
        } else if(method == "bayesprism") {
            message("Running BayesPrism...")
            deconv_res <- bayesprism_deconvolute(bulk_data, ref, cache_path = ".", ...)
        } else if(method == "bisque") {
            message("Running Bisque...")
            deconv_res <- bisque_deconvolute(bulk_data, ref)
        } else if(method == "autogenes") {
            message("Running AutoGeneS...")
            deconv_res <- autogenes_deconvolute(bulk_data, ref, ...)
        } else if(method == "scaden") {
            message("Running scaden...")
            deconv_res <- scaden_deconvolute(bulk_data, ref, ...)
        } else {
            stop(paste0("`", method, "` is not supported for this application."))
        }
        # deconv_res[sapply(deconv_res, is.na)] <- 0
        all_results[[level]] <- deconv_res
    }
    if(length(levels) > 1 & correct_finer) {
        for(i in 2:length(levels)) {
            finer <- paste0("l", i); coarser <- paste0("l", i-1)
            finer_deconv <- all_results[[finer]]
            coarser_deconv <- all_results[[coarser]]
            pop_hierarchy <- scref$hpop_table
            hierarchy_list <- pop_hierarchy %>% pull(!! finer, !! coarser) %>% split(., names(.))
            finer_deconv_norm <- .normalize_deconvolution_by_hierarchy(
                hierarchy_list, finer_deconv, coarser_deconv)

            all_results[[finer]] <- finer_deconv_norm
        }
    }
    class(all_results) <- c(class(all_results), "deconverse_results")
    return(all_results)

}

#' Deconvolute `Seurat` Spatial object using a chosen method
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay
#' @param scref an `screference` object containing single-cell reference data
#' and pre-computed references for methods that require them.
#' @param method a string with the name of the deconvolution method. For available
#' methods, consult `deconvolution_methods()`
#' @param bulk_norm the normalization to be applied in the bulk data. One of:
#' `"rp10k"` for reads per 10 thousand, `"none"` for raw counts, or `"proportional_fitting"`
#' for using the mean library size of the pseudobulk data to normalize. If already
#' normalized, make sure to choose `"none"` (default)
#' @param correct_finer normalize the resulting fractions from finer-grained
#' annotation using coarser-grained results
#' @param ... other parameters, passed to the method wrapper, enables the user
#' to change the method parameters. See: `deconvolute_{method}` where `{method}`
#' is the method name in lowercase for method-specific parameters
#'
#' @import tidyverse
#' @importFrom R.utils filePath
#'
#' @return an object of class `Seurat`, with deconverse results
#'
#' @rdname deconvolute
#' @export
deconvolute.Seurat <- function(spatial_obj, scref,
                               method = deconvolution_methods()[1],
                               normalization_method = c("none", "rp10k", "proportional_fitting")[1],
                               correct_finer = TRUE,
                               cache_path = "deconverse_cache",
                               ...) {
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(class(spatial_obj) == "Seurat")
    assert("Spatial" %in% names(spatial_obj@assays))
    assert(method %in% deconvolution_methods())

    if(class(scref) == "hscreference") {
        levels <- paste0("l", 1:scref$nlevels)
    } else {
        levels <- "l1"
    }
    all_results <- list()

    if(method %in% spatial_only_methods()) {
        for (level in levels) {
            if(class(scref) == "hscreference") {
                ref <- scref$screfs[[level]]
            } else {
                ref <- scref
            }
            if(method == "card") {
                message("Running CARD...")
                assert(!is.null(ref$cached_results[["card"]]))
                deconv_res <- card_deconvolute(spatial_obj, ref)
            } else if (method == "rctd") {
                message("Running RCTD")
                method_cache <- filePath(cache_path, "rctd", level)
                deconv_res <- rctd_deconvolute(spatial_obj, ref, cache_path = method_cache, ...)
            } else if (method == "spotlight") {
                message("Running SPOTlight")
                method_cache <- filePath(cache_path, "spotlight", level)
                deconv_res <- spotlight_deconvolute(spatial_obj, ref, cache_path = method_cache, ...)
            }
            all_results[[level]] <- deconv_res
        }
    } else {
        if(normalization_method == "rp10k") {
            bulk_data <- librarySizeNormalization(pb, 10^4)
        } else if (normalization_method == "proportional_fitting") {
            bulk_data <- librarySizeNormalization(pb, mean(colSums(pb)))
        }

        bulk_data <- spatial_obj@assays$Spatial@counts
        all_results <- deconvolute(as.matrix(bulk_data), scref,
                                   correct_finer = FALSE, bulk_norm = FALSE, ...)
    }
    if(length(levels) > 1 & correct_finer) {
        for(i in 2:length(levels)) {
            finer <- paste0("l", i); coarser <- paste0("l", i-1)
            finer_deconv <- all_results[[finer]]
            coarser_deconv <- all_results[[coarser]]
            pop_hierarchy <- scref$hpop_table
            hierarchy_list <- pop_hierarchy %>% pull(!! finer, !! coarser) %>% split(., names(.))
            finer_deconv_norm <- .normalize_deconvolution_by_hierarchy(
                hierarchy_list, finer_deconv, coarser_deconv)

            all_results[[finer]] <- finer_deconv_norm
        }
    }
    method_name <- all_results[[1]]$method[1]
    all_results_df <- lapply(levels, function(lv) {
        deconv_res <- all_results[[lv]]
        level_name <- if_else(length(levels) > 1, paste0("_", lv), "")
        deconv_df <- deconv_res %>%
            as.data.frame() %>%
            column_to_rownames("sample") %>%
            rename_with(.fn = ~ str_replace(., "frac_", paste0(method_name, level_name, "_")), .cols = starts_with("frac_")) %>%
            select(-method)
        return(deconv_df)
    })
    names(all_results_df) <- levels

    main_pop_df <- lapply(levels, function(lv) {
        deconv_df <- all_results_df[[lv]]
        level_name <- if_else(length(levels) > 1, paste0("_", lv), "")

        pop_names <- str_remove(colnames(deconv_df), paste0("^", method_name, level_name, "_"))
        main_pop_df <- data.frame(main_pop = pop_names[apply(deconv_df, 1, which.max)],
                                  row.names = rownames(deconv_df))
        colnames(main_pop_df) <- paste0(method_name, level_name, "_", "major_population")
        return(main_pop_df)
    }) %>% bind_cols()

    class(all_results) <- c(class(all_results), "deconverse_results")
    pop_vars <- lapply(all_results_df, colnames)

    metafeats_res <- left_df_join(main_pop_df, bind_cols(all_results_df))

    spatial_obj <- spatial_obj[,rownames(metafeats_res)]
    spatial_obj@tools[["deconverse"]][[method]] <- all_results
    spatial_obj@tools[["deconverse"]][["populations"]][[method]] <- pop_vars
    spatial_obj@tools[["deconverse"]][["major_population"]] <-
        c(spatial_obj@tools[["deconverse"]][["major_population"]], colnames(main_pop_df))
    spatial_obj@meta.data <- left_df_join(spatial_obj@meta.data, metafeats_res)

    return(spatial_obj)
}

#' Deconvolute `scbench` object using all methods with default settings
#'
#' @param scbench an `scbench` object already processed by `pseudobulks`
#' @param scref an `screference` object containing references for all methods that
#' require them
#' @param methods a vector of strings, which methods to use. For available
#' methods, consult `deconvolution_methods()`
#' @param ... passed to `deconvolute` for CIBERSORTx credentials
#'
#' @import tidyverse
#'
#' @return an object of class `scbench`
#' @rdname deconvolute_all
#'
#' @export
deconvolute_all.scbench <- function(scbench,
                                    scref,
                                    methods = deconvolution_methods(),
                                    ...){
    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(all(methods %in% deconvolution_methods()))

    #-- Get types
    types <- names(scbench$pseudobulk_counts)
    #-- Run deconvolutions
    for(type in types) {
        message("========= Deconvoluting pseudobulks for ", type, " analysis ==========")
        for(method in methods) {
            message("----------------> Using ", method, "...")
            if(method == "cibersortx") {
                scbench <- deconvolute(scbench, scref, method = method, type = type, ...)
            } else {
                scbench <- deconvolute(scbench, scref, method = method, type = type)
            }
        }
    }
    return(scbench)

}

#' Deconvolute bulk RNA-seq matrix using all methods with default settings
#'
#' @param bulk_data a count matrix of genes-by-samples
#' @param scref an `screference` object containing references for all methods that
#' require them
#' @param methods a vector of strings, which methods to use. For available
#' methods, consult `deconvolution_methods()`
#' @param ... passed to `deconvolute` for CIBERSORTx credentials
#'
#' @import tidyverse
#'
#' @return a list with each method's deconvolution results
#' @rdname deconvolute_all
#'
#' @export
deconvolute_all.matrix <- function(bulk_data,
                                   scref,
                                   methods = deconvolution_methods(),
                                   ...){
    #-- Error handling
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(all(methods %in% deconvolution_methods()))

    #-- Run deconvolutions
    all_deconv_res <- list()
    for(method in methods) {
        if(method == "cibersortx") {
            all_deconv_res[[method]] <- deconvolute(bulk_data, scref, method = method, ...)
        } else {
            all_deconv_res[[method]] <- deconvolute(bulk_data, scref, method = method)
        }
    }
    return(all_deconv_res)

}

#' Deconvolute Seurat spatial object using all methods with default settings
#'
#' @param bulk_data a count matrix of genes-by-samples
#' @param scref an `screference` object containing references for all methods that
#' require them
#' @param methods a vector of strings, which methods to use. For available
#' methods, consult `deconvolution_methods()`
#' @param ... passed to `deconvolute` for CIBERSORTx credentials
#'
#' @import tidyverse
#'
#' @return a list with each method's deconvolution results
#' @rdname deconvolute_all
#'
#' @export
deconvolute_all.Seurat <- function(spatial_obj,
                                   scref,
                                   methods = deconvolution_methods(),
                                   ...){
    #-- Error handling
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(all(methods %in% deconvolution_methods()))

    #-- Run deconvolutions
    for(method in methods) {
        if(method == "cibersortx") {
            spatial_obj <- deconvolute(spatial_obj, scref, method = method, ...)
        } else {
            spatial_obj <- deconvolute(spatial_obj, scref, method = method)
        }
    }
    return(spatial_obj)
}

# Gets -----------------------

#' Returns feature names of `deconverse` results that are stored in the Seurat object
#'
#' @param spatial_obj a Seurat object with a "Spatial" assay
#' @param method a string, corresponding to the method of deconvolution.
#' Must be one of `deconvolution_methods()`. If NULL, all available
#' deconvolution features are returned.
#' @param major_population if TRUE, the method returns a categorical
#' feature with the most abundant population predicted for each spot.
#' @param raw_results if TRUE, returns deconverse results as stored in the object
#' instead of feature names. Also available on `spatial_obj@tools$deconverse`
#'
#' @return a vector of strings correspoding to features in the Seurat object.
#' @rdname deconverse_results
#' @export
deconverse_results.Seurat <- function(spatial_obj, method = NULL, major_population = FALSE, raw_results = FALSE) {
    deconv_res <- spatial_obj@tools$deconverse

    if(raw_results) {
        return(deconv_res)
    }

    assert(method %in% names(deconv_res))
    assert(length(method) == 1)

    if(is.null(method)) {
        if(major_population) {
            results <- deconv_res$major_population
        } else {
            results <- deconv_res$populations
        }
        return(results)
    }

    method_name <- names(deconvolution_methods())[deconvolution_methods() == method]
    if(!is.null(deconv_res)) {
        if(major_population) {
            results <- str_subset(deconv_res$major_population, method_name)
        } else {
            results <- deconv_res$populations[[method]]
        }
    } else {
        stop("deconverse result for `", method, "` not found")
    }
    return(results)

}

# Generics -----------
#' @export
deconverse_results <- function(x, ...) {
    UseMethod("deconverse_results")
}
