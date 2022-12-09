# scbench ----------
#' Generate `scbench` object
#'
#' @param ref_scrna a reference Seurat object
#' @param pop_bounds list of population bounds, see example
#' @param annot_ids named vector with column names of annotations per level,
#' e.g. `c("l1" = "annot1", "l2" = "annot2")`
#' @param project_name a string indicating the project name, used for caching
#' @param cache_path path to the directory where results will be cached
#' @param batch_id NULL or a string of a column name with the batch information
#' for simulating pseudobulks by batch
#'
#' @import tidyverse
#' @return an object of class scbench
#'
#' @export
new_scbench <- function(ref_scrna,
                        pop_bounds,
                        annot_ids,
                        project_name = "project",
                        cache_path = "scbench_cache",
                        batch_id = NULL) {
    message("Generating new `scbench` object...")
    #-- Get number of levels
    levels <- str_remove(names(pop_bounds), "^l") %>% str_remove("_.*") %>% as.numeric()
    nl <- max(levels)

    #-- Checks
    assert(class(ref_scrna) == "Seurat")
    assert(str_detect(names(pop_bounds), "^l[0-9]+"))
    assert(all(sapply(pop_bounds, function(bounds) all(colnames(bounds) %in% c("population", "lower", "upper")))))
    if(nl > 1) {
        for(i in 2:nl) {
            subpops <- str_subset(names(pop_bounds), paste0("^l",i)) %>% str_remove(paste0("^l",i, "_"))
            lm1_bounds <- str_subset(names(pop_bounds), paste0("l", i-1))
            assert(any(sapply(lm1_bounds, function(lm1) subpops %in% pop_bounds[[lm1]]$population)))
        }
    }
    assert(all(annot_ids %in% colnames(ref_scrna@meta.data)))
    if(is.null(names(annot_ids))) {
        names(annot_ids) <- paste0("l", 1:nl)
    }
    colnames(ref_scrna@meta.data) <- plyr::mapvalues(colnames(ref_scrna@meta.data), annot_ids, names(annot_ids))
    if(!is.null(batch_id)) {
        assert(batch_id %in% colnames(ref_scrna@meta.data))
        colnames(ref_scrna@meta.data)[colnames(ref_scrna@meta.data) == batch_id] <- "batch_id"
    }
    pop_hierarchy <- .pop_hierarchy(pop_bounds, ref_scrna@meta.data)

    #-- Create object
    scbench <- list(ref_data = ref_scrna,
                    project_name = project_name,
                    pop_bounds = pop_bounds,
                    pop_hierarchy = pop_hierarchy,
                    nlevels = nl,
                    cache = cache_path,
                    ref_samps = NULL,
                    nsamps = NULL,
                    pop_props = NULL,
                    shrunk = FALSE,
                    deconvolution = NULL,
                    status = .update_status())
    class(scbench) <- "scbench"
    return(scbench)
}

# Mixtures ----------------------

#' Simulate sample population mixtures given bounds in an `scbench` object
#'
#' @param scbench an `scbench` object
#' @param nsamps an integer representing the number of samples to simulate
#' @param seed an integer representing a "seed", for reproducibility of mixtures.
#'
#' @import tidyverse testit
#' @importFrom hitandrun hitandrun
#' @importFrom fastDummies dummy_cols
#'
#' @return an object of class `scbench`
#'
#' @export
mixtures_population <- function(scbench, nsamps = 1000, seed = 0) {
    message("Simulating population mixtures given bounds...")
    assert(class(scbench) == "scbench")
    bounds <- scbench$pop_bounds
    #-- Hit-and-run MCMC
    samples <- lapply(bounds, .fit_hitandrun, n.samples = nsamps, seed)

    #-- Calculate higher level mixtures
    levels <- str_remove(names(bounds), "^l") %>% str_remove("_.*") %>% as.numeric()
    nlevs <- max(levels)
    pop_props <- list(l1 = samples[["l1"]])

    if(nlevs > 1) {
        for(i in 2:nlevs) {
            #-- Get information from higher level
            upper_lev <- paste0("l", i-1)
            lev_abbr <- paste0("l", i)
            coarse_tb <- pop_props[[upper_lev]]

            #-- Get subsetted populations
            ln_names <- str_subset(names(bounds), lev_abbr)
            rep_coarse <- str_remove(ln_names, "l[0-9]+_")
            missing_coarse <- setdiff(colnames(coarse_tb), rep_coarse)
            names(ln_names) <- rep_coarse

            #-- Calculate finer grained mixtures
            fine_pops <- lapply(1:length(ln_names), function(i) {
                coarse_tb[,names(ln_names)[i]] * samples[[ln_names[i]]]
            }) %>% bind_cols() %>% as.data.frame()
            fine_pops <- cbind(coarse_tb[,missing_coarse], fine_pops)
            colnames(fine_pops)[1:length(missing_coarse)] <- missing_coarse
            #-- Fix names if changed
            name_corresps <- split(pull(scbench$pop_hierarchy, !!lev_abbr), pull(scbench$pop_hierarchy, !!upper_lev)) %>%
                sapply(as.character)
            rep_coarse_new <- name_corresps[sapply(name_corresps, length) == 1] %>% unlist()
            colnames(fine_pops) <- str_replace_all(colnames(fine_pops), rep_coarse_new)
            pop_props[[lev_abbr]] <- fine_pops
        }
    }
    pop_props <- lapply(pop_props, as.data.frame)

    scbench[["mixtures"]][["population"]] <- pop_props
    scbench$status <- .update_status(scbench$status, "mixtures")

    return(scbench)
}

#' Simulate mixtures for measuring spillover between all pairs of populations
#' in an `scbench` object
#'
#' @param scbench an `scbench` object
#' @param step a float, defines the step in the sequence for creating mixtures
#' and simulating spillover effects
#'
#' @import tidyverse
#' @return an object of class `scbench`
#'
#' @export
mixtures_spillover <- function(scbench, step = 0.05) {
    message("Simulating spillover mixtures between population pairs...")
    #-- Checks
    assert(class(scbench) == "scbench")
    if(is.null(scbench[["mixtures"]][["population"]])) {
        stop("`population_mixtures` needs to be run before `spillover_mixtures`")
    }
    #-- Get combinations
    pop_props <- scbench[["mixtures"]][["population"]]
    base_props <- seq(0, 1, by = step)

    spillover_props <- lapply(pop_props, function(l_pop_props) {
        pops <- colnames(l_pop_props)
        l_combos <- t(combn(pops, 2))
        l_spillover_props <- lapply(1:nrow(l_combos), function(i) {
            combo <- l_combos[i,]
            combo_mat <- matrix(0, ncol = length(pops), nrow = length(base_props),
                                dimnames = list(NULL, pops))
            combo_mat[,combo[1]] <- base_props
            combo_mat[,combo[2]] <- 1-base_props
            combo_mat <- combo_mat %>%
                as.data.frame() %>%
                mutate(combo = paste(combo, collapse = "|")) %>%
                relocate(combo)
            return(combo_mat)
        }) %>% bind_rows()
    })
    scbench[["mixtures"]][["spillover"]] <- spillover_props
    scbench$status <- .update_status(scbench$status, "spillover")
    return(scbench)
}

#' Simulate mixtures to find the limit of detection of each population in
#' mixtures in an `scbench` object
#'
#' @param scbench an `scbench` object
#' @param max_prop a float, the maximum proportion of each target population
#' tested
#' @param step the step between 0 and `max_prop` in a sequence to be tested as
#' limit of detection
#' @param nsamp_per_step number of samples per step, used for statistics of the
#' results for the calibration curves
#' @param prop_noise_sd a value for added noise to the different samples generated
#' for each observation
#'
#' @return an object of class `scbench`
#' @export
mixtures_lod <- function(scbench, max_prop = 0.2, step = 0.01) {
    message("Simulating limits of detection for each population...")
    #-- Checks
    assert(class(scbench) == "scbench")
    if(is.null(scbench[["mixtures"]][["population"]])) {
        stop("`population_mixtures` needs to be run before `spillover_mixtures`")
    }
    #-- Get step
    pop_props <- scbench[["mixtures"]][["population"]]
    prop_steps <- seq(0, max_prop, by = step)
    nblanks <- 5
    blank_noise_sd <- 0.01
    #-- Get LoD mixtures
    lod_props <- lapply(pop_props, function(l_pop_props) {
        mean_prop <- colMeans(l_pop_props)
        pops <- colnames(l_pop_props)
        pop_refs <- mean_prop/sum(mean_prop)
        blank_refs <- lapply(1:(nblanks-1), function(i) {
            new_refs <- pop_refs + rnorm(length(pops), 0, blank_noise_sd)
            new_refs[new_refs < 0] <- 0
            new_refs/sum(new_refs)
            new_refs[new_refs < 0] <- 0
            return(new_refs)
            })
        pop_res <- lapply(pops, function(pop) {
            sapply(prop_steps, function(prop) {
                reduced_pops <- setdiff(pops, pop)
                step_props <- pop_refs
                step_props[reduced_pops] <- step_props[reduced_pops] + (pop_refs[pop]-prop)/length(reduced_pops)
                step_props[pop] <- prop
                return(step_props)
            }) %>%
                t() %>%
                as.data.frame() %>%
                tibble() %>%
                mutate(population = pop) %>%
                relocate(population)
        }) %>%
            bind_rows()

        extra_blanks <- lapply(pops, function(pop) {
            sapply(blank_refs, function(blank_ref) {
                blank_ref[pop] <- 0
                extra_blank <- blank_ref/sum(blank_ref)
                return(extra_blank)
            }) %>%
                t() %>%
                as.data.frame() %>%
                tibble() %>%
                mutate(population = pop) %>%
                relocate(population)
        }) %>%
            bind_rows()

        lod_mixes <- bind_rows(pop_res, extra_blanks) %>%
            mutate(sample = paste0("s", 1:n())) %>%
            relocate(sample)
        return(lod_mixes)
    })
    scbench[["mixtures"]][["lod"]] <- lod_props
    scbench$status <- .update_status(scbench$status, "lod")
    return(scbench)
}
# Pseudobulk ------------------------

#' Generate pseudobulk gene expression from single-cell RNA-seq given
#' population mixtures in an `scbench` object
#'
#' @param scbench an `scbench` object that has been evaluated by
#' `population_mixtures`
#' @param level a string character, the reference data level for pseudobulk pooling
#' (e.g. `"l1"`, `"l2"`)
#' @param ncells an integer, the number of cells to pool into pseudobulk.
#' @param shrink_obj a boolean, if TRUE, `scbench` object is reduced in size by
#' removing the single-cell expression data and keeping only its metadata after
#' estimating the pseudobulks
#' @param by_batch a boolean, if TRUE, pseudobulks are estimated by batch using
#' the `batch_id` provided when creating the object.
#' @param seed an integer representing a "seed", for reproducibility of mixtures.
#'
#' @import tidyverse
#'
#' @return an object of class `scbench`
#'
#' @export
pseudobulks <- function(scbench,
                        ncells = 2000,
                        ncores = 8,
                        shrink_obj = TRUE,
                        by_batch = FALSE,
                        seed = 0) {
    assert(class(scbench) == "scbench")
    assert("mixtures" %in% names(scbench))
    assert(!scbench$shrunk)

    #-- Get level (finest grained)
    levels <- names(scbench[["mixtures"]][["population"]])
    level <- levels[length(levels)]

    #-- Pseudobulks to generate
    types <- names(scbench[["mixtures"]])
    pseudobulks <- list()
    for (type in types) {
        if(type == "spillover") {
            ncells_run <- min(table(scbench$ref_data@meta.data[,level]))
        } else {
            ncells_run <- ncells
        }
        if(type == "population") {
            pb_res <- .get_pseudobulk(scbench, type,
                                      level = level, seed = seed,
                                      ncores = ncores, ncells = ncells_run,
                                      by_batch = by_batch)
            scbench[["mixtures"]][[type]][[level]] <- pb_res$props

            if(level != "l1") {
                scbench$pop_hierarchy
                finer_lvl <- str_remove(level, "l") %>% as.numeric()
                for(i in (finer_lvl-1):1) {
                    coarser_lvl <- paste0("l", i)
                    joins <- split(scbench$pop_hierarchy[[finer_lvl]], scbench$pop_hierarchy[[coarser_lvl]])
                    props_coarser <- sapply(joins, function(pops2join) {
                        rowSums(as.matrix(pb_res$props[,pops2join]))
                    })
                    rownames(props_coarser) <- rownames(pb_res$props)
                    scbench[["mixtures"]][[type]][[coarser_lvl]] <- props_coarser
                }
                scbench[["pseudobulk_counts"]][["population"]] <- pb_res$pseudobulk
                rm(pb_res)
            }
        } else {
            multi_pb_res <- lapply(levels, function(level) {
                .get_pseudobulk(scbench, type,
                                level = level, seed = seed,
                                ncores = ncores, ncells = ncells_run,
                                by_batch = by_batch)
            })
            pbs <- lapply(multi_pb_res, "[[", "pseudobulk")
            names(pbs) <- levels
            scbench[["pseudobulk_counts"]][[type]] <- pbs
            rm(multi_pb_res)

        }
    gc()
    }
    if(shrink_obj) {
        scbench$ref_data <- scbench$ref_data@meta.data
        gc()
    }
    scbench$shrunk <- shrink_obj
    scbench$status <- .update_status(scbench$status, "pseudobulks")

    return(scbench)
}


# Deconvolution ----------------------------------
#' @export
deconvolution_methods <- function() {
    c("ols", "dwls", "svr", "cibersortx", "music",  "bayesprism", "bisque", "autogenes")
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
#' @param ... other parameters, passed to the method wrapper, enables the user
#' to change the method parameters. See: `deconvolute_{method}` where `{method}`
#' is the method name in lowercase for method-specific parameters
#'
#' @import tidyverse
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
                        ...) {
    #-- Error handling
    # assert(class(scbench) == "scbench")
    assert(class(scref) %in% c("screference", "hscreference"))
    assert(!is.null(scbench$pseudobulk_counts[[type]]))

    #-- Reference picking
    if(class(scref) == "hscreference") {
        #-- Error handling
        assert(scref$nlevels == scbench$nlevels)
        levels <- paste0("l", 1:scbench$nlevels)
        same_populations <- sapply(1:scref$nlevels, function(lv) {
            bench_pops <- unique(scbench$pop_hierarchy[[lv]])
            all(bench_pops %in% intersect(bench_pops, unique(scref$hpop_table[[lv]])))
        })
        if(!all(same_populations)) {
            stop("Populations in `scref` and `scbench` are not the same at all levels.")
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
                deconv_res <- cibersortx_deconvolute(data, ref, cache_path = method_cache, ...)
            } else if (method == "music") {
                message("Running MuSiC...")
                deconv_res <- music_deconvolute(data, ref, ...)
            } else if (method == "dwls") {
                message("Running DWLS...")
                assert(!is.null(ref$cached_results[["dwls"]]))
                deconv_res <- dwls_deconvolute(data, ref, ...)
            } else if(method == "ols") {
                assert(!is.null(ref$cached_results[["dwls"]]))
                message("Running OLS using the DWLS signature matrix...")
                deconv_res <- ols_deconvolute(data, ref, ...)
            } else if(method == "svr") {
                assert(!is.null(ref$cached_results[["dwls"]]))
                message("Running SVR using the DWLS signature matrix...")
                deconv_res <- svr_deconvolute(data, ref, ...)
            } else if(method == "bayesprism") {
                message("Running BayesPrism...")
                deconv_res <- bayesprism_deconvolute(data, ref, cache_path = method_cache, ...)
            } else if(method == "bisque") {
                message("Running Bisque...")
                deconv_res <- bisque_deconvolute(data, ref)
            } else if(method == "autogenes") {
                message("Running AutoGeneS...")
                deconv_res <- autogenes_deconvolute(data, ref, ...)
            }
            #-- Save to cache
            saveRDS(deconv_res, file = filePath(method_cache, "deconv_res.RDS"))
            #-- Return
            scbench[["deconvolution"]][[level]][[type]][[method]] <- deconv_res
        }
    }
    if(length(levels) > 1 & type == "population") {
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
#' normalized, make sure to choose `"none"`.
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
                               bulk_norm = c("rpm", "none", "proportional_fitting")[3],
                               ...) {
    assert(class(scref) %in% c("screference", "hscreference"))

    if(class(scref) == "hscreference") {
        levels <- paste0("l", 1:scref$nlevels)
    } else {
        levels <- "l1"
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
            deconv_res <- cibersortx_deconvolute(data, scref, cache_path = cache_path, ...)
        } else if (method == "music") {
            message("Running MuSiC...")
            deconv_res <- music_deconvolute(data, ref, ...)
        } else if (method == "dwls") {
            message("Running DWLS...")
            assert(!is.null(scref$cached_results[["dwls"]]))
            deconv_res <- dwls_deconvolute(data, ref, ...)
        } else if(method == "ols") {
            assert(!is.null(scref$cached_results[["dwls"]]))
            message("Running OLS using the DWLS signature matrix...")
            deconv_res <- ols_deconvolute(data, ref, ...)
        } else if(method == "svr") {
            assert(!is.null(scref$cached_results[["dwls"]]))
            message("Running SVR using the DWLS signature matrix...")
            deconv_res <- svr_deconvolute(data, ref, ...)
        } else if(method == "bayesprism") {
            message("Running BayesPrism...")
            deconv_res <- bayesprism_deconvolute(data, ref, cache_path = cache_path, ...)
        } else if(method == "bisque") {
            message("Running Bisque...")
            deconv_res <- bisque_deconvolute(data, ref)
        } else if(method == "autogenes") {
            message("Running AutoGeneS...")
            deconv_res <- autogenes_deconvolute(data, ref, ...)
        }
        all_results[[level]] <- deconv_res
    }
    if(length(all_results) == 1) {
        return(unlist(all_results))
    } else {
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
        return(all_results)
    }
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
    assert(class(scref) == "screference")
    assert(all(methods %in% deconvolution_methods()))

    #-- Get types
    types <- names(scbench$pseudobulk_counts)
    #-- Run deconvolutions
    all_deconv_res <- list()
    for(method in methods) {
        if(method == "cibersortx") {
            all_deconv_res[[method]] <- deconvolute(bulk_data, scref, method = method, type = type, ...)
        } else {
            all_deconv_res[[method]] <- deconvolute(bulk_data, scref, method = method, type = type)
        }
    }
    return(all_deconv_res)

}

# Gets -------------------
# TODO: Improve document
#' Get benchmark results from `scbench` object
#'
#' @param scbench an `scbench` object already processed by `deconvolute` or
#' `deconvolute_all`
#' @param methods NULL, a string, or vector of strings. Results will be returned
#' for the methods indicated. See `deconvolution_methods()` for methods.
#' @param type type of deconvoluted mixtures. One of: `"population"`,
#' `"spillover"` or `"lod"`
#'
#' @import tidyverse
#'
#' @param a `tibble` with results.
#'
#' @export
get_benchmark_results <- function(scbench, methods = NULL, level = NULL, type = "population") {
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench[["mixtures"]]))
    assert(!is.null(scbench[["deconvolution"]]))

    if(is.null(methods)) {
        methods <- names(scbench[["deconvolution"]][[type]])
    }
    level <- .get_level(level, scbench)

    if(type == "population") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[[level]][["population"]][[method]]
            bench_tb <- .get_benchmark_fractions(scbench, deconv_res, level)
            results <- .get_bench_results(deconv_res, bench_tb) %>%
                mutate(method = method)
        }) %>%
            bind_rows()
    } else if(type == "spillover") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[[level]][["spillover"]][[method]]
            spillover_tb <- .get_spillover_mixtures(scbench, deconv_res, level)
            results <- .get_spillover_results(deconv_res, spillover_tb) %>%
                mutate(mixture = paste0(pop1, "|", pop2),
                       method = method) %>%
                select(sample, mixture, pop1, pop2, pop1_truth, pop1_pred = pop1_estimate, method)
        }) %>%
            bind_rows()
    } else if(type == "lod") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[[level]][["lod"]][[method]]
            lod_tb <- .get_lod_mixtures(scbench, deconv_res, level)
            results <- .get_lod_results(deconv_res, lod_tb)
            results <- results$lod_tb %>%
                select(sample, population, truth, pred = est, method) %>%
                left_join(select(results$lod_annot, population, method, lod = truth),
                          by = c("population", "method")) %>%
                relocate(method, .after = lod)
        }) %>%
            bind_rows()

    }
    return(results)
}

# Plotting -----------

#' Plot sample population mixtures from an `scbench` object
#' @param scbench an `scbench` object that has been evaluated by
#' `population_mixtures`
#' @param nshow an integer, the number of samples to plot
#' @param order_by a string of a name of a population to order the samples with
#' increasing mixtures of `order_by`
#'
#' @import tidyverse
#'
#' @return a ggplot object
#'
#' @export
plt_population_mixtures <- function(scbench, nshow = 50, order_by = NULL) {
    #-- Get data
    assert(class(scbench) == "scbench")
    assert("pop_props" %in% names(scbench))
    pop_props <- scbench[["mixtures"]][["population"]]

    pop_props <- lapply(pop_props, as_tibble)
    nlevs <- length(pop_props)
    #-- Plot examples
    snames <- paste0("S", 1:nshow)
    if(!is.null(order_by)) {
        sorder <- snames[order(as.data.frame(pop_props[[1]])[1:nshow, order_by])]
    } else {
        sorder <- snames
    }

    level_plots <- lapply(1:nlevs, function(i) {
        pop_tb <- pop_props[[paste0("l", i)]]
        pop_tb %>%
            slice(1:nshow) %>%
            mutate(sample = snames) %>%
            pivot_longer(cols = -sample, names_to = "population", values_to = "proportion") %>%
            mutate(population = factor(population, levels = colnames(pop_tb)),
                   sample = factor(sample, levels = sorder)) %>%
            ggplot(aes(sample, proportion, fill = population)) +
            geom_col() +
            labs(subtitle = paste0("Level ", i)) +
            theme_sparse()
    })
    return(level_plots)
}

#' Plot correlations by population between deconvolution results
#' for simulated populations and true mixtures used for the pseudobulk
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param method a string, one of `deconvolution_methods()`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#' @importFrom ggdensity geom_hdr_points
#' @importFrom ggheatmapper theme_scatter
#'
#' @return a ggplot object
#'
#' @export
plt_cors_scatter <- function(scbench, method, level = NULL) {
    assert(class(scbench) == "scbench")
    level <- .get_level(level, scbench)

    assert(!is.null(scbench$deconvolution[[level]]$population[[method]]))

    #-- Get deconv results
    deconv_res <- scbench$deconvolution[[level]][["population"]][[method]]

    #-- Get benchmark fractions
    bench_tb <- .get_benchmark_fractions(scbench, deconv_res, level)

    #-- Join data
    plot_tb <-.get_bench_results(deconv_res, bench_tb)

    #-- Plot
    plt_cor <- ggplot(plot_tb, aes(pred, truth)) +
        facet_wrap(~ population, scales = "free") +
        geom_hdr_points(probs = c(0.99,0.95,0.9,0.8,0.5,0.25,0.1), alpha = 0.8, size = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "black", lty = "dotted", size = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed", size = 0.5) +
        scale_color_brewer() +
        stat_cor(size = 3, method = "spearman") +
        labs(x = paste0("Fraction predicted by ", unique(deconv_res$method)), y = "True fraction") +
        guides(color = 'none') +
        theme_scatter()

    return(plt_cor)
}

#' Plot heatmap of correlation between deconvolution results for simulated
#' populations and true mixtures used for the pseudobulks
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `cor_table`, a
#' table summarizing results
#'
#' @importFrom ggheatmapper get_hmPlot update_hmPlot
#'
#' @export
plt_cor_heatmap <- function(scbench, level = NULL) {
    level <- .get_level(level, scbench)
    plot_tb <- .aggregate_deconvolution_results(scbench, level) %>%
        group_by(method, population) %>%
        summarize(pop_cors = cor(pred, truth, use = "pair")) %>%
        pivot_wider(names_from = method, values_from = pop_cors)

    mean_cors <- colMeans(plot_tb[,2:ncol(plot_tb)], na.rm = TRUE)
    methods <- names(mean_cors)[order(mean_cors)]

    plt <- ggheatmap(plot_tb,
              colv = "population",
              rowv = methods,
              cluster_rows = FALSE,
              hm_colors = viridis::viridis(100),
              hm_color_limits = c(0,1),
              colors_title = "Correlation",
              show_dend_row = FALSE,
              rows_title = "Method",
              column_title = "Populations")
    #-- Add numbers
    hm <- get_hmPlot(plt)
    annot_data <- hm$data %>%
        mutate(lev = ifelse(value < 0.5, "low", "high"))
    new_hm <- hm +
        geom_text(aes(x = observations, y = rows, label = signif(value, 3),
                      color = lev), size = 3, data = annot_data) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        guides(color = "none")
    plt <- update_hmPlot(plt, new_hm)

    return(list(heatmap = plt, cor_table = plot_tb))
}

#' Plot RMSE between deconvolution results for simulated populations and
#' true mixtures used for the pseudobulks
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `rmse_table`, a
#' table summarizing results
#'
#' @importFrom ggheatmapper ggheatmap align_to_hm theme_sparse2
#' @importFrom yardstick rmse_vec
#' @export
plt_rmse_heatmap <- function(scbench, level = NULL) {
    level <- .get_level(level, scbench)
    #-- Get plot data
    rmse_tb <- .aggregate_deconvolution_results(scbench, level) %>%
        group_by(method, population) %>%
        summarize(RMSE = yardstick::rmse_vec(truth, pred))
    plot_tb <- pivot_wider(rmse_tb, names_from = method, values_from = RMSE)

    #-- Overall RMSE
    annot <- rmse_tb %>%
        group_by(method) %>%
        summarize(RMSE = sum(RMSE)) %>%
        mutate(method = fct_reorder(method, RMSE, .desc = TRUE),
               metric = "RMSE",
               lev = ifelse(RMSE < max(RMSE)/2, "low", "high")) %>%
        arrange(method)

    #-- RMSE by population heatmap
    gghm <- ggheatmap(plot_tb,
              colv = "population",
              rowv = annot$method,
              cluster_rows = FALSE,
              hm_colors = viridis::magma(100),
              colors_title = "RMSE\nby population",
              show_dend_row = FALSE,
              rows_title = "Method",
              column_title = "Populations")


    annot_plt <- ggplot(annot, aes(metric, method, fill = RMSE)) +
        geom_tile() +
        geom_text(aes(label=signif(RMSE, 2), color = lev), size = 3) +
        scale_fill_viridis_c(option = "inferno", limits = c(0, max(annot$RMSE))) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        labs(x = "Metric", fill = "RMSE\ntotal") +
        guides(color = "none") +
        theme_sparse2()

    #-- Full heatmap
    plt <- align_to_hm(gghm, annot_plt, pos = "right", newplt_size_prop = 0.1, legend_action = "collect")

    return(list(heatmap = plt, rmse_table = plot_tb))
}

#' Plot correlations between deconvolution results and spillover mixtures
#' between each pair of populations
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param method a string, one of `deconvolution_methods()`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#'
#' @return a ggplot object
#'
#' @export
plt_spillover_scatter <- function(scbench, method, level = NULL) {
    level <- .get_level(level, scbench)

    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution[[level]]$spillover[[method]]))


    #-- Get deconv estimations
    deconv_res <- scbench$deconvolution[[level]]$spillover[[method]]
    #-- Get mixtures
    spillover_tb <- .get_spillover_mixtures(scbench, deconv_res, level)

    #-- Get estimations and join
    plot_tb <- .get_spillover_results(deconv_res, spillover_tb)

    #-- Plot
    plt <- ggplot(plot_tb, aes(pop1_estimate, pop1_truth)) +
        facet_grid(rows = vars(pop1), cols = vars(pop2)) +
        geom_point(color = "steelblue") +
        geom_smooth(method = "lm", se = FALSE, color = "black", lty = "dotted", size = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed", size = 0.5) +
        scale_color_brewer() +
        scale_x_continuous(breaks = c(0, 0.5, 1)) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        stat_cor(size = 3, method = "pearson") +
        labs(x = paste0("Population fraction predicted by ", unique(deconv_res$method)), y = "True population fraction") +
        guides(color = 'none') +
        theme_scatter()
    plt <- .remove_lowerti(plt, length(levels(plot_tb$pop1)), length(levels(plot_tb$pop2)))
    return(plt)
}

#' Plot limit of detection of deconvolution for each population as the proportion
#' of the pseudobulk where the method estimation is higher than at proportion=0 for
#' each population
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param method a string, one of `deconvolution_methods()`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#'
#' @return a ggplot object
#'
#' @export
plt_lod_scatter <- function(scbench, method, level = NULL) {
    #-- Limit of detection
    level <- .get_level(level, scbench)
    deconv_res <- scbench$deconvolution[[level]]$lod[[method]]

    #-- Get mixtures
    bench_tb <- .get_lod_mixtures(scbench, deconv_res, level)
    lod_res <- .get_lod_results(deconv_res, bench_tb)
    step <- lod_res$lod_tb$truth[2] - lod_res$lod_tb$truth[1]
    blanks <- lod_res$lod_tb %>%
        filter(truth == 0)

   plt <-
        ggplot(lod_res$lod_tb, aes(est, truth)) +
        facet_wrap(~ population) +
        geom_point(color = "steelblue") +
        geom_violin(aes(x = est), data = blanks, width = step, alpha = 0) +
        geom_segment(aes(x = blank_mean, xend = blank_mean, y = -step, yend = step),
                     data = lod_res$blank_metrics) +
        geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed", size = 0.5) +
        labs(x = paste0("Fraction predicted by ", unique(deconv_res$method)), y = "True fraction") +
        guides(color = 'none') +
        geom_hline(aes(yintercept = truth), data = lod_res$lod_annot, lty = "dotted") +
        geom_text(aes(x = 0.1, label = label), data = lod_res$lod_annot, size = 3, vjust = -0.3, hjust = 0.5) +
        theme_scatter()

    return(plt)
}

#' Plot heatmap summarizing spillover RMSE by population and method
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `rmse_table`, a
#' table summarizing results
#'
#' @importFrom data.table rbindlist
#' @export
plt_spillover_heatmap <- function(scbench, level = NULL) {
    level <- .get_level(level, scbench)

    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution[[level]]$spillover))
    methods <- names(scbench$deconvolution[[level]]$spillover)

    #-- Get deconv estimations
    all_deconv_res <- lapply(methods, function(method) {
        scbench$deconvolution[[level]]$spillover[[method]]
    })
    # level <- .match_level(scbench, all_deconv_res[[1]], "spillover")
    #-- Get mixtures
    spillover_tb <- .get_spillover_mixtures(scbench, all_deconv_res[[1]], level)

    #-- Get estimations and join
    all_res <- lapply(all_deconv_res, .get_spillover_results, spillover_tb)
    names(all_res) <- methods
    all_tb <- rbindlist(all_res, idcol = TRUE) %>%
        rename(method = .id) %>%
        tibble()

    #-- Get RMSE
    plot_tb <- all_tb %>%
        group_by(method, pop1, pop2) %>%
        summarize(rmse = yardstick::rmse_vec(pop1_truth, pop1_estimate)) %>%
        pivot_wider(names_from = method, values_from = rmse) %>%
        mutate(pop_pairs = paste0(pop1, "|", pop2)) %>%
        ungroup()
    #-- Total RMSE per method
    method_total_error <- all_tb %>%
        group_by(method, pop1, pop2) %>%
        summarize(rmse = yardstick::rmse_vec(pop1_truth, pop1_estimate)) %>%
        group_by(method) %>%
        summarize(sum_rmse = sum(rmse)) %>%
        mutate(metric = "RMSE") %>%
        arrange(sum_rmse) %>%
        mutate(method = factor(method, levels = method),
               lev = ifelse(sum_rmse < max(sum_rmse)/2, "low", "high"))

    #-- Heatmap
    spillover_hm <- ggheatmap(plot_tb,
                              colv = "pop_pairs",
                              rowv = rev(method_total_error$method),
                              cluster_rows = FALSE,
                              cluster_cols = TRUE,
                              hm_colors = viridis::magma(100),
                              colors_title = "RMSE",
                              show_dend_row = FALSE,
                              rows_title = "Method",
                              column_title = "Population pairs")
    #-- Add total track
    annot_plt <- ggplot(method_total_error, aes(metric, fct_rev(method), fill = sum_rmse)) +
        geom_tile() +
        geom_text(aes(label=signif(sum_rmse, 2), color = lev), size = 3) +
        scale_fill_viridis_c(option = "inferno", limits = c(0, max(method_total_error$sum_rmse))) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        labs(x = "Metric", fill = "RMSE\ntotal") +
        guides(color = "none") +
        theme_sparse2()

    spillover_plt <- align_to_hm(spillover_hm, annot_plt, pos = "right", newplt_size_prop = 0.1, legend_action = "collect")
    return(list(heatmap = spillover_plt, rmse_table = plot_tb))
}
#' Plot heatmap summarizing limit of detection by population and method
#' @param scbench an `scbench` object that has been evaluated by
#' `deconvolute`
#' @param level an annotation level to show. Finest-grained annotation (highest level)
#' by default.
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `lod_table`, a
#' table summarizing results
#'
#' @export
plt_lod_heatmap <- function(scbench, level = NULL) {
    #-- Error handling
    level <- .get_level(level, scbench)
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution[[level]]$lod))
    methods <- names(scbench$deconvolution[[level]]$lod)

    #-- Setup
    all_deconv_res <- lapply(methods, function(method) {
        deconv_res <- scbench$deconvolution[[level]]$lod[[method]]
        deconv_res[is.na(deconv_res)] <- 0
        return(deconv_res)
    })

    #-- Get mixtures and estimations, calculate limits of detection
    bench_tb <- .get_lod_mixtures(scbench, all_deconv_res[[1]], level)
    all_lod_res <- lapply(all_deconv_res, .get_lod_results, bench_tb)
    lod_tb <- lapply(all_lod_res, function(method_res) { method_res$lod_annot }) %>%
        bind_rows() %>%
        select(population, method, lod = truth)

    #-- Organize for plotting
    lods <- lapply(all_lod_res, function(method_res) { method_res$lod_annot }) %>%
        bind_rows() %>%
        select(population, method, lod = truth) %>%
        pivot_wider(names_from = method, values_from = lod) %>%
        ungroup()

    method_order <- lod_tb %>%
        group_by(method) %>%
        summarize(mean_lod = mean(lod)) %>%
        arrange(mean_lod) %>%
        pull(method) %>%
        rev()
    #-- Plot
    lod_hm <- ggheatmap(lods,
                        colv = "population",
                        rowv = method_order,
                        cluster_rows = FALSE,
                        cluster_cols = TRUE,
                        hm_colors = viridis::viridis(100),
                        hm_color_limits = c(0, max(bench_tb$truth)),
                        colors_title = "Limit of detection (fraction)",
                        show_dend_row = FALSE,
                        rows_title = "Method",
                        column_title = "Population")
    #-- Add numbers
    hm <- get_hmPlot(lod_hm)
    annot_data <- hm$data %>%
        mutate(lev = ifelse(value < max(bench_tb$truth)/2, "low", "high"))
    new_hm <- hm +
        geom_text(aes(x = observations, y = rows, label = value,
                      color = lev), size = 3, data = annot_data) +
        scale_color_manual(values = c("low" = "white", "high" = "black")) +
        guides(color = "none")
    lod_hm <- update_hmPlot(lod_hm, new_hm)

    return(list(heatmap = lod_hm, lod_table = lods))
}

# Generics ------------------
#' @export
deconvolute <- function(x, ...) {
    UseMethod("deconvolute")
}
#' @export
deconvolute_all <- function(x, ...) {
    UseMethod("deconvolute_all")
}

# Helpers ------------
#' @importFrom R.utils filePath
.scbench_cache_check <- function(method_cache) {
    method_fname <- filePath(method_cache, "deconv_res.RDS")
    if(file.exists(method_fname)) {
        message("Results found in cache, returning...")
        deconv_res <- readRDS(method_fname)
    } else {
        dir.create(method_cache, showWarnings = FALSE, recursive = TRUE)
        deconv_res <- NULL
    }
    return(deconv_res)
}

#' @import tidyverse
.match_level <- function(scbench, deconv_res, method) {
    pop_names <- str_subset(colnames(deconv_res), "frac") %>% str_remove("frac_")
    level_match <- sapply(scbench[["mixtures"]][[method]], function(props) {
        all(pop_names %in% colnames(props))
    }) %>% which()
    level <- names(level_match)
    return(level)
}

#' @importFrom Matrix colSums
librarySizeNormalization <- function(data, factor = 10^6) {
    data/(colSums(data)/factor)
}

#' @export
print.scbench <- function(x) {
    cat(str_glue("scbench object named {x$project_name} with {length(unique(unlist(x$pop_hierarchy)))} reference populations and {x$nlevels} levels of annotation"))
    cat("\n")
    cat(x$status)
    cat("\n")
    if(!is.null(x[["deconvolution"]])) {
        cat("deconvolution results: ")
        cat("\n")
        levels <- names(x[["deconvolution"]])
        for(level in levels) {
            types <- names(x[["deconvolution"]][[level]])
            cat(str_glue("   {level}: "))
            cat("\n")
            for (type in types) {
                cat(str_glue("      {type}:"), paste(names(x[["deconvolution"]][[level]][[type]]), collapse = ", "))
                cat("\n")
            }
        }
    }
}
#' @import hitandrun
.fit_hitandrun <- function(pop_bounds, seed = 0, ...) {
    if(!is.null(seed)) {
        set.seed(seed)
    }
    cnstr <- .constraints_from_bounds(pop_bounds)
    hr <- hitandrun(cnstr, ...)
    colnames(hr) <- pop_bounds$population
    rownames(hr) <- paste0("s", 1:nrow(hr))
    return(hr)
}

#' @import tidyverse
.pop_hierarchy <- function(bounds, ref_meta) {
    nlevels <- str_extract(names(bounds), "l.") %>% unique() %>% str_remove("^l") %>% as.numeric %>% max()
    pop_hierarchy <- tibble(l1 = bounds$l1$population)
    if(nlevels > 1) {
        for (i in 2:nlevels) {
            finer_lev <- paste0("l", i)
            coarser_lev <- paste0("l", i-1)
            upper_pops <- unique(ref_meta[,coarser_lev])
            split_pops <- unique(ref_meta[,finer_lev])

            new_hierarch <- table(ref_meta[,coarser_lev], ref_meta[,finer_lev]) %>%
                as.data.frame() %>%
                filter(Freq > 0) %>%
                select(-Freq)
            colnames(new_hierarch) <- c(coarser_lev, finer_lev)
            # Check splits
            split_pops <- upper_pops[sapply(split(new_hierarch[,2], new_hierarch[,1]), length) > 1]
            split_pops <- paste0(finer_lev, "_", split_pops)
            assert(split_pops %in% names(bounds))

            pop_hierarchy <- right_join(pop_hierarchy, new_hierarch)
        }
    }
    return(pop_hierarchy)
}

#' @import tidyverse
.normalize_deconvolution_by_hierarchy <- function(hierarchy_list, finer_deconv, coarser_deconv) {
    finer_deconv_norm <- lapply(1:length(hierarchy_list), function(j) {
        coarse_col <- paste0("frac_", names(hierarchy_list)[j])
        coarse_vec <- coarser_deconv[[coarse_col]]

        fine_cols <- paste0("frac_", hierarchy_list[[j]])
        fine_mat <- finer_deconv[,fine_cols]
        fine_mat_norm <- lapply(1:nrow(fine_mat), function(ii) {
            fine_mat[ii,] / (sum(fine_mat[ii,])/coarse_vec[ii])
        }) %>% bind_rows()
        return(fine_mat_norm)
    }) %>% bind_cols() %>%
        tibble() %>%
        mutate(sample = finer_deconv$sample,
               method = finer_deconv$method) %>%
        relocate(sample) %>%
        .[,colnames(finer_deconv)]

    return(finer_deconv_norm)
}

#' @import tidyverse
.constraints_from_bounds <- function(pop_bounds) {
    pops <- pop_bounds$population
    oh_pops <- dummy_cols(pop_bounds$population) %>%
        rename_with(~ str_replace(., "\\.data", "pop")) %>%
        select(-pop) %>%
        .[,paste0("pop_", pops)]
    A <- cbind(cbind(oh_pops, -oh_pops), total = rep(1,length(pops))) %>% t()
    rownames(A) <- NULL
    colnames(A) <- pops
    b <- c(pop_bounds$upper, -pop_bounds$lower, 1)
    d <- c(rep("<=", length(b)-1), "=")
    constr <- list(constr=A, rhs=b, dir=d)
    return(constr)
}
#' @import tidyverse
theme_sparse <- function (...)
{
    theme_minimal(...) +
        theme(panel.grid = element_blank(),
              axis.title.x = element_blank(), axis.ticks.x = element_blank(),
              axis.text.x = element_blank(), axis.ticks.y = element_line(color = "black"),
              axis.text.y = element_text(color = "black"))
}
#' @import tidyverse
.aggregate_deconvolution_results <- function(scbench, level) {
    pop_props <- scbench[["mixtures"]][["population"]]
    #-- Get deconv results
    deconv_res <- bind_rows(scbench$deconvolution[[level]]$population)

    #-- Get benchmark fractions
    # level <- .match_level(scbench, deconv_res, "population")
    bench_tb <- pop_props[[level]] %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        rename_with(~ paste0("truth_", .), .cols = -sample)

    deconv_all <- deconv_res %>%
        left_join(bench_tb, by = "sample")

    aggregated_data <- deconv_all %>%
        pivot_longer(cols = matches("frac|truth"), names_to = "population", values_to = "fraction") %>%
        mutate(type = ifelse(str_detect(population, "frac"), "pred", "truth"),
               population = str_remove(population, "frac_|truth_")) %>%
        pivot_wider(names_from = "type", values_from = "fraction")


    return(aggregated_data)
}
.update_status <- function(prev = NULL, status = NULL) {
    if(is.null(prev)) {
        status <- c("Bounds [x] | Mixtures [ ] | Spillover [ ] | Limit of Detection [ ] | Pseudobulks [ ]")
        return(status)
    }
    prev_vec <- str_split_fixed(prev, pattern = " \\| ", 5) %>% as.vector()
    if(status == "mixtures") {
        prev_vec[2] <- "Mixtures [x]"
    } else if(status == "spillover") {
        prev_vec[3] <- "Spillover [x]"
    } else if (status == "lod") {
        prev_vec[4] <- "Limit of Detection [x]"
    } else if (status == "pseudobulks") {
        prev_vec[5] <- "Pseudobulks [x]"
    }
    new_status <- paste(prev_vec, collapse = " | ")
    return(new_status)
}
#' @import tidyverse
.get_benchmark_fractions <- function(scbench, deconv_res, level) {
    # level <- .match_level(scbench, deconv_res, "population")
    bench_tb <- scbench[["mixtures"]][["population"]][[level]] %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        mutate(sample = ifelse(str_detect(sample, "s"), sample, paste0("s", sample))) %>%
        tibble() %>%
        rename_with(~ paste0("truth_", .), .cols = -sample)
    return(bench_tb)
}
#' @import tidyverse
.get_bench_results <- function(deconv_res, bench_tb) {
    deconv_res %>%
        select(-starts_with("model")) %>%
        left_join(bench_tb) %>%
        pivot_longer(cols = matches("frac|truth"), names_to = "population", values_to = "fraction") %>%
        mutate(type = ifelse(str_detect(population, "frac"), "pred", "truth"),
               population = str_remove(population, "frac_|truth_")) %>%
        pivot_wider(id_cols = c("sample", "population"), names_from = "type", values_from = "fraction")
}

#' @import tidyverse
#' @importFrom pbmcapply pbmclapply
#' @importFrom R.utils filePath
#' @importFrom Matrix rowSums
.get_pseudobulk <- function(scbench, type, level, seed, ncores, ncells, by_batch) {
    cache_dir <- filePath(scbench$cache, scbench$project_name, level, "pseudobulks", type)
    cache_loc <- filePath(cache_dir, "pseudobulks.RDS")
    if(file.exists(cache_loc)) {
        message("Pseudobulks for ", type, " found in cache. Skipping...")
        pb_res <- readRDS(cache_loc)
        return(pb_res)
    }
    message("Generating pseudobulks for ", type, " analysis...")
    pop_props <- scbench$mixtures[[type]]

    #-- Get metadata
    if(type == "population") {
        pb_props <- pop_props[[level]]
    } else if (type == "lod") {
        pb_props <- pop_props[[level]][,-c(1,2)]
    } else {
        pb_props <- pop_props[[level]][,-1]
    }
    pops <- colnames(pb_props)
    meta <- scbench$ref_data@meta.data
    nsamps <- dim(pb_props)[1]

    if(!is.null(seed)) set.seed(seed)

    if(by_batch & ("batch_id" %in% colnames(meta)) & type == "population") {
        message("-- Pooling cells by batch...")
        #-- Get reference samples for each pseudosample
        pids <- meta[,"batch_id"]
        cells_by_pid <- split(rownames(meta), pids)
        pops_pat_cells <- lapply(cells_by_pid, function(cell_ids) {
            split(rownames(meta[cell_ids,]), meta[cell_ids,level])[pops]
        })
        upid <- unique(pids)
        #------- Check ref samps
        ref_samp_missing <- lapply(pops_pat_cells, function(samp_cells) {
            colnames(pb_props)[!colnames(pb_props) %in% names(samp_cells)]
        })
        ref_incomplete <- sapply(ref_samp_missing, length) > 0
        incomplete_pids <- names(ref_incomplete)[ref_incomplete]
        if(length(incomplete_pids) > 0) {
            warning(str_glue("Batches {paste(incomplete_pids, collapse = ', ')} don't contain at least one cell with each annotation.
                             Adapting mixtures for these samples..."))
        }
        #-- Get cells per sample
        nsamp_id <- trunc(nsamps/length(upid))
        ref_samp <- rep(upid, nsamp_id)
        rmd <- nsamps%%length(upid)
        if(rmd > 0) ref_samp <- c(ref_samp, upid[1:rmd])
        scbench$ref_samps <- ref_samp

        pseudobulks <- pbmclapply(1:nsamps, function(i) {
            samp <- ref_samp[i]
            sprop <- pb_props[i,]
            #-- Adjust sprop if needed
            if(samp %in% incomplete_pids) {
                missing_pop <- ref_samp_missing[[samp]]
                prop_redistr <- sum(sprop[missing_pop])
                sprop <- sprop + prop_redistr/(length(sprop)-1)
                sprop[missing_pop] <- 0
                pb_props[i,] <- sprop
            }
            sample_ncells <- round(sprop * ncells)
            sample_cells <- lapply(1:length(sample_ncells), function(j) {
                ps_n <- as.integer(sample_ncells[j])
                pop_cells <- unlist(pops_pat_cells[[samp]][names(sample_ncells)[j]])
                if(length(pop_cells) < ps_n) {
                    warning("--- Some cell profiles had to be repeated... `ncells` might be set too high for this dataset.")
                }
                sample(pop_cells, ps_n, replace = TRUE)
            }) %>% reduce(c)
            rowSums(scbench$ref_data@assays$RNA@counts[,sample_cells])
        }, mc.cores = ncores)
    } else {
        #-- Get pseudobulks
        message("-- Pooling cells...")
        pops_cells <- split(rownames(meta), meta[,level])[pops]
        pseudobulks <- pbmclapply(1:nsamps, function(i) {
            sample_ncells <- round(pb_props[i,] * ncells)
            sample_ncells[sample_ncells < 0] <- 0
            sample_cells <- lapply(1:length(sample_ncells), function(j) {
                ps_n <- as.integer(sample_ncells[j])
                pop_cells <- unlist(pops_cells[names(sample_ncells)[j]])
                if(length(pop_cells) < ps_n) {
                    warning("--- Some cell profiles had to be repeated... `ncells` might be set too high for this dataset.")
                    sample(pop_cells, ps_n, replace = TRUE)
                } else {
                    sample(pop_cells, ps_n)

                }
            }) %>% reduce(c)
            rowSums(scbench$ref_data[,sample_cells]@assays$RNA@counts)
        }, mc.cores = ncores)
        message("-- Joining matrix...")
    }
    pb_mat <- reduce(pseudobulks, cbind)
    colnames(pb_mat) <- paste0("s", 1:ncol(pb_mat))
    pb_res <- list(pseudobulk = pb_mat, props = pb_props)
    gc()

    if(is.character(pb_res$pseudobulk[1,1])) {
        stop(str_glue("Pseudobulk for `{type}` failed..."))
    }
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(pb_res, file = cache_loc)
    return(pb_res)
}

#' @import tidyverse
#' @importFrom data.table as.data.table
.get_spillover_mixtures <- function(scbench, deconv_res, level) {
    level <- .match_level(scbench, deconv_res, method = "spillover")
    bench_tb <- scbench[["mixtures"]][["spillover"]][[level]] %>%
        tibble() %>%
        mutate(
            sample = paste0("s", 1:n()),
            pop1 = str_remove(combo, "\\|.*"),
            pop2 = str_remove(combo, ".*\\|")) %>%
        relocate(sample, pop1, pop2) %>%
        select(-combo) %>%
        as.data.table()

    bench_tb[,pop1_truth := get(pop1), by = pop1]
    bench_tb[,pop2_truth := get(pop2), by = pop2]
    spillover_tb <- bench_tb %>%
        select(sample, pop1, pop2, pop1_truth, pop2_truth) %>%
        as_tibble()
    return(spillover_tb)
}

#' @import tidyverse
.get_spillover_results <- function(deconv_res, spillover_tb) {
    spillover_deconv <- deconv_res %>%
        select(-starts_with("model")) %>%
        left_join(spillover_tb, by = "sample") %>%
        mutate(pop1_est = paste0("frac_", pop1),
               pop2_est = paste0("frac_", pop2)) %>%
        as.data.table()
    spillover_deconv[,pop1_estimate := get(pop1_est), by = pop1_est]
    spillover_deconv[,pop2_estimate := get(pop2_est), by = pop2_est]
    plot_tb <- spillover_deconv %>%
        as_tibble() %>%
        select(sample, pop1, pop2, pop1_truth, pop2_truth, pop1_estimate, pop2_estimate) %>%
        mutate(pop1 = factor(pop1, levels = unique(spillover_tb$pop1)),
               pop2 = factor(pop2, levels = unique(spillover_tb$pop2)))
}
#' @import tidyverse
#' @importFrom data.table as.data.table
.get_lod_mixtures <- function(scbench, deconv_res, level) {
    # level <- .match_level(scbench, deconv_res, "lod")
    bench_tb <- scbench$mixtures$lod[[level]] %>% as.data.table()
    bench_tb[,truth := get(population), by = population]
    bench_tb <- bench_tb %>%
        select(sample, population, truth) %>%
        as_tibble()
}
#' @import tidyverse
#' @importFrom data.table as.data.table
#' @importFrom broom tidy
.get_lod_results <- function(deconv_res, bench_tb) {
    lod_tb <- bench_tb %>%
        left_join(deconv_res, by = "sample") %>%
        mutate(pop_name = paste0("frac_", population)) %>%
        as.data.table()
    lod_tb[,est := get(pop_name), by = pop_name]
    lod_tb <- lod_tb %>%
        select(sample, population, method, truth, est) %>%
        as_tibble()

    #-- Get blank stats
    blank_metrics <- lod_tb %>%
        filter(truth == 0) %>%
        group_by(population, method) %>%
        summarize(blank_mean = mean(est, na.rm = TRUE), blank_sd = sd(est, na.rm = TRUE))

    lod_filt <- lod_tb %>%
        filter(truth != 0)

    models <- lod_filt %>%
        group_by(population, method) %>%
        do(model = lm(est ~ truth, data = .))

    reg_slopes <- sapply(models$model, function(model) {
        tidy(model)$estimate[1]
    })
    reg_intercepts <- sapply(models$model, function(model) {
        tidy(model)$estimate[2]
    })

    lod_metrics <- models %>%
        select(-model) %>%
        bind_cols(tibble(slopes = reg_slopes,
               intercepts = reg_intercepts)) %>%
        left_join(blank_metrics) %>%
        mutate(lod_min = (blank_mean + 2*blank_sd))

    lod_annot <- lod_tb %>%
        left_join(select(lod_metrics, population, method, lod_min),
                  by = c("population", "method")) %>%
        filter(truth > 0) %>%
        mutate(find_lod = est > lod_min) %>%
        filter(find_lod) %>%
        group_by(population) %>%
        slice_min(truth) %>%
        mutate(label = paste0("LoD = ", signif(truth, 3))) %>%
        select(population, method, truth, est, label)


    return(list(lod_tb = lod_tb, lod_annot = lod_annot, blank_metrics = blank_metrics))
}

.get_level <- function(level, scbench) {
    if(is.null(level)) {
        message("Using finest grained annotation...")
        level <- paste0("l", scbench$nlevels)
    } else{
        assert(level %in% names(scbench$deconvolution))
    }
    return(level)
}

# credit: https://stackoverflow.com/questions/57927458/get-rid-of-empty-panels-in-the-first-row-of-facet-grid
#' @import tidyverse
#' @importFrom cowplot plot_to_gtable
#' @importFrom ggpubr as_ggplot
.remove_lowerti <- function(plt, n_row_fct, n_col_fct) {
    #-- Get indices of panels to remove
    idx <- lower.tri(matrix(0, nrow = n_row_fct, ncol = n_col_fct))
    panels_remove <- lapply(1:ncol(idx), function(j) {
        i <- which(idx[,j])
        data.frame(i = i, j = rep(j, length(i)))
    } ) %>%
        bind_rows() %>%
        mutate(panel = paste0("panel-", i, "-", j)) %>%
        pull(panel)

    #-- Draw plot without removed panels
    gt <- plot_to_gtable(plt)
    to.delete <- which (gt$layout$name %in% panels_remove)

    gt$grobs[to.delete] <- NULL
    gt$layout <- gt$layout[-to.delete, ]
    as_ggplot(gt)
}
