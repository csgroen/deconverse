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
    if(nl > 1) {
        for(i in 2:nl) {
            subpops <- str_subset(names(pop_bounds), paste0("^l",i)) %>% str_remove(paste0("^l",i, "_"))
            assert(subpops %in% pop_bounds[[paste0("l", i-1)]]$population)
        }
    }
    assert(all(annot_ids %in% colnames(ref_scrna@meta.data)))
    colnames(ref_scrna@meta.data) <- plyr::mapvalues(colnames(ref_scrna@meta.data), annot_ids, names(annot_ids))
    if(!is.null(batch_id)) {
        assert(batch_id %in% colnames(ref_scrna@meta.data))
        colnames(ref_scrna@meta.data)[colnames(ref_scrna@meta.data) == batch_id] <- "batch_id"
    }

    #-- Create object
    scbench <- list(ref_data = ref_scrna,
                    project_name = project_name,
                    pop_bounds = pop_bounds,
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
mixtures_population <- function(scbench, nsamps = 1000, seed = NULL) {
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
            rep_coarse <- str_remove(ln_names, "l.*_")
            missing_coarse <- setdiff(colnames(coarse_tb), rep_coarse)
            names(ln_names) <- rep_coarse

            #-- Calculate finer grained mixtures
            fine_pops <- lapply(1:length(ln_names), function(i) {
                coarse_tb[,names(ln_names)[i]] * samples[[ln_names[i]]]
            }) %>% bind_cols() %>% as.data.frame()
            fine_pops <- cbind(coarse_tb[,missing_coarse], fine_pops)
            colnames(fine_pops)[1:length(missing_coarse)] <- missing_coarse

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
#'
#' @return an object of class `scbench`
#' @export
mixtures_lod <- function(scbench, max_prop = 0.1, step = 0.005) {
    message("Simulating limits of detection for each population...")
    #-- Checks
    assert(class(scbench) == "scbench")
    if(is.null(scbench[["mixtures"]][["population"]])) {
        stop("`population_mixtures` needs to be run before `spillover_mixtures`")
    }
    #-- Get step
    pop_props <- scbench[["mixtures"]][["population"]]
    prop_steps <- seq(0, max_prop, by = step)
    #-- Get LoD mixtures
    lod_props <- lapply(pop_props, function(l_pop_props) {
        mean_prop <- colMeans(l_pop_props)
        ref_props <- mean_prop/sum(mean_prop)

        pops <- colnames(l_pop_props)
        lapply(pops, function(pop) {
            sapply(prop_steps, function(prop) {
                reduced_pops <- setdiff(pops, pop)
                step_props <- ref_props
                step_props[reduced_pops] <- step_props[reduced_pops] + (ref_props[pop]-prop)/length(reduced_pops)
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
    } )
    scbench[["mixtures"]][["lod"]] <- lod_props
    scbench$status <- .update_status(scbench$status, "lod")
    return(scbench)
}
# Pseudobulk and deconvolution ------------------------
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
                        level = NULL,
                        ncells = 2000,
                        ncores = 8,
                        shrink_obj = TRUE,
                        by_batch = FALSE,
                        seed = NULL) {
    assert(class(scbench) == "scbench")
    assert("mixtures" %in% names(scbench))
    assert(!scbench$shrunk)

    #-- Get level
    if(is.null(level)) {
        level <- names(scbench[["mixtures"]][["population"]])[length(pop_props)]
    }

    #-- Pseudobulks to generate
    types <- names(scbench[["mixtures"]])
    pseudobulks <- list()
    for (type in types) {
        if(type == "spillover") {
            ncells_run <- min(table(scbench$ref_data@meta.data[,level]))
        } else {
            ncells_run <- ncells
        }
        pb_res <- .get_pseudobulk(scbench, type,
                                  level = level, seed = seed,
                                  ncores = ncores, ncells = ncells_run,
                                  by_batch = by_batch)
        if(type == "population") {
            scbench[["mixtures"]][[type]][[level]] <- pb_res$props
        }
        scbench[["pseudobulk_counts"]][[type]] <- pb_res$pseudobulk
        rm(pb_res)
        gc()

    }
    if(shrink_obj) {
        scbench$ref_data <- scbench$ref_data@meta.data
        gc()
    }
    scbench$shrunk <- shrink_obj
    scbench$status <- .update_status(scbench$status, "pseudobulks")
    scbench$level <- level

    return(scbench)
}
#' @export
deconvolution_methods <- function() {
    c("ols", "dwls", "svr", "cibersortx", "music",  "bayesprism", "bisque")
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
#'
#' @export
deconvolute_all <- function(scbench,
                            scref,
                            methods = deconvolution_methods(),
                            ...){
    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(class(scref) == "screference")
    assert(all(methods %in% deconvolution_methods()))

    #-- Get types
    types <- names(scbench$pseudobulk_counts)
    #-- Run deconvolutions
    for(type in types) {
        message("========= Deconvoluting pseudobulks for ", type, " analysis ==========")
        for(method in methods) {
            if(method == "cibersortx") {
                scbench <- deconvolute.scbench(scbench, scref, method = method, type = type, ...)
            } else {
                scbench <- deconvolute.scbench(scbench, scref, method = method, type = type)
            }
        }
    }
    return(scbench)

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
#' @export
deconvolute.scbench <- function(scbench,
                        scref,
                        method = deconvolution_methods()[1],
                        type = c("population", "spillover", "lod")[1],
                        pseudobulk_norm = c("rpm", "none", "proportional_fitting")[3],
                        ...) {
    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(class(scref) == "screference")
    assert(!is.null(scbench$pseudobulk_counts[[type]]))

    #-- Cache handling
    method_cache <- filePath(scbench$cache, scbench$project_name, scbench$level, method, type)
    deconv_res <- .scbench_cache_check(method_cache)
    if(!is.null(deconv_res)) {
        message("Found cached results. Returning...")
        scbench[["deconvolution"]][[type]][[method]] <- deconv_res
        return(scbench)
    }
    #-- Get pseudocounts and make linear transformation
    pb <- scbench$pseudobulk_counts[[type]]
    if(pseudobulk_norm == "rpm") {
        data <- librarySizeNormalization(pb, 10^6)
    } else if (pseudobulk_norm == "proportional_fitting") {
        data <- librarySizeNormalization(pb, mean(colSums(pb)))
    }
    #-- Compute
    if(method == "cibersortx") {
        message("Running CIBERSORTx...")
        assert(!is.null(scref$cached_results[["cibersortx"]]))
        deconv_res <- cibersortx_deconvolute(data, scref, cache_path = method_cache, ...)
    } else if (method == "music") {
        message("Running MuSiC...")
        deconv_res <- music_deconvolute(data, scref, ...)
    } else if (method == "dwls") {
        message("Running DWLS...")
        assert(!is.null(scref$cached_results[["dwls"]]))
        deconv_res <- dwls_deconvolute(data, scref, ...)
    } else if(method == "ols") {
        assert(!is.null(scref$cached_results[["dwls"]]))
        message("Running OLS using the DWLS signature matrix...")
        deconv_res <- ols_deconvolute(data, scref, ...)
    } else if(method == "svr") {
        assert(!is.null(scref$cached_results[["dwls"]]))
        message("Running SVR using the DWLS signature matrix...")
        deconv_res <- svr_deconvolute(data, scref, ...)
    } else if(method == "bayesprism") {
        message("Running BayesPrism...")
        deconv_res <- bayesprism_deconvolute(data, scref, cache_path = method_cache, ...)
    } else if(method == "bisque") {
        message("Running Bisque...")
        deconv_res <- bisque_deconvolute(data, scref)
    }

    #-- Save to cache
    saveRDS(deconv_res, file = filePath(method_cache, "deconv_res.RDS"))
    #-- Return
    scbench[["deconvolution"]][[type]][[method]] <- deconv_res
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
#' @return an object of class `scbench`
#'
#' @export
deconvolute.matrix <- function(bulk_data, scref,
                               method = deconvolution_methods()[1],
                               bulk_norm = c("rpm", "none", "proportional_fitting")[3],
                               ...) {
    assert(class(scref) == "screference")
    #-- Compute
    if(method == "cibersortx") {
        message("Running CIBERSORTx...")
        assert(!is.null(scref$cached_results[["cibersortx"]]))
        deconv_res <- cibersortx_deconvolute(data, scref, cache_path = method_cache, ...)
    } else if (method == "music") {
        message("Running MuSiC...")
        deconv_res <- music_deconvolute(data, scref, ...)
    } else if (method == "dwls") {
        message("Running DWLS...")
        assert(!is.null(scref$cached_results[["dwls"]]))
        deconv_res <- dwls_deconvolute(data, scref, ...)
    } else if(method == "ols") {
        assert(!is.null(scref$cached_results[["dwls"]]))
        message("Running OLS using the DWLS signature matrix...")
        deconv_res <- ols_deconvolute(data, scref, ...)
    } else if(method == "svr") {
        assert(!is.null(scref$cached_results[["dwls"]]))
        message("Running SVR using the DWLS signature matrix...")
        deconv_res <- svr_deconvolute(data, scref, ...)
    } else if(method == "bayesprism") {
        message("Running BayesPrism...")
        deconv_res <- bayesprism_deconvolute(data, scref, cache_path = method_cache, ...)
    } else if(method == "bisque") {
        message("Running Bisque...")
        deconv_res <- bisque_deconvolute(data, scref)
    }
    return(deconv_res)
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
get_benchmark_results <- function(scbench, methods = NULL, type = "population") {
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench[["mixtures"]]))
    assert(!is.null(scbench[["deconvolution"]]))

    if(is.null(methods)) {
        methods <- names(scbench[["deconvolution"]][[type]])
    }

    if(type == "population") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[["population"]][[method]]
            bench_tb <- .get_benchmark_fractions(scbench, deconv_res)
            results <- .get_bench_results(deconv_res, bench_tb) %>%
                mutate(method = method)
        }) %>%
            bind_rows()
    } else if(type == "spillover") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[["spillover"]][[method]]
            spillover_tb <- .get_spillover_mixtures(scbench, deconv_res)
            results <- .get_spillover_results(deconv_res, spillover_tb) %>%
                mutate(mixture = paste0(pop1, "|", pop2),
                       method = method) %>%
                select(sample, mixture, pop1, pop2, pop1_truth, pop1_pred = pop1_estimate, method)
        }) %>%
            bind_rows()
    } else if(type == "lod") {
        results <- lapply(methods, function(method) {
            deconv_res <- scbench$deconvolution[["lod"]][[method]]
            lod_tb <- .get_lod_mixtures(scbench, deconv_res)
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
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#' @importFrom ggdensity geom_hdr_points
#' @importFrom ggheatmapper theme_scatter
#'
#' @return a ggplot object
#'
#' @export
plt_cors_scatter <- function(scbench, method) {
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution$population[[method]]))
    #-- Get deconv results
    deconv_res <- scbench$deconvolution[["population"]][[method]]

    #-- Get benchmark fractions
    bench_tb <- .get_benchmark_fractions(scbench, deconv_res)

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
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `cor_table`, a
#' table summarizing results
#'
#' @importFrom ggheatmapper get_hmPlot update_hmPlot
#'
#' @export
plt_cor_heatmap <- function(scbench) {
    plot_tb <- .aggregate_deconvolution_results(scbench) %>%
        group_by(method, population) %>%
        summarize(pop_cors = cor(pred, truth)) %>%
        pivot_wider(names_from = method, values_from = pop_cors)

    mean_cors <- colMeans(plot_tb[,2:ncol(plot_tb)])
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
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `rmse_table`, a
#' table summarizing results
#'
#' @importFrom ggheatmapper ggheatmap align_to_hm theme_sparse2
#' @importFrom yardstick rmse_vec
#' @export
plt_rmse_heatmap <- function(scbench) {
    #-- Get plot data
    rmse_tb <- .aggregate_deconvolution_results(scbench) %>%
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
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#'
#' @return a ggplot object
#'
#' @export
plt_spillover_scatter <- function(scbench, method) {
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution$spillover[[method]]))

    #-- Get deconv estimations
    deconv_res <- scbench$deconvolution$spillover[[method]]
    #-- Get mixtures
    spillover_tb <- .get_spillover_mixtures(scbench, deconv_res)

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
        labs(x = paste0("Fraction predicted by ", unique(deconv_res$method)), y = "True fraction") +
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
#'
#' @import tidyverse
#' @importFrom ggpubr stat_cor
#'
#' @return a ggplot object
#'
#' @export
plt_lod_scatter <- function(scbench, method) {
    #-- Limit of detection
    deconv_res <- scbench$deconvolution$lod[[method]]
    level <- .match_level(scbench, deconv_res, "lod")

    #-- Get mixtures
    bench_tb <- .get_lod_mixtures(scbench, deconv_res)
    lod_res <- .get_lod_results(deconv_res, bench_tb)

    plt <- ggplot(lod_res$lod_tb, aes(est, truth)) +
        facet_wrap(~ population) +
        geom_point(color = "steelblue") +
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
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `rmse_table`, a
#' table summarizing results
#'
#' @importFrom data.table rbindlist
#' @export
plt_spillover_heatmap <- function(scbench) {
    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution$spillover))
    methods <- names(scbench$deconvolution$spillover)

    #-- Get deconv estimations
    all_deconv_res <- lapply(methods, function(method) {
        scbench$deconvolution$spillover[[method]]
    })
    level <- .match_level(scbench, all_deconv_res[[1]], "spillover")
    #-- Get mixtures
    spillover_tb <- .get_spillover_mixtures(scbench, all_deconv_res[[1]])

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
#'
#' @return a list, containing `heatmap`, a ggheatmap object, and `lod_table`, a
#' table summarizing results
#'
#' @export
plt_lod_heatmap <- function(scbench) {
    #-- Error handling
    assert(class(scbench) == "scbench")
    assert(!is.null(scbench$deconvolution$lod))
    methods <- names(scbench$deconvolution$lod)

    #-- Setup
    all_deconv_res <- lapply(methods, function(method) { scbench$deconvolution$lod[[method]]})
    level <- .match_level(scbench, all_deconv_res[[1]], "lod")

    #-- Get mixtures and estimations, calculate limits of detection
    bench_tb <- .get_lod_mixtures(scbench, all_deconv_res[[1]])
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

# Helpers ----
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
    cat(str_glue("scbench object named {x$project_name} with {ncol(x$ref_data)} reference cells and {x$nlevels} levels of annotation"))
    cat("\n")
    cat(x$status)
    cat("\n")
    if(!is.null(x[["deconvolution"]])) {
        cat("deconvolution results: ")
        cat(paste0(names(x[["deconvolution"]]), collapse = ", "))
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
.aggregate_deconvolution_results <- function(scbench) {
    pop_props <- scbench[["mixtures"]][["population"]]
    #-- Get deconv results
    deconv_res <- bind_rows(scbench$deconvolution$population)

    #-- Get benchmark fractions
    level <- .match_level(scbench, deconv_res, "population")
    bench_tb <- pop_props[[level]] %>%
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
.get_benchmark_fractions <- function(scbench, deconv_res) {
    level <- .match_level(scbench, deconv_res, "population")
    bench_tb <- scbench[["mixtures"]][["population"]][[level]] %>%
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
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
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
            sample_cells <- lapply(1:length(sample_ncells), function(j) {
                ps_n <- as.integer(sample_ncells[j])
                pop_cells <- unlist(pops_cells[names(sample_ncells)[j]])
                if(length(pop_cells) < ps_n) {
                    warning("--- Some cell profiles had to be repeated... `ncells` might be set too high for this dataset.")
                }
                sample(pop_cells, ps_n)
            }) %>% reduce(c)
            rowSums(scbench$ref_data[,sample_cells]@assays$RNA@counts)
        }, mc.cores = ncores)
        message("-- Joining matrix...")
    }
    pb_mat <- reduce(pseudobulks, cbind)
    colnames(pb_mat) <- paste0("s", 1:ncol(pb_mat))
    pb_res <- list(pseudobulk = pb_mat, props = pb_props)
    gc()
    saveRDS(pb_res, file = cache_loc)
    return(pb_res)
}

#' @import tidyverse
#' @importFrom data.table as.data.table
.get_spillover_mixtures <- function(scbench, deconv_res) {
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
.get_lod_mixtures <- function(scbench, deconv_res) {
    level <- .match_level(scbench, deconv_res, "lod")
    bench_tb <- scbench$mixtures$lod[[level]] %>% as.data.table()
    bench_tb[,truth := get(population), by = population]
    bench_tb <- bench_tb %>%
        mutate(sample = paste0("s", 1:n())) %>%
        select(sample, population, truth) %>%
        as_tibble()
}
#' @import tidyverse
#' @importFrom data.table as.data.table
.get_lod_results <- function(deconv_res, bench_tb) {
    lod_tb <- bench_tb %>%
        left_join(deconv_res, by = "sample") %>%
        mutate(pop_name = paste0("frac_", population)) %>%
        as.data.table()
    lod_tb[,est := get(pop_name), by = pop_name]
    lod_tb <- lod_tb %>%
        select(sample, population, method, truth, est) %>%
        as_tibble()
    est_at_zero <- lod_tb %>%
        group_by(population) %>%
        filter(truth == 0) %>%
        select(population, min_est = est)

    lod_annot <- lod_tb %>%
        group_by(population) %>%
        left_join(est_at_zero, by = "population") %>%
        filter(truth >= 0.005, est >= max(0.005, min_est)) %>%
        slice_min(truth) %>%
        mutate(label = paste0("LoD = ", signif(truth, 3)))

    return(list(lod_tb = lod_tb, lod_annot = lod_annot))
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
