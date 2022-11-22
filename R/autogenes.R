#' @import reticulate
.setup_deconv_conda <- function() {
    envs <- conda_list()
    if(!"deconverse" %in% envs$name) {
        message("Creating `deconverse` conda environment...")
        conda_create("deconverse",
                     python_version = "3.7",
                     packages = c("numpy", "pandas", "anndata",
                                  "scipy", "scikit-learn", "matplotlib"))
    } else {
        message("Using `deconverse` conda environment...")
    }
}
#' @import reticulate
.install_autogenes <- function() {
    .setup_deconv_conda()
    use_condaenv("deconverse")
    packages <- py_list_packages()
    if("autogenes" %in% packages$package) {
        message("-- AutoGeneS found")
    } else {
        message("-- Installing AutoGeneS...")
        conda_install("deconverse", "autogenes", pip = TRUE)
    }

}

#' Compute reference matrix from an `screference` object using
#' AutoGeneS
#'
#' @param scref an object of class `screference`
#' @param ngen an integer, number of generations for optimization
#' @param mode one of `"standard"` or `"fixed"`. If "fixed", `"nfeatures"` is used.
#' @param nfeatures integer, number of genes selected in fixed mode.
#' @param seed int, for reproducibility
#' @param population_size int, size of each generation (mu parameter)
#' @param offspring_size int, number of individuals per generation (lambda parameter)
#' @param crossover_pb float, crossover probability
#' @param mutation_pb float, mutation probability
#' @param crossover_thres int, crossover threshold (standard mode)
#' @param ind_standard_pb float, probability used to generate initial population in
#' standard mode
#' @param ... passed to [Seurat::FindVariableFeatures]
#'
#' @return an screference object updated with AutoGeneS reference centroids
#' @note Reference: Aliee, Hananeh, and Fabian J. Theis. 2021.
#' “AutoGeneS: Automatic Gene Selection Using Multi-Objective Optimization
#' for RNA-Seq Deconvolution.” Cell Systems 12 (7): 706-715.e4.
#' https://doi.org/10.1016/j.cels.2021.05.006.
#'
#' See also: https://github.com/theislab/AutoGeneS
#'
#' @import tidyverse reticulate
#' @importFrom Seurat SetIdent NormalizeData FindVariableFeatures VariableFeatures
#' @importFrom Matrix rowMeans
#' @export
autogenes_scref <- function(scref, ngen = 500, mode = "standard",
                            nfeatures = 2000,
                            seed = 0, population_size = 100,
                            offspring_size = 50,
                            crossover_pb = 0.7, mutation_pb = 0.3,
                            crossover_thresh = 1000,
                            int_standard_pb = 0.1, ...) {
    scdata <- scref$seurat_obj
    id_name <- "annot_id"

    #-- Pre-process
    scdata <- SetIdent(scdata, value = id_name)
    scdata <- NormalizeData(scdata)
    scdata <- FindVariableFeatures(scdata, ...)

    #-- Get centroids
    var_feats <- VariableFeatures(scdata)
    cells_in_pops <- split(colnames(scdata), scdata$annot_id)
    centroids <- sapply(cells_in_pops, function(cells) {
        rowMeans(scdata@assays$RNA@data[var_feats,cells])
    })

    #-- Centroid optimization
    .install_autogenes()
    message("Optimizing markers...")
    ag <- import("autogenes")
    centroid_df <- centroids %>% as.data.frame() %>% r_to_py()
    ag$init(centroid_df$T)
    if(mode == "standard") {
        ag$optimize(mode = mode,
                    seed = as.integer(seed), population_size = as.integer(population_size),
                    offspring_size = as.integer(offspring_size),
                    crossover_pb = crossover_pb, mutation_pb = mutation_pb,
                    crossover_thresh = as.integer(crossover_thresh),
                    int_standard_pb = int_standard_pb, verbose = FALSE)
    } else {
        ag$optimize(mode = mode,
                    nfeatures = as.integer(nfeatures),
                    seed = as.integer(seed), population_size = as.integer(population_size),
                    offspring_size = as.integer(offspring_size),
                    crossover_pb = crossover_pb, mutation_pb = mutation_pb,
                    crossover_thresh = as.integer(crossover_thresh),
                    int_standard_pb = int_standard_pb, verbose = FALSE)
    }
    index <- ag$select()
    centroids_sc_pareto <- centroid_df[index] %>% py_to_r()
    centroids_sc_pareto <- centroids_sc_pareto[,scref$populations]

    #-- Saving
    return(centroids_sc_pareto)
}

#' AutoGeneS deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param model a model to run AutoGeneS. One of: `'nusvr'`, `'nnls'` or `'linear'`
#'
#' @return  a tibble with deconvolution fractions
#' @note Reference: Aliee, Hananeh, and Fabian J. Theis. 2021.
#' “AutoGeneS: Automatic Gene Selection Using Multi-Objective Optimization
#' for RNA-Seq Deconvolution.” Cell Systems 12 (7): 706-715.e4.
#' https://doi.org/10.1016/j.cels.2021.05.006.
#'
#' See also: https://github.com/theislab/AutoGeneS
#'
#' @import tidyverse reticulate
autogenes_deconvolute <- function(bulk_data, scref,
                                  model = c("nusvr", "nnls", "linear")[1]) {
    assert(class(scref) == "screference")
    assert("autogenes" %in% names(scref$cached_results))
    assert(is.matrix(bulk_data))

    .install_autogenes()
    centroids <- scref$cached_results$autogenes %>% r_to_py()

    #-- Regression
    ag <- import("autogenes")
    ag$init(centroids$T)
    ag$optimize(mode = "fixed", nfeatures = as.integer(nrow(centroids)), ngen = as.integer(1))
    ag$select()
    message("Performing regression...")
    bulk_data_py <- bulk_data %>% t() %>% as.data.frame() %>% r_to_py()
    coefs <- ag$deconvolve(bulk_data_py, model = model) %>% as.matrix()

    #-- Normalize proportions
    props <- .autogenes_normalize_props(coefs)
    colnames(props) <- colnames(scref$cached_results$autogenes)
    rownames(props) <- colnames(bulk_data)

    autogene_prop <- props %>%
        as.data.frame() %>%
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        mutate(method = "AutoGeneS") %>%
        tibble()

    return(autogene_prop)
}

## Helpers ------------------------------
.autogenes_normalize_props <- function(data) {
    data[data < 0] <- 0
    norm_data <- apply(data, 1, function(x) { x/sum(x) }) %>% t()
    return(norm_data)
}

