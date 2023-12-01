#' @import reticulate
.install_scaden <- function() {
    .setup_deconv_conda()
    packages <- py_list_packages()
    if("scaden" %in% packages$package) {
        message("-- scaden found")
    } else {
        message("-- Installing scaden...")
        conda_install(envname = "deconverse", packages = "scaden", pip = TRUE)
    }

}
#' @import reticulate
.install_tensorflow <- function() {
    packages <- py_list_packages()
    if("tensorflow-gpu" %in% packages$package) {
        message("-- tensorflow-gpu found")
    } else {
        message("-- Installing tensorflow-gpu...")
        conda_install("deconverse", "tensorflow-gpu")
    }
}

#' Compute reference matrix from an `screference` object using scaden (Single-cell assisted deconvolutional network)
#'
#' @param scref an object of class `screference`
#' @param cache_path the path to the directory where the intermediate files
#' will be stored (default is "scaden")
#' @param bulk_mat a matrix containing the bulk gene expression data.
#' If NULL, a pseudobulk will be created from the reference.
#' @param force_retrain If TRUE, the model will be retrained even if a
#' trained model already exists (default is FALSE)
#' @param n_cells_sim the number of cells to simulate (default is 100)
#' @param n_samples_sim the number of samples to simulate (default is 1000)
#' @param batch_size the batch size used during training (default is 128)
#' @param learning_rate the learning rate used during training (default is 0.0001)
#' @param steps the number of training steps (default is 1000)
#' @param seed the random seed used during training (default is 0)
#' @param gpu if TRUE, installs Tensorflow and uses GPU acceleration (default is FALSE)
#'
#' @return The path to the directory containing the trained model
#' @note Reference: Menden, Kevin, Mohamed Marouf, Sergio Oller et al., 2020.
#' “Deep Learning–Based Cell Composition Analysis from Tissue Expression Profiles.”
#' Science Advances 6 (30): eaba2619. https://doi.org/10.1126/sciadv.aba2619.
#'
#' See also: https://github.com/theislab/AutoGeneS
#'
#'
#' @import tidyverse reticulate
#' @importFrom data.table fwrite
#' @export
scaden_scref <- function(scref,
                         cache_path = "scaden",
                         bulk_mat = NULL,
                         force_retrain = FALSE,
                         n_cells_sim = 100,
                         n_samples_sim = 1000,
                         batch_size = 128,
                         learning_rate = 0.0001,
                         steps = 1000,
                         seed = 0,
                         gpu = FALSE) {
    .install_reticulate()
    .install_scaden()
    if(gpu) .install_tensorflow()
    dir.create(cache_path, showWarnings = FALSE, recursive = TRUE)

    ## Make a fake "bulk" to run `scaden process` (ugly but works!) / or use the real bulk ----
    feats <- .write_bulk_data(scref, bulk_mat, cache_path = cache_path)

    # Write reference data -------
    data.table::fwrite(t(as.matrix(scref$seurat_obj@assays$RNA@counts))[,feats],
           file = str_glue("{cache_path}/data_counts.txt"),
           sep = "\t", row.names = FALSE)
    tibble(Celltype = scref$seurat_obj$annot_id) %>% write_tsv(str_glue("{cache_path}/data_celltypes.txt"))

    ## Run pipeline ------------
    model_path <- str_glue("{cache_path}/model")
    if(!dir.exists(model_path) | force_retrain) {
        if(dir.exists(model_path)) {
            system(str_glue("rm -r {model_path}"))
        }
        use_condaenv("deconverse", required = TRUE)
        sp <- import("subprocess")
        scaden_simulate <- str_glue('cd {cache_path}; scaden simulate --data . -n {n_samples_sim} -c {n_cells_sim} --pattern "*_counts.txt"')
        system(scaden_simulate)
        scaden_process <- str_glue("cd {cache_path}; scaden process data.h5ad bulk_ex.txt")
        system(scaden_process)
        scaden_train <- str_glue("cd {cache_path}; scaden train processed.h5ad --model_dir model --batch_size {batch_size} --learning_rate {learning_rate} --steps {steps} --seed {seed}")
        system(scaden_train)
    } else {
        message("-- Trained model found in cache")
    }

    return(model_path)
}

#' scaden (Single-cell assisted deconvolutional network) deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param cache_path path to cache the intermediate results
#' @param retrain if TRUE, model is re-processed with the `bulk_data` reference and retrained
#' @param ... passed to scaden_scref if retrain = TRUE
#'
#' @param a tibble with deconvolution fractions
#'
#' @import tidyverse reticulate
#' @importFrom data.table fwrite fread
#' @importFrom R.utils filePath getAbsolutePath
#' @note Reference: Menden, Kevin, Mohamed Marouf, Sergio Oller et al., 2020.
#' “Deep Learning–Based Cell Composition Analysis from Tissue Expression Profiles.”
#' Science Advances 6 (30): eaba2619. https://doi.org/10.1126/sciadv.aba2619.
#'
#' See also: https://github.com/theislab/AutoGeneS
#' @export
scaden_deconvolute <- function(bulk_data, scref,
                               cache_path = "scaden_results",
                               retrain = FALSE,
                               ...) {
    assert(class(scref) == "screference")
    is_cached <- "scaden" %in% names(scref$cached_results)
    assert(is.matrix(bulk_data))
    .install_scaden()

    # Train or retrain model -----------
    model_cache_path <- getAbsolutePath(filePath(scref$cache, scref$project_name, "scaden"))
    if(retrain | !is_cached) {
        message("-- Retraining model...")
        scaden_scref(scref,
                     cache_path = cache_path,
                     bulk_mat = bulk_data,
                     force_retrain = TRUE,
                     ...)
    }
    # Write bulk_data --------
    dir.create(cache_path, showWarnings = FALSE, recursive = TRUE)
    train_data <- fread(str_glue("{model_cache_path}/bulk_ex.txt"))
    df <- bulk_data %>% as.data.frame()
    df[,"symbol"] <- rownames(df)
    df %>%
        relocate("symbol") %>%
        filter(symbol %in% train_data$symbol) %>%
        fwrite(file = str_glue("{cache_path}/bulk_to_deconv.txt"), sep = "\t",
               row.names = FALSE, col.names = TRUE)

    # Predict ------------
    use_condaenv("deconverse", required = TRUE)
    sp <- import("subprocess")
    scaden_predict <- str_glue("cd {cache_path}; scaden predict --model_dir {model_cache_path}/model bulk_to_deconv.txt")
    system(scaden_predict)
    scaden_prop <- read_tsv(str_glue("{cache_path}/scaden_predictions.txt"),
                            show_col_types = FALSE) %>%
        dplyr::rename(sample = `...1`) %>%
        rename_with(~ paste0("frac_", .), .cols = -sample) %>%
        mutate(method = "scaden") %>%
        tibble()

    return(scaden_prop)
}

.write_bulk_data <- function(scref, bulk_mat, cache_path = cache_path) {
    if(is.null(bulk_mat)) {
        set.seed(0)
        sampled_cells <- list(sample(1:ncol(scref$seurat_obj), size = 100, replace = TRUE),
                              sample(1:ncol(scref$seurat_obj), size = 100, replace = TRUE))

        pb_mat <- sapply(sampled_cells, function(cells) {
            rowMeans(as.matrix(scref$seurat_obj@assays$RNA@counts)[,cells]) }) %>%
            as.data.frame()
        pb_mat[,"symbol"] <- rownames(pb_mat)
        pb_mat <- relocate(pb_mat, "symbol")
        fwrite(pb_mat,
               file = str_glue("{cache_path}/bulk_ex.txt"), sep = "\t",
               col.names = TRUE, row.names = FALSE)
        features = pb_mat$symbol
    } else {
        features <- intersect(rownames(scref$seurat_obj@assays$RNA@counts), rownames(bulk_mat))
        bulk_mat2 <- bulk_mat[features,] %>%
            as.data.frame()
        bulk_mat2[,"symbol"] <- rownames(bulk_mat2)
        relocate(bulk_mat2, "symbol") %>%
            fwrite(file = str_glue("{cache_path}/bulk_ex.txt"), sep = "\t",
                   row.names = FALSE, col.names = TRUE)
    }
    return(features)
}
