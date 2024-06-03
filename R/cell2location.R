#' Compute Reference Signature Matrix from `screference` Object Using Cell2location
#'
#' This function computes a reference signature matrix from a `screference` object using
#' the method implemented in the Python `cell2location` package.
#'
#' @param scref An object of class `screference` containing the reference single-cell RNA-seq data.
#' @param cache_path A string specifying the path to save the computed reference results. Default is `NULL`.
#' @param cell_count_cutoff An integer specifying the cutoff for the number of cells expressing a gene. Default is 5.
#' @param cell_percentage_cutoff A numeric value specifying the cutoff for the percentage of cells expressing a gene. Default is 0.03.
#' @param nonz_mean_cutoff A numeric value specifying the cutoff for the mean expression of non-zero genes. Default is 1.2.
#' @param max_epochs An integer specifying the maximum number of epochs for model training. Default is 1000.
#' @param categorical_covariates A character vector of categorical covariates. Default is `NULL`.
#' @param continuous_covariates A character vector of continuous covariates. Default is `NULL`.
#'
#' @return A matrix of inferred average expression values per cluster.
#' @import reticulate
#' @importFrom stringr str_remove
#' @export
cell2location_scref <- function(scref,
                                cache_path = NULL,
                                cell_count_cutoff = 5,
                                cell_percentage_cutoff = 0.03,
                                nonz_mean_cutoff = 1.2,
                                max_epochs=1000,
                                categorical_covariates = NULL,
                                continuous_covariates=NULL) {
    # Setup conda -----
    .install_cell2location()

    # From tutorial: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html
    # Get reference data ------
    adata_ref <- seurat_to_adata(scref$seurat_obj)
    c2l <- import("cell2location")
    selected <- c2l$utils$filtering$filter_genes(adata_ref,
                                                 cell_count_cutoff=cell_count_cutoff,
                                                 cell_percentage_cutoff2=cell_percentage_cutoff,
                                                 nonz_mean_cutoff=nonz_mean_cutoff)
    adata_ref <- adata_ref[,selected]$copy()
    # Fit model ----
    c2l$models$RegressionModel$setup_anndata(adata=adata_ref,
                                             batch_key="batch_id",
                                             labels_key="annot_id",
                                             continuous_covariate_keys = r_to_py(categorical_covariates),
                                             categorical_covariate_keys = r_to_py(continuous_covariates))
    mod <- c2l$models$RegressionModel(adata_ref)
    mod$train(max_epochs=r_to_py(as.integer(max_epochs)), log_every_n_steps=as.integer(1))
    adata_ref2 <- mod$export_posterior(adata_ref,
                                       sample_kwargs=r_to_py(list(num_samples=as.integer(1000),
                                                                  batch_size=as.integer(2500))))
    # Results for deconvolution -----
    inf_aver <- adata_ref2$varm["means_per_cluster_mu_fg"]
    colnames(inf_aver) <- str_remove(colnames(inf_aver), "means_per_clustre_mu_fg_")

    if(!is.null(cache_path)) {
        dir.create(cache_path, showWarnings = FALSE)
        saveRDS(inf_aver, "reference_res.RDS")
    }

    return(inf_aver)
}

#' Cell2location deconvolution of spatial transcriptomics data from an 'screference'
#'
#' This function performs deconvolution of spatial transcriptomics data using
#' the `cell2location` method, leveraging reference signatures from a `screference` object.
#'
#' @param spatial_obj A Seurat object containing the spatial transcriptomics data.
#' @param scref A `screference` object containing cached `cell2location` reference signatures.
#' @param N_cells_per_location An integer specifying the number of cells per location. Default is 30.
#' @param detection_alpha An integer specifying the detection alpha parameter. Default is 20.
#' @param max_epochs An integer specifying the maximum number of epochs for model training. Default is 5000.
#'
#' @return A tibble containing the deconvoluted cell type fractions for each sample location.
#'
#' @import reticulate
#' @importFrom dplyr rename_with mutate
#' @importFrom tibble rownames_to_column tibble
#'
#' @export
cell2location_deconvolute <- function(spatial_obj, scref,
                                      cache_path = "cell2location",
                                      N_cells_per_location=30, detection_alpha=20,
                                      max_epochs=5e3) {
    .install_cell2location()
    # Check cache ---
    if(!is.null(cache_path)) {
        dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
        res_file <- filePath(cache_path, "deconverse_results.RDS")
    }

    if (file.exists(res_file)) {
        message("Results found in cache. Returning...")
        deconv_res <- readRDS(res_file)
        return(deconv_res)
    }
    # Check reference -----
    if(is.null(scref$cached_results[["cell2location"]])) {
        stop("scref doesn't have cached `cell2location` reference signatures. Please run `compute_reference` on the scref object.")
    }
    c2l <- import("cell2location")

    # Preparing spatial obj and reference ------
    adata_vis <- seurat_to_adata(spatial_obj, assay = "Spatial")
    reference_res <- scref$cached_results$cell2location
    genes2keep <- intersect(rownames(reference_res), rownames(adata_vis$var))
    adata_vis <- adata_vis[,genes2keep]$copy()
    reference_res <- reference_res[genes2keep,]

    # Preparing model ------
    c2l$models$Cell2location$setup_anndata(adata=adata_vis)
    mod2 <- c2l$models$Cell2location(adata_vis,
                                     cell_state_df=r_to_py(reference_res),
                                     N_cells_per_location=r_to_py(as.integer(N_cells_per_location)),
                                     detection_alpha=r_to_py(as.integer(detection_alpha)))
    # Training model ------
    mod2$train(max_epochs=as.integer(max_epochs), train_size=1, log_every_n_steps=as.integer(1))

    # Getting results (quartile 0.05) ------
    adata_vis = mod2$export_posterior(
        adata_vis, sample_kwargs=list('num_samples' = as.integer(1000),
                                      'batch_size'= as.integer(mod2$adata$n_obs))
    )
    q05 <- as.matrix(adata_vis$obsm$as_dict()$q05_cell_abundance_w_sf)
    colnames(q05) <- str_remove(colnames(q05), "q05cell_abundance_w_sf_means_per_cluster_mu_fg_")


    q05 <- q05/apply(q05,1,sum) # normalizing into fractions

    deconv_res <- q05 |>
        as.data.frame() |>
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        mutate(method = "cell2Location")

    if(!is.null(cache_path)) saveRDS(deconv_res, file = res_file)

    return(deconv_res)
}

#' @import reticulate
.install_cell2location <- function() {
    .setup_deconv_conda()
    packages <- py_list_packages()
    if("cell2location" %in% packages$package) {
        message("-- cell2Location found")
    } else {
        conda_install("deconverse", packages = c("scvi-tools==1.1.2", "git+https://github.com/BayraktarLab/cell2location.git#egg=cell2location[tutorials]"), pip = TRUE, forge=FALSE)
    }

}
