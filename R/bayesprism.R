#' Compute reference matrix from an `screference` object using
#' BayesPrism
#'
#' @param scref an object of class `screference`
#' @param cache_path path to cache the results
#' @param data_type a string, one of: `"10X"` or `"Smart-seq"`
#' @param malignant_pop_id the ID in the screference `annot_id` corresponding
#' to the malignant cell population for deconvolution of cancer samples. If none
#' of the populations are malignant, leave as `NULL`.
#' @param pval_cutoff a float, p-value cutoff from differential expression
#' analysis for selection of markers
#' @param logFC_cutoff a float, cutoff for log fold-change from differential
#' expression analysis between populations for selection of markers
#' @param n_cores number of cores used for computation
#'
#' @return an screference object updated with BayesPrism reference
#' @note Reference: Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene
#' expression deconvolution with BayesPrism enables Bayesian integrative
#' analysis across bulk and single-cell RNA sequencing in oncology.
#' Nat Cancer 3, 505–517 (2022). https://doi-org.insb.bib.cnrs.fr/10.1038/s43018-022-00356-3
#' See also: https://github.com/Danko-Lab/BayesPrism
#'
#' @import tidyverse
#' @importFrom BayesPrism cleanup.genes select.gene.type get.exp.stat select.marker
#' @import snowfall
#' @export
bayesprism_scref <- function(scref,
                             cache_path = "bayes_prism",
                             data_type = c("10X", "Smart-seq")[1],
                             malignant_pop_id = NULL,
                             pval_cutoff = 0.01, logFC_cutoff = 0.1,
                             ncores = parallel::detectCores()/2
) {
    .install_bayesprism()
    assert(class(scref) == "screference")
    if(!"batch_id" %in% colnames(scref$seurat_obj@meta.data)) {
        stop("BayesPrism needs `batch_id` to be defined when creating the scref object for deconvolution.")
    }
    #-- Clean matrix with BayesPrism suggestions
    sc_mat_filt <- t(as.matrix(scref$seurat_obj@assays$RNA@counts)) |>
        cleanup.genes(input.type = "count.matrix",
                      species = "hs",
                      gene.group = c("Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY"),
                      exp.cells = 5) |>
        select.gene.type(gene.type = "protein_coding")

    #-- Get annotations
    type_labels <- scref$seurat_obj@meta.data[,"annot_id"]
    batch_annot <- scref$seurat_obj@meta.data[,"batch_id"]
    state_labels <- type_labels

    #-- Run expression
    pseudo_count_param <- ifelse(data_type == "10X", 0.1, 10)

    diff_exp_stat <- get.exp.stat(sc_mat_filt,
                                  cell.type.labels = type_labels,
                                  cell.state.labels = state_labels,
                                  psuedo.count = pseudo_count_param,
                                  cell.count.cutoff = 1,
                                  n.cores = ncores
    )
    gc()
    pc_sig <- select.marker(sc_mat_filt,
                            stat = diff_exp_stat,
                            pval.max = pval_cutoff,
                            lfc.min = logFC_cutoff)
    bayesprism_res <- list(pc_sig = pc_sig,
                           type_labels = type_labels,
                           state_labels = state_labels,
                           maglignant_pop_id = malignant_pop_id)
    gc()
    return(bayesprism_res)
}

#' BayesPrism deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param cache_path path to cache the results
#' @param outlier_cut,outlier.fraction two floats used to filter genes in `bulk_data`
#' whose expression fraction is greater than outlier.cut in more than outlier.fraction.
#' Typically for dataset with reasonable quality control, very few genes will be
#' filtered. Removal of outlier genes will ensure that the inference will not
#' be dominated by outliers, which sometimes may be resulted from poor QC in mapping.
#' See: [BayesPrism::new.prism()]
#' @param pseudo_min float, the desired minimum value to replace zero
#' after normalization. See: [BayesPrism::new.prism()].
#' @param n_cores number of cores used for computation
#'
#' @return  a tibble with deconvolution fractions
#' @note Reference: Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene
#' expression deconvolution with BayesPrism enables Bayesian integrative
#' analysis across bulk and single-cell RNA sequencing in oncology.
#' Nat Cancer 3, 505–517 (2022). https://doi-org.insb.bib.cnrs.fr/10.1038/s43018-022-00356-3
#' See also: https://github.com/Danko-Lab/BayesPrism
#'
#' @import tidyverse
#' @importFrom BayesPrism new.prism run.prism get.fraction
#' @importFrom snow setDefaultClusterOptions
#'
#' @export
bayesprism_deconvolute <- function(bulk_data,
                                   scref,
                                   cache_path = NULL,
                                   outlier_cut = 0.01,
                                   outlier_fraction = 0.1,
                                   pseudo_min = 1e-8,
                                   ncores = parallel::detectCores()/2) {
    .install_bayesprism()
    if(!is.null(cache_path)) {
        cache_fname <- ""
    }
    else {
        cache_fname <- filePath(cache_path, "bayes_prism_res.RData")
    }
    if(file.exists(cache_fname)) {
        load(cache_fname)
        return(bp_tb_res)
    }
    #-- Checks
    assert(class(scref) == "screference")
    assert("bayesprism" %in% names(scref$cached_results))
    assert(is.matrix(bulk_data))

    #-- Fit BayesPrism object
    prism_obj <- new.prism(
        reference = scref$cached_results$bayesprism$pc_sig,
        mixture = t(bulk_data),
        input.type = "count.matrix",
        cell.type.labels = scref$cached_results$bayesprism$type_labels,
        cell.state.labels = scref$cached_results$bayesprism$state_labels,
        key =  scref$cached_results$bayesprism$malignant_pop_id,
        outlier.cut = outlier_cut,
        outlier.fraction = outlier_fraction,
        pseudo.min = pseudo_min
    )
    #-- Run BayesPrism analysis
    bp_res <- run.prism(prism = prism_obj, n.cores = ncores)

    #-- Conform output
    bp_tb_res <- get.fraction(bp_res, "final", "type") %>%
        as.data.frame() %>%
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        mutate(method = "BayesPrism") %>%
        tibble()

    if(cache_fname != "") save(bp_tb_res, bp_res, file = cache_fname)

    return(bp_tb_res)
}

#' @importFrom remotes install_github
.install_bayesprism <- function() {
    if(!"BayesPrism" %in% installed.packages()) {
        message("R package BayesPrism not detected. Installing...")
        install.packages("snowfall")
        install_github("Danko-Lab/BayesPrism/BayesPrism")
    }
}
