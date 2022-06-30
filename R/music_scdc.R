#' MuSiC deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#'
#' @param a tibble with deconvolution fractions and model metrics
#'
#'
#' @import MuSiC
#' @importFrom Biobase ExpressionSet
#' @export
music_deconvolute <- function(bulk_data, scref) {
    if(!"batch_id" %in% colnames(scref$seurat_obj@meta.data)) {
        stop("MuSiC needs `batch_id` to be defined when creating the scref object for deconvolution.")
    }
    est_prop <- music_prop(bulk.eset = ExpressionSet(assayData = bulk_data),
                           sc.eset = ExpressionSet(as.matrix(scref$seurat_obj@assays$RNA@counts),
                                                   phenoData =AnnotatedDataFrame(scref$seurat_obj@meta.data)),
                           clusters = "annot_id",
                           samples = "batch_id")
    gc()
    deconv_res <- est_prop$Est.prop.weighted %>%
        as.data.frame() %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        rename_with(.fn = ~ paste0("frac_", .), .cols = -matches("sample")) %>%
        mutate(
            model_Rsquared = est_prop$r.squared.full,
            method = "MuSiC")

    return(deconv_res)
}

# SCDC ---------------------------------
# TODO: Not complete due to memory constraints
#
# library(SCDC)
# # Very memory intensive: re-implement?
# scdc_deconvolute <- function(bulk_data, scref, batch_id = NULL) {
#     message("Running QC on reference...")
#     scdc_qc <- SCDC_qc_ONE(Biobase::ExpressionSet(assayData = as.matrix(scref$seurat_obj@assays$RNA@counts),
#                                                   phenoData = AnnotatedDataFrame(scref$seurat_obj@meta.data)),
#                            ct.varname = "annot_id",
#                            ct.sub = scref$populations,
#                            sample = batch_id)
#     gc()
#     message("Running deconvolution...")
#     est_prop <- SCDC_prop(bulk_data = ExpressionSet(assayData = bulk_data),
#                            sc.eset = scdc_qc$sc.eset.qc,
#                            clusters = "annot_id",
#                            samples = batch_id)
#
#     est_prop$Var.prop %>% tibble()
#     est_prop$Est.prop.allgene %>% tibble()
#
#     deconv_res <- est_prop$Est.prop.weighted %>%
#         as.data.frame() %>%
#         rownames_to_column("sample") %>%
#         tibble() %>%
#         rename_with(.fn = ~ paste0("frac_", .), .cols = -matches("sample")) %>%
#         mutate(
#             model_Rsquared = est_prop$r.squared.full,
#             method = "MuSiC")
#
#     return(deconv_res)
# }
#
