#' Bisque deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of bulk RNA-seq data with rows representing genes
#' and columns representing samples.
#' @param scref an object of class `screference`
#'
#' @return a tibble containing the fraction of each cell type in each bulk sample.
#'
#' @note Reference:
#' Jew, B., Alvarez, M., Rahmani, E. et al. Accurate estimation of cell
#' composition in bulk expression through robust integration of single-cell
#' information. Nat Commun 11, 1971 (2020). https://doi.org/10.1038/s41467-020-15816-6
#' See also: https://github.com/cozygene/bisque
#'
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame
#' @importFrom BisqueRNA ReferenceBasedDecomposition
#' @export
bisque_deconvolute <- function(bulk_data, scref) {
    bulk_eset <- ExpressionSet(assayData = bulk_data)
    cell_annot <- scref$seurat_obj@meta.data
    cell_annot[,"cell_id"] <- rownames(cell_annot)
    sc_eset <- ExpressionSet(as.matrix(scref$seurat_obj@assays$RNA@counts),
                             phenoData = AnnotatedDataFrame(cell_annot))

    res <- ReferenceBasedDecomposition(bulk_eset, sc_eset,
                                       cell.types = "annot_id",
                                       subject.names = "cell_id",
                                       use.overlap = FALSE)
    deconv_res <- res$bulk.props %>%
        t() %>%
        as.data.frame() %>%
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        mutate(method = "Bisque")
               # residuals_norm = res$rnorm)
    return(deconv_res)
}
