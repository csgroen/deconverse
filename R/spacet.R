#' SpaCET deconvolution using an `screference`
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay
#' @param scref an `screference` object containing single-cell reference data
#' @param cancerType A character string specifying the type of cancer.
#' @param sc_lineage_tree An optional data frame representing the single-cell lineage tree. Default is NULL.
#' @param sc_nCellEachLineage An integer specifying the number of cells for each lineage. Default is 100.
#' @param ncores An integer specifying the number of cores to use for parallel processing. Default is \code{parallel::detectCores() - 1}.
#'
#' @details
#' This function utilizes the SpaCET package to perform deconvolution.
#' It first creates a SpaCET object with spatial transcriptomics data, and then
#' performs deconvolution using matched single-cell RNA-seq data.
#'
#' For detailed documentation on the underlying SpaCET functions, refer to the original SpaCET package documentation:
#' - \code{\link{create.SpaCET.object}}
#' - \code{\link{SpaCET.deconvolution.matched.scRNAseq}}
#'
#' @import Seurat tidyverse
#'
#' @export
spacet_deconvolute <- function(spatial_obj, scref,
                               cancerType,
                               sc_lineage_tree = NULL,
                               sc_nCellEachLineage = 100,
                               ncores = parallel::detectCores()-1) {
    .install_spacet()

    SpaCET_obj <- create.SpaCET.object(
        counts=as.matrix(GetAssayData(spatial_obj, "Spatial", "counts")),
        spotCoordinates=GetTissueCoordinates(spatial_obj)[,1:2],
        imagePath=NA,
        platform = "Visium"
    )

    SpaCET_obj <- SpaCET.deconvolution.matched.scRNAseq(
        SpaCET_obj = SpaCET_obj,
        cancerType = cancerType,
        sc_counts = as.matrix(GetAssayData(scref$seurat_obj, "RNA", "counts")),
        sc_annotation = data.frame(cellID = colnames(scref$seurat_obj),
                                   bio_celltype = scref$seurat_obj$annot_id),
        sc_lineageTree = sc_lineage_tree,
        sc_nCellEachLineage = sc_nCellEachLineage,
        coreNo = ncores
    )

    deconv_props <- SpaCET_obj@results$deconvolution$propMat
    deconv_props <- apply(deconv_props, 2, function(row) { row/sum(row) })

    deconv_res <- deconv_props |>
        t() |>
        as.data.frame() |>
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        mutate(method = "SpaCET") %>%
        tibble()

    return(deconv_res)
}

#' @importFrom remotes install_bioc install_github
.install_spacet <- function() {
    if(! "SpaCET" %in% installed.packages()) {
        message("R package SpaCET is not detected. Installing...")
        install_bioc(c("BiRewire", "limma", "UCell"))
        install_github(c("JEFworks/MUDAN", "data2intelligence/SpaCET"))
    }
    require(SpaCET)
}

