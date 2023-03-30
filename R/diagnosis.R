#' Evaluate annotation labels using Seurat markers
#'
#' Identifies differentially expressed markers between clusters and evaluates their performance.
#'
#' @param seurat_obj A Seurat object containing single-cell gene expression data
#' @param annotation A character string indicating the name of the column in \code{seurat_obj@meta.data} containing the annotation of interest.
#' @param min_FC Minimum fold change required for a marker to pass the diagnosis. The default value is 5.
#' @param ideal_FC Ideal fold change required for a marker to be considered outstanding. The default value is 20.
#' @param ncores Number of cores to use in parallelization for the \code{FindAllMarkers} function. The default value is 4.
#'
#' @return A list of three items:
#'  * top_mks: data.frame with top 10 markers for each cluster, with their
#'  average fold change and p-value.
#'  * full_mks: data.frame with all markers tested, with their average fold
#'  change, p-value, and other metrics.
#'  * mk_status: named vector with the status of each cluster based on the
#'  expression of their top markers.
#'
#'
#' @importFrom Seurat FindAllMarkers Idents
#' @importFrom dplyr group_by slice_max mutate relocate filter pull arrange
#' @importFrom future plan
#'
#' @export
test_markers <- function(seurat_obj, annotation,
                         min_FC = 5, ideal_FC = 20, ncores = 4) {

    ## Get markers for annotations ------
    plan("multisession", workers = ncores)
    Idents(seurat_obj) <- seurat_obj@meta.data[,annotation]
    mks <- FindAllMarkers(seurat_obj)

    top10_mks <- mks %>%
        group_by(cluster) %>%
        slice_max(avg_log2FC, n = 10) %>%
        mutate(avg_FC = 2^avg_log2FC,
               diff_mkExp = pct.1 - pct.2) %>%
        relocate(avg_FC, .before = avg_log2FC) %>%
        relocate(gene) %>%
        arrange(cluster, desc(avg_FC))

    ## Diagnosis ---------------
    labels <- as.character(unique(top10_mks$cluster))
    message("Running diagnosis...")
    mk_status <- c()
    for(label in labels) {
        label_mks <- top10_mks %>%
            filter(cluster == label)
        cl_mks <- pull(label_mks, avg_FC, gene)
        n_pass <- sum(cl_mks > min_FC)
        n_ideal <- sum(cl_mks > ideal_FC)

        cl_diffmkexp <- pull(label_mks, diff_mkExp, gene)[1:5]

        ## Annotation status --------
        if(n_pass == 10 & n_ideal > 1) {
            status <- "Outstanding"
        } else if(n_pass > 5 & n_ideal >= 1) {
            status <- "Good"
        } else if (n_pass > 3) {
            status <- "Pass"
        } else if(n_pass > 1) {
            status <- "Marginal pass"
        } else {
            status <- "Fail"
        }
        mk_status[label] <- status

        message("===================== ", label, " ======================")
        message("Top markers: ", paste0(names(cl_mks), collapse = ", "))
        message("------------ Fraction of cells expressing marker -----------")
        message("-- Interpretation: (Expressed equally in `", label, "` and other cells) 0 -> 1 (Only expressed in `", label, "`)")
        message(paste0(paste0(names(cl_diffmkexp), ": ", cl_diffmkexp), collapse = ", "))
        message("----------------------- Fold change ------------------------")
        message("Top fold-change: ", signif(max(cl_mks), 4))
        message("n(FC > ", min_FC, "): ", n_pass)
        message("n(FC > ", ideal_FC, "): ", n_ideal)
        message("Status: ", status)

    }
    return(list(top_mks = top10_mks %>% select(-diff_mkExp),
                full_mks = mks,
                mk_status = mk_status))
}
