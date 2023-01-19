#' Get Seurat PBMCs pre-processed example
.pp_PBMCs <- function() {
    #-- Get data
    download.file("https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                  "./pbmc3k_filtered_gene_bc_matrices.tar.gz")
    untar("pbmc3k_filtered_gene_bc_matrices.tar.gz")
    pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
    #-- Run Seurat pp
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    pbmc <- NormalizeData(pbmc)
    pbmc <- ScaleData(pbmc)
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    pbmc <- FindNeighbors(pbmc, dims = 1:10)
    pbmc <- FindClusters(pbmc, resolution = 0.5)
    pbmc <- RunUMAP(pbmc, dims = 1:10)
    # DimPlot(pbmc, reduction = "umap")
    #-- Annotate
    pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    pbmc.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)
    new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                         "NK", "DC", "Platelet")
    names(new.cluster.ids) <- levels(pbmc)
    pbmc <- RenameIdents(pbmc, new.cluster.ids)
    pbmc$Cell_minor_identities <- Idents(pbmc)
    table(pbmc$Cell_minor_identities)
    cell_identity_hierarchy <- c("TNK" = "Naive CD4 T", "TNK" = "Memory CD4 T", "TNK" = "CD8 T",
      "TNK" = "NK", "B" = "B", "Monocytic_lineage" = "CD14+ Mono",
      "Monocytic_lineage" = "FCGR3A+ Mono", "Monocytic_lineage" = "DC",
      "Platelet" = "Platelet")
    pbmc$Cell_major_identities <- fct_recode(pbmc$Cell_minor_identities,
               "TNK" = "Naive CD4 T", "TNK" = "Memory CD4 T", "TNK" = "CD8 T",
               "TNK" = "NK", "B" = "B", "Monocytic_lineage" = "CD14+ Mono",
               "Monocytic_lineage" = "FCGR3A+ Mono", "Monocytic_lineage" = "DC",
               "Platelet" = "Platelet")
    #-- Remove platelets (too few cells)
    pbmc <- pbmc[,!pbmc$Cell_major_identities == "Platelet"]
    pbmc$Cell_major_identities <- as.character(pbmc$Cell_major_identities)
    pbmc$Cell_minor_identities <- as.character(pbmc$Cell_minor_identities)
    dir.create("data", showWarnings = FALSE)
    usethis::use_data(pbmc, overwrite = TRUE)
    return(TRUE)
}

#' Seurat example PBCM dataset
#'
#' Cells annotated with the example pipeline from the Seurat package vignette,
#' already preprocessed, and with Platelet cluster removed (too few cells)
#'
#' @format `pbmc`
#' A Seurat Object.
#'
#' @source <https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>
"pbmc"

