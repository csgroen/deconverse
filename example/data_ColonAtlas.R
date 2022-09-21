library(tidyverse)
library(Seurat)
#-- CRC Atlas so
get_crc_atlas <- function() {
    if(!file.exists("data/CRC_Atlas/GSE178341_Seurat.RData")) {
        dir.create("data/CRC_Atlas", showWarnings = FALSE, recursive = FALSE)
        if(!file.exists("data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit.h5")) {
            message("Downloading files...")
            download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178341&format=file&file=GSE178341%5Fcrc10x%5Ffull%5Fc295v4%5Fsubmit%2Eh5",
                          "data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit.h5")
            download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178341&format=file&file=GSE178341%5Fcrc10x%5Ffull%5Fc295v4%5Fsubmit%5Fcluster%2Ecsv%2Egz",
                          "data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz")
            download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE178341&format=file&file=GSE178341%5Fcrc10x%5Ffull%5Fc295v4%5Fsubmit%5Fmetatables%2Ecsv%2Egz",
                          "data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
        }
        message("Creating Seurat object...")
        #-- Read raw data
        so <- Read10X_h5("data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit.h5")
        cell_meta <- read_csv("data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz")
        cell_meta2 <- read_csv("data/CRC_Atlas/GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz")
        
        cell_meta <- left_join(cell_meta, cell_meta2, by = c("sampleID" = "cellID")) %>%
            as.data.frame() %>%
            column_to_rownames("sampleID")
        #-- Create SeuratObject
        so <- CreateSeuratObject(so,
                                 project = "CRC_Atlas",
                                 meta.data = cell_meta)
        save(so, file = "data/CRC_Atlas/GSE178341_Seurat.RData")
    } else {
        message("Seurat object found in cache...")
        load("data/CRC_Atlas/GSE178341_Seurat.RData")
    }
    gc()
    invisible(so)
}

get_crc_subset <- function(pids,
                           version = "v2") {
    fname <- str_glue("data/CRC_Atlas/GSE178341_tsubset_Seurat_{version}.RData")
    if(!file.exists(fname)) {
        so <- get_crc_atlas()
        so_small <- so[,so$PatientTypeID %in% pids]
        save(so_small, file = fname)
    } else {
        load(fname)
    }
    return(so_small)
}