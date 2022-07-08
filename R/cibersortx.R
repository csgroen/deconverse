#' Installs CIBERSORTx container using docker or singularity
#'
#' @param container_type container software used. one of:
#' "docker" or "singularity"
#'
#' @importFrom purrr possibly
#' @export
install_cibersortx <- function(container_type = c("docker", "singularity")[1]) {
    if(container_type == "docker") {
        if(!.docker_check()) {
            stop("Error: docker is not installed or not available in PATH. See: https://docs.docker.com/get-docker/ or use singularity")
        }
        system("docker pull cibersortx/fractions")
    } else {
        #-- Check if singularity is available
        if(!.singularity_check()) {
            stop("Error: singularity is not installed or not available in PATH. See: https://docs.sylabs.io/guides/3.0/user-guide/installation.html or use docker.")
        }
        deconverse_path <- find.package("deconverse")[1]
        mv_cx <- paste0("mv fractions_latest.sif ", deconverse_path, "/cibersortx.sif")
        system("singularity pull docker://cibersortx/fractions")
        system(mv_cx)
    }
    message("CIBERSORTx successfully installed")
    invisible(TRUE)
}

.docker_check <- function() {
    purrr::possibly(function() {
        system("docker --help", intern = TRUE)
        return(TRUE)
    }, otherwise = FALSE)()
}

.singularity_check <- function() {
    purrr::possibly(function() {
        system("singularity --help", intern = TRUE)
        return(TRUE)
    }, otherwise = FALSE)()
}

#' Compute reference signature matrix from `screference` object using
#' method implemented in CIBERSORTx
#'
#' @param scref an object of `screference
#' @param cache_path path to cache results
#' @param username a username used to authenticate access to CIBERSORTx method.
#' See: https://cibersortx.stanford.edu/getoken.php
#' @param token a token used to authenticate access to CIBERSORTx method
#'
#' @return a path to cached signature matrix
#'
#' @importFrom R.utils getAbsolutePath
#' @import tidyverse
#' @export
cibersortx_scref <- function(scref,
                             cache_path = "cibersortx",
                             username = NULL,
                             token = NULL) {
    #-- Get reference data
    ref <- scref$seurat_obj
    ref_pops <- scref$populations

    #-- Create file structure
    dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
    path <- getAbsolutePath(cache_path)

    #-- Get paths
    data_path <- file.path(path, "data")
    out_path <- file.path(path, "out_dir")
    dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_path, showWarnings = FALSE)

    #-- Generate sigmat
    message("-- CIBERSORTx: Computing reference data matrix...")
    ref_fname <- file.path(path, "data/ref_data.txt")

    #----- Generate data matrix
    ref <- NormalizeData(ref, normalization.method = "RC", scale.factor = 10^5)
    if(is.null(ref_pops)) {
        ref_pops <- unique(ref@meta.data[,"annot_id"])
    }
    gc()
    cell_annots <- ref@meta.data %>%
        dplyr::filter(annot_id %in% ref_pops) %>%
        rownames_to_column("id") %>%
        pull(id, annot_id)
    ref_data <- ref@assays$RNA@data[,cell_annots] %>% as.matrix() %>% as.data.frame()
    colnames(ref_data) <- names(cell_annots)
    data.table::fwrite(ref_data, file = ref_fname, sep = "\t",
                           row.names = TRUE, col.names = TRUE)
    rm(ref_data)

    #-- Create docker call
    header <- case_when(
        .docker_check() ~ str_glue("docker run -v {data_path}:/src/data/ -v {out_path}:/src/outdir"),
        .singularity_check() ~ str_glue("singularity exec -B {data_path}:/src/data/ -B {out_path}:/src/outdir") ,
        TRUE ~ "")
    if(header == "") {
        stop("No CIBERSORTx container found. Please run `install_cibersortx()`")
    }
    container_address <- case_when(
        header == "docker run" ~ "cibersortx/fractions",
        TRUE ~ paste0(find.package("deconverse"), "/cibersortx.sif /src/CIBERSORTxFractions")
    )
    docker_call <- str_glue('{header} {container_address} --username {username} --token {token} --rmbatchSmode TRUE --verbose TRUE --QN FALSE --single_cell TRUE --refsample /src/data/ref_data.txt')
    rm(ref); gc()
    system(docker_call)

    system(str_glue("rm -r {data_path}"))
    return(out_path)

}

#' CIBERSORTx deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param cache_path path to cache deconvolution results
#' @param username a username used to authenticate access to CIBERSORTx method.
#' See: https://cibersortx.stanford.edu/getoken.php
#' @param token a token used to authenticate access to CIBERSORTx method
#'
#' @param a tibble with deconvolution fractions
#'
#' @importFrom R.utils getAbsolutePath filePath
#' @import tidyverse
#'
#' @export
cibersortx_deconvolute <- function(bulk_data,
                                   scref,
                                   cache_path = "cibersortx",
                                   username = NULL,
                                   token = NULL) {
    #-- Get reference data
    ref <- scref$cached_results$cibersortx

    #-- Check credentials
    if(is.null(username) | is.null(token)) {
        stop("Can't run CIBERSORTx without credentials. See: https://cibersortx.stanford.edu/getoken.php")
    }
    #-- Create file structure
    path <- getAbsolutePath(cache_path)
    data_path <- file.path(path, "data")
    out_path <- file.path(path, "out_dir")
    out_file <- file.path(path, "out_dir", "CIBERSORTx_Adjusted.txt")
    if(file.exists(out_file)) {
        message("-- Getting cached results...")
        deconv_res <- read_tsv(out_file)
        deconv_res <- deconv_res %>%
            rename(sample = Mixture, model_pvalue = `P-value`, model_cor = Correlation, model_RMSE = RMSE) %>%
            rename_with(.fn = ~ paste0("frac_", .), .cols = -matches("^model|sample")) %>%
            mutate(method = "CIBERSORTx")
        return(deconv_res)
    }

    dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
    dir.create(out_path, showWarnings = FALSE)

    #-- Writing data
    data_fname <- file.path(data_path, "data.txt")
    message("-- Writing data...")
    data.table::fwrite(as.data.frame(bulk_data), file = data_fname,
                       sep = "\t", row.names = TRUE, col.names = TRUE)
    #-- Grabbing reference
    assert(file.exists(ref))
    message("-- Using signature matrix...")
    file.copy(from = ref, to = filePath(data_path, "sig_mat.txt"), overwrite = TRUE)

    #-- Running CIBERSORTx in docker
    header <- case_when(
        .docker_check() ~ "docker run",
        .singularity_check() ~ "singularity run",
        TRUE ~ "")
    if(header == "") {
        stop("No CIBERSORTx container found. Please run `install_cibersortx()`")
    }
    container_address <- case_when(
        header == "docker run" ~ "cibersortx/fractions",
        TRUE ~ paste0(find.package("deconverse"), "/cibersortx.sif /src/CIBERSORTxFractions")
    )
    docker_call <- str_glue('{header} run -v {data_path}:/src/data/ -v {out_path}:/src/outdir {container_address} --username {username} --token {token} --rmbatchBmode TRUE --verbose TRUE --QN FALSE --mixture src/data/data.txt --sigmatrix /src/datasig_mat.txt')
    rm(ref); rm(bulk_data); gc()
    message("-- Running CIBERSORTx...")
    system(docker_call)

    system(str_glue("rm -r {data_path}"))
    deconv_res <- read_tsv(out_file)

    deconv_res <- deconv_res %>%
        rename(sample = Mixture, model_pvalue = `P-value`, model_cor = Correlation, model_RMSE = RMSE) %>%
        rename_with(.fn = ~ paste0("frac_", .), .cols = -matches("^model|sample")) %>%
        mutate(method = "CIBERSORTx")

    return(deconv_res)
}

# cibersortx_scRef <- function(data,
#                              ref = NULL,
#                              ref_annot = NULL,
#                              ref_pops = NULL,
#                              cache_path = "results/scref/CIBERSORTx/",
#                              username = "clarice.groeneveld@curie.fr",
#                              token = "9409049f84c88f699b5b541eb951466e") {
#     #-- Create file structure
#     if(is.null(cache_path)) {
#         dir.create("cibersort_tmp", showWarnings = FALSE)
#         path <- file.path(getwd(), "cibersort_tmp")
#     } else {
#         dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
#         path <- getAbsolutePath(cache_path)
#     }
#     data_path <- file.path(path, "data")
#     out_path <- file.path(path, "out_dir")
#     dir.create(data_path, recursive = TRUE, showWarnings = FALSE)
#     dir.create(out_path, showWarnings = FALSE)
#     out_file <- file.path(path, "out_dir", "CIBERSORTx_Adjusted.txt")
#
#     if(file.exists(out_file)) {
#         message("-- Getting cached results...")
#         deconv_res <- read_tsv(out_file)
#         return(deconv_res)
#     }
#
#     #--
#     message("-- Writing data...")
#     data_fname <- file.path(data_path, "data.txt")
#     if(!file.exists(data_fname)) {
#         data.table::fwrite(as.data.frame(data), file = data_fname,
#                            sep = "\t", row.names = TRUE, col.names = TRUE)
#     }
#
#     #-- Run from sigmat or not
#     if(is.character(ref)) {
#         message("-- Using signature matrix...")
#         file.copy(ref, filePath(data_path, "sig_mat.txt"), overwrite = TRUE)
#         docker_call <- str_glue('docker run -v {data_path}:/src/data/ -v {out_path}:/src/outdir cibersortx/fractions --username {username} --token {token} --rmbatchBmode TRUE --verbose TRUE --QN FALSE --mixture data.txt --sigmatrix sig_mat.txt')
#
#     } else {
#         message("-- Saving reference data matrix...")
#         ref_fname <- file.path(path, "data/ref_data.txt")
#
#         if(!file.exists(ref_fname)) {
#             if(!is.null(ref) & !is.null(ref_annot)) {
#                 #-- Generate data reference file
#                 ref <- NormalizeData(ref, normalization.method = "RC", scale.factor = 10^5)
#                 if(is.null(ref_pops)) {
#                     ref_pops <- unique(ref@meta.data[,ref_annot])
#                 }
#                 gc()
#                 cell_annots <- ref@meta.data %>%
#                     dplyr::rename(ref_annot = ref_annot) %>%
#                     dplyr::filter(ref_annot %in% ref_pops) %>%
#                     rownames_to_column("id") %>%
#                     pull(id, ref_annot)
#                 ref_data <- ref@assays$RNA@data[,cell_annots] %>% as.matrix() %>% as.data.frame()
#                 colnames(ref_data) <- names(cell_annots)
#                 data.table::fwrite(ref_data, file = ref_fname, sep = "\t",
#                                    row.names = TRUE, col.names = TRUE)
#                 rm(ref_data)
#             }
#         }
#         docker_call <- str_glue('docker run -v {data_path}:/src/data/ -v {out_path}:/src/outdir cibersortx/fractions --username {username} --token {token} --rmbatchSmode TRUE --rmbatchBmode TRUE --verbose TRUE --QN FALSE --single_cell TRUE --refsample ref_data.txt --mixture data.txt')
#     }
#     rm(ref); rm(data); gc()
#     message("-- Running CIBERSORTx...")
#     system(docker_call)
#
#     deconv_res <- read_tsv(out_file)
#     return(deconv_res)
#
# }


#-------------------------------------------------------------------------------
# cibersortx_default <- function(data,
#                            username = "clarice.groeneveld@curie.fr",
#                            token = "9409049f84c88f699b5b541eb951466e",
#                            refdata_path = "~/Documents/Tools/CIBERSORTx/data/",
#                            cache_path = "results/tests/CIBERSORTx/") {
#     #-- Create file structure
#     if(is.null(cache_path)) {
#         dir.create("cibersort_tmp", showWarnings = FALSE)
#         path <- file.path(getwd(), "cibersort_tmp")
#     } else {
#         dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
#         path <- getAbsolutePath(cache_path)
#     }
#
#
#     data_path <- file.path(path, "data")
#     dir.create(data_path, showWarnings = FALSE)
#     coarse_path <- file.path(path, "coarse")
#     fine_path <- file.path(path, "fine")
#
#     if(dir.exists(coarse_path) & dir.exists(fine_path)) {
#         message("-- Using cached deconvolution results...")
#     } else {
#         message("-- Writing data...")
#         #-- Write data
#         data %>%
#             as.data.frame() %>%
#             rownames_to_column("symbol") %>%
#             write_tsv(file.path(data_path, "data.txt"))
#
#         coarse_sig <- file.path(refdata_path, "FACS_signatures.txt")
#         fine_sig <- file.path(refdata_path, "LM22.txt")
#
#         copyFile(coarse_sig, file.path(data_path, "coarse_sig.txt"), overwrite=TRUE)
#         copyFile(fine_sig, file.path(data_path, "fine_sig.txt"), overwrite=TRUE)
#
#         docker_coarse <- str_glue('docker run -v {data_path}:/src/data/ -v {coarse_path}:/src/outdir cibersortx/fractions --username {username} --token {token} --rmbatchBmode TRUE --perm 1 --verbose TRUE --QN FALSE --mixture data.txt --sigmatrix coarse_sig.txt')
#         docker_fine <- str_glue('docker run -v {data_path}:/src/data/ -v {fine_path}:/src/outdir cibersortx/fractions --username {username} --token {token} --rmbatchBmode TRUE --perm 1 --verbose TRUE --QN FALSE --mixture data.txt --sigmatrix fine_sig.txt')
#
#         #-- Call dockers
#         message("-- Running coarse-grained deconvolution...")
#         system(docker_coarse)
#         message("-- Running fine-grained deconvolution...")
#         system(docker_fine)
#     }
#
#     #-- Read data
#     coarse_res <- read_tsv(filePath(coarse_path, "CIBERSORTx_Adjusted.txt"))
#     fine_res <- read_tsv(filePath(fine_path, "CIBERSORTx_Adjusted.txt"))
#
#     #-- Get scores
#     cx_res <- .calculate_cibersort_scores(coarse_res, fine_res)
#
#     return(cx_res)
# }
#
# .calculate_cibersort_scores <- function(coarse_res, fine_res) {
#     coarse_res %>%
#         select(-`P-value`, -Correlation, -RMSE) %>%
#         left_join(select(fine_res, -`P-value`, Correlation, RMSE), by = "Mixture") %>%
#         transmute(
#             sample = Mixture,
#             B_cells = CD45 * (`B cells naive` + `B cells memory`),
#             memory_B_cells = CD45 * `B cells memory`,
#             naive_B_cells = CD45 * `B cells naive`,
#             CD4_T_cells = CD45 * (`T cells CD4 naive` + `T cells CD4 memory resting` +
#                                       `T cells CD4 memory activated` + `T cells regulatory (Tregs)` +
#                                       `T cells follicular helper`),
#             memory_CD4_T_cells =  CD45 * (`T cells CD4 memory activated` + `T cells CD4 memory resting`),
#             naive_CD4_T_cells = CD45 * `T cells CD4 naive`,
#             regularoty_T_cells = CD45 * `T cells regulatory (Tregs)`,
#             CD8_T_cells = CD45 * (`T cells CD8` + `T cells gamma delta`),
#             NK_cells = CD45 * (`NK cells resting` + `NK cells activated`),
#             neutrophils = CD45 * Neutrophils,
#             monocytic_lineage = CD45 * (Monocytes + `Macrophages M0` +
#                                             `Macrophages M1` + `Macrophages M2` +
#                                             `Dendritic cells resting` +
#                                             `Dendritic cells activated`),
#             macrophages = CD45 * (`Macrophages M0` + `Macrophages M1` + `Macrophages M2`),
#             monocytes = CD45 * Monocytes,
#             myeloid_dendritic_cells = CD45 * (`Dendritic cells resting` +
#                                                   `Dendritic cells activated`),
#             endothelial_cells = CD31,
#             fibroblasts = CD10,
#             epithelial_cells = EPCAM)
#
# }
