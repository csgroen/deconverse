# source: https://bitbucket.org/yuanlab/dwls/
# reimplemented for performance improvements and code readability
# dwls -------------------------------------
#' Dampened weighted least squares (DWLS) deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param n_cores number of cores used for computation
#'
#' @return a tibble with deconvolution fractions
#'
#' @import tidyverse
#' @importFrom stats na.exclude
#' @importFrom pbmcapply pbmclapply
#'
#' @export
dwls_deconvolute <- function(bulk_data, scref, ncores = 4) {
    #-- Get signatures
    if(is.null(scref[["cached_results"]][["dwls"]])) {
        stop("scref doesn't have cached `dwls` reference signatures. Please run `compute_reference` on the scref object.")
    }
    sc_sigs <- scref$cached_results$dwls
    # colnames(sc_sigs) <- scref$populations
    #-- Pre-process
    data4dwls <- .trimData(sc_sigs, bulk_data)
    B_list <- lapply(1:ncol(data4dwls$bulk), function(i) data4dwls$bulk[,i])
    #-- Deconvolute
    # props <- pbmclapply(B_list, function(bulk_sample) {
    #     .solveDampenedWLS(data4dwls$sig, bulk_sample)
    # }, mc.cores = ncores)
    props <- lapply(B_list, function(bulk_sample) {
        .solveDampenedWLS(data4dwls$sig, bulk_sample)
    })
    #-- Prepare results
    deconv_res <- bind_rows(props) %>%
        dplyr::rename_with(~ paste0("frac_", .)) %>%
        mutate(sample = colnames(bulk_data),
               method = "DWLS") %>%
        relocate(sample)
    return(deconv_res)
}

#' Ordinary least squares (OLS) deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param n_cores number of cores used for computation
#'
#' @param a tibble with deconvolution fractions
#'
#' @import tidyverse
#' @importFrom pbmcapply pbmclapply
#'
#' @export

ols_deconvolute <- function(bulk_data, scref, ncores = 4) {
    #-- Get signatures
    if(is.null(scref[["cached_results"]][["dwls"]])) {
        stop("scref doesn't have cached `dwls` reference signatures. Please run `compute_reference` on the scref object.")
    }
    sc_sigs <- scref$cached_results$dwls
    #-- Pre-process
    data4dwls <- .trimData(sc_sigs, bulk_data)
    B_list <- lapply(seq_len(ncol(data4dwls$bulk)), function(i) data4dwls$bulk[,i])
    #-- Deconvolute
    props <- pbmclapply(B_list, function(bulk_sample) {
        .solveOLS(data4dwls$sig, bulk_sample)
    }, mc.cores = ncores)
    #-- Prepare results
    deconv_res <- bind_rows(props) %>%
        dplyr::rename_with(~ paste0("frac_", .)) %>%
        mutate(sample = colnames(bulk_data),
               method = "OLS") %>%
        relocate(sample)
    return(deconv_res)
}

#' Support vector regression (SVR) deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#' @param n_cores number of cores used for computation
#'
#' @import tidyverse
#' @importFrom pbmcapply pbmclapply
#'
#' @param a tibble with deconvolution fractions
#' @export
svr_deconvolute <- function(bulk_data, scref, ncores = 4) {
    #-- Get signatures
    if(is.null(scref[["cached_results"]][["dwls"]])) {
        stop("scref doesn't have cached `dwls` reference signatures. Please run `compute_reference` on the scref object.")
    }
    sc_sigs <- scref$cached_results$dwls
    #-- Pre-process
    data4dwls <- .trimData(sc_sigs, bulk_data)
    #-- Deconvolute
    B_list <- lapply(seq_len(ncol(data4dwls$bulk)), function(i) data4dwls$bulk[,i])
    props <- pbmclapply(B_list, function(bulk_sample) {
        .solveSVR(data4dwls$sig, bulk_sample)
    }, mc.cores = ncores)
    #-- Prepare results
    deconv_res <- bind_rows(props) %>%
        dplyr::rename_with(~ paste0("frac_", .)) %>%
        mutate(sample = colnames(bulk_data),
               method = "SVR") %>%
        relocate(sample)
    return(deconv_res)
}

#' Compute reference signature matrix from `screference` object using
#' method implemented in DWLS package for methods `dwls`, `svr` and `ols`
#'
#' @param scref an object of `screference`
#' @param cache_path path to cache results
#' @param logFC_cutoff a float, cutoff for log fold-change from differential
#' expression analysis between populations for selection of markers
#' @param pval_cutoff a float, p-value cutoff from differential expression
#' analysis for selection of markers
#'
#' @import tidyverse
#' @import Seurat
#' @importFrom pbmcapply pbmclapply
#' @importFrom purrr reduce
#'
#' @return a signature matrix
#' @export
dwls_scref <- function(scref, cache_path = "dwls",
                       logFC_cutoff = 0.5, pval_cutoff = 0.01) {
    scdata <- scref$seurat_obj
    id_name <- "annot_id"

#     #-- Get cached if exists
#     sig_fname <- filePath(cache_path, "Sig.RData")
#     if (file.exists(sig_fname)) {
#         message("-- Signature found in cache. Skipping...")
#         load(sig_fname)
#         return(Sig)
#     }
    #-- Run DEG analysis
    ids <- scdata@meta.data[, id_name]
    pops <- unique(ids)
    .dwls_DEAnalysis(scdata, id_name, pops, cache_path)

    #-- Load DEG results
    cluster_lr_tables <- lapply(pops, function(pop) {
        load(file = paste(cache_path, "/de_", pop, ".RData", sep = ""))
        DEGenes <- rownames(de_group)[de_group$p_val_adj < pval_cutoff & de_group$avg_log2FC > logFC_cutoff]
        DEGenes <- str_subset(DEGenes, "MIR|Mir", negate = TRUE)
        de_group[DEGenes, ] %>%
            arrange(desc(p_val_adj))
    })
    names(cluster_lr_tables) <- pops
    numberofGenes <- sapply(cluster_lr_tables, nrow)

    #-- Find best number of genes (minimize Kappa)
    message("-- Computing signature matrix...")
    ngenes <- seq(50, 200, by = 1)
    kappas <- sapply(ngenes, function(G) {
        Genes <- sapply(pops, function(pop) {
            deg_table <- cluster_lr_tables[[pop]]
            rownames(deg_table)[1:min(G, nrow(deg_table))]
        }) %>% reduce(c)
        Genes <- na.exclude(unique(Genes))
        ExprSubset <- scdata@assays$RNA@data[Genes, ]
        Sig <- sapply(pops, function(pop) {
            Matrix::rowMeans(ExprSubset[, ids == pop])
        })
        kappa(Sig)
    })
    gc()
    G <- which.min(kappas) + min(49, numberofGenes - 1)
    Genes <- sapply(pops, function(pop) {
        deg_table <- cluster_lr_tables[[pop]]
        rownames(deg_table)[1:min(G, nrow(deg_table))]
    }) %>% reduce(c)
    Genes <- unique(Genes)
    ExprSubset <- scdata@assays$RNA@data[Genes, ]
    Sig <- sapply(pops, function(pop) {
        Matrix::rowMeans(ExprSubset[, ids == pop])
    })
    save(Sig, file = paste(cache_path, "/Sig.RData", sep = ""))
    return(Sig)
}

# Helpers -------------------------------------------


# solve using OLS, constrained such that cell type numbers>0
#' @importFrom quadprog solve.QP
.solveOLS <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  solution <- solve.QP(D, d, A, bzero)$solution
  names(solution) <- colnames(S)
  print(round(solution / sum(solution), 5))
  return(solution / sum(solution))
}

# return cell number, not proportion
# do not print output
#' @importFrom quadprog solve.QP
.solveOLSInternal <- function(S, B) {
  D <- t(S) %*% S
  d <- t(S) %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  solution <- solve.QP(D, d, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# solve using WLS with weights dampened by a certain dampening constant
.solveDampenedWLS <- function(S, B) {
  # first solve OLS, use this solution to find a starting point for the weights
  solution <- .solveOLSInternal(S, B)
  # now use dampened WLS, iterate weights until convergence
  iterations <- 0
  changes <- c()
  # find dampening constant for weights using cross-validation
  j <- .findDampeningConstant(S, B, solution)
  change <- 1
  while (change > .01 & iterations < 1000) {
    newsolution <- .solveDampenedWLSj(S, B, solution, j)
    # decrease step size for convergence
    solutionAverage <- rowMeans(cbind(newsolution, matrix(solution, nrow = length(solution), ncol = 4)))
    change <- norm(as.matrix(solutionAverage - solution))
    solution <- solutionAverage
    iterations <- iterations + 1
    changes <- c(changes, change)
  }
  # print(round(solution / sum(solution), 5))
  return(solution / sum(solution))
}

# solve WLS given a dampening constant
#' @importFrom quadprog solve.QP
.solveDampenedWLSj <- function(S, B, goldStandard, j) {
  multiplier <- 1 * 2^(j - 1)
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsDampened <- wsScaled
  wsDampened[which(wsScaled > multiplier)] <- multiplier
  W <- diag(wsDampened)
  D <- t(S) %*% W %*% S
  d <- t(S) %*% W %*% B
  A <- cbind(diag(dim(S)[2]))
  bzero <- c(rep(0, dim(S)[2]))
  sc <- norm(D, "2")
  solution <- solve.QP(D / sc, d / sc, A, bzero)$solution
  names(solution) <- colnames(S)
  return(solution)
}

# find a dampening constant for the weights using cross-validation
#' @importFrom stats sd lm
.findDampeningConstant <- function(S, B, goldStandard) {
  solutionsSd <- NULL
  # goldStandard is used to define the weights
  sol <- goldStandard
  ws <- as.vector((1 / (S %*% sol))^2)
  wsScaled <- ws / min(ws)
  wsScaledMinusInf <- wsScaled
  # ignore infinite weights
  if (max(wsScaled) == "Inf") {
    wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
  }
  # try multiple values of the dampening constant (multiplier)
  # for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
    multiplier <- 1 * 2^(j - 1)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    solutions <- NULL
    seeds <- c(1:100)
    for (i in 1:100) {
      set.seed(seeds[i]) # make nondeterministic
      subset <- sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
      # solve dampened weighted least squares for subset
      fit <- lm(B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
      sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
      solutions <- cbind(solutions, sol)
    }
    solutionsSd <- cbind(solutionsSd, apply(solutions, 1, sd))
  }
  # choose dampening constant that results in least cross-validation variance
  j <- which.min(colMeans(solutionsSd^2))
  return(j)
}

#' @importFrom e1071 svm
.solveSVR <- function(S, B) {
  # scaling
  ub <- max(c(as.vector(S), B)) # upper bound
  lb <- min(c(as.vector(S), B)) # lower bound
  Bs <- ((B - lb) / ub) * 2 - 1
  Ss <- ((S - lb) / ub) * 2 - 1

  # perform SVR
  model <- svm(Ss, Bs, nu = 0.5, scale = TRUE, type = "nu-regression", kernel = "linear", cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[which(coef < 0)] <- 0
  coef <- as.vector(coef)
  names(coef) <- colnames(S)
  print(round(coef / sum(coef), 5))
  return(coef / sum(coef))
}
# trim bulk and single-cell data to contain the same genes
.trimData <- function(Signature, bulkData) {
    Genes <- intersect(rownames(Signature), rownames(bulkData))
    B <- bulkData[Genes,]
    S <- Signature[Genes, ]
    return(list("sig" = S, "bulk" = B))
}

#' @import Seurat
.dwls_DEAnalysis <- function(scdata, id_name, pops, cache_path) {
  scdata <- SetIdent(scdata, value = id_name)
  for (i in pops) {
    de_file <- paste(cache_path, "/de_", i, ".RData", sep = "")
    if (file.exists(de_file)) {
      message("-- DEGs for ", i, " found in cache...")
      next()
    }
    message("-- Calculating DEGs for ", i, "...")
    de_group <- FindMarkers(
      object = scdata,
      ident.1 = i,
      ident.2 = NULL,
      only.pos = TRUE,
      test.use = "bimod"
    )
    save(de_group, file = de_file)
  }
}
