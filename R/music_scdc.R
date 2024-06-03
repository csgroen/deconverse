#' MuSiC deconvolution of bulk data using an `screference`
#'
#' @param bulk_data a matrix of genes-by-samples with bulk mixtures
#' @param scref an object of class `screference`
#'
#' @return a tibble with deconvolution fractions and model metrics
#'
#' @importFrom Biobase ExpressionSet exprs AnnotatedDataFrame
#' @export
music_deconvolute <- function(bulk_data, scref) {
    if(!"batch_id" %in% colnames(scref$seurat_obj@meta.data)) {
        stop("MuSiC needs `batch_id` to be defined when creating the scref object for deconvolution.")
    }
    est_prop <- music_prop(bulk.eset = ExpressionSet(assayData = bulk_data),
                           sc.eset = ExpressionSet(as.matrix(GetAssayData(scref$seurat_obj, "RNA", "counts")),
                                                   phenoData = AnnotatedDataFrame(scref$seurat_obj[[]])),
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

    message("MuSiC deconvolution complete.")

    return(deconv_res)
}
# MuSiC internals ------
# TODO: optimize

#' @importFrom Biobase exprs pData phenoData
#' @importFrom Matrix rowMeans colMeans
music_prop <- function (bulk.eset, sc.eset, markers = NULL, clusters, samples,
                        select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE,
                        iter.max = 1000, nu = 1e-04, eps = 0.01, centered = FALSE,
                        normalize = FALSE, ...)
{
    bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) != 0]
    bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
    if (is.null(markers)) {
        sc.markers = bulk.gene
    }
    else {
        sc.markers = intersect(bulk.gene, unlist(markers))
    }
    sc.basis = music_basis(sc.eset, non.zero = TRUE, markers = sc.markers,
                           clusters = clusters, samples = samples, select.ct = select.ct,
                           cell_size = cell_size, ct.cov = ct.cov, verbose = verbose)
    cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
    if (is.null(markers)) {
        if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.eset)))
            stop("Too few common genes!")
    }
    else {
        if (length(cm.gene) < 0.2 * length(unlist(markers)))
            stop("Too few common genes!")
    }
    if (verbose) {
        message(paste("Used", length(cm.gene), "common genes..."))
    }
    m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
    m.bulk = match(cm.gene, bulk.gene)
    D1 = sc.basis$Disgn.mtx[m.sc, ]
    M.S = colMeans(sc.basis$S, na.rm = T)
    if (!is.null(cell_size)) {
        if (!is.data.frame(cell_size)) {
            stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
        }
        else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
            stop("Cell type names in cell_size must match clusters")
        }
        else if (any(is.na(as.numeric(cell_size[, 2])))) {
            stop("Cell sizes should all be numeric")
        }
        my_ms_names <- names(M.S)
        cell_size <- cell_size[my_ms_names %in% cell_size[, 1],
        ]
        M.S <- cell_size[match(my_ms_names, cell_size[, 1]),
        ]
        M.S <- M.S[, 2]
        names(M.S) <- my_ms_names
    }
    Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ])
    N.bulk = ncol(bulk.eset)
    if (ct.cov) {
        Sigma.ct = sc.basis$Sigma.ct[, m.sc]
        Est.prop.allgene = NULL
        Est.prop.weighted = NULL
        Weight.gene = NULL
        r.squared.full = NULL
        Var.prop = NULL
        for (i in 1:N.bulk) {
            if (sum(Yjg[, i] == 0) > 0) {
                D1.temp = D1[Yjg[, i] != 0, ]
                Yjg.temp = Yjg[Yjg[, i] != 0, i]
                Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
                if (verbose)
                    message(paste(colnames(Yjg)[i], "has common genes",
                                  sum(Yjg[, i] != 0), "..."))
            }
            else {
                D1.temp = D1
                Yjg.temp = Yjg[, i]
                Sigma.ct.temp = Sigma.ct
                if (verbose)
                    message(paste(colnames(Yjg)[i], "has common genes",
                                  sum(Yjg[, i] != 0), "..."))
            }
            lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp,
                                           M.S, Sigma.ct.temp, iter.max = iter.max, nu = nu,
                                           eps = eps, centered = centered, normalize = normalize)
            Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
            Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
            weight.gene.temp = rep(NA, nrow(Yjg))
            weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
            Weight.gene = cbind(Weight.gene, weight.gene.temp)
            r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
            Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
        }
    }
    else {
        Sigma = sc.basis$Sigma[m.sc, ]
        valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) ==
                                                       0) & (!is.na(M.S))
        if (sum(valid.ct) <= 1) {
            stop("Not enough valid cell type!")
        }
        if (verbose) {
            message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
        }
        D1 = D1[, valid.ct]
        M.S = M.S[valid.ct]
        Sigma = Sigma[, valid.ct]
        Est.prop.allgene = NULL
        Est.prop.weighted = NULL
        Weight.gene = NULL
        r.squared.full = NULL
        Var.prop = NULL
        for (i in 1:N.bulk) {
            if (sum(Yjg[, i] == 0) > 0) {
                D1.temp = D1[Yjg[, i] != 0, ]
                Yjg.temp = Yjg[Yjg[, i] != 0, i]
                Sigma.temp = Sigma[Yjg[, i] != 0, ]
                if (verbose)
                    message(paste(colnames(Yjg)[i], "has common genes",
                                  sum(Yjg[, i] != 0), "..."))
            }
            else {
                D1.temp = D1
                Yjg.temp = Yjg[, i]
                Sigma.temp = Sigma
                if (verbose)
                    message(paste(colnames(Yjg)[i], "has common genes",
                                  sum(Yjg[, i] != 0), "..."))
            }
            lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S,
                                        Sigma.temp, iter.max = iter.max, nu = nu, eps = eps,
                                        centered = centered, normalize = normalize)
            Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
            Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
            weight.gene.temp = rep(NA, nrow(Yjg))
            weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
            Weight.gene = cbind(Weight.gene, weight.gene.temp)
            r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
            Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
        }
    }
    colnames(Est.prop.weighted) = colnames(D1)
    rownames(Est.prop.weighted) = colnames(Yjg)
    colnames(Est.prop.allgene) = colnames(D1)
    rownames(Est.prop.allgene) = colnames(Yjg)
    names(r.squared.full) = colnames(Yjg)
    colnames(Weight.gene) = colnames(Yjg)
    rownames(Weight.gene) = cm.gene
    colnames(Var.prop) = colnames(D1)
    rownames(Var.prop) = colnames(Yjg)
    return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
                Weight.gene = Weight.gene, r.squared.full = r.squared.full,
                Var.prop = Var.prop))
}

music_basis <- function (x, non.zero = TRUE, markers = NULL, clusters, samples,
          select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE)
{
    if (!is.null(select.ct)) {
        s.ct = sampleNames(x)[as.character(pData(x)[, clusters]) %in%
                                  select.ct]
        x <- x[, s.ct, drop = FALSE]
    }
    if (non.zero) {
        nz.gene = rownames(x)[(rowSums(exprs(x)) != 0)]
        x <- x[nz.gene, , drop = FALSE]
    }
    clusters <- as.character(pData(x)[, clusters])
    samples <- as.character(pData(x)[, samples])
    M.theta <- sapply(unique(clusters), function(ct) {
        my.rowMeans(sapply(unique(samples), function(sid) {
            y = exprs(x)[, clusters %in% ct & samples %in% sid,
                         drop = FALSE]
            rowSums(y)/sum(y)
        }), na.rm = TRUE)
    })
    if (verbose) {
        message("Creating Relative Abundance Matrix...")
    }
    if (ct.cov) {
        nGenes = nrow(x)
        n.ct = length(unique(clusters))
        nSubs = length(unique(samples))
        Theta <- sapply(unique(clusters), function(ct) {
            sapply(unique(samples), function(sid) {
                y = exprs(x)[, clusters %in% ct & samples %in%
                                 sid, drop = FALSE]
                return(rowSums(y)/sum(y))
            })
        })
        if (!is.null(select.ct)) {
            m.ct = match(select.ct, colnames(Theta))
            Theta = Theta[, m.ct]
        }
        Sigma.ct = sapply(1:nGenes, function(g) {
            sigma.temp = Theta[nGenes * (0:(nSubs - 1)) + g,
            ]
            Cov.temp = cov(sigma.temp)
            Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes *
                                                               (0:(nSubs - 1)) + 1, ])) == 0, ])
            Cov.temp[which(colSums(is.na(sigma.temp)) > 0), ] = Cov.temp1[which(colSums(is.na(sigma.temp)) >
                                                                                    0), ]
            Cov.temp[, which(colSums(is.na(sigma.temp)) > 0)] = Cov.temp1[,
                                                                          which(colSums(is.na(sigma.temp)) > 0)]
            return(Cov.temp)
        })
        colnames(Sigma.ct) = rownames(x)
        if (!is.null(markers)) {
            ids <- intersect(unlist(markers), rownames(x))
            m.ids = match(ids, rownames(x))
            Sigma.ct <- Sigma.ct[, m.ids]
        }
        if (verbose) {
            message("Creating Covariance Matrix...")
        }
    } else {
        Sigma <- sapply(unique(clusters), function(ct) {
            apply(sapply(unique(samples), function(sid) {
                y = exprs(x)[, clusters %in% ct & samples %in%
                                 sid, drop = FALSE]
                rowSums(y)/sum(y)
            }), 1, var, na.rm = TRUE)
        })
        if (!is.null(select.ct)) {
            m.ct = match(select.ct, colnames(Sigma))
            Sigma = Sigma[, m.ct]
        }
        if (!is.null(markers)) {
            ids <- intersect(unlist(markers), rownames(x))
            m.ids = match(ids, rownames(x))
            Sigma <- Sigma[m.ids, ]
        }
        if (verbose) {
            message("Creating Variance Matrix...")
        }
    }
    S <- sapply(unique(clusters), function(ct) {
        my.rowMeans(sapply(unique(samples), function(sid) {
            y = exprs(x)[, clusters %in% ct & samples %in% sid,
                         drop = FALSE]
            sum(y)/ncol(y)
        }), na.rm = TRUE)
    })
    if(is.null(ncol(S))) {
        stop("`MuSiC requires more than one sample per cluster. See MuSiC documentation.")
    }
    if (verbose) {
        message("Creating Library Size Matrix...")
    }
    S[S == 0] = NA
    M.S = colMeans(S, na.rm = TRUE)
    if (!is.null(cell_size)) {
        if (!is.data.frame(cell_size)) {
            stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
        } else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
            stop("Cell type names in cell_size must match clusters")
        } else if (any(is.na(as.numeric(cell_size[, 2])))) {
            stop("Cell sizes should all be numeric")
        }
        my_ms_names <- names(M.S)
        cell_size <- cell_size[my_ms_names %in% cell_size[, 1],
        ]
        M.S <- cell_size[match(my_ms_names, cell_size[, 1]),
        ]
        M.S <- M.S[, 2]
        names(M.S) <- my_ms_names
    }
    D <- t(t(M.theta) * M.S)
    if (!is.null(select.ct)) {
        m.ct = match(select.ct, colnames(D))
        D = D[, m.ct]
        S = S[, m.ct]
        M.S = M.S[m.ct]
        M.theta = M.theta[, m.ct]
    }
    if (!is.null(markers)) {
        ids <- intersect(unlist(markers), rownames(x))
        m.ids = match(ids, rownames(x))
        D <- D[m.ids, ]
        M.theta <- M.theta[m.ids, ]
    }
    if (ct.cov) {
        return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta,
                    Sigma.ct = Sigma.ct))
    }
    else {
        return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta,
                    Sigma = Sigma))
    }
}
music.iter <- function (Y, D, S, Sigma, iter.max = 1000, nu = 1e-04, eps = 0.01,
                        centered = FALSE, normalize = FALSE)
{
    if (length(S) != ncol(D)) {
        common.cell.type = intersect(colnames(D), names(S))
        if (length(common.cell.type) <= 1) {
            stop("Not enough cell types!")
        }
        D = D[, match(common.cell.type, colnames(D))]
        S = S[match(common.cell.type, names(S))]
    }
    if (ncol(Sigma) != ncol(D)) {
        common.cell.type = intersect(colnames(D), colnames(Sigma))
        if (length(common.cell.type) <= 1) {
            stop("Not enough cell type!")
        }
        D = D[, match(common.cell.type, colnames(D))]
        Sigma = Sigma[, match(common.cell.type, colnames(Sigma))]
        S = S[match(common.cell.type, names(S))]
    }
    k = ncol(D)
    common.gene = intersect(names(Y), rownames(D))
    if (length(common.gene) < 0.1 * min(length(Y), nrow(D))) {
        stop("Not enough common genes!")
    }
    Y = Y[match(common.gene, names(Y))]
    D = D[match(common.gene, rownames(D)), ]
    Sigma = Sigma[match(common.gene, rownames(Sigma)), ]
    X = D
    if (centered) {
        X = X - mean(X)
        Y = Y - mean(Y)
    }
    if (normalize) {
        X = X/sd(as.vector(X))
        S = S * sd(as.vector(X))
        Y = Y/sd(Y)
    }
    else {
        Y = Y * 100
    }
    lm.D = music.basic(Y, X, S, Sigma, iter.max = iter.max, nu = nu,
                       eps = eps)
    return(lm.D)
}

music.iter.ct <- function (Y, D, S, Sigma.ct, iter.max = 1000, nu = 1e-04, eps = 0.01,
                           centered = FALSE, normalize = FALSE)
{
    if (length(S) != ncol(D)) {
        common.cell.type = intersect(colnames(D), names(S))
        if (length(common.cell.type) <= 1) {
            stop("Not enough cell types!")
        }
        D = D[, match(common.cell.type, colnames(D))]
        S = S[match(common.cell.type, names(S))]
    }
    k = ncol(D)
    common.gene = intersect(names(Y), rownames(D))
    common.gene = intersect(common.gene, colnames(Sigma.ct))
    if (length(common.gene) < 0.1 * min(length(Y), nrow(D), ncol(Sigma.ct))) {
        stop("Not enough common genes!")
    }
    Y = Y[match(common.gene, names(Y))]
    D = D[match(common.gene, rownames(D)), ]
    Sigma.ct = Sigma.ct[, match(common.gene, colnames(Sigma.ct))]
    X = D
    if (centered) {
        X = X - mean(X)
        Y = Y - mean(Y)
    }
    if (normalize) {
        X = X/sd(as.vector(X))
        S = S * sd(as.vector(X))
        Y = Y/sd(Y)
    }
    else {
        Y = Y * 100
    }
    lm.D = music.basic.ct(Y, X, S, Sigma.ct, iter.max = iter.max,
                          nu = nu, eps = eps)
    return(lm.D)
}

music.basic <- function (Y, X, S, Sigma, iter.max, nu, eps)
{
    k = ncol(X)
    lm.D = nnls(X, Y)
    r = resid(lm.D)
    weight.gene = 1/(nu + r^2 + colSums((lm.D$x * S)^2 * t(Sigma)))
    Y.weight = Y * sqrt(weight.gene)
    D.weight = sweep(X, 1, sqrt(weight.gene), "*")
    lm.D.weight = nnls(D.weight, Y.weight)
    p.weight = lm.D.weight$x/sum(lm.D.weight$x)
    p.weight.iter = p.weight
    r = resid(lm.D.weight)
    for (iter in 1:iter.max) {
        weight.gene = 1/(nu + r^2 + colSums((lm.D.weight$x *
                                                 S)^2 * t(Sigma)))
        Y.weight = Y * sqrt(weight.gene)
        D.weight = X * as.matrix(sqrt(weight.gene))[, rep(1,
                                                          k)]
        lm.D.weight = nnls(D.weight, Y.weight)
        p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
        r.new = resid(lm.D.weight)
        if (sum(abs(p.weight.new - p.weight)) < eps) {
            p.weight = p.weight.new
            r = r.new
            R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)
            fitted = X %*% as.matrix(lm.D.weight$x)
            var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
            return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x,
                        fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
                        p.weight = p.weight, q.weight = lm.D.weight$x,
                        fit.weight = fitted, resid.weight = Y - X %*%
                            as.matrix(lm.D.weight$x), weight.gene = weight.gene,
                        converge = paste0("Converge at ", iter), rsd = r,
                        R.squared = R.squared, var.p = var.p))
        }
        p.weight = p.weight.new
        r = r.new
    }
    fitted = X %*% as.matrix(lm.D.weight$x)
    R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)
    var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
    return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x,
                fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), p.weight = p.weight,
                q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight = Y -
                    X %*% as.matrix(lm.D.weight$x), weight.gene = weight.gene,
                converge = "Reach Maxiter", rsd = r, R.squared = R.squared,
                var.p = var.p))
}

#' @importFrom nnls nnls
#' @importFrom stats resid
music.basic.ct <- function (Y, X, S, Sigma.ct, iter.max, nu, eps)
{
    k = ncol(X)
    lm.D = nnls(X, Y)
    r = resid(lm.D)
    weight.gene = 1/(nu + r^2 + weight.cal.ct(lm.D$x * S, Sigma.ct))
    Y.weight = Y * sqrt(weight.gene)
    D.weight = sweep(X, 1, sqrt(weight.gene), "*")
    lm.D.weight = nnls(D.weight, Y.weight)
    p.weight = lm.D.weight$x/sum(lm.D.weight$x)
    p.weight.iter = p.weight
    r = resid(lm.D.weight)
    for (iter in 1:iter.max) {
        weight.gene = 1/(nu + r^2 + weight.cal.ct(lm.D$x * S,
                                                  Sigma.ct))
        Y.weight = Y * sqrt(weight.gene)
        D.weight = X * as.matrix(sqrt(weight.gene))[, rep(1,
                                                          k)]
        lm.D.weight = nnls(D.weight, Y.weight)
        p.weight.new = lm.D.weight$x/sum(lm.D.weight$x)
        r.new = resid(lm.D.weight)
        if (sum(abs(p.weight.new - p.weight)) < eps) {
            p.weight = p.weight.new
            r = r.new
            R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)
            fitted = X %*% as.matrix(lm.D.weight$x)
            var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
            return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x,
                        fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D),
                        p.weight = p.weight, q.weight = lm.D.weight$x,
                        fit.weight = fitted, resid.weight = Y - X %*%
                            as.matrix(lm.D.weight$x), weight.gene = weight.gene,
                        converge = paste0("Converge at ", iter), rsd = r,
                        R.squared = R.squared, var.p = var.p))
            break
        }
        p.weight = p.weight.new
        r = r.new
    }
    fitted = X %*% as.matrix(lm.D.weight$x)
    R.squared = 1 - var(Y - X %*% as.matrix(lm.D.weight$x))/var(Y)
    var.p = diag(solve(t(D.weight) %*% D.weight)) * mean(r^2)/sum(lm.D.weight$x)^2
    return(list(p.nnls = (lm.D$x)/sum(lm.D$x), q.nnls = lm.D$x,
                fit.nnls = fitted(lm.D), resid.nnls = resid(lm.D), p.weight = p.weight,
                q.weight = lm.D.weight$x, fit.weight = fitted, resid.weight = Y -
                    X %*% as.matrix(lm.D.weight$x), weight.gene = weight.gene,
                converge = "Reach Maxiter", rsd = r, R.squared = R.squared,
                var.p = var.p))
}

my.rowMeans <- function (x, na.rm = FALSE, dims = 1L)
{
    if (is.data.frame(x))
        x <- as.matrix(x)
    if (length(dn <- dim(x)) < 2L) {
        return(x)
    }
    if (!is.array(x) || length(dn <- dim(x)) < 2L)
        stop("'x' must be an array of at least two dimensions")
    if (dims < 1L || dims > length(dn) - 1L)
        stop("invalid 'dims'")
    p <- prod(dn[-(id <- seq_len(dims))])
    dn <- dn[id]
    z <- if (is.complex(x))
        .Internal(rowMeans(Re(x), prod(dn), p, na.rm)) + (0 +
                                                              (0+1i)) * .Internal(rowMeans(Im(x), prod(dn), p,
                                                                                           na.rm))
    else .Internal(rowMeans(x, prod(dn), p, na.rm))
    if (length(dn) > 1L) {
        dim(z) <- dn
        dimnames(z) <- dimnames(x)[id]
    }
    else names(z) <- dimnames(x)[[1L]]
    z
}
relative.ab <- function (X, by.col = TRUE)
{
    if (sum(X < 0) > 0) {
        stop("Negative entry appears!")
    }
    if (by.col == T) {
        RX = sweep(X, 2, colSums(X), "/")
    }
    else {
        RX = sweep(X, 1, rowSums(X), "/")
    }
    return(RX)
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
