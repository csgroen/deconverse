#' Construct the mean gene expression basis matrix (B) from `screference`
#' object using the method implemented in the CARD package
#'
#' @param scref an object of `screference`
#'
#' @importFrom Matrix colSums
#'
#' @return a Basis matrix
#' @export
card_scref <- function(scref) {
    # Adapt to scref input
    countMat <- as(scref$seurat_obj@assays$RNA@counts,"sparseMatrix")
    ct.select <- scref$populations
    sample.id <- scref$seurat_obj$batch_id
    ct.id <- scref$seurat_obj$annot_id

    # Run CARD internal
    ct_sample.id <- paste(ct.id, sample.id, sep = "$*$")
    colSums_countMat <- colSums(countMat)
    colSums_countMat_Ct = aggregate(colSums_countMat ~ ct.id + sample.id, FUN = 'sum')
    colSums_countMat_Ct_wide = reshape(colSums_countMat_Ct, idvar = "sample.id", timevar = "ct.id", direction = "wide")
    colnames(colSums_countMat_Ct_wide) = gsub("colSums_countMat.","",colnames(colSums_countMat_Ct_wide))
    rownames(colSums_countMat_Ct_wide) = colSums_countMat_Ct_wide$sample.id
    colSums_countMat_Ct_wide$sample.id <- NULL
    tbl <- table(sample.id,ct.id)
    colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[,match(colnames(tbl),colnames(colSums_countMat_Ct_wide))]
    colSums_countMat_Ct_wide = colSums_countMat_Ct_wide[match(rownames(tbl),rownames(colSums_countMat_Ct_wide)),]
    S_JK <- colSums_countMat_Ct_wide / tbl
    S_JK <- as.matrix(S_JK)
    S_JK[S_JK == 0] = NA
    S_JK[!is.finite(S_JK)] = NA
    S = colMeans(S_JK, na.rm = TRUE)
    S = S[match(unique(ct.id),names(S))]
    if(nrow(countMat) > 10000 & ncol(countMat) > 50000){ ### to save memory
        seqID = seq(1,nrow(countMat),by = 10000)
        Theta_S_rowMean = NULL
        for(igs in seqID){
            if(igs != seqID[length(seqID)]){
                Theta_S_rowMean_Tmp <- .rowGrpMeans(as.matrix(countMat[(igs:(igs+9999)),]), grp = ct_sample.id, na.rm = TRUE)
            }else{
                Theta_S_rowMean_Tmp <- .rowGrpMeans(as.matrix(countMat[igs:nrow(countMat),]), grp = ct_sample.id, na.rm = TRUE)

            }
            Theta_S_rowMean <- rbind(Theta_S_rowMean,Theta_S_rowMean_Tmp)

        }
    } else {
        Theta_S_rowMean <- .rowGrpMeans(as.matrix(countMat), grp = ct_sample.id, na.rm = TRUE)
    }
    tbl_sample = table(ct_sample.id)
    tbl_sample = tbl_sample[match(colnames(Theta_S_rowMean),names(tbl_sample))]
    Theta_S_rowSums <- sweep(Theta_S_rowMean,2,tbl_sample,"*")
    Theta_S <- sweep(Theta_S_rowSums,2,colSums(Theta_S_rowSums),"/")
    grp <- sapply(strsplit(colnames(Theta_S),split="$*$",fixed = TRUE),"[",1)
    Theta = .rowGrpMeans(Theta_S, grp = grp, na.rm = TRUE)
    Theta = Theta[,match(unique(ct.id),colnames(Theta))]
    S = S[match(colnames(Theta),names(S))]
    basis = sweep(Theta,2,S,"*")
    colnames(basis) = colnames(Theta)
    rownames(basis) = rownames(Theta)

    # Selecting informative genes (moved from CARD_deconvolution)
    basis <- .select_informative(basis, ct.select, countMat, ct.id)

    return(basis)
}

#' CARD deconvolution using an `screference`
#'
#' @param spatial_obj an `Seurat` object with a Spatial assay
#' @param scref an `screference` object containing single-cell reference data
#' and pre-computed references
#' @importFrom stringr str_subset
#' @importFrom rdist rdist
#' @importFrom gtools rdirichlet
#' @importFrom CARD CARDref
#' @export
card_deconvolute <- function(spatial_obj, scref) {
    # Adapting to spatial_obj and scref input
    assert(class(scref) == "screference")
    assert(class(spatial_obj) == "Seurat")
    assert("Spatial" %in% names(spatial_obj@assays))

    B <- scref$cached_results$card[,scref$populations]
    spatial_count <- spatial_obj@assays$Spatial@counts
    ct.select <- scref$populations

    commonGene <- intersect(rownames(spatial_count),rownames(B))
    commonGene <- str_subset(commonGene, "MT-|mt-", negate = TRUE)

    # Card internals
    ##### match the common gene names
    spatial_count = spatial_count[order(rownames(spatial_count)),]
    B = B[order(rownames(B)),]
    B = B[rownames(B) %in% commonGene,]
    spatial_count = spatial_count[rownames(spatial_count) %in% commonGene,]
    ##### filter out non expressed genes or cells again
    spatial_count = spatial_count[rowSums(spatial_count) > 0,]
    spatial_count = spatial_count[,colSums(spatial_count) > 0]
    ##### normalize count data
    colsumvec = colSums(spatial_count)
    spatial_count_norm = sweep(spatial_count,2,colsumvec,"/")
    B = B[rownames(B) %in% rownames(spatial_count_norm),]
    B = B[match(rownames(spatial_count_norm),rownames(B)),]
    #### spatial location
    spatial_location = spatial_obj@images$slice1@coordinates[,2:3]
    #--- adapt
    colnames(spatial_location) <- c("x","y")
    #-- end
    spatial_location = spatial_location[rownames(spatial_location) %in% colnames(spatial_count_norm),]
    spatial_location = spatial_location[match(colnames(spatial_count_norm),rownames(spatial_location)),]
    ##### normalize the coordinates without changing the shape and relative position
    norm_cords = spatial_location[ ,c("x","y")]
    norm_cords$x = norm_cords$x - min(norm_cords$x)
    norm_cords$y = norm_cords$y - min(norm_cords$y)
    scaleFactor = max(norm_cords$x,norm_cords$y)
    norm_cords$x = norm_cords$x / scaleFactor
    norm_cords$y = norm_cords$y / scaleFactor
    ##### initialize the proportion matrix
    ED <- rdist(as.matrix(norm_cords))##Euclidean distance matrix
    message("-- CARD: Computing deconvolution...")
    set.seed(20200107)
    Vint1 = as.matrix(rdirichlet(ncol(spatial_count_norm), rep(10,ncol(B))))
    colnames(Vint1) = colnames(B)
    rownames(Vint1) = colnames(spatial_count_norm)
    b = rep(0,length(ct.select))
    ###### parameters that need to be set
    isigma = 0.1 ####construct Gaussian kernel with the default scale /length parameter to be 0.1
    epsilon = 1e-04  #### convergence epsion
    phi = c(0.01,0.1,0.3,0.5,0.7,0.9,0.99) #### grided values for phi
    kernel_mat <- exp(-ED^2 / (2 * isigma^2))
    diag(kernel_mat) <- 0
    rm(ED)
    rm(spatial_count)
    rm(norm_cords)
    gc()
    ###### scale the Xinput_norm and B to speed up the convergence.
    mean_X = mean(spatial_count_norm)
    mean_B = mean(B)
    spatial_count_norm = spatial_count_norm * 1e-01 / mean_X
    B = B * 1e-01 / mean_B
    gc()
    ResList = list()
    Obj = c()
    for(iphi in 1:length(phi)){
        res = CARDref(
            XinputIn = as.matrix(spatial_count_norm),
            UIn = as.matrix(B),
            WIn = as.matrix(kernel_mat),
            phiIn = phi[iphi],
            max_iterIn =1000,
            epsilonIn = epsilon,
            initV = Vint1,
            initb = rep(0,ncol(B)),
            initSigma_e2 = 0.1,
            initLambda = rep(10,length(ct.select)))
        rownames(res$V) = colnames(spatial_count_norm)
        colnames(res$V) = colnames(B)
        ResList[[iphi]] = res
        Obj = c(Obj,res$Obj)
    }
    Optimal = which(Obj == max(Obj))
    Optimal = Optimal[length(Optimal)] #### just in case if there are two equal objective function values
    OptimalPhi = phi[Optimal]
    OptimalRes = ResList[[Optimal]]
    message("-- CARD: End of deconvolution.")
    proportion_CARD <- sweep(OptimalRes$V,1,rowSums(OptimalRes$V),"/")
    # Adapt output --
    deconv_res <- proportion_CARD %>%
        as.data.frame() %>%
        rename_with(~ paste0("frac_", .)) %>%
        rownames_to_column("sample") %>%
        tibble() %>%
        mutate(method = "CARD")
    return(deconv_res)
}

.rowGrpMeans <- function(x, grp, na.rm=TRUE){
    if(!is.matrix(x) && !is.data.frame(x)) stop(" 'x' should be data.frame or matrix")
    if(length(dim(x)) !=2) stop(" 'x' should be data.frame or matrix of 2 dimensions")
    if(length(grp) != ncol(x)) stop(" 'grp' should be of length of number of cols in 'x'")
    if(length(grp) <1 || all(is.na(grp))) stop(" 'grp' appears to be empty or all NAs")
    if(!is.factor(grp)) grp <- as.factor(grp)
    if(!is.matrix(x)) x <- matrix(as.matrix(x), nrow=nrow(x), dimnames=if(length(dim(x)) >1) dimnames(x) else list(names(x),NULL))
    if(length(na.rm) !=1 || !is.logical(na.rm)) na.rm <- TRUE
    ## main
    out <- ..rowGrpMeans(x, grp, na.rm=na.rm)
    chNan <- is.nan(out)
    if(any(chNan)) out[which(chNan)] <- NA
    out }

..rowGrpMeans <- function(x, grp, na.replVa=NULL, na.rm=TRUE){
    ## determine means of rows conditional as (multiple) groups
    ## NAs (eg from counting data) can be replaced by specified value 'na.replVa', eg 0)
    ## 'grp' expected as factor !!
    if(!is.null(na.replVa)) x[is.na(x)] <- na.replVa
    grNa <- unique(..naOmit(as.character(grp)))
    out <- matrix(nrow=nrow(x), ncol=length(grNa), dimnames=list(rownames(x),grNa))
    for(i in 1:length(grNa)) {
        useC <- which(as.character(grp)==grNa[i])
        out[,i] <- if(length(useC) >1) base::rowMeans(x[,useC], na.rm=na.rm) else x[,useC] }
    out }

..naOmit <- function(x) {chNa <- is.na(x); if(all(chNa)) NULL else x[which(!chNa)]}


.select_informative <- function(Basis, ct.select, countMat, ct.id){
    #### log2 mean fold change >0.5
    gene1 = lapply(ct.select,function(ict) {
        rest = rowMeans(Basis[,colnames(Basis) != ict])
        FC = log((Basis[,ict] + 1e-06)) - log((rest + 1e-06))
        rownames(Basis)[FC > 1.25 & Basis[,ict] > 0]
    })
    gene1 = unique(unlist(gene1))
    # gene1 = intersect(gene1,commonGene)
    countMat = countMat[rownames(countMat) %in% gene1,]
    ##### only check the cell type that contains at least 2 cells
    ct.select = names(table(ct.id)[table(ct.id) > 1])
    sd_within = sapply(ct.select,function(ict){
        temp = countMat[, (ct.id == ict)]
        apply(temp,1,var) / apply(temp,1,mean)
    })
    ##### remove the outliers that have high dispersion across cell types
    gene2 = rownames(sd_within)[apply(sd_within,1,mean,na.rm = T) < quantile(apply(sd_within,1,mean,na.rm = T),prob = 0.99,na.rm = T)]
    Basis[gene2,]

    return(Basis)
}
