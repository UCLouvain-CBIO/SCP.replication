##' @title ComBat v3.34 batch correction
##'
##' Older version of ComBat (sva version 3.34.0). This is for
##' compatibility with the replication of SCoPE2.
##' 
##' ComBat allows users to adjust for batch effects in datasets where
##' the batch covariate is known, using methodology described in 
##' Johnson et al. 2007. It uses either parametric or non-parametric 
##' empirical Bayes frameworks for adjusting data for batch effects.
##' Users are returned an expression matrix that has been corrected
##' for batch effects. The input data are assumed to be cleaned and
##' normalized before batch effect removal.
##' 
##' @param dat Genomic measure matrix (dimensions probe x sample) - 
##'     for example, expression matrix
##'
##' @param batch Batch covariate (only one batch allowed)
##'
##' @param mod Model matrix for outcome of interest and other 
##'     covariates besides batch
##'
##' @param par.prior (Optional) TRUE indicates parametric adjustments
##'     will be used, FALSE indicates non-parametric adjustments will
##'     be used
##'
##' @param prior.plots (Optional) TRUE give prior plots with black as
##'      a kernel estimate of the empirical batch effect density and 
##'      red as the parametric
##'
##' @param mean.only (Optional) FALSE If TRUE ComBat only corrects the
##'     mean of the batch effect (no scale adjustment)
##'
##' @param ref.batch (Optional) NULL If given, will use the selected
##'     batch as a reference for batch adjustment.
##'
##' @param BPPARAM (Optional) BiocParallelParam for parallel operation
##'
##' @export
##'
##' @import sva
##'
##' @return data A probe x sample genomic measure matrix, adjusted for batch effects.
##' 
##' @examples 
##' library(bladderbatch)
##' data(bladderdata)
##' dat <- bladderEset[1:50,]
##' 
##' pheno = pData(dat)
##' edata = exprs(dat)
##' batch = pheno$batch
##' mod = model.matrix(~as.factor(cancer), data=pheno)
##' 
##' # parametric adjustment
##' combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
##' 
##' # non-parametric adjustment, mean-only version
##' combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)
##' 
##' # reference-batch version, with covariates
##' combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
##' 
ComBatv3.34 <- function (dat, batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE,
                         mean.only = FALSE, ref.batch = NULL, BPPARAM = bpparam("SerialParam")) {
    stopifnot(require(sva)) ## works for svq_3.37.0
    if (mean.only) {
        message("Using the 'mean only' version of ComBat")
    }
    if (length(dim(batch)) > 1) {
        stop("This version of ComBat only allows one batch variable")
    }
    batch <- as.factor(batch)
    batchmod <- model.matrix(~-1 + batch)
    if (!is.null(ref.batch)) {
        if (!(ref.batch %in% levels(batch))) {
            stop("reference level ref.batch is not one of the levels of the batch variable")
        }
        cat("Using batch =", ref.batch, "as a reference batch (this batch won't change)\n")
        ref <- which(levels(as.factor(batch)) == ref.batch)
        batchmod[, ref] <- 1
    }
    else {
        ref <- NULL
    }
    message("Found", nlevels(batch), "batches")
    n.batch <- nlevels(batch)
    batches <- list()
    for (i in 1:n.batch) {
        batches[[i]] <- which(batch == levels(batch)[i])
    }
    n.batches <- sapply(batches, length)
    if (any(n.batches == 1)) {
        mean.only = TRUE
        message("Note: one batch has only one sample, setting mean.only=TRUE")
    }
    n.array <- sum(n.batches)
    design <- cbind(batchmod, mod)
    check <- apply(design, 2, function(x) all(x == 1))
    if (!is.null(ref)) {
        check[ref] <- FALSE
    }
    design <- as.matrix(design[, !check])
    message("Adjusting for", ncol(design) - ncol(batchmod),
            "covariate(s) or covariate level(s)")
    if (qr(design)$rank < ncol(design)) {
        if (ncol(design) == (n.batch + 1)) {
            stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
        }
        if (ncol(design) > (n.batch + 1)) {
            if ((qr(design[, -c(1:n.batch)])$rank < ncol(design[,
                                                                -c(1:n.batch)]))) {
                stop("The covariates are confounded! Please remove one or more of the covariates so the design is not confounded")
            }
            else {
                stop("At least one covariate is confounded with batch! Please remove confounded covariates and rerun ComBat")
            }
        }
    }
    NAs <- any(is.na(dat))
    if (NAs) {
        message(c("Found", sum(is.na(dat)), "Missing Data Values"),
                sep = " ")
    }
    cat("Standardizing Data across genes\n")
    if (!NAs) {
        B.hat <- solve(crossprod(design), tcrossprod(t(design),
                                                     as.matrix(dat)))
    }
    else {
        B.hat <- apply(dat, 1, Beta.NA, design)
    }
    if (!is.null(ref.batch)) {
        grand.mean <- t(B.hat[ref, ])
    }
    else {
        grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,
                                                         ])
    }
    if (!NAs) {
        if (!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- ((ref.dat - t(design[batches[[ref]],
                                               ] %*% B.hat))^2) %*% rep(1/n.batches[ref], n.batches[ref])
        }
        else {
            var.pooled <- ((dat - t(design %*% B.hat))^2) %*%
                rep(1/n.array, n.array)
        }
    }
    else {
        if (!is.null(ref.batch)) {
            ref.dat <- dat[, batches[[ref]]]
            var.pooled <- rowVars(ref.dat - t(design[batches[[ref]],
                                                     ] %*% B.hat), na.rm = TRUE)
        }
        else {
            var.pooled <- rowVars(dat - t(design %*% B.hat),
                                  na.rm = TRUE)
        }
    }
    stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
    if (!is.null(design)) {
        tmp <- design
        tmp[, c(1:n.batch)] <- 0
        stand.mean <- stand.mean + t(tmp %*% B.hat)
    }
    s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1,
                                                             n.array)))
    message("Fitting L/S model and finding priors")
    batch.design <- design[, 1:n.batch]
    if (!NAs) {
        gamma.hat <- solve(crossprod(batch.design), tcrossprod(t(batch.design),
                                                               as.matrix(s.data)))
    }
    else {
        gamma.hat <- apply(s.data, 1, Beta.NA, batch.design)
    }
    delta.hat <- NULL
    for (i in batches) {
        if (mean.only == TRUE) {
            delta.hat <- rbind(delta.hat, rep(1, nrow(s.data)))
        }
        else {
            delta.hat <- rbind(delta.hat, rowVars(s.data[, i],
                                                  na.rm = TRUE))
        }
    }
    gamma.bar <- rowMeans(gamma.hat)
    t2 <- rowVars(gamma.hat)
    a.prior <- apply(delta.hat, 1, sva:::aprior)
    b.prior <- apply(delta.hat, 1, sva:::bprior)
    if (prior.plots && par.prior) {
        par(mfrow = c(2, 2))
        tmp <- density(gamma.hat[1, ])
        plot(tmp, type = "l", main = expression(paste("Density Plot of First Batch ",
                                                      hat(gamma))))
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        lines(xx, dnorm(xx, gamma.bar[1], sqrt(t2[1])), col = 2)
        qqnorm(gamma.hat[1, ], main = expression(paste("Normal Q-Q Plot of First Batch ",
                                                       hat(gamma))))
        qqline(gamma.hat[1, ], col = 2)
        tmp <- density(delta.hat[1, ])
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        tmp1 <- list(x = xx, y = dinvgamma(xx, a.prior[1], b.prior[1]))
        plot(tmp, typ = "l", ylim = c(0, max(tmp$y, tmp1$y)),
             main = expression(paste("Density Plot of First Batch ",
                                     hat(delta))))
        lines(tmp1, col = 2)
        invgam <- 1/qgamma(1 - ppoints(ncol(delta.hat)), a.prior[1],
                           b.prior[1])
        qqplot(invgam, delta.hat[1, ], main = expression(paste("Inverse Gamma Q-Q Plot of First Batch ",
                                                               hat(delta))), ylab = "Sample Quantiles", xlab = "Theoretical Quantiles")
        lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
    }
    gamma.star <- delta.star <- matrix(NA, nrow = n.batch, ncol = nrow(s.data))
    if (par.prior) {
        message("Finding parametric adjustments")
        results <- bplapply(1:n.batch, function(i) {
            if (mean.only) {
                gamma.star <- postmean(gamma.hat[i, ], gamma.bar[i],
                                       1, 1, t2[i])
                delta.star <- rep(1, nrow(s.data))
            }
            else {
                temp <- sva:::it.sol(s.data[, batches[[i]]], gamma.hat[i,
                                                                       ], delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i],
                                     b.prior[i])
                gamma.star <- temp[1, ]
                delta.star <- temp[2, ]
            }
            list(gamma.star = gamma.star, delta.star = delta.star)
        }, BPPARAM = BPPARAM)
        for (i in 1:n.batch) {
            gamma.star[i, ] <- results[[i]]$gamma.star
            delta.star[i, ] <- results[[i]]$delta.star
        }
    }
    else {
        message("Finding nonparametric adjustments")
        results <- bplapply(1:n.batch, function(i) {
            if (mean.only) {
                delta.hat[i, ] = 1
            }
            temp <- int.eprior(as.matrix(s.data[, batches[[i]]]),
                               gamma.hat[i, ], delta.hat[i, ])
            list(gamma.star = temp[1, ], delta.star = temp[2,
                                                           ])
        }, BPPARAM = BPPARAM)
        for (i in 1:n.batch) {
            gamma.star[i, ] <- results[[i]]$gamma.star
            delta.star[i, ] <- results[[i]]$delta.star
        }
    }
    if (!is.null(ref.batch)) {
        gamma.star[ref, ] <- 0
        delta.star[ref, ] <- 1
    }
    message("Adjusting the Data\n")
    bayesdata <- s.data
    j <- 1
    for (i in batches) {
        bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i,
                                                           ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1,
                                                                                                               n.batches[j])))
        j <- j + 1
    }
    bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1,
                                                          n.array)))) + stand.mean
    if (!is.null(ref.batch)) {
        bayesdata[, batches[[ref]]] <- dat[, batches[[ref]]]
    }
    return(bayesdata)
}

##' @title Count Unique Features
##' 
##' This function counts the number of unique features per sample. A 
##' grouping structure can be provided to count higher level features
##' from assays, for example counting the number of unique proteins 
##' from PSM data. 
##'
##' @param object An object of class `QFeatures`.
##'
##' @param i  A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##'
##' @param groupBy A `character(1)` indicating the variable name in 
##'     the `rowData` that contains the grouping variable, for 
##'     instance to count the unique number of peptides or proteins 
##'     expressed in each samples (column). If `groupBy` is missing, 
##'     the number of non zero elements per sample will be stored.
##'     
##' @param colDataName A `character(1)` giving the name of the new 
##'     variable in the `colData` where the number of unique features
##'     will be stored. The name cannot already exist in the 
##'     `colData`.
##'
##' @export
##'
##' @return An object of class `QFeatures`.
##' 
countUniqueFeatures <- function(object,
                                i,
                                groupBy = NULL,
                                colDataName = "count") {
    ## Check the colData does not already contain the name
    if (colDataName %in% colnames(colData(object)))
        stop("'", colDataName, "' is already present in the colData.")

    snames <- unlist(colnames(object)[i])
    ## Avoid that a sample is contained in 2 different assays
    if (anyDuplicated(snames))
        stop("The same sample is present in multupl")

    ## Initialize the vector containing the feature counts
    fcounts <- vector(length = length(snames), mode = "integer")
    names(fcounts) <- snames

    if (is.null(groupBy)) {
        ## If no  grouping is supplied, count the non-missing
        ## features per sample
        for (ii in i) {
            fcount <- apply(assay(object[[ii]]), 2, function(x) {
                sum(!is.na(x))
            })
            fcounts[names(fcount)] <- fcount
        }
    } else {
        ## Count the number of unique entries of groupBy
        for (ii in i) {
            fcount <- apply(assay(object[[ii]]), 2, function(x) {
                length(unique(rowData(object[[ii]])[!is.na(x), groupBy]))
            })
            fcounts[names(fcount)] <- fcount
        }
    }

    ## Store the counts in the colData
    colData(object)[names(fcounts), colDataName] <- fcounts
    object
}

##' @title KNN imputation (for replication)
##'
##' The function imputes missing data using the k-nearest neighbours 
##' (KNN) approach using Euclidean distance as a similarity measure 
##' between the cells. This function was provided as part of the 
##' SCoPE2 publication (Specht et al. 2021). 
##' 
##' 
##' @param object A `QFeatures` object
##'
##' @param i  A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##'
##' @param name A `character(1)` naming the new assay name. Default is
##' `KNNimputedAssay`. 
##'
##' @param k An `integer(1)` giving the number of neighbours to use.
##'
##' @export
##'
##' @return An object of class `QFeatures` containing an extra assay
##'     with imputed values. 
##' 
##' @references Specht, Harrison, Edward Emmott, Aleksandra A. Petelski,
##'     R. Gray Huffman, David H. Perlman, Marco Serra, Peter Kharchenko, 
##'     Antonius Koller, and Nikolai Slavov. 2021. “Single-Cell Proteomic
##'     and Transcriptomic Analysis of Macrophage Heterogeneity Using 
##'     SCoPE2.” Genome Biology 22 (1): 50.
##' 
##' @source The implementation of the imputation was retrieved from
##'     https://github.com/SlavovLab/SCoPE2
##' 
scp_imputeKNN <- function(object, i, name = "KNNimputedAssay", k = 3){

    oldi <- i
    exp <- object[[i]]
    dat <- assay(exp)

    # Create a copy of the data, NA values to be filled in later
    dat.imp<-dat

    # Calculate similarity metrics for all column pairs (default is Euclidean distance)
    dist.mat<-as.matrix( dist(t(dat)) )
    #dist.mat<-as.matrix(as.dist( dist.cosine(t(dat)) ))

    # Column names of the similarity matrix, same as data matrix
    cnames<-colnames(dist.mat)

    # For each column in the data...
    for(X in cnames){

        # Find the distances of all other columns to that column
        distances<-dist.mat[, X]

        # Reorder the distances, smallest to largest (this will reorder the column names as well)
        distances.ordered<-distances[order(distances, decreasing = F)]

        # Reorder the data matrix columns, smallest distance to largest from the column of interest
        # Obviously, first column will be the column of interest, column X
        dat.reordered<-dat[ , names(distances.ordered ) ]

        # Take the values in the column of interest
        vec<-dat[, X]

        # Which entries are missing and need to be imputed...
        na.index<-which( is.na(vec) )

        # For each of the missing entries (rows) in column X...
        for(i in na.index){

            # Find the most similar columns that have a non-NA value in this row
            closest.columns<-names( which( !is.na(dat.reordered[i, ])  ) )

            # If there are more than k such columns, take the first k most similar
            if( length(closest.columns)>k ){
                # Replace NA in column X with the mean the same row in k of the most similar columns
                vec[i]<-mean( dat[ i, closest.columns[1:k] ] )
            }

            # If there are less that or equal to k columns, take all the columns
            if( length(closest.columns)<=k ){
                # Replace NA in column X with the mean the same row in all of the most similar columns
                vec[i]<-mean( dat[ i, closest.columns ])
            }
        }
        # Populate a the matrix with the new, imputed values
        dat.imp[,X]<-vec
    }

    assay(exp) <- dat.imp
    object <- addAssay(object, exp, name = name)
    addAssayLinkOneToOne(object, from = oldi, to = name)
}


##' Replace columns in the rowData
##' 
##' The function replaces columns in the rowData of assays from a 
##' `QFeatures` object. This function can also be used to remove column
##' in the rowData (see examples).
##'
##' @param object An object of class `QFeatures`
##' 
##' @param i  A `numeric()` or `character()` vector indicating from 
##'     which assays the `rowData` should be taken.
##'     
##' @param col A `numeric()` or `character()` vector indicating the 
##'     column(s) to replace.
##'     
##' @param value A `list()` of same length as `i`. If `i` is provided 
##'     as a character vector, the elements of `value` must be named
##'     after `i`. Each element should contain the replacement values to insert in the rowData. If
##'     `length(col) > 1`, each element should be a matrix with the 
##'     same rows as its corresponding assay and the same number of 
##'     columns as the length of col.
##'
##' @return An object of class `QFeatures`
##' 
##' @export
##'
##' @examples
##' data("feat2")
##' rowData(feat2)[[1]]
##' ## Add a new column to all assay
##' values <- lapply(rowData(feat2), function(x) rep("bar", nrow(x)))
##' feat2 <- replaceRowDataCols(feat2, 1:3, "foo", values)
##' rowData(feat2)[[1]]
##' ## Update column
##' values <- lapply(rowData(feat2), function(x) rep("barNew", nrow(x)))
##' feat2 <- replaceRowDataCols(feat2, 1:3, "foo", values)
##' rowData(feat2)[[1]]
##' ## Remove column
##' values <- lapply(1:3, function(x) NULL)
##' feat2 <- replaceRowDataCols(feat2, 1:3, "foo", values)
##' rowData(feat2)[[1]]
replaceRowDataCols <- function(object, i, col, value) {
    stopifnot(inherits(object, "QFeatures"))
    stopifnot(is.list(value))
    stopifnot(length(value) == length(i))
    if(is.character(i)) stopifnot(all(i %in% names(value)))
    el <- experiments(object)
    for (ii in i)
        rowData(el[[ii]])[, col] <- value[[ii]]
    BiocGenerics:::replaceSlots(object,
                                ExperimentList = el,
                                check = FALSE)
}


##' PCA plot suggested by SCoPE2
##'
##' The function performs a weighted principal component analysis (PCA)
##' as suggested by Specht et al. The PCA is performed on the 
##' correlation matrix where the rows (features) are weighted according
##' to the sum of the correlation with the other rows.
##' 
##' @param object A `SingleCellExperiment` object
##'
##' @param center A `logical(1)` indicating whether the columns of the
##'     weighted input matrix should be mean centered.
##'     
##' @param scale A `logical(1)` indicating whether the columns of the
##'     weighted input matrix should be scaled to unit variance.
##'     
##' @return An object of class `eigen` containing the computed 
##'     eigenvector and eigenvalues.
##' 
##' @export 
##'
##' @references 
##' 
##' Specht, Harrison, Edward Emmott, Aleksandra A. Petelski, R. Gray 
##' Huffman, David H. Perlman, Marco Serra, Peter Kharchenko, Antonius
##' Koller, and Nikolai Slavov. 2021. "Single-Cell Proteomic and
##' Transcriptomic Analysis of Macrophage Heterogeneity Using SCoPE2.”
##' Genome Biology 22 (1): 50.
##' [link to article](http://dx.doi.org/10.1186/s13059-021-02267-5),
##' [link to preprint](https://www.biorxiv.org/content/10.1101/665307v5)
##' 
##' @examples
##' data("feat2")
##' sce <- as(feat2[[3]], "SingleCellExperiment")
##' pcaSCoPE2(sce)
##' 
pcaSCoPE2 <- function(object, scale = FALSE, center = FALSE) {
    ## Extract the protein expression matrix
    X <- assay(object)
    ## Compute the weights
    w <- rowSums(cor(t(X))^2)
    ## Compute the PCs (code taken from the SCoPE2 script)
    Xw <- diag(w) %*%  X
    ## Optionally center and/or scale
    if (center)
        Xw <- sweep(Xw, 2, colMeans(Xw), "-")
    if (scale)
        Xw <- sweep(Xw, 2, colSds(Xw), "/")
    Xcor <- cor(Xw)
    eigen(Xcor)
}

