

scale_cols <- function(A,s) {
    assert_that(ncol(A) == length(s))
    t(t(A)*s)
}

# Least squares
# Where multiple solutions exist, one should be chosen arbitrarily
least_squares <- function(A,w,b) {
    #sw <- sqrt(w)
    # qr.solve(A*sw,b*sw)

    #decomp <- svd(A*sw)
    #as.vector( decomp$v %*% ((t(decomp$u) %*% (b*sw))/decomp$d) )

    WA <- A*w
    good <- apply(WA != 0, 2, any)
    WA <- WA[,good,drop=F]
    A <- A[,good,drop=F]
    tAWA <- crossprod(WA, A)
    tAWb <- as.vector(crossprod(WA,b))
    L <- chol(tAWA)

    result <- rep(0, length(good))
    result[good] <- backsolve(L, forwardsolve(L, tAWb, upper.tri=TRUE, transpose=TRUE))
    result
}

fit_all_cols <- function(x,y,w) {
    result <- lapply(seq_len(ncol(y)), function(i) {
        wi <- w[,i]
        present <- wi > 0
        least_squares(x[present,,drop=F],wi[present],y[present,i])
    })
    
    do.call(rbind, result)
}

fit_all_rows <- function(x,y,w) {
    result <- lapply(seq_len(nrow(y)), function(i) {
        wi <- w[i,]
        present <- wi > 0
        least_squares(x[present,,drop=F],wi[present],y[i,present])
    })
    
    do.call(rbind, result)
}

# Turn decomposition rows %*% t(cols) into an orthogonal version
# Furthermore, the columns of rows will have unit variance
orthogonalize_decomp <- function(rows,cols) {
    n <- nrow(rows)
    decomp <- svd(rows)
    rows <- decomp$u * sqrt(n)
    cols <- scale_cols(cols %*% decomp$v, decomp$d / sqrt(n))
    decomp <- svd(cols)
    cols <- scale_cols(decomp$u, decomp$d)
    rows <- rows %*% decomp$v
    list(rows=rows,cols=cols)
}

## Test:
# n <- 10
# m <- 7
# p <- 3
# A <- matrix(rnorm(n*p),nrow=n)
# B <- matrix(rnorm(m*p),nrow=m)
# decomp <- orthogonalize_decomp(A,B)
# (A %*% t(B)) - (decomp$rows %*% t(decomp$cols))


#' @export
elist_factors <- function(elist, p=2, design=NULL, use_varimax=TRUE, max_iter=100, tol=1e-5, verbose=TRUE, initial=NULL) {
    y <- elist$E
    weights <- elist$weights
    rownames(y) <- NULL
    colnames(y) <- NULL
    rownames(weights) <- NULL
    colnames(weights) <- NULL
    # Standardize missingness, no NAs
    weights[is.na(y)] <- 0
    y[weights == 0] <- 0

    n <- nrow(y)
    m <- ncol(y)
    assert_that(n >= p, m >= p)

    if (is.null(design))
        design <- cbind(intercept=rep(1,m))

    p_design <- ncol(design)
    p_total <- p_design+p
    ind_design <- seq_len(p_design)
    ind_factors <- p_design+seq_len(p)

    if (is.null(colnames(design)))
        colnames(design) <- paste0("design", seq_len(p_design))

    col_mat <- design
    if (!is.null(initial))
        col_mat <- cbind(col_mat, initial)
    col_mat <- cbind(col_mat, 
        matrix(rnorm(m*(p_total-ncol(col_mat))), nrow=m))

    # Centering to compare against
    row_mat_center <- fit_all_rows(design, y, weights)
    centered <- y - row_mat_center %*% t(design)
    ss_total <- sum(centered^2*weights)

    R2 <- -Inf

    for(i in seq_len(max_iter)) {
        start <- proc.time()["elapsed"]

        #gc()

        # Esure col_mat[,ind_factors] are orthogonal col_mat[,ind_design]
        #col_mat[,ind_factors] <- qr.Q(qr(col_mat))[,ind_factors,drop=F]

        # Update row_mat
        row_mat <- fit_all_rows(col_mat, y, weights)

        # Update col_mat
        centered <- y - row_mat[,seq_len(p_design),drop=F] %*% t(design)
        col_mat <- fit_all_cols(row_mat[,p_design+seq_len(p),drop=F], centered, weights)
        col_mat <- cbind(design, col_mat)

        # Make decomposition matrices orthogonal
        decomp <- orthogonalize_decomp(row_mat[,ind_factors,drop=F], col_mat[,ind_factors,drop=F])
        row_mat[,ind_factors] <- decomp$rows
        col_mat[,ind_factors] <- decomp$cols

        # Check R^2
        resid <- y - row_mat %*% t(col_mat)
        ss_resid <- sum(resid^2*weights) 
        #ss_total <- sum(centered^2*weights)
        ratio <- ss_resid/ss_total
        last_R2 <- R2
        R2 <- 1-ratio

        end <- proc.time()["elapsed"]
        if (verbose) {
            cat("Iteration",i,"R^2:",R2,"seconds:",end-start,"\n")
        }

        if (R2 - last_R2 <= tol) 
            break
    }

    # Apply varimax rotation
    if (p > 1 && use_varimax) {
        scaling <- sqrt(colMeans(col_mat[,ind_factors,drop=F]^2))
        rotation <- varimax(scale_cols(row_mat[,ind_factors,drop=FALSE],scaling), normalize=FALSE)
        row_mat[,ind_factors] <- row_mat[,ind_factors] %*% rotation$rotmat
        col_mat[,ind_factors] <- col_mat[,ind_factors] %*% rotation$rotmat
    }
    
    # Ensure largest factor first
    scaling <- sqrt(colMeans(col_mat[,ind_factors,drop=F]^2))
    reordering <- ind_factors[order(scaling, decreasing=TRUE)]
    row_mat[,ind_factors] <- row_mat[,reordering,drop=F]
    col_mat[,ind_factors] <- col_mat[,reordering,drop=F]

    # Ensure positive skew (outliers positive)
    flips <- ind_factors[colSums(col_mat[,ind_factors,drop=F]^3) < 0]
    row_mat[,flips] <- -row_mat[,flips,drop=F]
    col_mat[,flips] <- -col_mat[,flips,drop=F]

    # Use original row and column names
    # Give factors in the decomposition meaningful names
    rownames(row_mat) <- rownames(elist)
    rownames(col_mat) <- colnames(elist)
    colnames(row_mat) <- c(colnames(design), paste0("factor",seq_len(p)))
    colnames(col_mat) <- colnames(row_mat)

    list(row=row_mat, col=col_mat, R2=R2, iters=i)
}


#' Permute measurements weighted residuals
#'
#' @export
elist_permute <- function(elist, design=NULL) {
    if (is.null(design))
        design <- cbind(intercept=rep(1,ncol(elist)))

    result <- map(seq_len(nrow(elist)), function(i) {
        row <- rep(NA_real_,ncol(elist))
        present <- which(elist$weight[i,]>0)
        if (length(present) > 0) {
            permute <- sample.int(length(present))
            y <- elist$E[i,present]
            w <- elist$weights[i,present]
            design_present <- design[present,,drop=F]
            coef <- least_squares(design_present, w, y)
            m <- as.vector(design_present %*% coef)    
            row[present[permute]] <- (y-m)*sqrt(w/w[permute])+m[permute]
        }
        row
    })

    elist$E <- do.call(rbind, result)
    elist 
}

#' Generate a random normally distributed version of an EList object.
#' Weights are used to choose the standard deviation of values generated.
#'
#' @export
elist_randomize <- function(elist) {
    elist$E <- 
        matrix(rnorm(nrow(elist)*ncol(elist)), nrow=nrow(elist)) /
        sqrt(elist$weights)
    elist
}


#' Progressively produce factor analyses with increasing numbers of factors
#'
#' @export
elist_factors_seq <- function(elist, p=2, design=NULL, use_varimax=TRUE, max_iter=100, tol=1e-6, verbose=TRUE) {
    result <- vector("list",p)

    if (verbose)
        cat("Finding 1 factor\n")

    result[[1]] <- elist_factors(elist, p=1, 
        design=design, use_varimax=use_varimax, max_iter=max_iter, tol=tol, verbose=verbose)
    
    for(i in seq_len(p-1)+1) {
        if (verbose)
            cat("Finding",i,"factors\n")

        n <- ncol(result[[i-1]]$col)
        result[[i]] <- elist_factors(elist, p=i, 
            design=design, use_varimax=use_varimax, max_iter=max_iter, tol=tol, verbose=verbose) #,
            #This may actually produce worse results:
            #initial=result[[i-1]]$col[,n-(i-1)+seq_len(i-1),drop=F])
    }

    result
}




