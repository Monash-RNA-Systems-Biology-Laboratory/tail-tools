

fit_all <- function(x,y,w) {
    n <- ncol(y)

    result <- lapply(seq_len(n), function(i) {
        wi <- w[,i]
        present <- wi > 0
        swpi <- sqrt(wi[present])
        qr.solve(x[present,,drop=F]*swpi,y[present,i]*swpi)
    })

    do.call(rbind, result)
}


#' @export
elist_factors <- function(elist, p=2, design=NULL, use_varimax=TRUE, max_iter=100, tol=1e-5, verbose=TRUE) {
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

    ind_design <- seq_len(p_design)
    ind_factors <- p_design+seq_len(p)

    df_null <- sum(weights>0)-n*p_design
    df <- df_null-n*p-m*p

    col_mat <- matrix(rnorm(m*p), ncol=p)
    col_mat <- cbind(design, col_mat)

    R2 <- -Inf

    for(i in seq_len(max_iter)) {
        # Update row_mat
        row_mat <- fit_all(col_mat, t(y), t(weights))

        # Make sure factors are orthogonal, sensibly ordered
        part <- row_mat[,ind_factors,drop=F]
        scaling <- sqrt(colMeans(col_mat[,ind_factors,drop=F]^2))
        part <- t(t(part)*scaling)
        decomp <- svd(part)
        row_mat[,ind_factors] <- decomp$u
        
        # Update col_mat
        centered <- y - row_mat[,seq_len(p_design),drop=F] %*% t(design)
        col_mat <- fit_all(row_mat[,p_design+seq_len(p),drop=F], centered, weights)
        col_mat <- cbind(design, col_mat)

        resid <- y - row_mat %*% t(col_mat)
        ss_resid <- sum(resid^2*weights) 
        ss_total <- sum(centered^2*weights)
        ratio <- ss_resid/ss_total
        last_R2 <- R2
        R2 <- 1-ratio
        R2adj <- 1-ratio*(df_null)/(df)
        if (verbose) {
            cat("Iteration",i,"R^2:",R2,"Adjusted R^2:",R2adj,"\n")
        }

        if (i >= 5 && R2 - last_R2 <= tol) 
            break
    }

    if (p > 1 && use_varimax) {
        scaling <- sqrt(colMeans(col_mat[,ind_factors,drop=F]^2))
        rotation <- varimax(t(t(row_mat[,ind_factors,drop=FALSE])*scaling), normalize=FALSE)
        row_mat[,ind_factors] <- row_mat[,ind_factors] %*% rotation$rotmat
        col_mat[,ind_factors] <- col_mat[,ind_factors] %*% rotation$rotmat
    }
    
    # Ensure largest factor first
    scaling <- sqrt(colMeans(col_mat[,ind_factors,drop=F]^2))
    reordering <- ind_factors[order(scaling, decreasing=TRUE)]
    row_mat[,ind_factors] <- row_mat[,reordering,drop=F]
    col_mat[,ind_factors] <- col_mat[,reordering,drop=F]   

    rownames(row_mat) <- rownames(elist)
    rownames(col_mat) <- colnames(elist)
    if (is.null(colnames(design)))
        colnames(design) <- paste0("design", seq_len(p_design))
    colnames(row_mat) <- c(colnames(design), paste0("factor",seq_len(p)))
    colnames(col_mat) <- colnames(row_mat)

    list(row=row_mat, col=col_mat, R2=R2, R2adj=R2adj, iters=i)
}
