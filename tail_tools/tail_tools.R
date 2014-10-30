
library(nesoni)
library(limma)
library(edgeR)
library(MASS)
library(Matrix)
library(parallel)

variance.scorer <- function(design,data,counts, read.var,bio.var,bio.var.tail,bio.var.tail2) {
    vars <- read.var/counts + bio.var + bio.var.tail*data + bio.var.tail2*data*data

    n <- Null(design)
    tn <- t(n)
    
    x <- data %*% n

    total <- 0
    for(i in 1:nrow(data)) {
        covar <- tn %*% (vars[i,] * n)
        total <- total - log(det(covar)) - sum(x[i,] * solve(covar,x[i,]))
    }
    #cat(sprintf("%.1f %.5f %.5f %.5f -> %.1f\n", read.var,bio.var,bio.var.tail,bio.var.tail2,total))
    total
}

my.optim <- function(initial, func) {
    result <- optim(rep(0,length(initial)), function(x) func(exp(x)*initial))
    exp(result$par) * initial
}


elist.tails.for.fitnoise <- function(tails,counts,design,  genes=NULL,reads.cutoff=10, use.fitnoise=FALSE) {
# Returns an elist suitable for use with fitnoise models model.t.patseq or model.normal.patseq
# - filter low count genes
    tails <- as.matrix(tails)
    counts <- as.matrix(counts)

    good <- c()
    for(i in 1:nrow(tails))
        good[i] <- any(counts[i,]>=reads.cutoff) && rankMatrix(design[counts[i,]>=reads.cutoff,, drop=FALSE ]) == ncol(design)

    tail.elist <- new("EList")
    tail.elist$E <- tails[good,]
    tail.elist$other$counts <- counts[good,]
    if (!is.null(genes))
        tail.elist$genes <- genes[good,]

    tail.elist$info <- sprintf(
        paste(
            '%d of %d features kept after filtering\n',
            '(required sufficient samples with %d poly(A) reads to fit linear model)\n',
            sep=''), 
        sum(good),length(good),reads.cutoff)

    tail.elist
}


elist.tails <- function(tails,counts,design,  genes=NULL,reads.cutoff=10, use.fitnoise=FALSE) {
# Returns an elist
# - filter low count genes
# - give appropriate weights
    tails <- as.matrix(tails)
    counts <- as.matrix(counts)

    good <- c()
    for(i in 1:nrow(tails))
        good[i] <- any(counts[i,]>=reads.cutoff) && rankMatrix(design[counts[i,]>=reads.cutoff,, drop=FALSE ]) == ncol(design)

    allgood <- good & (row.apply(counts,min) > 0)
    
    cat('\nFitting variance model\n')
    opt <- my.optim(c(250.0,0.001),function (v) -variance.scorer(design,tails[allgood,],counts[allgood,],v[1],0.0,0.0,v[2]) )
    read.var <- opt[1]
    bio.var <- 0.0
    bio.var.tail <- 0.0
    bio.var.tail2 <- opt[2]
    
    weights <- 1.0 / (read.var/counts+bio.var+bio.var.tail*tails+bio.var.tail2*tails*tails)    
        
    tail.elist <- new("EList")
    tail.elist$E <- tails[good,]
    tail.elist$weights <- weights[good,]
    if (!is.null(genes))
        tail.elist$genes <- genes[good,]
    
    tail.elist$read.var <- read.var
    tail.elist$bio.var.tail2 <- bio.var.tail2
    tail.elist$info <- sprintf(
        paste(
            '%d of %d features kept after filtering\n',
            '(required sufficient samples with %d poly(A) reads to fit linear model)\n',
            'variance = %.1f^2 / polya_read_count + %.6f^2 * tail_length^2\n',
            sep=''), 
        sum(good),length(good),reads.cutoff,sqrt(read.var), sqrt(bio.var.tail2))
    
    tail.elist
}





