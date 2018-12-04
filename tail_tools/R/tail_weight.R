
#' Weight tail lengths
#'
#' Tail lengths are log2 transformed, and given weights based on read counts.
#'
#' Weights are calculated based on a variance model including a technical variation component inversely proportional to the number of reads plus a constant biological variance component.
#' 
#' Weighting is fine-tuned based on genes with at least min_reads in each sample.
#'
#' Features are filtered so that only features that can be fitted using samples with min_reads reads are used. (Samples with less than min_reads reads will still be used in the calculation, so this is stricter than absolutely necessary.)
#'
#' Returns an EList of weighted log2 tail lengths.
#'
#' @export
weighted_log2_tails <- function(tails, tail_counts, design, genes=NULL, min_reads=10, max_weight=10000, biovar=TRUE) {
    tails <- as.matrix(tails)
    tail_counts <- as.matrix(tail_counts)
    
    log2_tails <- log2(tails)
    log2_tails[tail_counts == 0] <- 0.0 #Remove NAs, will be weighted zero
    
    #all_present <- apply(tail_counts >= min_reads, 1, all)
    #if (!any(all_present)) 
    #    stop("No features with all samples exceeding minimum read count, while weighting tail lengths.")
    
    #ap_log2_tails <- log2_tails[all_present,,drop=F]
    #ap_tail_counts <- tail_counts[all_present,,drop=F]
    
    # Calculate Ordinary Least Squares residuals
    #residulator <- diag(nrow(design)) - design %*% MASS::ginv(design)
    #resids2 <- t(residulator %*% t(ap_log2_tails)) ^ 2
    #n <- length(resids2)

    # Choose optimium weight for technical variance component
    #
    # variance model for log2 tails is: overall_variance * ( technical_var_weight / read_count + 1)
    #
    # This is based on Maximum Likelihood for the residuals (assumed normally distributed).
    # The ML overall_variance given the weights can be directly found and substituted in, yielding this optimization: 
    #score_weights <- function(param) {
    #    weights <- ap_tail_counts / (param + ap_tail_counts)
    #    n*log(mean(resids2*weights)) - sum(log(weights))
    #}
    #param <- optimize(score_weights, c(0, max_weight))$minimum
    
    
    # Filter for features where the design can actually be fitted using samples with min_reads reads
    full_rank <- function(selection) {
        if (sum(selection) < ncol(design)) return(FALSE)
        sum(abs(svd(design[selection,,drop=F])$d)>=1e-10) == ncol(design)
    }
    good <- apply(tail_counts>=min_reads, 1, full_rank)

    good_log2_tails <- log2_tails[good,,drop=F]
    good_tail_counts <- tail_counts[good,,drop=F]
    good_genes <- genes[good,,drop=F]
    
    #new("EList", list(
    #    E=good_log2_tails, 
    #    weights=good_tail_counts/(param+good_tail_counts),
    #    genes=good_genes,
    #    technical_var_weight=param
    #))
    
    elist <- new("EList", list(
        E=good_log2_tails, 
        weights=good_tail_counts,
        genes=good_genes
    ))
    
    if (biovar)
        elist <- biovar_reweight(elist, design)
        
    elist
}


