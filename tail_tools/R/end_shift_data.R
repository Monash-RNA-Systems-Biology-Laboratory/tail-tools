
# Utility functions for processing peak information into easily usable data.

#' @export
get_grouped_peaks <- function(
        pipeline_dir, samples=NULL,  
        antisense=F, colliders=F, non_utr=F, collapse_utr=F,
        min_reads=10, min_group=2) {

    gene_counts_filename <- paste0(pipeline_dir, "/expression/genewise/counts.csv")
    gene_dat <- read_grouped_table(gene_counts_filename)

    counts_filename <- paste0(pipeline_dir, "/expression/peakwise/counts.csv")
    dat <- read_grouped_table(counts_filename)
    
    counts <- as.matrix(dat$Count)
    peak_info <- dplyr::as_data_frame(dat$Annotation)
    peak_info$id <- rownames(counts)

    if (is.null(samples))
        samples <- colnames(counts)

    assert_that(all(samples %in% colnames(counts)))
    
    for(name in colnames(peak_info))
        if (is.factor(peak_info[[name]]))
            peak_info[[name]] <- as.character(peak_info[[name]])

    # A peak may be sense to one gene and antisense to another.
    #   ( Tail Tools does not consider situations any more complex than this. )
    # Duplicate peaks antisense to a gene so they can be included in both genes
    #   ( Peaks not assigned a sense gene may already be labelled antisense to another gene,
    #     these don't need to be duplicated. )
    if (antisense && colliders && "antisense_parent" %in% colnames(peak_info)) {
        anti <- (peak_info$antisense_parent != "") & 
                (peak_info$relation != "Antisense")
        anti_counts <- counts[anti,,drop=F]
        anti_info <- peak_info[anti,,drop=F]
        
        # Incorporate antisense peaks
        anti_info <- anti_info %>% 
            dplyr::transmute_(
                id =~ paste0(id,"-collider"),
                start =~ start,
                end =~ end,
                strand =~ strand,
                relation =~ "Antisense",
                gene =~ antisense_gene,
                product =~ antisense_product,
                biotype =~ antisense_biotype,
                parent =~ antisense_parent
            )
        rownames(anti_counts) <- anti_info$id
        
        peak_info <- peak_info %>% 
            dplyr::select_(~id,~start,~end,~strand,~relation,~gene,~product,~biotype,~parent)
                
        counts <- rbind(counts, anti_counts)
        peak_info <- dplyr::bind_rows(peak_info, anti_info)
    }

    peak_info$product <- stringr::str_match(peak_info$product, "^[^ ]+ (.*)$")[,2]
    
    # Filter by relation to gene
    keep <- peak_info$parent != ""
    
    if (!antisense)
        keep <- keep & peak_info$relation != "Antisense"
        
    if (!non_utr)
        keep <- keep & peak_info$relation == "3'UTR"
        
    counts <- counts[keep,samples,drop=F]
    peak_info <- peak_info[keep,,drop=F]


    # Minimum read count filter
    keep2 <- rowSums(counts) >= min_reads
    counts <- counts[keep2,,drop=F]
    peak_info <- peak_info[keep2,,drop=F]


    # Order by position within gene
    position <- ifelse(peak_info$strand>0, peak_info$end, peak_info$start)
    strand <- peak_info$strand
    anti <- peak_info$relation == "Antisense"
    strand[anti] <- strand[anti] * -1
    
    ord <- order(peak_info$parent, strand*position)
    counts <- counts[ord,,drop=F]
    peak_info <- peak_info[ord,,drop=F]

    parent <- peak_info$parent

    display_members <- split(sub("-collider$","",rownames(counts)), parent)

    if (collapse_utr) {
        prev_parent <- ""
        in_utr <- FALSE
        mapping <- rep(NA, nrow(counts))
        j <- 0
        for(i in seq_len(nrow(counts))) {
            if (parent[i] != prev_parent) {
                prev_parent <- parent[i]
                in_utr <- FALSE
            }

            if (!in_utr)
                j <- j + 1
            mapping[i] <- j

            if (peak_info$relation[i] %in% c("3'UTR","Downstrand"))
                in_utr <- TRUE
        }


        n <- mapping[length(mapping)]
        new_parent <- rep("",n)
        new_counts <- matrix(0, nrow=n, ncol=ncol(counts))
        for(i in seq_along(mapping)) {
            j <- mapping[i]
            new_parent[j] <- parent[i]
            new_counts[j,] <- new_counts[j,] + counts[i,]
        }
        parent <- new_parent
        counts <- new_counts
    }
    
    grouping <- dplyr::data_frame(group=parent, name=rownames(counts)) %>%
        dplyr::group_by_(~group) %>%
        dplyr::filter_(~dplyr::n() >= min_group) %>%
        dplyr::ungroup()
    
    list(
        counts=counts, 
        grouping=grouping, 
        genes=gene_dat$Annotation[unique(grouping$group),,drop=F])
}


#
# Converts a matrix of read counts into a vector of shifts relative to the average
# Also provides appropriate weights (1/variance) for technical variation
#
# Samples with less than min_read reads are given shift=0, weight=0
#
weighted_shift <- function(mat, min_reads=1) {
    n <- nrow(mat)
    totals <- colSums(mat)
    good <- totals >= min_reads
    good_mat <- mat[,good,drop=F]
    good_totals <- totals[good]
    
    props <- t(t(good_mat)/good_totals)
    
    mid <- rowMeans(props)
    cummid <- cumsum(mid[-n])
    befores <- c(0,cummid)
    afters <- c(1-cummid, 0)
    pos_score <- befores-afters
    
    # Note: sum(pos_score * mid) == 0, sum(mid) == 1
    per_read_var <- sum(pos_score^2 * mid)

    shifts <- rep(0, length(good))
    shifts[good] <- colSums(props * pos_score)
    
    weights <- rep(0, length(good))
    weights[good] <- good_totals/per_read_var
    
    list(shifts=shifts, weights=weights, totals=totals, per_read_var=per_read_var)
}

#
# Given an EList with weights appropriate to technical variances,
# produces new EList with weights based on technical variance plus biological variance
# Biological variance assumed constant across genes
#
# Weights may be taken as 1/variance (although this doesn't allow for different variabilities between genes). 
#
# If design is not given, there is a hack to estimate the residual variance based on the median singular value.
#
#' @export
biovar_reweight <- function(elist, design=NULL, bio_weights=1) {
    E <- elist$E
    tv_weights <- elist$weights
    
    if (length(bio_weights) == 1)
        bio_weights <- rep(bio_weights, nrow(E))
    stopifnot(length(bio_weights) == nrow(E))
    
    all_present <- apply(tv_weights > 0, 1, all)
    if (!any(all_present)) 
        stop("No features with all weights present, while applying biological variation reweighting.")
    
    ap_E <- E[all_present,,drop=F]
    ap_tv_weights <- tv_weights[all_present,,drop=F]
    ap_bio_weights <- bio_weights[all_present]
    
    if (is.null(design)) {
        # Blind, robust residual variance estimation based on median singular value     
        ap_E_centered <- ap_E - rowMeans(ap_E)
        calc_var <- function(ap_weights) {
            d <- svd(sqrt(ap_weights) * ap_E_centered)$d
            median(d^2)/max(nrow(ap_E),ncol(ap_E))
        }
    } else {
        # Conventional residual variance estimation
        stopifnot(nrow(design) == ncol(E))
    
        # Calculate Ordinary Least Squares residuals
        residulator <- diag(nrow(design)) - design %*% MASS::ginv(design)
        resids2 <- t(residulator %*% t(ap_E)) ^ 2
        calc_var <- function(weights) {
            sum(resids2*ap_weights) / (nrow(resids2)*(ncol(resids2)-ncol(design)))
        }
    }

    # Choose optimium weight for technical variance component
    #
    # variance model is: overall_variance * ((1-p)*tv + p*bv)
    # hence weight is: 1/((1-p)/tw+p/bw)
    #
    # This is based on Maximum Likelihood for the residuals (assumed normally distributed).
    # The ML overall_variance given the weights can be directly found and substituted in, yielding this optimization: 
    score_weights <- function(param) {
        weights <- ap_tv_weights/(1-param + param/ap_bio_weights*ap_tv_weights)
        #n*log(mean(resids2*weights)) - sum(log(weights))
        length(weights)*log(calc_var(weights)) - sum(log(weights))
    }
    param <- optimize(score_weights, c(0, 1))$minimum

    elist$weights <- tv_weights / (1-param+param/bio_weights*tv_weights)
    
    # Allow weights to be used directly as precisions
    ap_weights <- elist$weights[all_present,,drop=F]
    residual_var <- calc_var(ap_weights)
    elist$techvar <- residual_var*(1-param)
    elist$biovar <- residual_var*param
    elist$weights <- elist$weights / residual_var
    #fit <- lmFit(elist, design) %>% eBayes()
    #elist$weights <- elist$weights / fit$s2.prior
    
    elist
}


#
# Shift values relative to overall, with associated weight
# Missing values are given value 0 and weight 0
#
# grouping data frame should have columns "group" and "name",
# "name" corresponding to rownames of counts
#
# Biological variance estimation:
# If design is not given, there is a hack to esimtate residual variance based on the median singular value.
# If biovar is FALSE, this step is skipped, and weights represent 1/(technical variance).
#
# Rows containing less than two samples with enough reads are discarded.
#
#' @export
weighted_shifts <- function(counts, grouping, design=NULL, biovar=TRUE, min_reads=1, genes=NULL) {
    groups <- split(grouping$name, grouping$group)

    # Only use groups of 2 or more peaks
    good <- map_int(groups, length) >= 2
    groups <- groups[good]
    
    results <- purrr::map(groups, function(members) 
        weighted_shift(counts[members,,drop=FALSE], min_reads=min_reads))
    
    shifts <- do.call(rbind, purrr::map(results, "shifts"))
    weights <- do.call(rbind, purrr::map(results, "weights"))
    totals <- do.call(rbind, purrr::map(results, "totals"))
    
    bio_weights <- 1/purrr::map_dbl(results, "per_read_var")
    
    rownames(shifts) <- names(groups)
    rownames(weights) <- names(groups)
    rownames(totals) <- names(groups)
    colnames(shifts) <- colnames(counts)
    colnames(weights) <- colnames(counts)
    colnames(totals) <- colnames(counts)
    
    good <- rowSums(weights > 0) >= 2
    shifts <- shifts[good,,drop=F]
    weights <- weights[good,,drop=F]
    totals <- totals[good,,drop=F]
    bio_weights <- bio_weights[good]
    
    elist <- new("EList", list(
        E=shifts, 
        weights=weights,
        total_counts=totals,
        genes=genes[rownames(shifts),,drop=F],
        bio_weights=bio_weights
    ))
    
    if (biovar)
        elist <- biovar_reweight(elist, design, bio_weights=bio_weights)
    
    elist
}