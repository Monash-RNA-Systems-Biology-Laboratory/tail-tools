


# function to build results table

# Mann-Whitney r value
# a and b are equal length vectors giving read counts from two groups for an ordered list of peaks
mann_whitney_r <- function(a,b) {
    # Avoid integer overflow
    a <- as.numeric(a)
    b <- as.numeric(b)
    
    na <- sum(a)
    nb <- sum(b)
    
    if (na == 0 || nb == 0)
        return(data_frame(r=0.0, var=Inf))
    
    ca <- cumsum(a)
    cb <- cumsum(b)
    
    # Twice the Mann-Whitney-Wilcoxon U statistics, but centered at zero
    U2 <- sum(ca*b) - sum(cb*a)
    
    # Normal approximation to the null for U, accounting for ties
    # (but no continuity correction)
    ties <- a+b
    var_U <- na*nb/12 * (na+nb+1 - sum(ties*ties*ties-ties)/((na+nb)*(na+nb-1)))
    var_U2 <- 4*var_U
    
    # An "r" statistic between -1 and 1
    n <- na*nb
    r <- U2 / n
    var_r <- var_U2 / (n*n)
    
    data_frame(
        r = r,
        var = var_r)
}

# mat - int matrix - read counts for peak x sample
# condition - bool - division of samples into control and experimental
# group - factor - grouping of samples (eg by patient or batch or time) if such groups exist
combined_r <- function(mat, condition, group) {
    n_samples <- ncol(mat)
    
    group_result <-
        split(seq_len(n_samples), group) %>%
        map_df(function(cols) {
            submat <- mat[,cols,drop=F]
            subcondition <- condition[cols]
            a <- rowSums(submat[,!subcondition,drop=F])
            b <- rowSums(submat[,subcondition,drop=F])
            mann_whitney_r(a,b)
        }) %>%
        mutate(weight = 1/var)
    
    total_weight <- sum(group_result$weight)
    data_frame(
        r = sum(group_result$r*group_result$weight) / total_weight,
        var = 1/total_weight,
        mean_reads = sum(mat) / ncol(mat)
    )
}

#
# Permute condition, respecting group
#
grouped_permutations <- function(condition, group) {
    assert_that(is.logical(condition))
    assert_that(length(group) == length(condition))
    
    groups <- split(seq_along(group), group)
    
    permutations <- list( rep(F,length(group)) )
    for(group in groups) {
        n_true <- sum(condition[group])
        combs <- combn(group, n_true)
        new_permutations <- vector("list", length(permutations) * ncol(combs))
        i <- 1
        for(old_perm in permutations)
            for(j in seq_len(ncol(combs))) {
                new_perm <- old_perm
                new_perm[combs[,j]] <- T
                new_permutations[[i]] <- new_perm
                i <- i + 1
            }
        permutations <- new_permutations
    }
    
    permutations
}

#grouped_permutations(c(T,F,T,T,F),c(1,1,1,2,2))


#
# peak_info should contain: id, parent
#
edger_end_shift <- function(counts, peak_info, condition, group) {
    group <- as.factor(group)    
    if (length(unique(group)) <= 1)
        design <- model.matrix(~ condition)
    else
        design <- model.matrix(~ condition + group)
    
    # Terms are: intercept, condition, group terms (other than first level)
    # Second term is to be tested.    
    
    dge <- 
      DGEList(counts, genes=peak_info) %>%
      calcNormFactors() %>%
      estimateDisp(design) #Maybe use robust=T
    
    #dge$prior.df
    #dge$common.disp    
    #plotBCV(dge)
    
    fit <- glmQLFit(dge, design) #Maybe use robust=T
    
    # A normal edgeR analysis would then go on to:
    # qlf <- glmQLFTest(fit, coef=2)
    # topTags(qlf)
    # summary(decideTestsDGE(qlf, p.value=0.01))
    
    sp <- diffSpliceDGE(fit, coef=2, geneid="parent", exonid="id")
    top <- topSpliceDGE(sp, test="gene", n=Inf)
    
    list(
       top = top,
       dge = dge,
       design = design
    )
}


#
# peak_info should contain: id, parent
#
limma_end_shift <- function(counts, peak_info, condition, group) {
    group <- as.factor(group)    
    if (length(unique(group)) <= 1)
        design <- model.matrix(~ condition)
    else
        design <- model.matrix(~ condition + group)
    
    # Terms are: intercept, condition, group terms (other than first level)
    # Second term is to be tested.    
    
    dge <- 
      DGEList(counts, genes=peak_info) %>%
      calcNormFactors()

    fit <-
        voom(dge, design) %>% 
        lmFit(design) %>% 
        eBayes

    lsp <- diffSplice(fit, geneid="parent")
    top <- topSplice(lsp, coef=2, test="F", n=Inf)

    list(
        top = top,
        fit = fit,
        design = design
    )
}


#'
#' End shift statistics, generic function
#'
#' peak_info should contain columns: id, start, end, strand (+/-1), parent
#' strand should be the strand of the *gene*, if including antisense features
#'
#' min_reads is minimum mean-reads-per-sample for each gene
#'
#' @export
end_shift <- function(counts, peak_info, condition, group=NULL, 
                      gene_info_columns=c("gene","product","biotype"), 
                      ci=0.95, fdr=T, edger=T, limma=T, min_reads=10,
                      title="End-shift test") {
    counts <- as.matrix(counts)    
    if (is.null(group)) group <- rep(1,length(condition))
    
    assert_that(nrow(counts) > 0)
    assert_that(ncol(counts) == length(condition))
    assert_that(ncol(counts) == length(group))
    
    peak_info <- as_data_frame(peak_info) %>%
        mutate(parent = as.character(parent)) %>%
        mutate(parent = ifelse(parent == "", NA, parent))
    assert_that(nrow(peak_info) == nrow(counts))
    
    ci_sds <- qnorm((1+ci)/2)    
    
    cat("Split\n")
    splitter <-
        seq_len(nrow(counts)) %>% 
        split(peak_info$parent) %>%
        keep(function(item) {
            length(item) > 1 &&
            sum(counts[item,,drop=F])/ncol(counts) >= min_reads
        }) %>%
        lapply(function(item) {
            item[ order(peak_info$start[item] * peak_info$strand[item]) ]
        })
    
    keep_peaks <- peak_info$parent %in% names(splitter)
    
    if (edger) {
        cat("Run edgeR\n")
        edger_result <- edger_end_shift(
            counts[keep_peaks,,drop=F], 
            peak_info[keep_peaks,], 
            condition, group)
    } else {
        edger_result <- NULL
    }

    if (limma) {
        cat("Run limma\n")
        limma_result <- limma_end_shift(
            counts[keep_peaks,,drop=F], 
            peak_info[keep_peaks,], 
            condition, group)
    } else {
        limma_result <- NULL
    }
        
    cat("Score\n")
    scores <- map(splitter, ~ combined_r(counts[.,,drop=F],condition,group) )
    
    if (fdr) {
        cat("Permute\n")
        # Get the null distribution of "interest"
        null_interest <-
            grouped_permutations(condition, group) %>%
            map(function(condition_perm) {
                cat(paste0(ifelse(condition_perm,"+","-")),"\n")
                
                splitter %>%
                map_dbl(function(peaks) {
                    item <- combined_r(counts[peaks,,drop=F], condition_perm, group)
                    abs(item$r) - sqrt(item$var)*ci_sds
                })
            }) %>% 
            unlist %>%
            replace(.,is.na(.),-Inf) %>%
            sort(decreasing=T)
    }
        
    cat("Annotate\n")
    gene_info <- peak_info %>% 
        { .[,c("parent",gene_info_columns)] } %>%
        group_by(parent) %>%
        summarise_each(funs(paste(unique(as.character(.)),collapse="/") ))

    result <- 
        data_frame(parent=names(scores), scores=scores) %>% unnest(scores) %>%
        mutate(
            sd = sqrt(var),
            r_low = r - sd*ci_sds,
            r_high = r + sd*ci_sds,
            interest = ( abs(r) - sd*ci_sds ) %>% replace(.,is.na(.),-Inf)
        ) %>%
        arrange(desc(interest)) %>%
        mutate( rank=seq_len(n()) ) %>%
        select(rank, parent, r, r_low, r_high, interest, mean_reads)
        
    if (fdr) {
        cat("FDR\n")
        fdr <- numeric(nrow(result))
        j <- 0
        for(i in seq_along(fdr)) {
            interest <- result$interest[i]
            while(j < length(null_interest) && null_interest[j+1] >= interest) j <- j+1
            
            fdr[i] <- j/length(null_interest)*length(fdr) / i
        }
        fdr <- fdr%>%rev%>%cummin%>%rev
        
        result$fdr <- fdr
    }
    
    if (edger) {
        result <- left_join(
            result,
            edger_result$top %>% transmute(
                parent=parent, 
                edger_rank=seq_len(n()), 
                edger_fdr=FDR
            ),
            "parent")
    }
    
    if (limma) {
        result <- left_join(
            result,
            limma_result$top %>% transmute(
                parent=parent, 
                limma_rank=seq_len(n()), 
                limma_fdr=FDR
            ),
            "parent")
    }
    
    list(
        table = left_join(result, gene_info, "parent"),
        splitter = splitter,
        peak_info = peak_info,
        counts = counts,
        condition = condition,
        group = group,
        ci = ci,
        edger = edger_result,
        limma = limma_result,
        min_reads = min_reads,
        title = title
    )
}




#'
#' End shift statistics from tail-tools pipeline output.
#'
#' @param fdr Produce permutation based q values, can be very slow.
#'
#' Samples with NA in condition will be omitted.
#'
end_shift_pipeline <- function(path, condition, group=NULL, ci=0.95, fdr=T, edger=T, limma=T, antisense=T, non_utr=T, min_reads=10, title="End-shift test") {
    dat <- read.grouped.table(paste0(path,"/expression/peakwise/counts.csv"))
    
    counts <- as.matrix(dat$Count)
    peak_info <- dat$Annotation
    peak_info$id <- rownames(peak_info)
    peak_info$product <- str_match(peak_info$product, "^[^ ]+ (.*)$")[,2]
    
    anti <- peak_info$relation == "Antisense"
    peak_info$strand[anti] <- peak_info$strand[anti] * -1
    
    # Filter
    keep <- rep(T,nrow(peak_info))
    
    sample_keep <- !is.na(condition)
    
    if (!antisense)
        keep <- keep & peak_info$relation != "Antisense"
        
    if (!non_utr)
        keep <- keep & peak_info$relation == "3'UTR"
        
    counts <- counts[keep,sample_keep,drop=F]
    peak_info <- peak_info[keep,,drop=F]
    
    end_shift(counts, peak_info, condition[sample_keep], group[sample_keep], 
              ci=ci, fdr=fdr, edger=edger, limma=limma, min_reads=min_reads, title=title)
}



