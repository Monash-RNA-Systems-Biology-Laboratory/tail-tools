
library(knitr)
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(assertthat)
library(ggplot2)
library(nesoni)
library(edgeR)
library(limma)
library(seriation)
library(gsubfn)
library(xtable)

library(shiny)
library(DT) # github version devtools::install_github("rstudio/DT")
library(dplyr)

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
        var = 1/total_weight
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



#'
#' End shift statistics, generic function
#'
#' peak_info should contain columns: start, end, strand (+/-1), parent
#' strand should be the strand of the *gene*, if including antisense features
#'
#' @export
end_shift <- function(counts, peak_info, condition, group=NULL, 
                      gene_info_columns=c("gene","product","biotype"), 
                      ci=0.95, fdr=T) {
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
    
    cat("split\n")
    splitter <- 
        split(seq_len(nrow(counts)), peak_info$parent) %>%
        keep(~ length(.) > 1) %>%
        map(~ .[ order(peak_info$start[.] * peak_info$strand[.]) ])
    
    cat("score\n")
    scores <- map(splitter, ~ combined_r(counts[.,,drop=F],condition,group) )
    
    if (fdr) {
        cat("permute\n")
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
        
    cat("annotate\n")
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
        select(rank, parent, r, r_low, r_high, interest)
        
    if (fdr) {
        cat("fdr\n")
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
    
    list(
        table = left_join(result, gene_info, "parent"),
        splitter = splitter,
        peak_info = peak_info,
        counts = counts,
        condition = condition,
        group = group,
        ci = ci
    )
}




#'
#' End shift statistics from tail-tools pipeline output.
#'
#' @param fdr Produce permutation based q values, can be very slow.
#'
end_shift_pipeline <- function(path, condition, group=NULL, ci=0.95, fdr=T, antisense=T, non_utr=T) {
    dat <- read.grouped.table(paste0(path,"/expression/peakwise/counts.csv"))
    
    counts <- as.matrix(dat$Count)
    peak_info <- dat$Annotation
    peak_info$peak <- rownames(peak_info)
    peak_info$product <- str_match(peak_info$product, "^[^ ]+ (.*)$")[,2]
    
    anti <- peak_info$relation == "Antisense"
    peak_info$strand[anti] <- peak_info$strand[anti] * -1
    
    # Filter
    keep <- rep(T,nrow(peak_info))
    
    if (!antisense)
        keep <- keep & peak_info$relation != "Antisense"
        
    if (!non_utr)
        keep <- keep & peak_info$relation == "3'UTR"
        
    counts <- counts[keep,,drop=F]
    peak_info <- peak_info[keep,,drop=F]
    
    end_shift(counts, peak_info, condition, group, ci=ci, fdr=fdr)
}



