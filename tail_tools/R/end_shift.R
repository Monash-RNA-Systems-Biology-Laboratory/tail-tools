


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
    
    list(
        r = r,
        var = var_r
    )
}

# mat - int matrix - read counts for peak x sample
# condition - bool - division of samples into control and experimental
# group - factor - grouping of samples (eg by patient or batch or time) if such groups exist
combined_r <- function(mat, condition, sample_splitter) {
    total_r <- 0.0
    total_weight <- 0.0
    for(cols in sample_splitter) {
        submat <- mat[,cols,drop=F]
        subcondition <- condition[cols]
        a <- rowSums(submat[,!subcondition,drop=F])
        b <- rowSums(submat[,subcondition,drop=F])
        result <- mann_whitney_r(a,b)
        weight <- 1/result$var
        total_weight <- total_weight + weight
        total_r <- total_r + result$r * weight
    }
    
    list(
        r = total_r / total_weight,
        var = 1 / total_weight
    )
}

#
# Permute condition, respecting group
#

random_permutations <- function(condition, group, n) {
    groups <- split(seq_along(group), group)
    
    #Definitely include true permutation
    permutations <- list(condition)
    
    while(length(permutations) < n) {
        permutation <- rep(F,length(group))
        for(group in groups)
            permutation[group] <- condition[group][sample(length(group))]
        
        # Inefficient!
        seen <- map_lgl(permutations, function(item) identical(permutation,item)) %>% any()
        if (!seen) permutations[[length(permutations)+1]] <- permutation
    }
    
    permutations
}


grouped_permutations <- function(condition, group, max_n=Inf) {
    assert_that(is.logical(condition))
    assert_that(length(group) == length(condition))
    
    groups <- split(seq_along(group), group)
    
    # Switch to random if too many.
    n <- 1
    for(item in groups) {
        n <- n * choose(length(item), sum(condition[item]))
        if (n > max_n)
            return(random_permutations(condition,group,max_n))
    }
    
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
# gene information should contain: id, parent
#
edger_end_shift <- function(dge, condition, group) {
    group <- as.factor(group)    
    if (length(unique(group)) <= 1)
        design <- model.matrix(~ condition)
    else
        design <- model.matrix(~ condition + group)
    
    # Terms are: intercept, condition, group terms (other than first level)
    # Second term is to be tested.    
    
    dge <- estimateDisp(dge, design) #Maybe use robust=T
    
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
# gene information should contain: id, parent
#
limma_end_shift <- function(dge, condition, group) {
    group <- as.factor(group)    
    if (length(unique(group)) <= 1)
        design <- model.matrix(~ condition)
    else
        design <- model.matrix(~ condition + group)
    
    # Terms are: intercept, condition, group terms (other than first level)
    # Second term is to be tested.    
    
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
#' @param peak_info should be a data frame containing columns: id, position, strand (+/-1), parent. position is position in chromosome of transcription stop site. strand should be the strand of the *gene*, if including antisense features.
#'
#' @param condition is a logical vector splitting samples into control and experimental groups.
#'
#' @param group if given is a factor splitting samples into batches, if there is a batch effect. 
#'
#' @param gene_info_columns Columns of peak_info to retain in per-gene output.
#'
#' @param ci Confidence interval. A value closer to 1 will more heavily demote genes with low read counts.
#'
#' @param fdr Perform permutation based FDR?
#'
#' @param fdr_max_permute If there are more than this many possible permutations, just use this many randomly sampled distinct permutations (but certainly including the true permutation).
#'
#' @param edger Perform edgeR differential exon usage test. May need to disable if too few replicates.
#'
#' @param limma Perform limma differential exon usage test. May need to disable if no replicates.
#'
#' @param min_reads Discard genes with lower than this average number of reads per sample.
#'
#' 
#'
#' @export
end_shift <- function(counts, peak_info, condition, group=NULL,
                      gene_info_columns=c("gene","product","biotype"), 
                      ci=0.95, fdr=T, edger=T, limma=T, min_reads=10,
                      title="End-shift test",
                      fdr_max_permute=1000) {
    counts <- as.matrix(counts)    
    if (is.null(group)) group <- rep(1,length(condition))
    
    assert_that(nrow(counts) > 0)
    assert_that(ncol(counts) == length(condition))
    assert_that(ncol(counts) == length(group))
    
    peak_info <- dplyr::as_data_frame(peak_info) %>%
        dplyr::mutate(parent = as.character(parent)) %>%
        dplyr::mutate(parent = ifelse(parent == "", NA, parent))
    assert_that(nrow(peak_info) == nrow(counts))
    
    ci_sds <- qnorm((1+ci)/2)    


    cat("Lib size\n")
    
    dge <- 
      DGEList(counts, genes=peak_info) %>%
      calcNormFactors()
    
    lib_size <- dge$samples$lib.size * dge$samples$norm.factors
    
    
    cat("Split\n")

    orderer <- peak_info$position * peak_info$strand 

    splitter <-
        seq_len(nrow(counts)) %>% 
        split(peak_info$parent) %>%
        purrr::keep(function(item) {
            length(item) > 1 &&
            sum(counts[item,,drop=F])/ncol(counts) >= min_reads
        }) %>%
        lapply(function(item) {
            item[ order(orderer[item]) ]
        })
    
    keep_peaks <- peak_info$parent %in% names(splitter)

    
    if (edger) {
        cat("Run edgeR\n")
        edger_result <- edger_end_shift(
            dge[keep_peaks,],
            condition, group)
    } else {
        edger_result <- NULL
    }

    if (limma) {
        cat("Run limma\n")
        limma_result <- limma_end_shift(
            dge[keep_peaks,], 
            condition, group)
    } else {
        limma_result <- NULL
    }
        
    cat("Score\n")
    sample_splitter <- split(seq_along(group), group)
    
    scores <- splitter %>% 
        lapply(function(item) {
            combined_r(counts[item,,drop=F],condition,sample_splitter) 
        })
    
    if (fdr) {
        cat("Permute\n")

        # Get the null distribution of "interest"
        perms <- grouped_permutations(condition, group, fdr_max_permute)
        cat(length(perms),"permutations\n")         
        null_interest <-
            lapply(perms, function(condition_perm) {
                cat(paste0(ifelse(condition_perm,"+","-")),"\n")
                
                map_dbl(splitter, function(peaks) {
                    item <- combined_r(counts[peaks,,drop=F], condition_perm, sample_splitter)
                    abs(item$r) - sqrt(item$var)*ci_sds
                })
            }) %>% 
            unlist %>%
            replace(.,is.na(.),-Inf) %>%
            sort(decreasing=T)
    }
        
    cat("Annotate\n")
    gene_info <- peak_info[,c("parent",gene_info_columns)] %>%
        dplyr::group_by(parent) %>%
        dplyr::summarise_each(funs(paste(unique(as.character(.)),collapse="/") ))

    result <- 
        dplyr::data_frame(parent=names(scores), score=scores) %>% 
        dplyr::rowwise() %>%
        dplyr::mutate(
            r = score$"r",
            var = score$"var",
            score = NULL,
            mean_reads = sum(counts[splitter[[parent]],,drop=F]) / ncol(counts)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
            sd = sqrt(var),
            r_low = r - sd*ci_sds,
            r_high = r + sd*ci_sds,
            interest = ( abs(r) - sd*ci_sds ) %>% replace(.,is.na(.),-Inf)
        ) %>%
        dplyr::arrange(desc(interest)) %>%
        dplyr::mutate( rank=seq_len(n()) ) %>%
        dplyr::select(rank, parent, r, r_low, r_high, interest, mean_reads)
        
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
        lib_size = lib_size,
        condition = condition,
        group = group,
        ci = ci,
        edger = edger_result,
        limma = limma_result,
        min_reads = min_reads,
        title = title,
        fdr_max_permute = fdr_max_permute,
        fdr_actual_permute = length(perms)
    )
}




#'
#' End shift statistics from tail-tools pipeline output.
#'
#' @param fdr Produce permutation based q values, can be very slow.
#'
#' Samples with NA in condition will be omitted.
#' Samples can also be selected using select=, which should contain a logical vector of samples to retain.
#'
end_shift_pipeline <- function(path, condition, group=NULL, select=NULL, ci=0.95, fdr=T, edger=T, limma=T, antisense=T, colliders=T, non_utr=T, min_reads=10, title="End-shift test", fdr_max_permute=1000) {
    dat <- read.grouped.table(paste0(path,"/expression/peakwise/counts.csv"))
    
    counts <- as.matrix(dat$Count)
    peak_info <- dplyr::as_data_frame(dat$Annotation)
    peak_info$id <- rownames(counts)
    
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
                (peak_info$relation != "Antisense") &
                (peak_info$antisense_parent != "")
        anti_counts <- counts[anti,,drop=F]
        anti_info <- peak_info[anti,,drop=F]
        
        # Incorporate antisense peaks
        anti_info <- anti_info %>% 
            dplyr::transmute(
                id = paste0(id,"-collider"),
                start = start,
                end = end,
                strand = strand,
                relation = "Antisense",
                gene = antisense_gene,
                product = antisense_product,
                biotype = antisense_biotype,
                parent = antisense_parent
            )
        rownames(anti_counts) <- anti_info$id
        
        peak_info <- peak_info %>% 
            dplyr::select(id,start,end,strand,relation,gene,product,biotype,parent)
                
        counts <- rbind(counts, anti_counts)
        peak_info <- dplyr::bind_rows(peak_info, anti_info)
    }


    peak_info$product <- str_match(peak_info$product, "^[^ ]+ (.*)$")[,2]
    
    peak_info$position <- ifelse(peak_info$strand>0, peak_info$end, peak_info$start)
    anti <- peak_info$relation == "Antisense"
    peak_info$strand[anti] <- peak_info$strand[anti] * -1
    
    # Filter
    keep <- rep(T,nrow(peak_info))
    
    sample_keep <- !is.na(condition) 
    if (!is.null(select))
        sample_keep <- sample_keep & select
    
    if (!antisense)
        keep <- keep & peak_info$relation != "Antisense"
        
    if (!non_utr)
        keep <- keep & peak_info$relation == "3'UTR"
        
    counts <- counts[keep,sample_keep,drop=F]
    peak_info <- peak_info[keep,,drop=F]
    
    end_shift(counts, peak_info, condition[sample_keep], group[sample_keep], 
              ci=ci, fdr=fdr, edger=edger, limma=limma, min_reads=min_reads, title=title, fdr_max_permute=fdr_max_permute)
}



