
# Differential testing, topconfects based

test_variants = list(
    "test_vs" = c(
       "End shift"="test_end_shift", 
       "End shift, old quasi-likelihood method"="test_end_shift_ql",
       "Differential expression"="test_diff_exp")
)
# For legacy code
test_variants$test_end_shift <- test_variants$test_vs


#' Test for differential expression 
#'
#' @param min_reads There must be this many reads for an item to be included in the test, summing over all relevant samples.
#'
#' @export
test_diff_exp <- function(pipeline_dir, design, contrast=NULL, coef1=NULL, coef2=NULL, min_reads=10, samples=NULL, title=NULL, step=0.001, fdr=0.05) {
    tc <- read_tail_counts(paste0(pipeline_dir, "/expression/genewise/counts.csv"))

    if (is.null(samples))
        samples <- tc$samples$sample
        
    if (is.null(contrast)) {
        contrast <- rep(0, ncol(design))
        names(contrast) <- colnames(design)
        contrast[coef1] <- -1
        contrast[coef2] <- 1
    }

    tc <- tail_counts_subset_samples(tc, samples)

    mat <- tail_counts_get_matrix(tc, "count")
    keep <- rowSums(mat) >= min_reads

    voomed <- 
        edgeR::DGEList(
            mat[keep,,drop=FALSE],
            genes=select_(tc$features[keep,,drop=FALSE],~-feature)) %>%
        edgeR::calcNormFactors() %>%
        limma::voom(design)
    
    effect <- topconfects::effect_contrast(contrast)

    result <- topconfects::limma_nonlinear_confects(voomed, design, effect, step=step, fdr=fdr)
    result$pipeline_dir <- pipeline_dir
    result$title <- paste0(title, " - log2 fold change in expression")

    result
}


#test_diff_tail <- function(counts_filename, min_reads=10, ...) {
#
#}



#test_shiftexp <- function(pipeline_dir, subset="", min_reads=10, design, coef1, coef2) {
#    # Use effect_shift_log2 or ... hmm ... thingything
#}


#test_shift_tail <- function(pipeline_dir, subset="", min_reads=10, design, coef1, coef2) {
#    # Use effect_rss
#}



#' Test for 3' end shift
#'
#' @param pipeline_dir Pipeline directory.
#'
#' @param design Design matrix. Should have a column for expression level in condition 1, and a column for expression level in condition 2.
#'
#' @param coef1 Column number in the design matrix for expression level in condition 1.
#'
#' @param coef2 Column number in the design matrix for expression level in condition 2.
#'
#' @param min_reads Minimum total reads in order to use a peak.
#'
#' @param samples Sample names to be used. Rows of the design matrix should correspond to these samples. Defaults to all samples in the order they were given when the pipeline was run.
#'
#' @param fdr False Discovery Rate to maintain.
#'
#' @param step Accuracy of "confect" value. Should not need to be changed.
#'
#' @param antisense Include most antisense peaks.
#'
#' @param colliders In combination with antisense=TRUE, include antisense peaks, even those that are more likely sense peaks for another gene.
#'
#' @param non_utr Include non-3'UTR peaks.
#'
#' @param collapse_utr Use in combination with non_utr=TRUE. Collapse 3'UTR peaks into a "single peak", to look for shifts to/from non-3'UTR peaks.
#'
#' @param title A title for this test.
#'
#' @export
test_end_shift <- function(
        pipeline_dir, design, coef1, coef2, min_reads=10, samples=NULL,  
        fdr=0.05, step=0.001,
        antisense=F, colliders=F, non_utr=F, collapse_utr=F,
        title=NULL) {
    assert_that(length(coef1) == 1)
    assert_that(length(coef2) == 1)

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
    assert_that(nrow(design) == length(samples))
    
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

    display_members <- split(rownames(counts), parent)

    if (collapse_utr) {
        #mapping <- cumsum(!(peak_info$relation %in% c("3'UTR","Downstrand")) | parent != dplyr::lag(parent,default=""))
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
        dplyr::filter_(~n() >= 2) %>%
        dplyr::ungroup()

    # Perform test
    voomed <- 
        edgeR::DGEList(counts) %>%
        edgeR::calcNormFactors() %>%
        limma::voom(design)
    
    group_effect <- topconfects::group_effect_shift_unlog2(coef1, coef2)

    result <- topconfects::limma_group_confects(
        voomed, design, grouping, group_effect, step=step, fdr=fdr)

    result$table <- cbind(result$table, gene_dat$Annotation[result$table$name,,drop=F])

    result$pipeline_dir <- pipeline_dir
    result$title <- paste0(title, " - end shift")
    result$display_members <- display_members

    result
}


#' @export
test_end_shift_ql <- function(
        pipeline_dir, design, coef1, coef2, min_reads=10, samples=NULL,  
        fdr=0.05, step=0.01,
        antisense=F, colliders=F, non_utr=F, collapse_utr=F,
        title=NULL) {
    assert_that(length(coef1) == 1)
    assert_that(length(coef2) == 1)

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
    assert_that(nrow(design) == length(samples))
    
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

    display_members <- split(rownames(counts), parent)

    if (collapse_utr) {
        #mapping <- cumsum(!(peak_info$relation %in% c("3'UTR","Downstrand")) | parent != dplyr::lag(parent,default=""))
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

    # Perform test
    fit <- 
        edgeR::DGEList(counts) %>%
        edgeR::calcNormFactors() %>%
        edgeR::estimateDisp(design) %>%
        edgeR::glmQLFit(design)
    
    group_effect <- topconfectsql::group_effect_shift_log2(design, coef1, coef2)

    result <- topconfectsql::edger_group_confects(
        fit, parent, group_effect, step=step, fdr=fdr)

    result$table <- cbind(result$table, gene_dat$Annotation[result$table$name,,drop=F])

    result$pipeline_dir <- pipeline_dir
    result$title <- paste0(title, " - end shift (old quasi-likelihood method)")
    result$display_members <- display_members

    result
}

