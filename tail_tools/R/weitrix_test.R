
as_two_coef_test <- function(func) {
    func

    function(pipeline_dir, design, coef1=NULL, coef2=NULL, ...) {
        contrast <- rep(0, ncol(design))
        names(contrast) <- colnames(design)
        contrast[coef1] <- -1
        contrast[coef2] <- 1

        func(pipeline_dir, design, contrast, ...)
    }
}

as_contrast <- function(contrast, design) {
   if (is.character(contrast)) {
       stopifnot( contrast %in% colnames(design) )
       contrast <- as.numeric( colnames(design) == contrast )
   }

   contrast <- as.matrix(contrast)
   stopifnot(nrow(contrast) == ncol(design))
   stopifnot(ncol(contrast) == 1)
   contrast
}


#' Test for differential expression 
#'
#' Doesn't actually use weitrix, but method is similar.
#'
#' @param min_reads There must be this many reads for an item to be included in the test, summing over all relevant samples.
#'
#' @export
test_diff_exp <- function(pipeline_dir, design, contrast=NULL, coef1=NULL, coef2=NULL, min_reads=10, samples=NULL, title=NULL, step=0.001, fdr=0.05, what="genewise") {

    contrast <- as_contrast(contrast, design)
    design_coef <- limma::contrastAsCoef(design, contrast)

    tc <- read_tail_counts(paste0(pipeline_dir, "/expression/", what, "/counts.csv"))

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
            genes=select(tc$features[keep,,drop=FALSE], -feature)) %>%
        edgeR::calcNormFactors() %>%
        limma::voom(design)
    
    # limma's use of weights is not exact for contrasts (old solution)
    #effect <- topconfectswald::effect_contrast(contrast)
    #result <- topconfectswald::limma_nonlinear_confects(voomed, design, effect, step=step, fdr=fdr, full=TRUE)

    # limma's use of weights is not exact for contrasts
    # so alter design matrix to make contrast a coefficient
    fit <- limma::lmFit(voomed, design_coef$design)
    result <- topconfects::limma_confects(fit, design_coef$coef, full=TRUE, step=step, fdr=fdr)

    result$pipeline_dir <- pipeline_dir
    result$effect_desc <- "log2 fold change in expression"
    result$title <- paste0(title, " - log2 fold change in expression - ", what, " - at least ", min_reads, " reads")

    text <- capture.output({
        cat("Samples\n")
        print(colnames(mat))
        cat("\nContrast\n")
        print(contrast)
        cat("\nDesign matrix\n")
        print(design)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=weitrix::weitrix_calplot(voomed, design, cat=col))

    result$diagnostics[["Calibration vs expression"]] <- list(
        plot=weitrix::weitrix_calplot(voomed, design, covar=mu))

    result
}

test_diff_exp_20 <- function(...) test_diff_exp(..., min_reads=20)
test_diff_exp_50 <- function(...) test_diff_exp(..., min_reads=50)



#' @export
test_end_shift_weitrix <- function(
        pipeline_dir, design, contrast, samples=NULL, fdr=0.05,
        antisense=F, colliders=F, non_utr=F, collapse_utr=F, min_reads=10,
        title="End shift test (weitrix based)") {

    contrast <- as_contrast(contrast, design)
    design_coef <- limma::contrastAsCoef(design, contrast)

    wei <- pipeline_weitrix_shift(
        pipeline_dir=pipeline_dir, samples=samples,  
        antisense=antisense, colliders=colliders, non_utr=non_utr, collapse_utr=collapse_utr,
        min_reads=min_reads)
    
    cal <- weitrix::weitrix_calibrate_all(wei, design_coef$design)

    fit <- limma::lmFit(weitrix::weitrix_elist(cal), design_coef$design)

    result <- topconfects::limma_confects(fit, design_coef$coef, full=TRUE, fdr=fdr)

    # AveExpr should actually be about expression (log2 RPM)
    total_reads <- SummarizedExperiment::rowData(cal)$total_reads
    ave_expr <- log2( total_reads*1e6 / sum(total_reads) )
    result$table$AveExpr <- ave_expr[result$table$index]

    result$limits <- c(-1,1)
    result$pipeline_dir <- pipeline_dir
    result$title <- paste0(title, " - end shift (weitrix based) - at least ", min_reads, " reads")
    result$effect_desc <- "shift"
    result$display_members <- S4Vectors::metadata(wei)$display_members

    text <- capture.output({
        cat("Samples\n")
        print(colnames(wei))
        cat("\nContrast\n")
        print(contrast)
        cat("\nDesign matrix\n")
        print(design)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, cat=col))

    result$diagnostics[["Calibration vs shift"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, covar=mu))

    result$diagnostics[["Calibration vs per_read_var"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, covar=per_read_var))

    result$diagnostics[["Calibration vs read count"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, covar=log2(weitrix::weitrix_weights(wei))) + labs(x="log2 reads"))

    result
}

test_end_shift_weitrix_twocoef <- as_two_coef_test(test_end_shift_weitrix)

test_end_shift_weitrix_nonutr <- function(...) test_end_shift_weitrix(..., non_utr=TRUE)
test_end_shift_weitrix_nonutr_twocoef <- as_two_coef_test(test_end_shift_weitrix_nonutr)



#' @export
test_diff_tail_weitrix <- function(
        pipeline_dir, design, contrast, samples=NULL, fdr=0.05,
        what="genewise", min_reads=50,
        title="Differential tail length (weitrix based)") {

    contrast <- as_contrast(contrast, design)
    design_coef <- limma::contrastAsCoef(design, contrast)

    wei <- pipeline_weitrix_tail(
        pipeline_dir=pipeline_dir, samples=samples, what=what,
        min_reads=min_reads, design=design_coef$design)

    cal <- weitrix::weitrix_calibrate_all(wei, design_coef$design)

    fit <- limma::lmFit(weitrix::weitrix_elist(cal), design_coef$design)
    
    result <- topconfects::limma_confects(fit, design_coef$coef, full=TRUE, fdr=fdr)

    total_reads <- rowSums(weitrix::weitrix_weights(wei))
    ave_expr <- log2( total_reads*1e6 / sum(total_reads) )
    result$table$AveExpr <- ave_expr[result$table$index]

    ave_tail <- weitrix::weitrix_components(cal, design=~1, verbose=F)$row[,1]
    result$table$AveTail <- ave_tail[result$table$index]

    result$pipeline_dir <- pipeline_dir
    result$title <- paste0(title, " - tail length (weitrix based) - at least ", min_reads, " reads in sufficient samples")
    result$effect_desc <- "change in tail length"
    result$magnitude_column <- "AveTail"
    result$magnitude_desc <- "Average tail length"

    text <- capture.output({
        cat("Samples\n")
        print(colnames(wei))
        cat("\nContrast\n")
        print(contrast)
        cat("\nDesign matrix\n")
        print(design)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, cat=col))

    result$diagnostics[["Calibration vs tail length"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, covar=mu))

    result$diagnostics[["Calibration vs read count"]] <- list(
        plot=weitrix::weitrix_calplot(cal, design, covar=log2(weitrix::weitrix_weights(wei))) + labs(x="log2 reads with tail"))

    result  
}

test_diff_tail_weitrix_twocoef <- as_two_coef_test(test_diff_tail_weitrix)

test_diff_tail_weitrix_100 <- function(...) test_diff_tail_weitrix(..., min_reads=100)
test_diff_tail_weitrix_100_twocoef <- as_two_coef_test(test_diff_tail_weitrix_100)

test_diff_tail_weitrix_200 <- function(...) test_diff_tail_weitrix(..., min_reads=200)
test_diff_tail_weitrix_200_twocoef <- as_two_coef_test(test_diff_tail_weitrix_200)


test_diff_tail_weitrix_primary_peak_50 <- function(...) test_diff_tail_weitrix(..., min_reads=50, what="primarypeakwise")
test_diff_tail_weitrix_primary_peak_50_twocoef <- as_two_coef_test(test_diff_tail_weitrix_primary_peak_50)

test_diff_tail_weitrix_primary_peak_100 <- function(...) test_diff_tail_weitrix(..., min_reads=100, what="primarypeakwise")
test_diff_tail_weitrix_primary_peak_100_twocoef <- as_two_coef_test(test_diff_tail_weitrix_primary_peak_100)

test_diff_tail_weitrix_primary_peak_200 <- function(...) test_diff_tail_weitrix(..., min_reads=200, what="primarypeakwise")
test_diff_tail_weitrix_primary_peak_200_twocoef <- as_two_coef_test(test_diff_tail_weitrix_primary_peak_200)
