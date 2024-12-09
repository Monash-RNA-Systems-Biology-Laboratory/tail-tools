
as_two_coef_test <- function(func) {
    func

    function(pipeline_dir, design, coef1=NULL, coef2=NULL, ...) {
        contrast <- rep(0, ncol(design))
        names(contrast) <- colnames(design)
        contrast[coef1] <- -1
        contrast[coef2] <- 1

        func(pipeline_dir=pipeline_dir, design=design, contrast=contrast, ...)
    }
}

as_contrast <- function(contrast, design) {
   if (is.character(contrast)) {
       stopifnot( contrast %in% colnames(design) )
       contrast <- as.numeric( colnames(design) == contrast )
   }

   contrast <- as.matrix(contrast)
   stopifnot(nrow(contrast) == ncol(design))
   stopifnot(ncol(contrast) >= 1)
   contrast
}


# Detrend each sample vs overall mean, calibrate
detrend_samples_and_calibrate <- function(wei, design) {
    cal <- weitrix::weitrix_calibrate_all(wei, design)

    long <- weitrix::matrix_long(weitrix::weitrix_x(cal))
    long$weight <- as.vector(weitrix::weitrix_weights(cal))
    long <- long %>%
        dplyr::group_by(name) %>%
        dplyr::mutate(mu=weighted.mean(value,weight)) %>%
        dplyr::ungroup()
    fit <- lm( value ~ weitrix::well_knotted_spline(mu,3) * col, data=long, weight=weight)

    long$pred <- predict(fit, newdata=long)

    mat <- matrix(long$value - long$pred + long$mu, nrow=nrow(wei))
    rownames(mat) <- rownames(wei)
    colnames(mat) <- colnames(wei)

    fixed <- wei
    weitrix::weitrix_x(fixed) <- mat

    weitrix::weitrix_calibrate_all(fixed, design)
}


perform_weitrix_test <- function(weitrix, design, coef, contrast, fdr, step, dispersion_est) {
    # Backwards compatability, should not be needed
    if (!is.null(contrast)) {
        contrast <- as_contrast(contrast, design)
    }

    if (is.null(coef) && is.null(contrast)) {
        weitrix::weitrix_sd_confects(weitrix, design, fdr=fdr, step=step)
    } else {
        weitrix::weitrix_confects(weitrix, design, coef=coef, contrasts=contrast, fdr=fdr, step=step, dispersion_est=dispersion_est)
    }
}

# Save disk space and loading time for diagnostic plots
compact_plot <- function(p) patchwork::wrap_elements(cowplot::as_grob(p))


#' Legacy test for differential expression 
#' @export
test_diff_exp <- function(pipeline_dir, design, contrast=NULL, coef1=NULL, coef2=NULL, min_reads=10, samples=NULL, title=NULL, step=0.001, fdr=0.05, what="genewise") {
    tc <- read_tail_counts(paste0(pipeline_dir, "/expression/", what, "/counts.csv"))

    if (is.null(samples))
        samples <- tc$samples$sample
        
    if (is.null(contrast)) {
        contrast <- rep(0, ncol(design))
        names(contrast) <- colnames(design)
        contrast[coef1] <- -1
        contrast[coef2] <- 1
    }

    test_diff_exp_weitrix(pipeline_dir=pipeline_dir, design=design, contrast=contrast, min_reads=min_reads, samples=samples, title=title, step=step, fdr=fdr, what=what)
}

test_diff_exp_20 <- function(...) test_diff_exp(..., min_reads=20)
test_diff_exp_50 <- function(...) test_diff_exp(..., min_reads=50)


#' @export
test_diff_exp_weitrix <- function(
        pipeline_dir, design, coef=NULL, contrast=NULL, min_reads=10, 
        samples=NULL, title=NULL, step=NULL, fdr=0.05, 
        calibration_design=design, dispersion_est="ebayes_limma",
        what="genewise") {
    tc <- read_tail_counts(paste0(pipeline_dir, "/expression/", what, "/counts.csv"))
    tc <- tail_counts_subset_samples(tc, samples)

    mat <- tail_counts_get_matrix(tc, "count")
    keep <- rowSums(mat) >= min_reads

    voomed <- 
        edgeR::DGEList(
            mat[keep,,drop=FALSE],
            genes=select(tc$features[keep,,drop=FALSE], -feature)) %>%
        edgeR::calcNormFactors() %>%
        limma::voom(calibration_design)
    
    result <- perform_weitrix_test(weitrix::as_weitrix(voomed), design=design, coef=coef, contrast=contrast, fdr=fdr, step=step, dispersion_est=dispersion_est)

    result$pipeline_dir <- pipeline_dir
    result$effect_desc <- paste0(result$effect_desc, " of log2 RPM")
    result$title <- paste0(title, " - ", result$effect_desc, " - ", what, " - at least ", min_reads, " reads")

    text <- capture.output({
        cat("\nDesign matrix\n")
        design <- result$design
        rownames(design) <- colnames(voomed)
        print(design)
        cat("\nContrasts\n")
        contrasts <- result$contrasts
        if (!is.null(contrasts))
            rownames(contrasts) <- colnames(design)
        print(contrasts)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(voomed, design, cat=col)))

    result$diagnostics[["Calibration vs expression"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(voomed, design, covar=mu)))

    result
}

test_diff_exp_weitrix_20 <- function(...) test_diff_exp_weitrix(..., min_reads=20)
test_diff_exp_weitrix_50 <- function(...) test_diff_exp_weitrix(..., min_reads=50)

test_diff_exp_weitrix_peak_10 <- function(...) test_diff_exp_weitrix(..., min_reads=10, what="peakwise")
test_diff_exp_weitrix_peak_20 <- function(...) test_diff_exp_weitrix(..., min_reads=20, what="peakwise")
test_diff_exp_weitrix_peak_50 <- function(...) test_diff_exp_weitrix(..., min_reads=50, what="peakwise")

#' @export
test_end_shift_weitrix <- function(
        pipeline_dir, design, coef=NULL, contrast=NULL, samples=NULL, fdr=0.05, step=NULL,
        antisense=F, colliders=F, non_utr=F, collapse_utr=F, min_reads=10,
        calibration_design=design, dispersion_est="ebayes_limma",
        title="End shift test (weitrix based)") {

    # Take smaller steps than default for sd effect size
    if (is.null(step) && (length(coef) >= 2 || (is.matrix(contrast) && ncol(contrast) >= 2)))
        step <- 0.01

    wei <- pipeline_weitrix_shift(
        pipeline_dir=pipeline_dir, samples=samples,  
        antisense=antisense, colliders=colliders, non_utr=non_utr, collapse_utr=collapse_utr,
        min_reads=min_reads, design=calibration_design)
    
    cal <- weitrix::weitrix_calibrate_all(wei, calibration_design)
    
    result <- perform_weitrix_test(cal, design=design, coef=coef, contrast=contrast, fdr=fdr, step=step, dispersion_est=dispersion_est)

    # AveExpr should actually be about expression (log2 RPM)
    total_reads <- SummarizedExperiment::rowData(cal)$total_reads
    ave_expr <- log2( total_reads*1e6 / sum(total_reads) )
    result$table$AveExpr <- ave_expr[result$table$index]
    result$table$row_mean <- NULL
    
    if (is.null(result$limits))
        result$limits <- c(-1,1)
    result$pipeline_dir <- pipeline_dir
    result$effect_desc <- paste0(result$effect_desc, " of APA")
    result$title <- paste0(
        title, 
        " - ", result$effect_desc, " - ", 
        if (antisense || colliders || collapse_utr) "custom peaks" 
        else if (non_utr) "all sense peaks" 
        else "3' UTR peaks",
        " - at least ", min_reads, " reads")
    result$display_members <- S4Vectors::metadata(wei)$display_members

    text <- capture.output({
        cat("\nDesign matrix\n")
        design <- result$design
        rownames(design) <- colnames(wei)
        print(design)
        cat("\nContrasts\n")
        contrasts <- result$contrasts
        if (!is.null(contrasts))
            rownames(contrasts) <- colnames(design)
        print(contrasts)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, cat=col)))

    result$diagnostics[["Calibration vs shift"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=mu)))

    result$diagnostics[["Calibration vs per_read_var"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=per_read_var)))

    result$diagnostics[["Calibration vs read count"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=log2(weitrix::weitrix_weights(wei))) + labs(x="log2 reads")))

    result
}

test_end_shift_weitrix_twocoef <- as_two_coef_test(test_end_shift_weitrix)

test_end_shift_weitrix_nonutr <- function(...) test_end_shift_weitrix(..., non_utr=TRUE)
test_end_shift_weitrix_nonutr_twocoef <- as_two_coef_test(test_end_shift_weitrix_nonutr)



#' @export
test_diff_tail_weitrix <- function(
        pipeline_dir, design, coef=NULL, contrast=NULL, samples=NULL, fdr=0.05, step=NULL, 
        what="genewise", detrend_samples=FALSE,
        min_reads=50, calibration_design=design, dispersion_est="ebayes_limma",
        title="Differential tail length (weitrix based)") {

    # Take bigger steps than default for sd effect size
    if (is.null(step) && (length(coef) >= 2 || (is.matrix(contrast) && ncol(contrast) >= 2)))
        step <- 1.0

    wei <- pipeline_weitrix_tail(
        pipeline_dir=pipeline_dir, samples=samples, what=what,
        min_reads=min_reads, design=calibration_design)

    if (detrend_samples)
        cal <- detrend_samples_and_calibrate(wei, calibration_design)
    else
        cal <- weitrix::weitrix_calibrate_all(wei, calibration_design)
    
    result <- perform_weitrix_test(cal, design=design, coef=coef, contrast=contrast, fdr=fdr, step=step, dispersion_est=dispersion_est)

    total_reads <- rowSums(weitrix::weitrix_weights(wei))
    ave_expr <- log2( total_reads*1e6 / sum(total_reads) )
    result$table$AveExpr <- ave_expr[result$table$index]

    #ave_tail <- weitrix::weitrix_components(cal, design=~1, verbose=F)$row[,1]
    #result$table$AveTail <- ave_tail[result$table$index]

    result$pipeline_dir <- pipeline_dir
    result$effect_desc <- paste0(result$effect_desc, " of tail length")
    result$title <- paste0(
        title, " - ", 
        result$effect_desc, " - ",
        what, " - ",
        if (detrend_samples) "detrended samples - ", 
        "at least ", min_reads, " reads in sufficient samples")
    result$magnitude_column <- "row_mean"
    result$magnitude_desc <- "Average tail length"

    text <- capture.output({
        cat("\nDesign matrix\n")
        design <- result$design
        rownames(design) <- colnames(wei)
        print(design)
        cat("\nContrasts\n")
        contrasts <- result$contrasts
        if (!is.null(contrasts))
            rownames(contrasts) <- colnames(design)
        print(contrasts)
    })

    result$diagnostics <- list()
    result$diagnostics[["Details"]] <- list(
        text=text)

    result$diagnostics[["Calibration vs sample"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, cat=col)))

    result$diagnostics[["Calibration vs tail length"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=mu)))

    result$diagnostics[["Calibration vs sample and tail length"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=mu, cat=col)))

    result$diagnostics[["Calibration vs read count"]] <- list(
        plot=compact_plot(weitrix::weitrix_calplot(cal, design, covar=log2(weitrix::weitrix_weights(wei))) + labs(x="log2 reads with tail")))

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


test_diff_tail_weitrix_peak_50 <- function(...) test_diff_tail_weitrix(..., min_reads=50, what="peakwise")
test_diff_tail_weitrix_peak_50_two_coef <- as_two_coef_test(test_diff_tail_weitrix_peak_50)

test_diff_tail_weitrix_peak_100 <- function(...) test_diff_tail_weitrix(..., min_reads=100, what="peakwise")
test_diff_tail_weitrix_peak_100_two_coef <- as_two_coef_test(test_diff_tail_weitrix_peak_100)

test_diff_tail_weitrix_peak_200 <- function(...) test_diff_tail_weitrix(..., min_reads=200, what="peakwise")
test_diff_tail_weitrix_peak_200_two_coef <- as_two_coef_test(test_diff_tail_weitrix_peak_200)



test_diff_tail_weitrix_50_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=50)
test_diff_tail_weitrix_100_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=100)
test_diff_tail_weitrix_200_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=200)
test_diff_tail_weitrix_primary_peak_50_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=50, what="primarypeakwise")
test_diff_tail_weitrix_primary_peak_100_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=100, what="primarypeakwise")
test_diff_tail_weitrix_primary_peak_200_detrend <- function(...) test_diff_tail_weitrix(..., detrend_samples=TRUE, min_reads=200, what="primarypeakwise")
