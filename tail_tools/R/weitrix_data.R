
# Filter for features where the design can actually be fitted using samples with min_reads reads
weitrix_filter_full_rank <- function(wei, design, min_reads) {
    full_rank <- function(selection) {
        if (sum(selection) < ncol(design)) return(FALSE)
        sum(abs(svd(design[selection,,drop=F])$d)>=1e-10) == ncol(design)
    }
    good <- apply(weitrix::weitrix_weights(wei)>=min_reads, 1, full_rank)
    wei[good,]
}

#' @export
pipeline_weitrix_shift <- function(
        pipeline_dir, samples=NULL,  
        antisense=F, colliders=F, non_utr=T, collapse_utr=F,
        min_reads=0, design=NULL) {
    data <- get_grouped_peaks(
        pipeline_dir,samples=samples,antisense=antisense,colliders=colliders,
        non_utr=non_utr,collapse_utr=collapse_utr,min_reads=min_reads,min_group=2)
    
    wei <- weitrix::counts_shift(data$counts, data$grouping)
    SummarizedExperiment::rowData(wei) <- cbind(
        SummarizedExperiment::rowData(wei), 
        data$gene[match(rownames(data$gene),rownames(wei)),,drop=FALSE])
    
    S4Vectors::metadata(wei)$display_members <-
        split(sub("-collider$","",data$grouping$name), data$grouping$group)

    # Remove any rows that will produce NA coefficients.
    # (Note: main filtering is earlier by min_reads against the total number of reads in the row)
    if (!is.null(design))
        wei <- weitrix_filter_full_rank(wei, design, 1)

    wei
}


#' @export
pipeline_weitrix_tail <- function(
        pipeline_dir, what="genewise", samples=NULL, min_reads=0, design=NULL) {
    
    data <- read_grouped_table(paste0(pipeline_dir, "/expression/", what, "/counts.csv"))

    wei <- weitrix::as_weitrix(as.matrix(data$Tail), as.matrix(data$Tail_count))

    SummarizedExperiment::rowData(wei) <- data$Annotation
    S4Vectors::metadata(wei)$weitrix$calibrate_all_formula <- "~splines::ns(mu,4)+splines::ns(log(weight),4)"

    if (!is.null(samples))
        wei <- wei[,samples]

    if (min_reads > 0)
        wei <- weitrix_filter_full_rank(wei, design, min_reads)

    wei
}


