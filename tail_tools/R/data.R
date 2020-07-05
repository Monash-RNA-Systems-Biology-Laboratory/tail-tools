

#
# Loading and manipulating data
#

factor_retaining_order <- function(chr) factor(chr, unique(chr))



#' Get BAM filenames from a Tail Tools pipeline output directory
#'
#' @export
pipeline_bams <- function(pipeline_dir) {
    meta <- jsonlite::fromJSON(file.path(pipeline_dir,"plotter-config.json"))
    bam_filenames <- meta$samples$bam
    names(bam_filenames) <- meta$samples$name
    bam_filenames
}


#' Get sample information from a Tail Tools pipeline output directory
#'
#' @export
pipeline_samples <- function(pipeline_dir) {
    meta <- jsonlite::fromJSON(file.path(pipeline_dir,"plotter-config.json"))
    samples <- meta$samples

    tags_lists <- samples$tags
    samples$tags <- NULL
    tags <- tags_lists %>% unlist %>% unique
    for(tag in tags)
        samples[[tag]] <- purrr::map_lgl(tags_lists, ~ tag %in% .)

    samples
}

#' Augment sample data frame with a column grouping several tags (logical columns) into a factor column
#'
#' If treatment=FALSE, contrasts is set to "contr.sum".
#'
#' @export
samples_group_tags <- function(samples, group_name, tags, treatment=TRUE) {
    result <- factor(rep(NA, nrow(samples)), tags)
    for(tag in tags) {
        vec <- samples[[tag]]
        !is.null(vec) || stop(paste0(tag," not present"))
        all(is.na(result[vec])) || stop("Tags not mutually exclusive.")
        result[vec] <- tag
    }

    if (!treatment)
        contrasts(result) <- "contr.sum"

    samples[[group_name]] <- result
    samples
}


#' @export
read_grouped_table <- function(filename, require=c(), default.group='All') {
    groups <- c()
    tab.separated <- FALSE
    
    skip <- 0
    
    f <- file(filename,'r')
    repeat {
        line <- readLines(f,1)
        if (length(line) == 0) break;
        if (substr(line,1,1) != '#') {
            tab.separated <- (length(grep('\t',line)) > 0)
            break;
        }
        
        skip <- skip + 1
        
        parts <- strsplit(line,',')[[1]]
        if (parts[1] == '#Groups') {
            groups <- parts[seq_len(length(parts)-1)+1]
        }
    }
    close(f)

    if (tab.separated)
        data <- read.delim(filename, skip=skip, check.names=FALSE)    
    else    
        data <- read.csv(filename, skip=skip, check.names=FALSE)
    
    rows <- data[,1]
    cols <- colnames(data)[seq_len(ncol(data)-1)+1]    
    data <- data[ ,seq_len(ncol(data)-1)+1, drop=FALSE]
    rownames(data) <- rows
    colnames(data) <- cols

    # === Fallbacks if groups not given ===
        
    if (!length(groups)) {
        rpkms <- grep('^RPKM', colnames(data))
        if (length(rpkms)) {
            # === Legacy count file ===
            n_samples <- rpkms[1] - 1
            groups <- c(
                rep('Count', n_samples),
                rep('RPKM', n_samples),
                rep('Annotation', ncol(data)-n_samples*2)
            )
        }
    }
    if (!length(groups)) {
        groups <- c(default.group)
    }

    
    i <- 2
    while(i <= ncol(data)) {
        if (is.null(groups[i]) || is.na(groups[i]) || groups[i] == '')
            groups[i] <- groups[i-1]
        i <- i + 1
    }
        
    groups <- factor(groups)
    
    result <- list()
    for(name in levels(groups)) {
        result[[name]] <- data[,groups == name,drop=FALSE]
    }
    
    for(item in require)
        if (is.null( result[[item]] ))
            stop('Table ',filename,' has no group ',item,'. Is it in the right format?')
    
    result
}





#' @export
read_tail_counts <- function(filename) {
   tab <- read_grouped_table(filename)
   
   features <- 
       tibble(
           feature = factor_retaining_order(rownames(tab$Count))) %>%
       cbind(tab$Annotation)
   
   samples <-
       tibble(
           sample = factor_retaining_order(colnames(tab$Count)))
   
   # observations
   obs <- 
       tibble(
           feature = rep(features$feature, nrow(samples)),
           sample = rep(samples$sample, each=nrow(features)),
           count = do.call(c, tab$Count),
           tail_count = do.call(c, tab$Tail_count),
           tail = do.call(c, tab$Tail))
   
   list(
       obs = obs, 
       features = features,
       samples = samples
   )
}


#' Subset features in a tail_counts
#'
#' @export
tail_counts_subset_features <- function(tc, features) {
    stopifnot(!any(duplicated(features)))   
    
    tc$features <- tc$features %>%
        dplyr::filter(feature %in% features) %>%
        dplyr::mutate(feature = factor(feature, features)) %>%
        dplyr::arrange(feature)
    
    tc$obs <- tc$obs %>%
        dplyr::filter(feature %in% features) %>%
        dplyr::mutate(feature = factor(feature, features))
    
    tc
}


#' Subset samples in a tail_counts
#'
#' @export
tail_counts_subset_samples <- function(tc, samples) {
    stopifnot(!any(duplicated(samples)))
    
    tc$samples <- tc$samples %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(sample = factor(sample, samples)) %>%
        dplyr::arrange(sample)
    
    tc$obs <- tc$obs %>%
        dplyr::filter(sample %in% samples) %>%
        dplyr::mutate(sample = factor(sample, samples))
    
    tc
}


#' Extract a column from tail_counts$obs as a matrix
#'
#' @export
tail_counts_get_matrix <- function(tc, column_name) {
    tapply(tc$obs[[column_name]], list(tc$obs$feature, tc$obs$sample), identity)
}


#' Augment a tail_counts by performing varistran's vst on the counts
#'
#' @export
tail_counts_vst <- function(tc) {
    mat <- tail_counts_get_matrix(tc, "count")
    
    vmat <- varistran::vst(mat)
    
    lib_size <- attr(vmat,"lib.size")
    
    tc$obs$norm_count <- c(t(t(mat) / (lib_size/mean(lib_size))))
    tc$obs$log2_norm_count <- c(vmat)
    tc$samples$true_lib_size <- attr(vmat, "true.lib.size")
    tc$samples$effective_lib_size <- attr(vmat, "lib.size")
    tc$samples$normalizer <- tc$samples$effective_lib_size / mean(tc$samples$effective_lib_size)
    tc$vst_dispersion <- attr(vmat, "dispersion") 
    tc$vst_method <- attr(vmat, "method")
    tc$vst_lib_size_method <- attr(vmat, "lib.size.method")
    tc$vst_cpm <- attr(vmat, "cpm")
    
    tc
}


#' Extract the vst values from a tail_counts
#'
#' @export
tail_counts_get_vst <- function(tc) {
    mat <- tail_counts_get_matrix(tc, "log2_norm_count")
    attr(mat, "true.lib.size") <- tc$sample$true_lib_size
    attr(mat, "lib.size") <- tc$sample$effective_lib_size
    attr(mat, "dispersion") <- tc$vst_dispersion
    attr(mat, "method") <- tc$vst_method
    attr(mat, "lib.size.method") <- tc$vst_lib_size_method
    attr(mat, "cpm") <- tc$vst_cpm
    
    mat
}





