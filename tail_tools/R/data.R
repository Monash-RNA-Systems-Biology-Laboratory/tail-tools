

#
# Loading and manipulating data
#

factor_retaining_order <- function(chr) factor(chr, unique(chr))


read_tail_counts <- function(filename) {
   tab <- read.grouped.table(filename)
   
   features <- 
       data_frame(
           feature = factor_retaining_order(rownames(tab$Count))) %>%
       cbind(tab$Annotation)
   
   samples <-
       data_frame(
           sample = factor_retaining_order(colnames(tab$Count)))
   
   # observations
   obs <- 
       data_frame(
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
tail_counts_subset_features <- function(tc, features) {
    stopifnot(!any(duplicated(features)))   
    
    tc$features <- tc$features %>%
        dplyr::filter_(~ feature %in% features) %>%
        dplyr::mutate_(feature =~ factor(feature, features)) %>%
        dplyr::arrange_(~ feature)
    
    tc$obs <- tc$obs %>%
        dplyr::filter_(~ feature %in% features) %>%
        dplyr::mutate_(feature =~ factor(feature, features))
    
    tc
}


#' Subset samples in a tail_counts
#'
tail_counts_subset_samples <- function(tc, samples) {
    stopifnot(!any(duplicated(samples)))
    
    tc$samples <-
        filter_(~ sample %in% samples) %>%
        mutate_(sample =~ factor(sample, samples)) %>%
        arrange_(~ sample)
    
    tc$obs <-
        filter(~ sample %in% samples) %>%
        mutate_(sample =~ factor(sample, samples))
    
    tc
}


#' Extract a column from tail_counts$obs as a matrix
#'
tail_counts_get_matrix <- function(tc, column_name) {
    tapply(tc$obs[[column_name]], list(tc$obs$feature, tc$obs$sample), identity)
}


#' Augment a tail_counts by performing varistran's vst on the counts
#'
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





