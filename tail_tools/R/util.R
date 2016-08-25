
#' @export
is_reverse <- function(granges) as.logical(strand(granges) == "-")

#' @export
is_forward <- function(granges) !is_reverse(granges)

#' @export
cached <- function(prefix, func, args_var=list(), args_const=list(), version=NULL) {
    filename <- paste0(prefix,".Rds")
    key <- list(args=args_var, version=version)
    
    if (file.exists(filename)) {
        cached_value <- readRDS(filename)
        if (identical(cached_value$key, key))
            return(cached_value$result)        
    }
    
    result <- do.call(func, c(args_var, args_const))
    
    saveRDS(list(key=key,result=result), filename)    
    
    result
}
