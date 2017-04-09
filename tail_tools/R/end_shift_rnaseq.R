
#
# end-shift test adapted to RNA-Seq data
#

#' @export
bigwig_get <- function(bigwig_fwd, bigwig_rev, pos) {
   splitter <- paste0(GenomeInfoDb::seqnames(pos), BiocGenerics::strand(pos))
   assert_that(noNA(splitter))

   split_pos <- split(pos, splitter, drop=TRUE)
   
   result <- list()
   for(item in names(split_pos)) {
       this_pos <- split_pos[[item]]
       
       this_is_reverse <- as.logical(BiocGenerics::strand(this_pos)[1] == "-")
       
       bigwig <- if (this_is_reverse == "TRUE") bigwig_rev else bigwig_fwd
       
       #s <- summary(bigwig, pos, size=width(pos), type="mean")
       #vec <- mcols(s[[1]])$score
       
       #got <- import.bw(bigwig,selection=BigWigSelection(pos))
       #vec <- rep(got$score, width(got))
       
       # import.bw shuffles results by seqname!
       vecs <- rtracklayer::import.bw(bigwig,"bigWig",which=this_pos,as="NumericList") %>% as.list
       
       assert_that(length(vecs) == length(this_pos))
       
       vecs <- purrr::map(vecs, function(item) { item[is.na(item)] <- 0 ; item })
       
       if (this_is_reverse == "TRUE")
           vecs <- purrr::map(vecs, rev)
       
       result[[item]] <- vecs
   }
   
   result <- unsplit(result, splitter)
   
   assert_that(all( sapply(result,length) == BiocGenerics::width(pos) ))
   
   result
}



end_shift_rnaseq_batch <- function(samples, utrs, extended_utrs, exons) {
    condition <- samples$condition
    keep <- !is.na(condition)
    
    samples <- samples[keep,,drop=F]
    condition <- condition[keep]

    if (is.null(samples$group))
        samples$group <- rep(1,nrow(samples))
    
    sample_splitter <- split(seq_along(samples$group), samples$group)
    #todo: fix extension
        
    cat("get\n")
    all_covers <- purrr::map(seq_len(nrow(samples)), function(i) {
        bigwig_get(samples$cover_fwd[i],samples$cover_rev[i], extended_utrs)
    })
    all_ends <- purrr::map(seq_len(nrow(samples)), function(i) {
        bigwig_get(samples$end_fwd[i],samples$end_rev[i], extended_utrs)
    })
    cat("think\n")

    all_forward <- is_forward(utrs)

    result <- purrr::map_df(seq_along(utrs), function(i) {
        forward <- all_forward[i]
        utr_width <- BiocGenerics::width(utrs)[i]
        extended_width <- BiocGenerics::width(extended_utrs)[i]
        
        exonic_depth <- rep(0L, extended_width)
        if (utr_width < extended_width)
            exonic_depth[utr_width+1L] <- 1L
        
        this_exons <- exons[[ utrs[i]$ID ]]
        exon_starts <- 
            if (forward) 
                start(this_exons) - start(utrs)[i]
            else
                end(utrs)[i] - end(this_exons)
        exon_widths <- BiocGenerics::width(this_exons)
        for(j in seq_along(exon_starts)) {
            this_start <- exon_starts[j]
            this_width <- exon_widths[j]
            if (this_start < 0L) {
                this_width <- this_width + this_start
                this_start <- 0L
            }
            if (this_width > 0L) {
                exonic_depth[this_start+1L] <- exonic_depth[this_start+1L] + 1L
                this_end <- this_start + this_width
                if (this_end < extended_width)
                    exonic_depth[this_end+1L] <- exonic_depth[this_end+1L] - 1L
            }
        }
        exonic_depth <- cumsum(exonic_depth)
        exonic <- exonic_depth > 0L
    
        covers <- lapply(seq_len(nrow(samples)), function(j) all_covers[[j]][[i]])
        
        total_cover <- purrr::reduce(covers, `+`)    
        
        cutoff <- 0.01 * max(total_cover[seq_len(utr_width)])
        distal <- min(length(total_cover), utr_width) + 1
        while(distal <= length(total_cover) && total_cover[distal] > cutoff) 
            distal <- distal + 1
        distal <- distal - 1
        
        range <- seq_len(distal)
        exonic_range <- which(exonic[range])
        
        mat <- lapply(seq_len(nrow(samples)), function(j) all_ends[[j]][[i]][exonic_range]) 
        mat <- do.call(cbind,mat)
        
        result <- combined_r(mat, condition, sample_splitter)
        
        #steps <- function(cond) {
        #    span <- lapply(covers[condition == cond], function(item) item[range]) %>%  #!!!!!
        #        { purrr::reduce(.,`+`) }
        #    #reg <- sort(span, decreasing=T)
        #    #delayed <- c(reg[-1],reg[length(reg)])
        #    #reg - delayed
        #    span
        #}
        #
        #a <- steps(F)
        #b <- steps(T)
        #
        #if (max(a) < 20 || max(b) < 20)
        #    return(NULL)
        #
        #na <- sum(a)
        #nb <- sum(b)
        #
        ##if (na < 30 || nb < 30)
        #
        #ca <- cumsum(a)
        #cb <- cumsum(b)
        #U2 <- sum(ca*b) - sum(cb*a)
        #r <- U2/(na*nb)
        
        data_frame(
            r = result$r,
            var = result$var,
            mean_reads = mean(apply(mat,2,sum)),
            min_reads = min(apply(mat,2,sum)),
            transcript_id = utrs$ID[i],
            transcript_name = utrs$Name[i],
            gene_id = utrs$gene_id[i],
            gene_name = utrs$gene[i],
            description = utrs$description[i],
            chromosome = GenomeInfoDb::seqnames(utrs)[i] %>% as.character,
            strand = BiocGenerics::strand(utrs)[i] %>% as.character,
            five_prime = 
                ( if (forward) BiocGenerics::start(utrs)[i] else BiocGenerics::end(utrs)[i] ),
            three_prime_given = 
                ( if (forward) BiocGenerics::end(utrs)[i] else BiocGenerics::start(utrs)[i] ),
            three_prime_extended = 
                ( if (forward) BiocGenerics::start(utrs)[i]+distal-1 else BiocGenerics::end(utrs)[i]-distal+1 )
        )
    })
    
    result    
}


#' End-shift test for RNA-Seq data
#'
#' @param samples should be a data frame with columns: name = name of sample, condition = logical vector with NA for samples to ignore, bigwig filenames giving depth of coverage cover_fwd, cover_rev. 
#'
#' @param utrs a GRanges of 3' UTRS. Should have metadata columns: ID, Name, gene_id, gene, description.
#'
#' @param extended_utrs if given tries extending 3' UTRS beyond the annotated end point. A GRanges of validly extended 3' UTRS, in the same order as utrs. Suggestion is UTRs extended up to 10kb but not into the following gene on the same strand.
#'
#' @param exons if given, restricts testing to exonic regions. A GRanges of transcript exons, which should have metadata column Parent corresponding to the ID in the utrs parameter.
#'
#' @export
end_shift_rnaseq <- function(samples, utrs, extended_utrs=NULL, exons=NULL, ci=0.95, min_min_reads=20) {
    if (is.null(extended_utrs))
        extended_utrs <- utrs
    
    if (is.null(exons)) {
        exons <- utrs
        S4Vectors::mcols(exons)$Parent <- exons$ID
    }
    
    # import() likes to make Parent a CharacterList
    unlisted <- BiocGenerics::unlist(exons$Parent)
    assert_that(length(unlisted) == length(exons))
    S4Vectors::mcols(exons)$Parent <- unlisted 

    exons <- GenomicRanges::split(exons, exons$Parent)

    reorder <- order(
        GenomeInfoDb::seqnames(utrs) %>% as.character, 
        BiocGenerics::strand(utrs) %>% as.character, 
        BiocGenerics::start(utrs) %>% as.integer)
    utrs <- utrs[ reorder ]
    extended_utrs <- extended_utrs[ reorder ]
    
    range <- seq_along(utrs)
    batches <- range %/% 1000
    split_ranges <- split(range, batches)
    result <- purrr::map_df(seq_along(split_ranges), function(i) {
        cat(i,"of",length(split_ranges),"\n")
        batch <- split_ranges[[i]]
        end_shift_rnaseq_batch(samples, utrs[batch], extended_utrs[batch], exons) 
    })
    
    ci_sds <- qnorm((1+ci)/2)    
    
    result <- result %>% 
        dplyr::filter_(~min_reads >= min_min_reads) %>%
        dplyr::mutate_(
            sd =~ sqrt(var),
            r_low =~ r - sd*ci_sds,
            r_high =~ r + sd*ci_sds,
            interest =~ ( abs(r) - sd*ci_sds ) %>% na_replace(-Inf)
        ) %>%
        dplyr::arrange_(~desc(interest)) %>%
        dplyr::mutate_( rank=~seq_len(n()) )
    
    #arrange(-abs(r))
    
    gene_result <- result %>% 
        group_by_(~gene_id) %>%
        filter_(~seq_len(n()) == which.min(rank)) %>%
        ungroup() %>%
        mutate_( rank=~seq_len(n()) )

    list(
        samples = samples,
        transcripts = result,
        genes = gene_result,
        min_min_reads = min_min_reads
    )
}





