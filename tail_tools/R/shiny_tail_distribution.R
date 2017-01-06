
#'
#' Get BAM filenames from a Tail Tools pipeline output directory
#'
#' @export
pipeline_bams <- function(pipeline_dir) {
    meta <- jsonlite::fromJSON(file.path(pipeline_dir,"plotter-config.json"))
    bam_filenames <- meta$samples$bam
    names(bam_filenames) <- meta$samples$name
    bam_filenames
}


read_info <- function(bam_filename, query) {
    library(dplyr)

    sbp <- Rsamtools::ScanBamParam(
        what=c("qname","pos","cigar","strand"),
        tag=c("AN"), 
        flag=Rsamtools::scanBamFlag(isMinusStrand=is_reverse(query)),
        which=query)
    
    result <- Rsamtools::scanBam(bam_filename, param=sbp)[[1]]
    AN <- as.integer(result$tag$AN) # Convert NULL to integer(0)
    AN[is.na(AN)] <- 0L
    
    if (is_reverse(query)) {
        three_prime <- result$pos
        width <- end(query) - three_prime + 1L
    } else {
        three_prime <- result$pos + 
            GenomicAlignments::cigarWidthAlongReferenceSpace(result$cigar) - 1L
        width <- three_prime - start(query) + 1L
    }
    good <- three_prime > query$three_prime_min & three_prime < query$three_prime_max
    
    tibble(
        length=AN[good],
        width=width[good]
    ) %>% 
    count(length, width)
}


meltdown <- function(mat, rows, columns, values) {
    result <- reshape2::melt(t(as.matrix(mat)))
    colnames(result) <- c(rows,columns,values)
    result
}


# sum n, binned by length, respecting existing group_by
bin_lengths <- function(df, stride) {
    df %>%
        mutate_(bin =~ floor((length+0.5)/stride)) %>%
        group_by_(~bin, add=TRUE) %>%
        summarize_(n =~ sum(n)) %>%
        mutate_(
            length_low =~ bin*stride-0.5,
            length_mid =~ (bin+0.5)*stride-0.5,
            length_high =~ (bin+1)*stride-0.5
        )
}



#' Display information about a genomic region
#'
shiny_tail_distribution <- function(
        peak_names=NULL,
        sample_normalizer=NULL,
        pipeline_dir,
        peaks=NULL,
        max_tail=300,
        grouping=NULL,        
        prefix="") {
    ns <- NS(prefix)

    if (is.null(peaks)) {
        count_filename <- paste0(pipeline_dir,"/expression/peakwise/counts.csv")    
        tc <- read_tail_counts(count_filename)
        peaks <- tc$features
    }
    
    have_grouping <- !is.null(grouping)
    
    bam_filenames <- pipeline_bams(pipeline_dir)
    samples <- names(bam_filenames)


    # Create relevant genomic ranges
    ranges <- GenomicRanges::GRanges(
        peaks$chromosome, 
        IRanges::IRanges(peaks$start, peaks$end),
        ifelse(peaks$strand < 0, "-","+"))
    names(ranges) <- peaks$feature
    
    ranges$three_prime <- ifelse(BiocGenerics::strand(ranges) == "+", end(ranges), start(ranges))
    
    bounds <- ranges %>% as.data.frame
    bounds$name <- names(ranges)
    bounds <- bounds %>% 
        arrange_(~seqnames,~strand,~three_prime) %>% 
        group_by_(~seqnames,~strand) %>%
        mutate_(
    	    three_prime_min =~ pmax(-Inf, (three_prime + lag(three_prime))*0.5, na.rm=T),
    	    three_prime_max =~ pmin(Inf, (three_prime + lead(three_prime))*0.5, na.rm=T)
        ) %>%
        ungroup()
    matching <- match(names(ranges), bounds$name)
    ranges$three_prime_min <- bounds$three_prime_min[matching]
    ranges$three_prime_max <- bounds$three_prime_max[matching]        
    



    # Tail lengths
    tail_distribution <- shiny_plot(prefix=ns("tail_distribution"), dlname="tail-distribution", width=800)
    
    
    ## 3' end positions
    #end_distribution <- shiny_plot(
    #    prefix=ns("end_distribution"), dlname="end-distribution", width=800, 
    #    function(env) withProgress(message="Plotting 3' end distribution", {                
    #        this_samples <- env$read_info()$this_samples
    #        samples_called <- env$read_info()$samples_called
    #        normalizer <- env$read_info()$normalizer
    #        describer <- env$read_info()$describer
    #        labeller <- env$read_info()$labeller
    #        transformer <- env$read_info()$transformer
    #        widths <- env$read_info()$read_info %>% 
    #            group_by(sample, width) %>% summarize(n=sum(n)) %>% ungroup()
    #        
    #        widths <- widths %>%
    #            left_join(normalizer, "sample") %>%
    #            mutate(n = n / normalizer)
    #        
    #        print(
    #                widths %>%
    #                arrange(sample, desc(width)) %>% 
    #                group_by(sample) %>%
    #                mutate(
    #                    cumn = cumsum(n),
    #                    cumn_lag = dplyr::lag(cumn,1,0),
    #                    width_lead = dplyr::lead(width,1,0)
    #                ) %>%
    #                ungroup() %>%
    #                ggplot(aes(color=sample)) +
    #                #geom_point(aes(x=length,y=cumn_lag)) +
    #                ggplot2::geom_segment(aes(x=width,xend=width,y=transformer(cumn),yend=transformer(cumn_lag))) +
    #                ggplot2::geom_segment(aes(x=width_lead,xend=width,y=transformer(cumn),yend=transformer(cumn))) +
    #                #scale_x_continuous(limits=c(0,tail_max), oob=function(a,b)a) +
    #                scale_y_continuous(labels = labeller) +
    #                labs(x="Position relative to primer start", y=describer, color=samples_called) +
    #                theme_bw()
    #        )
    #    }))

        
    ui <- function(request) div(
        uiOutput(ns("description")),
        shiny::selectInput(ns("samples"), "Select samples", 
                    selected=samples, choices=samples, multiple=TRUE),
        fluidRow(
            column(6,
                checkboxInput(ns("tail_tail"), "Only show reads with poly(A) tail", value=TRUE),
                checkboxInput(ns("tail_percent"), "Show as percent within each sample", value=TRUE),
                checkboxInput(ns("tail_log"), "Log scale", value=FALSE)),
            column(3, conditionalPanel(paste0("input.",ns("tail_tail")),
                numericInput(ns("tail_min"), "Minimum tail length to include in plots", 4)))),
        br(),
        #h3("Tail length distribution"),
        fluidRow(
            column(3, radioButtons(
                ns("tail_style"), "Display", 
                choices=c("Cumulative","Density","Heatmap"), selected="Cumulative")),
            column(3, numericInput(ns("tail_max"), "Maximum tail length", max_tail, min=1)),
            column(3, conditionalPanel(
                paste0("input['",ns("tail_style"),"' != 'Cumulative'"), 
                numericInput(ns("tail_bin"), "Tail length bin size", 1, min=1)))),
        br(),
        tail_distribution$component_ui(request))
        #br(),
        #h3("Templated sequence"),
        #end_distribution$component_ui(request))

    
    server <- function(env) { 
        e <- function(name) env[[ns(name)]]()
        i <- function(name) env$input[[ns(name)]]
        
        peak_names <- ensure_reactive(peak_names, ns("peak_names"), env)
        sample_normalizer <- ensure_reactive(sample_normalizer, ns("sample_normalizer"), env)
    
        env$output[[ns("description")]] <- renderUI({
            shiny::h2(paste(peak_names(), collapse=", "))
        })    
        

        read_info_basic <- reactive(withProgress(message="Reading tail lengths", {
            this_samples <- i("samples")
            
            lapply(peak_names(), function(peak)
                 tibble(
                     sample=factor(this_samples, this_samples),
                     info=lapply(bam_filenames[this_samples], read_info, ranges[peak])
                 ) %>%
                 unnest_(~info)
            ) %>%
            bind_rows
        }))


        read_info_full <- reactive({
            this_samples <- i("samples")
            read_info <- read_info_basic()
            
            if (i("tail_tail"))
                read_info <- read_info %>%
                    filter_(~length >= i("tail_min"))
            
            if (i("tail_percent")) {
                normalizer <- read_info %>% group_by_(~sample) %>% summarize_(normalizer=~sum(n))
                describer <- "Percent reads"
                labeller <- function(x) {
                    x <- x*100
                    x <- round(x,max(0,2-log10(x)))
                    paste0(sapply(x,scales::comma),"%")
                }
            } else {
                normalizer <- sample_normalizer()
                describer <- "Normalized count"
                labeller <- function(x) {
                    x <- round(x,max(0,2-log10(x)))
                    sapply(x,scales::comma)
                }
            }
            
            samples_called <- "Sample"
            #if (have_grouping && env$input$group_samples) {
            #    samples_called <- "Group"
            #    respectful_grouping <- env$normed()$grouping
            #    this_samples <- levels(respectful_grouping$group)                    
            #    read_info <- read_info %>%
            #        left_join(respectful_grouping, "sample") %>%
            #        group_by(group, length, width) %>%
            #        summarize(n=sum(n)) %>%
            #        ungroup() %>%
            #        select(sample=group, length, width, n)
            #}
            #
            #if (have_grouping && env$input$group_samples) {
            #    normalizer <- normalizer %>%
            #        left_join(respectful_grouping, "sample") %>%
            #        group_by(group) %>%
            #        summarize(normalizer = sum(normalizer)) %>%
            #        select(sample=group, normalizer)
            #}
            
            transformer <- identity
            if (i("tail_log")) {
                #transformer <- log2
                #if (env$input$tail_percent) describer <- "Proportion of reads"
                #describer <- paste("log2", describer)
                #labeller <- waiver()
                transformer <- log10
                old_labeller <- labeller
                labeller <- function(x) old_labeller(10**x)
            }
            
            list(
                read_info = read_info,
                normalizer = normalizer,
                describer = describer,
                labeller = labeller,
                transformer = transformer,
                this_samples = this_samples,
                samples_called = samples_called
            )            
        })
        
        
        env[[ns("tail_distribution-callback")]] <- function(env) withProgress(message="Plotting tail distribution", {                
            tail_bin <- max(1,i("tail_bin"))
            tail_max <- max(1,i("tail_max"))
            
            this_samples <- read_info_full()$this_samples
            samples_called <- read_info_full()$samples_called
            normalizer <- read_info_full()$normalizer
            describer <- read_info_full()$describer
            labeller <- read_info_full()$labeller
            transformer <- read_info_full()$transformer
            tail_lengths <- read_info_full()$read_info %>% 
                group_by_(~sample, ~length) %>% summarize_(n=~sum(n)) %>% ungroup()
            
            tail_lengths <- tail_lengths %>%
                left_join(normalizer, "sample") %>%
                mutate_(n =~ n / normalizer)
            
            if (i("tail_style") == "Cumulative") {
                print(
                    tail_lengths %>% 
                    arrange_(~sample, ~desc(length)) %>% 
                    group_by_(~sample) %>%
                    mutate_(
                        cumn =~ cumsum(n),
                        cumn_lag =~ dplyr::lag(cumn,1,0),
                        length_lead =~ dplyr::lead(length,1,0)
                    ) %>%
                    ungroup() %>%
                    ggplot(aes_(color=~sample)) +
                    #geom_point(aes(x=length,y=cumn_lag)) +
                    ggplot2::geom_segment(aes_(x=~length,xend=~length,y=~transformer(cumn),yend=~transformer(cumn_lag))) +
                    ggplot2::geom_segment(aes_(x=~length_lead,xend=~length,y=~transformer(cumn),yend=~transformer(cumn))) +
                    scale_x_continuous(limits=c(0,tail_max), oob=function(a,b)a) +
                    scale_y_continuous(labels = labeller) +
                    labs(x="poly(A) tail length", y=describer, color=samples_called) +
                    theme_bw()
                )
            } else if (i("tail_style") == "Density") {
                print(
                    tail_lengths %>%
                    complete(sample=factor(this_samples,this_samples), length=seq(0,tail_max), fill=list(n=0)) %>%
                    group_by_(~sample) %>% bin_lengths(tail_bin) %>% ungroup() %>%
                    ggplot(aes_(color=~sample,group=~sample,x=~length_mid,y=~transformer(n))) + 
                    geom_line() +
                    scale_x_continuous(limits=c(0,tail_max), oob=function(a,b)a) +
                    scale_y_continuous(labels = labeller) +
                    labs(x="poly(A) tail length", y=describer, color=samples_called) +
                    theme_bw()   
                )
            } else {
                print(
                    tail_lengths %>%
                    #group_by(sample) %>%
                    #mutate(n = n/max(n)) %>%
                    #ungroup() %>%
                    complete(sample=factor(this_samples,this_samples), length=seq(0,tail_max), fill=list(n=0)) %>%
                    group_by_(~sample) %>% bin_lengths(tail_bin) %>% ungroup() %>%
                    ggplot(aes_(x=~sample,y=~length_mid,fill=~transformer(n))) + 
                    geom_tile(height=tail_bin) +
                    scale_y_continuous(limits=c(0,tail_max), oob=function(a,b)a) +
                    scale_fill_viridis(guide=FALSE) +
                    labs(x="", y="poly(A) tail length") +
                    theme_minimal() +
                    theme(panel.grid=element_blank(), 
                          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
                )
            }
        })
    
        tail_distribution$component_server(env)
        #end_distribution$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}





