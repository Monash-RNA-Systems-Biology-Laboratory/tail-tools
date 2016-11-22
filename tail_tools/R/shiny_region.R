
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
        mutate(bin = floor((length+0.5)/stride)) %>%
        group_by(bin, add=TRUE) %>%
        summarize(n = sum(n)) %>%
        mutate(
            length_low = bin*stride-0.5,
            length_mid = (bin+0.5)*stride-0.5,
            length_high = (bin+1)*stride-0.5
        )
}



#' Display information about a genomic region
#'
shiny_region <- function(
        region,         
        pipeline_dir,
        max_tail=300,
        grouping=NULL,        
        prefix="") {
    ns <- NS(prefix)
    
    # A GRanges
    region <- ensure_reactable(region)
    
    have_grouping <- !is.null(grouping)
    
    bam_filenames <- pipeline_bams(pipeline_dir)

    # Tail lengths
    tail_distribution <- shiny_plot(
        prefix="tail_distribution", dlname="tail-distribution", width=800, 
        function(env) withProgress(message="Plotting tail distribution", {                
            tail_bin <- max(1,env$input$tail_bin)
            tail_max <- max(1,env$input$tail_max)
            
            this_samples <- env$read_info()$this_samples
            samples_called <- env$read_info()$samples_called
            normalizer <- env$read_info()$normalizer
            describer <- env$read_info()$describer
            labeller <- env$read_info()$labeller
            transformer <- env$read_info()$transformer
            tail_lengths <- env$read_info()$read_info %>% 
                group_by(sample, length) %>% summarize(n=sum(n)) %>% ungroup()
            
            tail_lengths <- tail_lengths %>%
                left_join(normalizer, "sample") %>%
                mutate(n = n / normalizer)
            
            if (env$input$tail_style == "Cumulative") {
                print(
                    tail_lengths %>% 
                    arrange(sample, desc(length)) %>% 
                    group_by(sample) %>%
                    mutate(
                        cumn = cumsum(n),
                        cumn_lag = dplyr::lag(cumn,1,0),
                        length_lead = dplyr::lead(length,1,0)
                    ) %>%
                    ungroup() %>%
                    ggplot(aes(color=sample)) +
                    #geom_point(aes(x=length,y=cumn_lag)) +
                    ggplot2::geom_segment(aes(x=length,xend=length,y=transformer(cumn),yend=transformer(cumn_lag))) +
                    ggplot2::geom_segment(aes(x=length_lead,xend=length,y=transformer(cumn),yend=transformer(cumn))) +
                    scale_x_continuous(lim=c(0,tail_max), oob=function(a,b)a) +
                    scale_y_continuous(labels = labeller) +
                    labs(x="poly(A) tail length", y=describer, color=samples_called) +
                    theme_bw()
                )
            } else if (env$input$tail_style == "Density") {
                print(
                    tail_lengths %>%
                    complete(sample=factor(this_samples,this_samples), length=seq(0,tail_max), fill=list(n=0)) %>%
                    group_by(sample) %>% bin_lengths(tail_bin) %>% ungroup() %>%
                    ggplot(aes(color=sample,group=sample,x=length_mid,y=transformer(n))) + 
                    geom_line() +
                    scale_x_continuous(lim=c(0,tail_max), oob=function(a,b)a) +
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
                    group_by(sample) %>% bin_lengths(tail_bin) %>% ungroup() %>%
                    ggplot(aes(x=sample,y=length_mid,fill=transformer(n))) + 
                    geom_tile(height=tail_bin) +
                    scale_y_continuous(lim=c(0,tail_max), oob=function(a,b)a) +
                    scale_fill_viridis(guide=FALSE) +
                    labs(x="", y="poly(A) tail length") +
                    theme_minimal() +
                    theme(panel.grid=element_blank(), 
                          axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
                )
            }
        }))
    
    
    # 3' end positions
    end_distribution <- shiny_plot(
        prefix="end_distribution", dlname="end-distribution", width=800, 
        function(env) withProgress(message="Plotting 3' end distribution", {                
            this_samples <- env$read_info()$this_samples
            samples_called <- env$read_info()$samples_called
            normalizer <- env$read_info()$normalizer
            describer <- env$read_info()$describer
            labeller <- env$read_info()$labeller
            transformer <- env$read_info()$transformer
            widths <- env$read_info()$read_info %>% 
                group_by(sample, width) %>% summarize(n=sum(n)) %>% ungroup()
            
            widths <- widths %>%
                left_join(normalizer, "sample") %>%
                mutate(n = n / normalizer)
            
            print(
                    widths %>%
                    arrange(sample, desc(width)) %>% 
                    group_by(sample) %>%
                    mutate(
                        cumn = cumsum(n),
                        cumn_lag = dplyr::lag(cumn,1,0),
                        width_lead = dplyr::lead(width,1,0)
                    ) %>%
                    ungroup() %>%
                    ggplot(aes(color=sample)) +
                    #geom_point(aes(x=length,y=cumn_lag)) +
                    ggplot2::geom_segment(aes(x=width,xend=width,y=transformer(cumn),yend=transformer(cumn_lag))) +
                    ggplot2::geom_segment(aes(x=width_lead,xend=width,y=transformer(cumn),yend=transformer(cumn))) +
                    #scale_x_continuous(lim=c(0,tail_max), oob=function(a,b)a) +
                    scale_y_continuous(labels = labeller) +
                    labs(x="Position relative to primer start", y=describer, color=samples_called) +
                    theme_bw()
            )
        }))

        
    ui <- div(
        textInput(ns("region"), "Region"),
        fluidRow(
            column(6,
                checkboxInput(ns("tail_tail"), "Only show reads with poly(A) tail", value=TRUE),
                checkboxInput(ns("tail_percent"), "Show as percent within each sample", value=TRUE),
                checkboxInput(ns("tail_log"), "Log scale", value=FALSE)),
            column(6, conditionalPanel("input.tail_tail",
                numericInput(ns("tail_min"), "Minimum tail length to include in plots", 4)))),
        br(),
        h3("Tail length distribution"),
        fluidRow(
            column(2, radioButtons(
                ns("tail_style"), "Display", 
                choices=c("Cumulative","Density","Heatmap"), selected="Cumulative")),
            column(3, numericInput(ns("tail_max"), "Maximum tail length", max_tail, min=1)),
            column(3, conditionalPanel(
                paste0("input['",ns("tail_style"),"' != 'Cumulative'"), 
                numericInput("tail_bin", "Tail length bin size", 1, min=1)))),
        tail_distribution$component_ui,                
        br(),
        h3("Templated sequence"),
        end_distribution$component_ui)

    
    server <- function(env) {
        # Push acceptimator
        observe({
            new_region <- region(env)
            if (length(new_region) != 1) return()
            
            updateTextInput(ns("region"), value=as.character(new_region))
        })
        
        
    
        env$output[[ns("description")]] <- renderUI({
            div("Hello", as.character(region(env)) )
        })    
    
        tail_distribution$component_server(env)
        end_distribution$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}





