
            
default_plot_count <- function(df, is_grouped, norm_count_name, ...) { 
    print(
        ggplot(df, aes(x=sample, y=norm_count)) +
        geom_bar(stat="identity", color="#000000", fill="#cccccc") +
        ylab(norm_count_name) +
        xlab("") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
}


default_plot_tail <- function(df, is_grouped, ...) {
    df_bad <- filter_(df, ~is_low_tail_count)
    
    print(
        ggplot(df, aes(x=sample, y=tail)) +
        geom_bar(stat="identity", color="#000000", fill="#ccffcc") +
        geom_bar(data=df_bad, 
                 stat="identity", color="#000000", fill="#ff0000") +
        ylab("poly(A) tail length") +
        xlab("") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
    )
}


#' Shiny report for mPAT results
#'
#' @param filename A counts.csv file as produced by "tail-tools analyse-tail-lengths:".
#'
#' @param normalizing_gene ID of housekeeping gene to normalize to. A vector of IDs may be given. NULL means perform no normalization.
#'
#' @param title Title of Shiny report.
#'
#' @param pipeline_dir Pipeline output directory, if you want to view detailed poly(A) tail length distributions extracted from BAM files.
#'
#' @param max_tail Default limit on tail length for tail length distribution plot.
#'
#' @grouping If you want to be able to group samples, a data frame with a "group" and a "sample" column.
#'
#' @return A shiny appobj. Printing it runs the app.
#'
#' @export
shiny_mpat <- function(
        filename, 
        normalizing_gene=NULL, 
        title="mPAT results",
        pipeline_dir=NULL,
        max_tail=300,
        grouping=NULL,
        
        plot_count=NULL,
        plot_tail=NULL) {
    library(shiny)
    library(nesoni)
    library(reshape2)
    library(varistran)
    library(ggplot2)
    library(gridExtra)
    library(rtracklayer)
    library(dplyr)
    library(tidyr)
    library(viridis)

    tables <- read.grouped.table(filename)
    genes <- rownames(tables$Count)
    samples <- colnames(tables$Count)
    
    raw_data <-
        meltdown(tables$Count, "sample","gene","count") %>%
        left_join(meltdown(tables$Tail, "sample","gene","tail"), 
                  c("sample","gene")) %>%
        left_join(meltdown(tables$Tail_count, "sample", "gene", "tail_count"), 
                  c("sample","gene"))
    
    
    normalizers <- c("Reads Per Million", genes)    
    
    
    have_plotters <- !is.null(plot_count) | !is.null(plot_tail)
    if (is.null(plot_count)) plot_count <- default_plot_count
    if (is.null(plot_tail)) plot_tail <- default_plot_tail
    
    have_grouping <- !is.null(grouping)
    
    have_bams <- !is.null(pipeline_dir)
    if (have_bams) {
        bam_filenames <- pipeline_bams(pipeline_dir)
    
        # Create relevant genomic ranges
        ranges <- GRanges(
            tables$Annotation$chromosome, 
            IRanges(tables$Annotation$start, tables$Annotation$end),
            ifelse(tables$Annotation$strand < 0, "-","+"))
        names(ranges) <- rownames(tables$Annotation)
        
        ranges$three_prime <- ifelse(strand(ranges) == "+", end(ranges), start(ranges))
        
        bounds <- ranges %>% as.data.frame
        bounds$name <- names(ranges)
        bounds <- bounds %>% 
            arrange(seqnames,strand,three_prime) %>% 
            group_by(seqnames,strand) %>%
            mutate(
                three_prime_min = pmax(-Inf, (three_prime + lag(three_prime))*0.5, na.rm=T),
                three_prime_max = pmin(Inf, (three_prime + lead(three_prime))*0.5, na.rm=T)
            ) %>%
            ungroup()
        matching <- match(names(ranges), bounds$name)
        ranges$three_prime_min <- bounds$three_prime_min[matching]
        ranges$three_prime_max <- bounds$three_prime_max[matching]        
    }


    # Plots

    overview <- shiny_plot(prefix="overview", dlname="overview", width=800,height=800, function(env) withProgress(message="Plotting overview", {
        normed <- env$normed()
        
        print( 
            ggplot(normed$norm_data,aes(x=sample,y=log2_fold_norm_count)) + 
            geom_point() + 
            facet_wrap(~ gene, drop=F) + 
            ylab(ifelse(normed$is_differential,"log2 fold change",paste0("log2 ",normed$norm_count_name))) +
            xlab("") +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )        
    }))
    
    overview_tail <- shiny_plot(prefix="overview_tail", dlname="overview_tail", width=800,height=800, function(env) withProgress(message="Plotting overview", {
        normed <- env$normed()
        df <- filter_(normed$norm_data, ~tail_count > 0)
        df_bad <- filter_(df, ~is_low_tail_count)
        print(
            ggplot(df,aes(x=sample,y=relative_tail)) + 
            geom_point() + 
            geom_point(data=df_bad,color="red") +
            facet_wrap(~ gene, drop=F) +
            ylab(ifelse(normed$is_differential,"change in poly(A) tail length","poly(A) tail length")) +
            xlab("") + 
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )
    }))
    
    individual_count <- shiny_plot(prefix="individual_count", dlname="individual-count", width=400,height=400, function(env) {
        normed <- env$normed()
        individual_gene <- env$input$individual_gene
        df <- filter_(normed$norm_data, ~gene == individual_gene)
        
        if (have_plotters && env$input$plotter_custom)
           plotter <- plot_count
        else
           plotter <- default_plot_count
        
        plotter(df=df, is_grouped=normed$is_grouped, norm_count_name=normed$norm_count_name)
    })

    individual_tail <- shiny_plot(prefix="individual_tail", dlname="individual-tail", width=400,height=400, function(env) {
        normed <- env$normed()
        individual_gene <- env$input$individual_gene
        df <- filter_(normed$norm_data, ~gene == individual_gene)
        
        if (have_plotters && env$input$plotter_custom)
           plotter <- plot_tail
        else
           plotter <- default_plot_tail
        
        plotter(df, is_grouped=normed$is_grouped)        
    })
    
    if (have_bams) {
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
    }
    
    

    # App

    ui <- div(
        titlePanel(title),
        navlistPanel(
            widths=c(2,10),
            well=FALSE,
            tabPanel("Select",
                selectInput("normalizing_gene","Normalize to",choices=normalizers,selected=normalizing_gene,multiple=TRUE),
                p("If multiple genes are selected, the geometric mean will be used."),
                br(),
                
                selectInput("genes","Genes to display",choices=genes,selected=genes,multiple=TRUE,width="100%"),
                p(
                    actionButton("genes_all", "All")
                ),
                br(),
                
                selectInput("samples","Samples to display",choices=samples,selected=samples,multiple=TRUE,width="100%"),
                p(
                    actionButton("samples_all", "All"),
                    if (have_grouping) checkboxInput("group_samples", "Group samples", value=FALSE)
                ),
                br(),
                
                selectInput("normalizing_samples","Samples to normalize to in overview",choices=samples,selected=c(),multiple=TRUE,width="100%"),
                p(
                    actionButton("normalizing_samples_all", "All"),
                    actionButton("normalizing_samples_none", "None")
                ),
                p("Tail lengths and log2 read counts are shown relative to the average of these samples, in the overview."), 
                p("Select none to show absolute log2 read counts and tail lengths in overview."),
                br(),
                
                numericInput("highlight_low", "Highlight in red tail lengths based on less than this many reads", 50),
                br(),
                
                h1("Normalization"),
                p("Read counts are divided by this to obtain normalized read counts."),
                htmlOutput("normalizer_table")
            ),
            tabPanel("Overview",
                h1("Expression levels"),
                p("Note these are shown on a log scale."),
                overview$component_ui,
                h1("Tail lengths"),
                p("Where there are few poly(A) reads, the length is shown in red."),
                overview_tail$component_ui
            ),
            tabPanel("Individual genes",
                selectInput("individual_gene", "Gene", choices=sort(genes)),
                if (have_plotters) checkboxInput("plotter_custom", "Use customized plots", value=TRUE),
                p("Note expression levels are *not* log transformed in this plot."),
                fluidRow(
                    column(6, individual_count$component_ui),
                    column(6, individual_tail$component_ui)),
                br(),
                br(),
                if (have_bams) h2("Detail"),
                if (have_bams) fluidRow(
                    column(6,
                        checkboxInput("tail_tail", "Only show reads with poly(A) tail", value=TRUE),
                        checkboxInput("tail_percent", "Show as percent within each sample", value=TRUE),
                        checkboxInput("tail_log", "Log scale", value=FALSE)),
                    column(6, conditionalPanel("input.tail_tail",
                        numericInput("tail_min", "Minimum tail length to include in plots", 4)))),
                if (have_bams) br(),
                if (have_bams) h3("Tail length distribution"),
                if (have_bams) fluidRow(
                    column(2, radioButtons(
                        "tail_style", "Display", 
                        choices=c("Cumulative","Density","Heatmap"), selected="Cumulative")),
                    column(3, numericInput("tail_max", "Maximum tail length", max_tail, min=1)),
                    column(3, conditionalPanel("input.tail_style != 'Cumulative'", 
                        numericInput("tail_bin", "Tail length bin size", 1, min=1)))),
                if (have_bams) tail_distribution$component_ui,                
                if (have_bams) br(),
                if (have_bams) h3("Templated sequence"),
                if (have_bams) end_distribution$component_ui
            )
        )
    )
    
    server <- function(env) {
        input <- env$input
        output <- env$output
        
        observeEvent(input$genes_all, {
            updateSelectInput(env$session, "genes", selected=genes)
        })
    
        observeEvent(input$samples_all, {
            updateSelectInput(env$session, "samples", selected=samples)
        })
    
        observeEvent(input$normalizing_samples_all, {
            updateSelectInput(env$session, "normalizing_samples", selected=samples)
        })
    
        observeEvent(input$normalizing_samples_none, {
            updateSelectInput(env$session, "normalizing_samples", selected=character(0))
        })
    
        env$normed <- reactive({
             normalizing_gene <- input$normalizing_gene
             selected_samples <- input$samples
             selected_genes <- input$genes
             highlight_low <- input$highlight_low
                          
             if (length(normalizing_gene) == 0) {
                 norm_count_name = "Count"
                 normalizer <- tibble(sample=factor(samples,samples), normalizer=1)
             
             } else if (identical(normalizing_gene, "Reads Per Million")) {
                 norm_count_name = "RPM"
                 normalizer <- raw_data %>%
                     group_by_(~sample) %>%
                     summarize_(normalizer =~ sum(count)/1000000)
                     
             } else {
                 if ("Reads Per Million" %in% normalizing_gene)
                     stop("Can't mix RPM and normalizing genes.")
             
                 norm_count_name = "Normalized count"
                 normalizer <- raw_data %>%
                     filter_(~ gene %in% normalizing_gene) %>%
                     group_by_(~ sample) %>%
                     summarize_(count =~ exp(mean(log(count)))) %>%
                     mutate_(normalizer =~ count / mean(count)) %>%
                     select_(~sample, ~normalizer)
             }
             
             norm_data <- raw_data %>%
                 left_join(normalizer, "sample") %>%
                 mutate_(
                     norm_count =~ count / normalizer,
                     log2_norm_count =~ log2(norm_count + 0.5),
                     log2_fold_norm_count =~ log2_norm_count,
                     relative_tail =~ tail,
                     is_low_tail_count =~ tail_count < highlight_low
                 )
             
             is_differential <- length(env$input$normalizing_samples) > 0
             
             if (is_differential) {
                 baseline <- norm_data %>%
                     filter_(~sample %in% env$input$normalizing_samples) %>%
                     group_by_(~gene) %>%
                     summarize_(
                         baseline =~ mean(log2_fold_norm_count),
                         baseline_tail =~ mean(tail, na.rm=TRUE)
                     ) %>%
                     select_(~gene, ~baseline, ~baseline_tail)
             
                 norm_data <- norm_data %>%
                     left_join(baseline, "gene") %>%
                     mutate_(
                         log2_fold_norm_count =~ log2_fold_norm_count - baseline,
                         relative_tail =~ tail - baseline_tail
                     )
             }     
             
             norm_data <- norm_data %>%
                 filter_(
                     ~sample %in% selected_samples,
                     ~gene %in% selected_genes
                 ) %>%
                 mutate_(
                     sample =~ factor(sample, selected_samples),
                     gene =~ factor(gene, selected_genes)
                 )
             
             is_grouped <- have_grouping && env$input$group_samples
             respectful_grouping <- NULL
             if (is_grouped) {
                 respectful_grouping <- grouping %>% 
                      filter_(~sample %in% selected_samples) %>%
                      mutate_(sample =~ factor(sample, selected_samples)) %>%
                      arrange_(~sample) %>%
                      mutate_(group =~ factor(group, unique(group)))
                 
                 norm_data <- norm_data %>%
                      left_join(respectful_grouping, "sample") %>%
                      group_by_(~group, ~gene) %>%
                      summarize_(
                          norm_count =~ mean(norm_count),
                          log2_norm_count =~ mean(log2_norm_count), #Hmm
                          log2_fold_norm_count =~ mean(log2_fold_norm_count),
                          tail =~ mean(tail, na.rm=TRUE),
                          tail_count =~ mean(tail_count),
                          relative_tail =~ mean(relative_tail, na.rm=TRUE),
                          is_low_tail_count =~ any(is_low_tail_count)
                      ) %>%
                      ungroup() %>%
                      select_(sample=~group, ~gene, ~norm_count, ~log2_norm_count, ~log2_fold_norm_count, ~tail, ~tail_count, ~relative_tail, ~is_low_tail_count) 
             }
             
             list(
                 norm_data=norm_data,
                 normalizer=normalizer,
                 norm_count_name=norm_count_name,
                 is_differential=is_differential,
                 is_grouped=is_grouped,
                 grouping=respectful_grouping
             )
        })
    
        output$normalizer_table <- renderTable({
             env$normed()$normalizer
        }, digits=6)


        if (have_bams) {
            env$read_info_basic <- reactive(withProgress(message="Reading tail lengths", {
                this_samples <- env$input$samples
                query <- ranges[ env$input$individual_gene ]
                
                tibble(
                    sample=factor(this_samples, this_samples),
                    info=lapply(bam_filenames[this_samples], read_info, query)
                ) %>%
                unnest(info)
            }))

            env$read_info <- reactive({
                this_samples <- env$input$samples
                read_info <- env$read_info_basic()
                
                if (env$input$tail_tail)
                    read_info <- read_info %>%
                        filter(length >= env$input$tail_min)
                
                if (env$input$tail_percent) {
                    normalizer <- read_info %>% group_by(sample) %>% summarize(normalizer=sum(n))
                    describer <- "Percent reads"
                    labeller <- function(x) paste0(sapply(x*100,scales::comma),"%")
                } else {
                    normalizer <- env$normed()$normalizer
                    describer <- if (length(env$input$normalizing_gene) == 0) "Count" else env$normed()$norm_count_name
                    labeller <- function(x) sapply(x,scales::comma)
                }
                
                samples_called <- "Sample"
                if (have_grouping && env$input$group_samples) {
                    samples_called <- "Group"
                    respectful_grouping <- env$normed()$grouping
                    this_samples <- levels(respectful_grouping$group)                    
                    read_info <- read_info %>%
                        left_join(respectful_grouping, "sample") %>%
                        group_by(group, length, width) %>%
                        summarize(n=sum(n)) %>%
                        ungroup() %>%
                        select(sample=group, length, width, n)
                }

                if (have_grouping && env$input$group_samples) {
                    normalizer <- normalizer %>%
                        left_join(respectful_grouping, "sample") %>%
                        group_by(group) %>%
                        summarize(normalizer = sum(normalizer)) %>%
                        select(sample=group, normalizer)
                }
                
                transformer <- identity
                if (env$input$tail_log) {
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
        }
            
        overview$component_server(env)
        overview_tail$component_server(env)
        individual_count$component_server(env)
        individual_tail$component_server(env)
        if (have_bams) tail_distribution$component_server(env)
        if (have_bams) end_distribution$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}



#' @export
shiny_repat <- function(...) shiny_mpat(...)

