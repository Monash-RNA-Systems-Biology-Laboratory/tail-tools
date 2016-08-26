
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


tail_lengths <- function(bam_filename, query) {
    library(dplyr)

    sbp <- Rsamtools::ScanBamParam(
        what=c("qname","pos","cigar","strand"),
        tag=c("AN"), 
        flag=Rsamtools::scanBamFlag(isMinusStrand=is_reverse(query)),
        which=query)
    
    result <- Rsamtools::scanBam(bam_filename, param=sbp)[[1]]
    if (is_reverse(query)) {
        three_prime <- result$pos
    } else {
        three_prime <- result$pos + 
            GenomicAlignments::cigarWidthAlongReferenceSpace(result$cigar) - 1L
    }
    good <- three_prime > query$three_prime_min & three_prime < query$three_prime_max
    
    tibble(length=na.omit(as.integer(result$tag$AN[good]))) %>% 
        count(length)
}


meltdown <- function(mat, rows, columns, values) {
    result <- reshape2::melt(t(as.matrix(mat)))
    colnames(result) <- c(rows,columns,values)
    result
}

#'
#' Shiny report for mPAT results
#'
#' @export
shiny_mpat <- function(
        filename, 
        normalizing_gene=NULL, 
        title="mPAT results",
        pipeline_dir=NULL,
        max_tail=300) {
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
    
    
    possible_normalizers <- c("None", sort(genes))
    if (is.null(normalizing_gene))
        normalizing_gene <- "None"
    
    
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
            ylab(ifelse(normed$is_differential,"log2 fold change","log2 normalized count")) +
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
    
    individual <- shiny_plot(prefix="individual", dlname="individual", width=800, function(env) {
        normed <- env$normed()
        individual_gene <- env$input$individual_gene
        df <- filter_(normed$norm_data, ~gene == individual_gene)
        df_bad <- filter_(df, ~is_low_tail_count)        
        
        plot.new()
        grid.arrange(
            ggplot(df, aes(x=sample, y=norm_count)) +
            geom_bar(stat="identity", color="#000000", fill="#cccccc") +
            ylab(if (env$input$normalizing_gene == "None") "Count" else "Normalized count") +
            xlab("") +
            title(paste(env$input$individual_gene,"expression")) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
            
            ggplot(df, aes(x=sample, y=tail)) +
            geom_bar(stat="identity", color="#000000", fill="#ccffcc") +
            geom_bar(data=df_bad, 
                     stat="identity", color="#000000", fill="#ff0000") +
            ylab("poly(A) tail length") +
            xlab("") +
            title(paste(env$input$individual_gene,"poly(A) tail"))+
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
            
            ncol=2, newpage=F
        )    
    })
    
    if (have_bams)
        tail_distribution <- shiny_plot(
            prefix="distribution", dlname="distribution", width=800, 
            function(env) withProgress(message="Plotting tail distribution", {
                this_samples <- env$input$samples
                query <- ranges[ env$input$individual_gene ]
                
                tail_lengths <- tibble(
                       sample=factor(this_samples, this_samples),
                       lengths=lapply(bam_filenames[this_samples], tail_lengths, query)
                    ) %>%
                    unnest(lengths)
                
                tail_lengths <- tail_lengths %>%
                    group_by(sample) %>%
                    mutate(n = n/sum(n)) %>%
                    ungroup()
                
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
                        ggplot2::geom_segment(aes(x=length,xend=length,y=cumn,yend=cumn_lag)) +
                        ggplot2::geom_segment(aes(x=length_lead,xend=length,y=cumn,yend=cumn)) +
                        scale_x_continuous(lim=c(0,max_tail), oob=function(a,b)a) +
                        scale_y_continuous(labels = scales::percent) +
                        labs(x="poly(A) tail length", y="Cumulative distribution", color="Sample") +
                        theme_bw()
                    )
                } else if (env$input$tail_style == "Density") {
                    print(
                        tail_lengths %>%
                        complete(sample=factor(this_samples,this_samples), length=seq(0,max_tail), fill=list(n=0)) %>%
                        ggplot(aes(color=sample,group=sample,x=length,y=n)) + 
                        geom_line() +
                        scale_x_continuous(lim=c(0,max_tail), oob=function(a,b)a) +
                        scale_y_continuous(labels = scales::percent) +
                        labs(x="poly(A) tail length", y="Percent reads", color="Sample") +
                        theme_bw()   
                    )
                } else {
                    print(
                        tail_lengths %>%
                        #group_by(sample) %>%
                        #mutate(n = n/max(n)) %>%
                        #ungroup() %>%
                        complete(sample=factor(this_samples,this_samples), length=seq(0,max_tail), fill=list(n=0)) %>%
                        ggplot(aes(x=sample,y=length,fill=n)) + 
                        geom_tile() +
                        scale_y_continuous(lim=c(0,max_tail), oob=function(a,b)a) +
                        scale_fill_viridis(guide=FALSE) +
                        labs(x="", y="poly(A) tail length") +
                        theme_minimal() +
                        theme(panel.grid=element_blank(), 
                              axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
                    )
                }
        }))


    # App

    ui <- div(
        titlePanel(title),
        navlistPanel(
            widths=c(2,10),
            well=FALSE,
            tabPanel("Select",
                selectInput("normalizing_gene","Normalizing gene",choices=possible_normalizers,selected=normalizing_gene),
                br(),
                
                selectInput("genes","Genes to display",choices=genes,selected=genes,multiple=TRUE,width="100%"),
                br(),
                
                selectInput("samples","Samples to display",choices=samples,selected=samples,multiple=TRUE,width="100%"),
                p(
                    actionButton("samples_all", "All")
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
                p("Note expression levels are *not* log transformed in this plot."),
                individual$component_ui,
                if (have_bams) h2("Tail length distribution"),
                if (have_bams) selectInput("tail_style", "Display", choices=c("Cumulative","Density","Heatmap"), selected="Cumulative"),
                if (have_bams) tail_distribution$component_ui
            )
        )
    )
    
    server <- function(env) {
        input <- env$input
        output <- env$output
        
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
             
             if (normalizing_gene == "None") {
                 normalizer <- tibble(sample=factor(samples,samples), normalizer=1)
             } else {
                 normalizer <- raw_data %>%
                     filter_(~ gene == normalizing_gene) %>%
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
                 baseline <- norm_counts %>%
                     select_(~sample %in% env$input$normalizing_samples) %>%
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
                 
             
             list(
                 norm_data=norm_data,
                 normalizer=normalizer, 
                 is_differential=is_differential
             )
        })
    
        output$normalizer_table <- renderTable({
             env$normed()$normalizer
        }, digits=6)
        
        overview$component_server(env)
        overview_tail$component_server(env)
        individual$component_server(env)
        if (have_bams) tail_distribution$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}



#' @export
shiny_repat <- shiny_mpat

