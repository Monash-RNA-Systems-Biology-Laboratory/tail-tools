
#
# Shiny modules associated with the RNA-Seq end shift test 
#

get_rnaseq_depth <- function(samples) {
    mean_depth <- function(filename) {
        result <- 
            rtracklayer::summary(BigWigFile(filename), type="mean") %>% 
            BiocGenerics::unlist()
        sum(result$score * width(result)) / sum(as.numeric(width(result)))
    }

    samples %>%
    rowwise() %>%
    mutate(
        mean_cover = 0.5*(mean_depth(cover_fwd)+mean_depth(cover_rev))
    ) %>%
    ungroup() %>%
    mutate(
        depth_normalizer = mean_cover / mean(mean_cover)
    )
}


#' @export
shiny_end_shift_rnaseq <- function(result, samples, features, prefix="") {
    result <- ensure_reactable(result)
    samples <- ensure_reactable(samples)
    features <- ensure_reactable(features)
    
    composable_shiny_panels_app(
        shiny_end_shift_rnaseq_panels(prefix),
        
        function(env) {
            callModule(shiny_end_shift_rnaseq_server, prefix, 
                reactive( result(env) ), 
                reactive( samples(env) ), 
                reactive( features(env) ),
                outer_session=env$session,
                session=env$session)
        }
    )
}

#' @export
shiny_end_shift_rnaseq_panels <- function(id) {
    ns <- NS(id)
    
    list(
        tabPanel("Overview",
            shiny_plot_ui(ns("mr_plot"), width=800, 
                brush=brushOpts(id=ns("mr_plot_brush"), delay=600000)),
            uiOutput(ns("blurb"))
        ),
    
        tabPanel("Results",
            DT::dataTableOutput(ns("table")),
            downloadButton(ns("table_download"), "Download CSV file")
        ),
    
        tabPanel("Gene viewer",
            shiny_genome_browser_ui(ns("browser"), "")        
        )
    )
}

#' @export
shiny_end_shift_rnaseq_server <- function(input, output, session, result, samples, features, outer_session=NULL) {
    result <- ensure_reactive(result)
    samples <- ensure_reactive(samples)
    features <- ensure_reactive(features)

    col_selection <- c(
            "rank", "r", "r_low", "r_high", "transcript_name", "gene_name", "description",
            "transcript_id", "gene_id", "chromosome", "strand", 
            "five_prime", "three_prime_given", "three_prime_extended",
            "mean_reads", "min_reads"
    )
    cols <- as.list(seq_along(col_selection)-1)
    names(cols) <- col_selection

    df <- reactive({
        result()$gene[,col_selection,drop=FALSE]
    })
    
    callModule(shiny_plot_server, "mr_plot", session=session, callback=function() {
        ggplot(df(),aes(x=min_reads, y=r)) +
            scale_x_log10() +
            coord_cartesian(ylim=c(-1,1)) +
            ggplot2::geom_segment(aes(xend=min_reads, y=r_low, yend=r_high), color="#bbbbbb") +
            geom_point() +
            theme_bw() +
            labs(x = "Min reads per sample (log scale)",
                 y = "r")        
    })
    
    output$blurb <- renderUI({
       tags$div(
           tags$p("Transcripts with a sample having less than", sprintf("%d",result()$min_min_reads), "reads were filtered.") 
       )
    })

    selected_df <- reactive({
        result <- df()
        if (!is.null(input$mr_plot_brush)) {
            result <- result %>%
                filter(!is.na(r)) %>%
                brushedPoints(input$mr_plot_brush, "min_reads", "r")
        }
        
        result
    })
    
    observeEvent(input$mr_plot_brush, {
        if (!is.null(input$mr_plot_brush))
            updateTabsetPanel(outer_session, "tabset", selected="Results")
    })

    output$table <- DT::renderDataTable(
        selection=list(mode="single",selection=integer(0)),
        rownames=F, 
        options=list(
            pageLength=25,
            columnDefs=list(
                list(
                    targets=c(cols$r_low, cols$r_high),
                    visible=FALSE
                ),
                r_coldef(cols$r, cols$r_low, cols$r_high)
            )
        ),
        selected_df()
    )

    output$table_download <- downloadHandler(
        filename="end-shift.csv",
        content=function(file) {
            write_csv(selected_df(), file)
        }
    )

    context <- reactive({
        row <- input$table_rows_selected
        if (length(row) != 1 || row > nrow(selected_df()))
            return()
        
        selected_df()[row,]
    })
    
    observeEvent(input$table_rows_selected, {
        if (!is.null(outer_session) && length(input$table_rows_selected) > 0) {
            updateTabsetPanel(outer_session, "tabset", selected="Gene viewer")
        }
    })
    
    observe({
        con <- context()
        if (is.null(con)) return()
    
        a <- con$five_prime
        b <- con$three_prime_extended
        w <- abs(a-b) %/% 2
        pos <- GRanges( 
             con$chromosome,
             IRanges(min(a,b)-w,max(a,b)+w),
             strand=con$strand)       
         
        updateTextInput(session, "browser-location", value=as.character(pos))
    })
    
    plot_it <- function(pos) {
        con <- isolate( context() )
        if (is.null(con))
            p <- plot_genome(features(), samples(), pos, c())
        else
            p <- plot_genome(features(), samples(), pos, c( con$five_prime, con$three_prime_extended ))
        
        print(p)
    }
    
    callModule(shiny_genome_browser_server, "browser", plot_it)
}







#' Interface to multiple RNA-Seq end-shift tests.
#'
#' @param samples A data frame of samples. Should have columns: name = name of sample, bigwig files: cover_fwd, cover_rev, span_fwd, span_rev.
#'
#' @param conditions A list of condition vectors.
#'
#' @param references A named character vector of reference directories created with "tail-tools make-rnaseq-reference:".
#'
#' @export
shiny_end_shift_rnaseq_multiple <- function(
        samples, tests, references, cache_prefix="cache", prefix="", title="RNA-Seq end-shift test") {

    composable_shiny_panels_app(
        shiny_end_shift_rnaseq_multiple_panels(prefix, samples, tests, references),
        
        function(env) {
            callModule(shiny_end_shift_rnaseq_multiple_server, prefix, 
                samples, tests, references, cache_prefix, 
                outer_session=env$session, session=env$session)
        },
        
        title=title
    )    
}

#' @export
shiny_end_shift_rnaseq_multiple_panels <- function(id, samples, tests, references) {
    ns <- NS(id)
    
    c(
        list(
            tabPanel("Select test",
                selectizeInput(ns("test"), "Test", choices=names(tests), 
                    selected=names(tests)[1], width="100%"),
                selectizeInput(ns("reference"), "Features to use", choices=names(references),
                    selected=names(references)[1]),
                selectizeInput(ns("samples"), "Samples to view", choices=samples$name,
                    selected=samples$name, multiple=TRUE),
                actionButton(ns("cache_all"), "Ensure all tests are cached")                     
            )
        ),
        shiny_end_shift_rnaseq_panels(ns("end_shift_rnaseq"))
    )
}

#' @export
shiny_end_shift_rnaseq_multiple_server <- function(
        input,output,session,
        samples, tests, references, cache_prefix="cache",
        outer_session=NULL) {

    samples <- cached(paste0(cache_prefix,"-samples"), 
        get_rnaseq_depth, list(samples=samples))

    get <- function(test, reference) {
        if (is.null(test) || is.null(reference))
            return(NULL)
    
        name <- paste0(cache_prefix,"-test-",test,"--",reference)
        
        this_samples <- samples
        this_samples$condition <- tests[[test]]$condition
        this_samples$group <- tests[[test]]$group
        this_samples <- dplyr::filter(this_samples, !is.na(condition))
    
        reference_dir <- references[reference]
        utrs <- import(file.path(reference_dir,"utr.gff"))
        extended_utrs <- import(file.path(reference_dir,"utr_extended.gff"))
        exons <- import(file.path(reference_dir,"utr_part.gff"))
        
        withProgress(message=paste0(test," ",reference),
            cached(name, end_shift_rnaseq,
                list(samples=this_samples),
                list(utrs=utrs, extended_utrs=extended_utrs, exons=exons),
                version=list(reference_dir, 13)
            )
        )
    }
    
    observeEvent(input$cache_all, {
        for(reference in names(references))
            for(test in names(tests))
                get(test, reference)
    })
    
    result <- reactive({
        get(input$test,input$reference)
    })
    
    selected_samples <- reactive({
        samples[ match(input$samples, samples$name),,drop=F ]
    })
    
    features <- reactive({
        if (is.null(input$reference))
            return(NULL)
        load_genome_features(references[input$reference])
    })
    
    callModule(shiny_end_shift_rnaseq_server, "end_shift_rnaseq", 
        result, selected_samples, features,
        outer_session=outer_session, session=session)
}














