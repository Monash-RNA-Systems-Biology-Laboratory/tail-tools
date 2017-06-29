

shiny_counts_report <- function(tc, pipeline_dir=NULL, species=NULL, title="Tail Tools count report", is_peaks=FALSE, peak_tc=NULL, prefix="") {
    ns <- NS(prefix)
    
    what <- if (is_peaks) "Peak" else "Gene"
    
    feature_labels <- tc$features$gene
    if (is_peaks)
        feature_labels <- paste(tc$features$feature, feature_labels)

    varistran_app <- shiny_report(
        y=tail_counts_get_vst(tc),
        counts=tail_counts_get_matrix(tc, "count"),
        feature_labels=feature_labels,
        prefix=ns("varistran")
        )

    # Replace me    
    dat <- list()
    dat$Count <- as.data.frame(tail_counts_get_matrix(tc, "count"))
    dat$Tail <- as.data.frame(tail_counts_get_matrix(tc, "tail"))
    dat$Tail_count <- as.data.frame(tail_counts_get_matrix(tc, "tail_count"))
    dat$Annotation <- as.data.frame(tc$features)
    rownames(dat$Annotation) <- dat$Annotation$feature
    dat$Annotation$feature <- NULL
    heatmap_app <- shiny_patseq_heatmap(dat, prefix=ns("heatmap"))
    
    #featureset_get <- function(env) {
    #    set <- env[[ns("varistran-rows-selected")]]()
    #    universe <- env[[ns("varistran-rows")]]()
    #    if (length(set) == 0)
    #        set <- universe
    #    list(
    #        set = tc$features$feature[set],
    #        universe = tc$features$feature[universe]
    #    )
    #}
    
    featureset <- shiny_featureset(tc, species=species, is_peaks=is_peaks, prefix=ns("featureset"))
        
    #feature_get <- function(env) {
    #    featureset_get(env)[ env$input[[ns("featureset-table_rows_selected")]] ]
    #}
    
    feature <- shiny_feature(tc, is_peak=is_peaks, peak_tc=peak_tc, prefix=ns("feature"))
    
    tail_distribution <- shiny_tail_distribution(
        pipeline_dir = pipeline_dir,
        peaks = if (is_peaks) tc$features else peak_tc$features,
        prefix=ns("tail_distribution"))
    
    panels <- c(
        list("Expression"),
        varistran_app$component_panels,
        list("Tail length"),
        heatmap_app$component_panels[c(1,2)],
        list(
            paste0(what," details"),
            function(request) 
                shiny::tabPanel(paste(what,"set"), 
                    shiny::radioButtons(ns("setsource"), "Source of feature set",
                        c("All features" = "all",
                          "Expression heatmap selection" = "heatmap",
                          "Tail length heatmap selection" = "tailmap"),
                        "all"),                           
                    featureset$component_ui(request) ),
            function(request)
                shiny::tabPanel(what, 
                    shiny::textInput(ns("feature_name"), what, ""),
                    feature$component_ui(request))),
            function(request)
                shiny::tabPanel("- poly(A) tail distribution",
                    tail_distribution$component_ui(request)))
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        i <- function(name) env$input[[ns(name)]]
        
        env[[ns("featureset-fs")]] <- reactive({
            if (env$input[[ns("setsource")]] == "all") {
                universe <- tc$features$feature
                set <- universe
            } else if (env$input[[ns("setsource")]] == "heatmap") {
                set <- tc$features$feature[ e("varistran-rows-selected") ]
                universe <- tc$features$feature[ e("varistran-rows") ]
            } else if (env$input[[ns("setsource")]] == "tailmap") {
                set <- e("heatmap-plot-rows-selected")
                universe <- e("heatmap-rows")
            } else {
                set <- c()
                universe <- c()
            }
            
            list(
                set = set,
                universe = universe
            )
        })
        
        #env[[ns("feature-feature")]] <- reactive({
        #    env[[ns("featureset-feature-selected")]]()
        #})
        
        env[[ns("feature-feature")]] <- reactive({
            if (i("feature_name") %in% tc$features$feature)
                i("feature_name")
            else
                NULL
        })
        
        env[[ns("tail_distribution-peak_names")]] <- reactive({
            feature <- e("feature-feature")
            if (is_peaks)
                return(feature)
            else
                return(peak_tc$features$feature[
                    peak_tc$features$parent %in% feature &
                    peak_tc$features$relation != "Antisense"])
        })
        
        env[[ns("tail_distribution-sample_normalizer")]] <- reactive({
            select_(tc$samples, ~sample, ~normalizer)
        })
        
        
        # Sub-components

        varistran_app$component_server(env)
        heatmap_app$component_server(env)
        featureset$component_server(env)
        feature$component_server(env)
        tail_distribution$component_server(env)
        

        # Non-reactive event handling
        
        observeEvent({ e("varistran-rows-selected") }, {
            if (!length(e("varistran-rows-selected"))) return()
            
            updateRadioButtons(env$session, ns("setsource"), selected="heatmap")
            updateTabsetPanel(env$session, ns("tabset"), selected=paste(what,"set"))
        })
        
        observeEvent({ e("heatmap-plot-rows-selected") }, {
            if (!length(e("heatmap-plot-rows-selected"))) return()
            
            updateRadioButtons(env$session, ns("setsource"), selected="tailmap")
            updateTabsetPanel(env$session, ns("tabset"), selected=paste(what,"set"))
        })
        
        
        # Normalize case, resolve nomeclature gene names
        observeEvent({ i("feature_name") }, {
            if (!is.null(e("feature-feature"))) return()
            
            hits <- tolower(tc$features$feature) == tolower(i("feature_name"))
            if (!is_peaks && sum(hits) == 0)
                hits <- tolower(tc$features$gene) == tolower(i("feature_name"))
            
            if (sum(hits) == 1) 
                updateTextInput(env$session, ns("feature_name"), value=tc$features$feature[hits])
        })
        
        # Go to feature when selected in featureset
        observeEvent({ e("featureset-feature-selected") }, {
            selected <- e("featureset-feature-selected")
            if (length(selected) == 1) {
                updateTextInput(env$session, ns("feature_name"), value=selected)
                updateTabsetPanel(env$session, ns("tabset"), selected=what)
            }
        })
    }
    
    composable_shiny_panels_app(panels, server, prefix=prefix, title=title)
}


#'
#' Shiny report based on pipeline output
#'
#' @param path Directory containing pipeline output.
#'
#' @param species Species for tail heatmap. Currently supports Human ("Hs"), Saccharomyces cerevisiae ("Sc"), Caenorhabditis elegans ("Ce"), Mus musculus ("Mm")
#'
#' @param title Title for report.
#'
#' @export
shiny_tailtools_report <- function(path, species=NULL, title="Tail Tools report", prefix="") {
    ns <- NS(prefix)

    genewise_filename <- paste0(path,"/expression/genewise/counts.csv")
    genewise_tc <- 
        read_tail_counts(genewise_filename) %>%
        tail_counts_vst()

    peakwise_filename <- paste0(path,"/expression/peakwise/counts.csv")
    peakwise_tc <- 
        read_tail_counts(peakwise_filename) %>%
        tail_counts_vst()

    genewise_report <- shiny_counts_report(genewise_tc, pipeline_dir=path, species=species, prefix=ns("genewise"), is_peaks=FALSE, peak_tc=peakwise_tc, title=NULL)
    genewise_panel <- function(request) tabPanel("Genes", genewise_report$component_ui(request))

    peakwise_report <- shiny_counts_report(peakwise_tc, pipeline_dir=path, species=species, prefix=ns("peakwise"), is_peaks=TRUE, title=NULL)
    peakwise_panel <- function(request) tabPanel("Peaks", peakwise_report$component_ui(request))
    
    panels <- list(
        genewise_panel,
        peakwise_panel
    )
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        i <- function(name) env$input[[ns(name)]]

        genewise_report$component_server(env)
        peakwise_report$component_server(env)
                          
        observeEvent({ i("genewise-feature-peak_table-table_rows_selected") }, {
            rows <- i("genewise-feature-peak_table-table_rows_selected")
            cat("rows");print(rows)
            selected <- e("genewise-feature-peak_table-df")[rows,1]
            cat("selected");print(selected)
            if (length(selected) == 1) {
                updateTextInput(env$session, ns("peakwise-feature_name"), value=selected)
                updateTabsetPanel(env$session, ns("peakwise-tabset"), selected="Peak")
                updateTabsetPanel(env$session, ns("tabset"), selected="Peaks")
            }
        })
        
        #observe({
        #    print(names(env$input))
        #})
    }
    
    composable_shiny_panels_app(panels, server, prefix=prefix, title=title, top=TRUE)
}

