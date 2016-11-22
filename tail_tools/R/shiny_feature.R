

shiny_feature <- function(tc, feature=NULL, is_peak=FALSE, peak_tc=NULL, prefix="") {
    ns <- NS(prefix)
    
    table <- shiny_table(title="Data", prefix=ns("table"))
    
    if (!is_peak)
       peak_table <- shiny_table(title="Peaks", prefix=ns("peak_table"))
    
    ui <- function(request)
        shiny::div(
            shiny::uiOutput(ns("ui")),
            if (!is_peak) peak_table$component_ui(request),
            table$component_ui(request))
    
    server <- function(env) {
        feature <- ensure_reactive(feature, ns("feature"), env)
        
        subtc <- reactive({ 
            tail_counts_subset_features(tc, feature()) 
        })
    
        env$output[[ns("ui")]] <- renderUI({
            if (is.null(feature()))
                return(shiny::div())
            
            shiny::div(
                shiny::h1( feature(), subtc()$features$gene ),
                shiny::p(subtc()$features$product))
        })
        
        if (!is_peak) env[[ns("peak_table-df")]] <- reactive({
            hits <- peak_tc$features$parent == feature() | peak_tc$features$antisense_parent == feature()
            peak_tc$features[hits,c("feature","relation","parent","antisense_parent")]
        })
        
        env[[ns("table-df")]] <- reactive({
            subtc()$obs
        })
        
        if (!is_peak) peak_table$component_server(env)
        table$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}

