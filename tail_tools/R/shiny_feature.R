

shiny_feature <- function(tc, feature=NULL, is_peak=FALSE, peak_tc=NULL, prefix="") {
    ns <- NS(prefix)
    
    table <- shiny_table(title="Data", prefix=ns("table"))
    
    if (!is_peak)
       peak_table <- shiny_table(title="Peaks", prefix=ns("peak_table"))
    
    ui <- function(request)
        shiny::div(
            shiny::uiOutput(ns("ui")),
            if (!is_peak) peak_table$component_ui(request),
            table$component_ui(request),
            shiny::p("The norm_count column contains counts adjusted to equalize the (TMM adjusted) library size. Low values of the log2_norm_count column are moderated upward using a Variance Stabilizing Transformation. See the \"Transform\" page."))
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        
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
            df <- peak_tc$features[hits,c("feature","relation","parent","antisense_parent")]
            clear <- df$parent == feature() & df$antisense_parent == feature()
            df$parent[clear] <- ""
            df$relation[clear] <- ""
            dplyr::rename_(df, peak="feature", gene="parent", antisense_gene="antisense_parent") 
        })
        
        env[[ns("table-df")]] <- reactive({
            subtc()$obs
        })
        
        env[[ns("table-options")]] <- reactive({
            col_names <- colnames(e("table-df"))
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names
            
            list(
                pageLength=20,
                columnDefs=list(
                    fixed_coldef(cols$tail, 1),
                    fixed_coldef(cols$norm_count, 1),
                    fixed_coldef(cols$log2_norm_count, 1)))
        })
        
        if (!is_peak) peak_table$component_server(env)
        table$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}

