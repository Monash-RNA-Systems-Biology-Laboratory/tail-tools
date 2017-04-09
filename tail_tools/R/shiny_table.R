

#' Put most important columns first in a data frame, for display
prioritize_columns <- function(df, columns) {
    reordering <- order(match(colnames(df), columns))
    df[,reordering,drop=FALSE]
}


#' Shiny dataTable with download button
#'
shiny_table <- function(df=NULL, options=NULL, title="Table", filename="table.csv", prefix="") {
    ns <- shiny::NS(prefix)

    ui <- function(request)
        shiny::div(
            shiny::fluidRow(
                shiny::column(8,
                    shiny::h3(title)),
                shiny::column(4, 
                    style="text-align: right",
                    shiny::br(),
                    shiny::downloadButton(ns("download"), "CSV"))),
            DT::dataTableOutput(ns("table")))
    
    server <- function(env) {
        df <- ensure_reactive(df, ns("df"), env)
        options <- ensure_reactive(options, ns("options"), env, default=function() list())
        
        options_plus <- reactive({
            result <- list(
                pageLength=50, 
                selected=c())
            opts <- options()
            for(name in names(opts))
                result[[name]] <- opts[[name]]
            result
        })
        
        env$output[[ns("table")]] <- DT::renderDataTable(
            selection = "single",
            options = options,
            rownames = FALSE,
            { df() }
        )
    
        env$output[[ns("download")]] <- shiny::downloadHandler(
            filename=filename,
            content=function(file) {
                readr::write_csv(df(), file)
            }
        )
    }
    
    composable_shiny_app(ui, server, title=title)
}

