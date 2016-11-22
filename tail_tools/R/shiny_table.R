
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
                    shiny::downloadButton(ns("download"), "Download CSV file"))),
            DT::dataTableOutput(ns("table")))
    
    server <- function(env) {
        df <- ensure_reactive(df, ns("df"), env)
        options <- ensure_reactive(options, ns("options"), env, default=list(pageLength=50))
        
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

