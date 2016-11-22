

peaks_to_genes <- function(tc, peaks) {
    genes <- tc$features$parent[ tc$features$feature %in% peaks & tc$features$relation != "Antisense" ]
    genes <- genes[ genes != "" ]
    unique(genes)
}


#' Display information about a subset of rows in a counts file
#'
shiny_featureset <- function(tc, fs=NULL, species=NULL, is_peaks=FALSE, prefix="") {
    ns <- NS(prefix)
    
    table <- shiny_table(
        title="Features",
        filename="featureset.csv", 
        prefix=ns("table"))
    
    have_species <- !is.null(species)
    if (have_species)
        enrich_table <- shiny_table(
            title="Enriched gene-sets",
            filename="enrichment.csv",
            prefix=ns("enrich"))
    
    ui <- function(requst)
        shiny::div(
            table$component_ui(request),            
            if (have_species) shiny::div(
                shiny::h2("Gene-set enrichment"),
                shiny::fluidRow(
                    shiny::column(4,
                        shiny::checkboxInput(ns("BP"), "GO Biological Pathways", TRUE),
                        shiny::checkboxInput(ns("MF"), "GO Molecular Functions", TRUE),
                        shiny::checkboxInput(ns("CC"), "GO Cellular Components", TRUE),
                        shiny::checkboxInput(ns("KEGG"), "KEGG Pathways", TRUE)),
                    shiny::column(8,
                        shiny::numericInput(ns("min_set"), "Minimum set size", 2))),
                enrich_table$component_ui(request)))
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()

        fs <- ensure_reactive(fs, ns("fs"), env)
    
        env[[ns("table-df")]] <- reactive({
            dplyr::left_join(
                dplyr::data_frame(feature=fs()$set),
                tc$features,
                "feature") %>%
            dplyr::select_(~-Length, ~-product)
        })
        
        env[[ns("feature-selected")]] <- reactive({
            if (nrow(e("table-df")) == 1)
                e("table-df")$feature
            else
                e("table-df")$feature[ env$input[[ns("table-table_rows_selected")]] ]
        })
        
        
        if (have_species) env[[ns("enrich-df")]] <- reactive(withProgress(
            message="Testing gene sets", {
            if (length(fs()$set) == 0 || all(fs()$universe %in% fs()$set))
                return(NULL)
            
            geneset <- fs()
            if (is_peaks) {
                geneset$set <- peaks_to_genes(tc, geneset$set)
                geneset$universe <- peaks_to_genes(tc, geneset$universe)
            }
            
            featureset_test_genesets(geneset, species, 
                minimum_set_size=env$input[[ns("min_set")]],
                ontologies = c(
                    if (env$input[[ns("BP")]]) "BP",
                    if (env$input[[ns("MF")]]) "MF",
                    if (env$input[[ns("CC")]]) "CC",
                    if (env$input[[ns("KEGG")]]) "KEGG"))
        }))
        
        if (have_species) env[[ns("enrich-options")]] <- reactive({
            col_names <- colnames(e("enrich-df"))
            cols <- seq_len(col_names)
            colnames(cols) <- col_names
            
            list(
                pageLength=20,
                columnDefs=list(
                    
                )
            )
        })
        
        table$component_server(env)
        if (have_species) 
            enrich_table$component_server(env)
    }

    composable_shiny_app(ui, server)
}