

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
    
    ui <- function(request)
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
                shiny::uiOutput(ns("description")),
                enrich_table$component_ui(request),
                shiny::p("p values and FDR are from Fisher's Exact Test. They may be overly liberal as this test does not account for inter-gene correlation.")))
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()

        fs <- ensure_reactive(fs, ns("fs"), env)
    
        env[[ns("table-df")]] <- reactive({
            dplyr::left_join(
                dplyr::data_frame(feature=fs()$set),
                tc$features,
                "feature") %>%
            dplyr::select_(~-one_of(c("Length","product","antisense_product"))) %>%
            dplyr::rename_(chr="chromosome")
        })
        
        env[[ns("table-options")]] <- reactive({
            col_names <- colnames(e("table-df"))
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names
            
            list(
                pageLength=20,
                columnDefs=list(
                    fixed_coldef(cols$"mean-tail", 1),
                    fixed_coldef(cols$"proportion-with-tail", 2)))
        })
        
        env[[ns("feature-selected")]] <- reactive({
            table_df <- isolate(e("table-df"))
            rows_selected <- env$input[[ns("table-table_rows_selected")]]
            
            #if (nrow(table_df) == 1)
            #    table_df$feature
            #else
            table_df$feature[rows_selected]
        })
        
        if (have_species) geneset <- reactive({
            geneset <- fs()
            if (is_peaks) {
                geneset$set <- peaks_to_genes(tc, geneset$set)
                geneset$universe <- peaks_to_genes(tc, geneset$universe)
            }
            geneset            
        })
        
        if (have_species) enrich_result <- reactive(withProgress(
            message="Testing gene sets", {
            if (length(geneset()$set) == 0 || all(geneset()$universe %in% geneset()$set))
                return(NULL)
            
            featureset_test_genesets(geneset(), species, 
                minimum_set_size=env$input[[ns("min_set")]],
                ontologies = c(
                    if (env$input[[ns("BP")]]) "BP",
                    if (env$input[[ns("MF")]]) "MF",
                    if (env$input[[ns("CC")]]) "CC",
                    if (env$input[[ns("KEGG")]]) "KEGG"))
        }))

        if (have_species) env$output[[ns("description")]] <- renderUI({
            shiny::div(
                shiny::p( 
                    paste0(
                        length(geneset()$set), " genes from a universe of ", 
                        length(geneset()$universe), "."),
                    if (!is.null(enrich_result())) shiny::p(
                        "Present in the gene-sets tested: ",
                        length(enrich_result()$set), " genes from a universe of ",
                        length(enrich_result()$universe), ".")))
        })
        
        
        if (have_species) env[[ns("enrich-df")]] <- reactive( enrich_result()$table )
        
        if (have_species) env[[ns("enrich-options")]] <- reactive({
            col_names <- colnames(e("enrich-df"))
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names
            
            list(
                pageLength=20,
                columnDefs=list(
                    precision_coldef(cols$FDR, 3),
                    precision_coldef(cols$p_value, 3)))
        })
        
        table$component_server(env)
        if (have_species) 
            enrich_table$component_server(env)
    }

    composable_shiny_app(ui, server)
}