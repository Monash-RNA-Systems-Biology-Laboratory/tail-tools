
#' @export
shiny_test <- function(confects=NULL, prefix="") {
    ns <- shiny::NS(prefix)

    me_plot <- shiny_plot(prefix=ns("me_plot"))

    results_table <- shiny_table(
        title="Results", 
        filename="results.csv", 
        prefix=ns("results"))

    gene_view <- shiny_feature(is_peak=FALSE, prefix=ns("gene"))

    overview_panel <- function(request)
        shiny::tabPanel("Overview", 
            shiny::uiOutput(ns("description")),
            me_plot$component_ui(request),
            shiny::p("Grey dots show estimated effect size -- these estimates tend to be noisy when there are few reads. Black dots show confident effect size, a confident smallest bound on the effect size at the specified FDR."))

    results_panel <- function(request)
        shiny::tabPanel("Results",
            shiny::uiOutput(ns("result_description")),
            results_table$component_ui(request),
            gene_view$component_ui(request)
            )

    panels <- list(overview_panel, results_panel)

    server <- function(env) {
        confects <- ensure_reactive(confects, ns("confects"), env)

        title <- reactive({
            if (!is.null(confects()$title))
                confects()$title
            else
                "Differential test"
        })

        env[[ns("me_plot-callback")]] <- function() {
            topconfects::confects_plot_me(confects()) %>% 
            print
        }

        env$output[[ns("description")]] <- renderUI({
            shiny::div(
                shiny::h2( title() ),
                shiny::pre( topconfects:::confects_description(confects()) ))
        })

        env$output[[ns("result_description")]] <- renderUI({
            shiny::h2( title() )
        })

        env[[ns("results-df")]] <- reactive({ 
            confects()$table %>%
            select_(~-index) %>%
            prioritize_columns(
                c("rank", "confect", "effect", "logCPM", "AveExpr", 
                  "gene", "biotype", "name", "product"))
        })

        env[[ns("results-options")]] <- reactive({ 
            col_names <- colnames(env[[ns("results-df")]]())
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names
            
            if (is.null(confects()$limits)) {
                min_effect <- min(0, confects()$table$effect, na.rm=TRUE)
                max_effect <- max(0, confects()$table$effect, na.rm=TRUE)
                if (!is.finite(max_effect)) {
                    low <- -1.05
                    high <- 1.05
                } else if (min_effect >= 0) {
                    low <- 0.0
                    high <- max_effect*1.05
                } else {
                    high <- max(max_effect, -min_effect)*1.05
                    low <- -high
                }
            } else {
                low <- confects()$limits[1]
                high <- confects()$limits[2]
            }

            list(
                pageLength=20,

                columnDefs=list(
                    list(targets=cols$confect, visible=FALSE),
                    confect_coldef(cols$effect, cols$confect, low, high),
                    fixed_coldef(cols$logCPM,1))) 
        })

        env[[ns("gene-tc")]] <- reactive(withProgress(message="Loading", {
            read_tail_counts(paste0(confects()$pipeline_dir, "/expression/genewise/counts.csv")) %>%
            tail_counts_vst
        }))

        env[[ns("gene-peak_tc")]] <- reactive(withProgress(message="Loading", {
            read_tail_counts(paste0(confects()$pipeline_dir, "/expression/peakwise/counts.csv")) %>%
            tail_counts_vst
        }))

        env[[ns("gene-feature")]] <- reactive({
            rows_selected <- env$input[[ns("results-table_rows_selected")]]
            if (length(rows_selected) != 1)
                NULL
            else 
                confects()$table$name[rows_selected]
        })
        
        # Awfully hacky
        env[[ns("gene-peak_names")]] <- reactive({
            feature <- env[[ns("gene-feature")]]()
            if (is.null(feature))
                return(NULL)

            if (!is.null(confects()$display_members))
                return(confects()$display_members[[ feature ]])

            if (is.null(confects()$members))
                return(NULL)

            # Awfully hacky fallback
            rownames(confects()$edger_fit)[ confects()$members[[feature]] ]
        })

        me_plot$component_server(env)
        results_table$component_server(env)
        gene_view$component_server(env)
    }

    composable_shiny_panels_app(panels, server, prefix=prefix)
}


#' Present a collection of differential tests
#'
#' @param tests A named list of tests. Each item is a list("function_name", ...arguments...).
#'
#' @export
shiny_tests <- function(tests, cache_prefix="cache_", title="Differential tests", prefix="") {
    ns <- shiny::NS(prefix)

    titles <- names(tests)
    for(i in seq_along(tests))
        if (!is.null(tests[[i]]$title))
            titles[i] <- tests[[i]]$title

    version <- 8
    get <- function(name, func, fdr) {
        filename <- paste0(cache_prefix,name,"_",func,"_fdr",fdr,".rds")
        call <- c(tests[[name]][-1], list(fdr=fdr))

        if (file.exists(filename)) {
            cached_value <- readRDS(filename)
            if (identical(cached_value$func, func) &&
                identical(cached_value$call, call) &&
                identical(cached_value$version, version))
                return(cached_value$result)
        }

        result <- withProgress(
            message=paste0("Calculating ",filename), value=1,
            do.call(func, call))
        saveRDS(list(func=func, call=call, version=version, result=result), filename)
        result
    }


    subapp <- shiny_test(prefix=ns("test"))

    choices <- names(tests)
    names(choices) <- titles
    
    panels <- list(
        function(request) tabPanel("Select test",
            h2("Select test"),
            selectizeInput(ns("test"), "Test", choices=choices, selected=NULL, width="100%"),
            uiOutput(ns("test_variant_control")),
            
            numericInput(ns("fdr"), "False Discovery Rate", 0.05),
            div(style="height: 4em") #,
            #actionButton(ns("cache_all"), "Ensure all tests are cached") 
        )
    )
    
    panels <- c(panels, subapp$component_panels)
        
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        
        #observe({
        #    if (env$input[[ns("cache_all")]] == 0) return()
        #    
        #    for(name in names(tests))
        #        for(fdr in c(0.05, 0.01))
        #            get(name, fdr)
        #})
        
        env$output[[ns("test_variant_control")]] <- renderUI({
            variants <- test_variants[[ tests[[env$input[[ns("test")]]]][[1]] ]]
            
            selectizeInput(ns("variant"), "Test variant", choices=variants, width="100%")
        })
    
        env[[ns("test-confects")]] <- reactive({
            get(env$input[[ns("test")]], env$input[[ns("variant")]], env$input[[ns("fdr")]])
        })
    
        subapp$component_server(env)
    }
    
    composable_shiny_panels_app(panels, server, title=title)
}