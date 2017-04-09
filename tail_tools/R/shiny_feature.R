

shiny_feature <- function(tc=NULL, feature=NULL, is_peak=FALSE, peak_tc=NULL, peak_names=NULL, prefix="") {
    ns <- NS(prefix)
    
    table <- shiny_table(title="Data", prefix=ns("table"))
    
    if (!is_peak) {
       sample_peak_table <- shiny_table(title="Samples and peaks", prefix=ns("sample_peak_table"))
       peak_table <- shiny_table(title="Peaks", prefix=ns("peak_table"))
    }

    ui <- function(request)
        shiny::div(
            shiny::uiOutput(ns("ui")),
            if (!is_peak) sample_peak_table$component_ui(request),
            if (!is_peak) shiny::p("These are read counts adjusted to equalize the (TMM adjusted) library size."),
            if (!is_peak) peak_table$component_ui(request),
            table$component_ui(request),
            shiny::p("The norm_count column contains read counts adjusted to equalize the (TMM adjusted) library size. Low values of the log2_norm_count column are moderated upward using a Variance Stabilizing Transformation."))
    
    server <- function(env) {
        e <- function(name) env[[ns(name)]]()
        
        tc <- ensure_reactive(tc, ns("tc"), env)
        peak_tc <- ensure_reactive(peak_tc, ns("peak_tc"), env)
        peak_names <- ensure_reactive(peak_names, ns("peak_names"), env, default=function() NULL)
        feature <- ensure_reactive(feature, ns("feature"), env)
        
        subtc <- reactive({ 
            tail_counts_subset_features(tc(), feature()) 
        })
    
        env$output[[ns("ui")]] <- renderUI({
            if (is.null(feature()))
                return(shiny::div())
            
            shiny::div(
                shiny::h1( feature(), subtc()$features$gene ),
                shiny::p(subtc()$features$product))
        })

        if (!is_peak) {
            relevant_peaks <- reactive({
                # Basically, just shiny_test overrides this
                fixed_peak_names <- peak_names()

                # Default if no peaks given
                if (is.null(fixed_peak_names)) {
                    if (is.null(feature()))
                        return(NULL)

                    hits <- peak_tc()$features$parent == feature() | peak_tc()$features$antisense_parent == feature()

                    if (!any(hits))
                        return(NULL)
                
                    fixed_peak_names <- peak_tc()$features$feature[hits]
                }

                if (length(fixed_peak_names) == 0)
                    return(NULL)

                result <- tail_counts_subset_features(peak_tc(), fixed_peak_names)
                ordering <- order(ifelse(result$features$strand >= 0,result$features$end,-result$features$start))
                tail_counts_subset_features(result, result$features$feature[ordering])
            })

            env[[ns("sample_peak_table-df")]] <- reactive({
                if (is.null(relevant_peaks()))
                    return(NULL)

                relevant_peaks()$obs %>%
                select_(~sample, ~feature, ~norm_count) %>%
                tidyr::spread("feature", "norm_count") 
                #tidyr::gather_("stat","value", c("norm_count","tail_count","tail")) %>%
                #mutate_(column =~ paste0(stat,"__",feature)) %>%
                #select_(~sample, ~column, ~value) %>%
                #tidyr::spread("column", "value")
            })

            env[[ns("sample_peak_table-options")]] <- reactive({
                if (is.null(relevant_peaks()))
                    return(list())

                col_names <- colnames(e("sample_peak_table-df"))

                maximum <- max(relevant_peaks()$obs$norm_count)
                norm_count_cols <- seq_len(length(col_names)-1)

                list(
                    dom="t", 
                    paging=FALSE,
                    columnDefs=list(
                        n_coldef(norm_count_cols, maximum, 1)
                        ))
            })
        
            env[[ns("peak_table-df")]] <- reactive({
                if (is.null(relevant_peaks()))
                    return(NULL)

                df <- relevant_peaks()$features[,c("feature","relation","parent","antisense_parent")]
                clear <- df$parent == feature() & df$antisense_parent == feature()
                df$parent[clear] <- ""
                df$relation[clear] <- ""
                dplyr::rename_(df, peak="feature", gene="parent", antisense_gene="antisense_parent") 
            })

            env[[ns("peak_table-options")]] <- reactive({
                list(
                    dom="t", 
                    paging=FALSE)
            })
        }


        env[[ns("table-df")]] <- reactive({
            subtc()$obs %>%
            select(-feature)
        })
        
        env[[ns("table-options")]] <- reactive({
            col_names <- colnames(e("table-df"))
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names

            maximum <- max(e("table-df")$norm_count)
            
            list(
                dom="t",
                paging=FALSE,
                columnDefs=list(
                    fixed_coldef(cols$tail, 1),
                    n_coldef(cols$norm_count, maximum, 1),
                    fixed_coldef(cols$log2_norm_count, 1)))
        })
        
        if (!is_peak) sample_peak_table$component_server(env)
        if (!is_peak) peak_table$component_server(env)
        table$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}

