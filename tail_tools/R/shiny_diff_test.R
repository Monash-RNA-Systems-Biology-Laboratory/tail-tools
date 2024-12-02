

# Confident sum of squares of x (with measurement variances v)
# 2024: I think this assumes v is known, which is not in general true, so I've disabled use of this
ss_confidence_bound <- function(x,v,confidence=0.95) {
    good <- !is.na(x) & is.finite(v)
    x <- x[good]
    v <- v[good]

    # We want to explain excess mean square, given variances
    ms <- sum(x^2) / sum(v)

    # Welch-Satterthwaite degrees of freedom approximation
    df <- sum(v)^2 / sum(v^2)

    score <- function(ncp)
        pchisq(ms*df,df,ncp=ncp*df) - confidence

    if (score(0) <= 0) return(0)

    uniroot(score, c(0,ms), extendInt="downX")$root * sum(v)
}


#' @export
shiny_test <- function(confects=NULL, prefix="") {
    ns <- shiny::NS(prefix)
    
    # confect and effect shown separately on same plot (turned out to be confusing)
    me_plot <- shiny_plot(prefix=ns("me_plot"))
    
    # More conventional mean-difference style plot, color by confect
    me2_plot <- shiny_plot(prefix=ns("me2_plot"), width=700)
    
    results_table <- shiny_table(
        title="Results", 
        filename="results.csv", 
        prefix=ns("results"))
    
    gene_view <- shiny_feature(is_peak=FALSE, prefix=ns("gene"))
    
    overview_panel <- function(request)
        shiny::tabPanel("Overview", 
            shiny::uiOutput(ns("description")),
            shiny::fluidRow(
                shiny::column(2, shiny::numericInput(ns("ymin"), "y-axis min", value=NA)),
                shiny::column(2, shiny::numericInput(ns("ymax"), "y-axis max", value=NA)),
                shiny::column(2, shiny::textInput(ns("ylabel"), "y-axis label", value="")),
                shiny::column(2, shiny::textInput(ns("breaks"), "Confect breaks", value=""))),
            me2_plot$component_ui(request),
            shiny::p("This plot shows \"effect\" values on the y axis. \"Confect\" values are indicated by color."),
            shiny::p("Confect breaks can be specified as a space separated list of numbers, eg: 0 0.5 1 1.5"),
            shiny::h2("Old style plot"),
            me_plot$component_ui(request),
            shiny::p("Grey dots show estimated effect size -- these estimates tend to be noisy when there are few reads. Black dots show confident effect size, a confident smallest bound on the effect size at the specified FDR. Each black dot has a corresponding grey dot. Grey dots without a corresponding black dot are not significantly different from zero."))

    results_panel <- function(request)
        shiny::tabPanel("Results",
            shiny::uiOutput(ns("result_description")),
            results_table$component_ui(request),
            gene_view$component_ui(request)
            )

    enrichment_panel <- function(request)
        shiny::tabPanel("Enrichment lists",
            shiny::selectInput(ns("enrichment_ranking"),
                label="Put at the top of the list",
                selected=c("abs"),
                choices=c(
                    "top confects"="abs",
                    "top confects positive"="up",
                    "top confects negative"="down",
                    "most \"significant\""="fdr_zero",
                    "most \"significant\" positive"="fdr_zero_up",
                    "most \"significant\" negative"="fdr_zero_down")),
            shiny::p(
                "These lists can be pasted into web-based enrichment tools, for example", 
                shiny::a("gProfiler",href="https://biit.cs.ut.ee/gprofiler/gost", target="_blank")),
            shiny::p(
                "Gene expression and poly(A)-tail length tests, signed or unsigned ranking: recommend using the all-genes list as an ordered query and also as the background set of genes."),
            shiny::p(
                "Shift tests, unsigned ranking: recommend using the all-genes list as an ordered query, and not providing a background set of genes."),
            shiny::p(
                "Shift tests, signed ranking: recommend using the significant genes list as an ordered query, and not providing a background set of genes."),
            shiny::p(
                "Note that confect ranking falls back to ranking by p-value for non-significant results."),
            shiny::uiOutput(ns("enrichment_out")))

    diagnostic_plot <- shiny_plot(prefix=ns("diagnostic_plot"))

    diagnostics_panel <- function(request)
        shiny::tabPanel("Diagnostics",
            shiny::uiOutput(ns("diagnostic_ui")),
            shiny::uiOutput(ns("diagnostic_text")),
            diagnostic_plot$component_ui(request))

    panels <- list(overview_panel, results_panel, enrichment_panel, diagnostics_panel)

    counts_reader <- memoise::memoise(function(path) tail_counts_vst(read_tail_counts(path)))

    server <- function(env) {
        confects <- ensure_reactive(confects, ns("confects"), env)

        title <- reactive({
            if (!is.null(confects()$title))
                confects()$title
            else
                "Differential test"
        })
        
        
        ylim <- reactive({
            ymin <- env$input[[ns("ymin")]]
            if (is.na(ymin))
                ymin <- min(confects()$table$effect,na.rm=T)
            ymax <- env$input[[ns("ymax")]]
            if (is.na(ymax))
                ymax <- max(confects()$table$effect,na.rm=T)
            c(ymin, ymax)
        })
        
        confects_with_desc <- reactive({
            result <- confects()
            ylabel <- env$input[[ns("ylabel")]]
            if (ylabel != "")
                result$effect_desc <- ylabel
            result
        })
        
        env[[ns("me_plot-callback")]] <- function() {
            result <- topconfects::confects_plot_me(confects_with_desc()) +
                coord_cartesian(ylim=ylim())
            print(result)
        }
        
        env[[ns("me2_plot-callback")]] <- function() {
            breaks <- env$input[[ns("breaks")]]
            breaks <- strsplit(breaks, "(,|\\s)+")[[1]]
            breaks <- na.omit(as.numeric(breaks))
            
            # TODO: move function to topconfects
            result <- topconfects::confects_plot_me2(confects_with_desc(), breaks=breaks) +
                coord_cartesian(ylim=ylim())
            print(result)
        }
        
        env$output[[ns("description")]] <- renderUI({
            desc <- topconfects:::confects_description(confects())
            if (!is.null(confects()$technical_var))
                desc <- sprintf("%s\nPer-read variance is %.1f^2 times per-sample variance", 
                    desc, sqrt(confects()$technical_var))
            
            if (!is.null(confects()$table$se)) {
                ss_effect <- sum(confects()$table$effect ^ 2, na.rm=TRUE)
                ss_se <- sum(confects()$table$se^2, na.rm=TRUE)
                ss_true <- ss_effect - ss_se
                #ss_lower <- ss_confidence_bound(confects()$table$effect, confects()$table$se^2, 0.975)
                #ss_upper <- ss_confidence_bound(confects()$table$effect, confects()$table$se^2, 0.025)
                ss_confect <- sum(confects()$table$confect ^ 2, na.rm=TRUE)
                #desc <- sprintf("%s\nEstimated Sum of Squares of real effects %.2f\n95%% CI [%.2f,%.2f]\n(SS effects %.2f minus SS standard errors %.2f)\n(SS confects %.2f)",
                #    desc, ss_true, ss_lower,ss_upper, ss_effect, ss_se, ss_confect)
                desc <- sprintf("%s\nEstimated Sum of Squares of real effects %.2f\n(SS effects %.2f minus SS standard errors %.2f)\n(SS confects %.2f)",
                    desc, ss_true, ss_effect, ss_se, ss_confect)
            }
            
            shiny::div(
                shiny::h2( title() ),
                shiny::pre( desc ))
        })

        env$output[[ns("result_description")]] <- renderUI({
            shiny::h2( title() )
        })

        env[[ns("results-df")]] <- reactive({
            contrast_cols <- colnames(confects()$contrasts)

            confects()$table %>%
            select(-index) %>%
            prioritize_columns(
                c("rank", "confect", "effect", contrast_cols, 
                  "AveTail", "logCPM", "AveExpr", "row_mean",
                  "fdr_zero", "gene", "biotype", "name", "product"))
        })

        env[[ns("results-options")]] <- reactive({ 
            col_names <- colnames(env[[ns("results-df")]]())
            cols <- as.list(seq_along(col_names)-1)
            names(cols) <- col_names
            
            if (is.null(confects()$limits) || any(is.na(confects()$limits))) {
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

            contrast_cols <- colnames(confects()$contrasts)

            list(
                pageLength=20,

                columnDefs=c(
                    list(
                        list(targets=cols$confect, visible=FALSE),
                        confect_coldef(cols$effect, cols$confect, low, high),
                        fixed_coldef(cols$logCPM,1),
                        fixed_coldef(cols$AveExpr,1),
                        fixed_coldef(cols$AveTail,1),
                        fixed_coldef(cols$row_mean,1),
                        fixed_coldef(cols$typical_obs_err,2),
                        precision_coldef(cols$fdr_zero,2),
                        precision_coldef(cols$se,2),
                        fixed_coldef(cols$df,1),
                        fixed_coldef(cols$`mean-tail`,1),
                        fixed_coldef(cols$`proportion-with-tail`,2)),
                    purrr::map(contrast_cols, ~fixed_coldef(cols[[.]], 2))
                    )) 
        })

        env[[ns("pipeline_dir")]] <- reactive({
            confects()$pipeline_dir
        })

        env[[ns("gene-tc")]] <- reactive(withProgress(message="Loading genes", {
            counts_reader(paste0( confects()$pipeline_dir, "/expression/genewise/counts.csv"))
        }))

        env[[ns("gene-peak_tc")]] <- reactive(withProgress(message="Loading peaks", {
            counts_reader(paste0( confects()$pipeline_dir, "/expression/peakwise/counts.csv"))
        }))

        env[[ns("gene-feature")]] <- reactive({
            rows_selected <- env$input[[ns("results-table_rows_selected")]]
            if (length(rows_selected) != 1)
                NULL
            else if ("parent" %in% colnames(confects()$table))
                as.character(confects()$table$parent[rows_selected])
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

        env$output[[ns("enrichment_out")]] <- renderUI({
            df <- confects()$table
            df$rev_rank <- nrow(df) - df$rank

            # For primarypeakwise results
            if ("parent" %in% names(df))
                df$name <- df$parent

            ranking <- env$input[[ns("enrichment_ranking")]]
            if (ranking == "abs") {
            } else if (ranking == "up") {
                df <- arrange(df, -sign(.data$effect)*.data$rev_rank)
            } else if (ranking == "down") {
                df <- arrange(df, sign(.data$effect)*.data$rev_rank)
            } else if (ranking == "fdr_zero") {
                df <- arrange(df, .data$fdr_zero)
            } else if (ranking == "fdr_zero_up") {
                df <- arrange(df, -sign(.data$effect)*-log(.data$fdr_zero), -sign(.data$effect)*.data$rev_rank)
            } else if (ranking == "fdr_zero_down") {
                df <- arrange(df, sign(.data$effect)*-log(.data$fdr_zero), sign(.data$effect)*.data$rev_rank)
            } else {
                stop("Unknown ranking")
            }
            sig <- df[!is.na(df$confect),]

            if (ranking == "up" || ranking == "fdr_zero_up")
                sig <- sig[sig$effect>0,]
            else if (ranking == "down" || ranking == "fdr_zero_down")
                sig <- sig[sig$effect<0,]

            shiny::verticalLayout(
                shiny::h3("All ", nrow(df), " genes"),
                shiny::tags$textarea( 
                    onfocus="this.select();this.scrollTop=0;", rows="10", cols="30",
                    paste(df$name, collapse="\n") ),
                shiny::h3(nrow(sig), "significant genes"),
                shiny::tags$textarea( 
                    onfocus="this.select();this.scrollTop=0;", rows="10", cols="30",
                    paste(sig$name, collapse="\n")))
        })

        env$output[[ns("diagnostic_ui")]] <- renderUI({
            diagnostics <- names(confects()$diagnostics)            
            selectizeInput(ns("diagnostic_wanted"), "Diagnostic", choices=diagnostics, width="100%")
        })

        env$output[[ns("diagnostic_text")]] <- renderUI({
            wanted <- env$input[[ns("diagnostic_wanted")]]
            shiny::req(wanted)
            text <- confects()$diagnostics[[wanted]]$text
            shiny::req(text)
            shiny::pre(paste(text, collapse="\n"))
        })

        env[[ns("diagnostic_plot-callback")]] <- function() {
            wanted <- env$input[[ns("diagnostic_wanted")]]
            shiny::req(wanted)

            diagnostic <- confects()$diagnostics[[wanted]]
            diagnostic$plot
        }
        
        me_plot$component_server(env)
        me2_plot$component_server(env)
        results_table$component_server(env)
        gene_view$component_server(env)
        diagnostic_plot$component_server(env)
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

    version <- 9
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
            div(style="height: 4em"), #,
            #actionButton(ns("cache_all"), "Ensure all tests are cached") 

            shiny::p("Results are ranked by FDR-adjusted confident effect size (confect) (see Bioconductor package topconfects). If you prefer ranking by \"significance\", rank results by the fdr_zero column, which is the traditional FDR adjusted p-value for the null hypothesis of zero effect."),
            shiny::p("Paul says: As you can see, there are now a lot of variants of the different tests. I would like to focus on using the weitrix-based tests going forward, so please use these. Weitrix methods are supported by the weitrix package which is made available in Bioconductor, and match the methods documented in the weitrix package vignettes (weitrix version 1.1.2 and higher)."),
            shiny::p("Also note the weitrix-based poly(A) tail length test is no longer log2 tail length, it's untransformed tail length."),
            shiny::p("Detrending of samples in tail length tests aims to remove systematic differences between samples. This might be samples having overall longer tails, or samples with a trend versus the mean tail length for each gene. Use the diagnostic plots to decide if this is necessary. This should only be used if a small number of genes are expected to change tail lengths! If the experiment causes a global change in tail lengths, using this may erase the global effect or even flip the sign of the effect for some genes!")
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
            selected <- isolate( env$input[[ns("variant")]] )
            
            selectizeInput(ns("variant"), "Test variant", choices=variants, selected=selected, width="100%")
        })
    
        env[[ns("test-confects")]] <- reactive({
            get(env$input[[ns("test")]], env$input[[ns("variant")]], env$input[[ns("fdr")]])
        })
    
        subapp$component_server(env)
    }
    
    composable_shiny_panels_app(panels, server, title=title)
}