

r_coldef <- function(target, r_low, r_high) {
    list(
        targets=target,
        className="dt-center",
        render=JS('
        function(data,type,row,meta) {
            var r = parseFloat(data);
            var r_low = parseFloat(row[',r_low,']);
            var r_high = parseFloat(row[',r_high,']);
            function tr(x) { return (x+1)*100+1.5; }
            return "<tt style=\\"white-space: pre\\">" + 
                (" "+r.toFixed(2)).slice(-5) + 
                "<span style=\\"font-size: 66%;\\"> [" + (" "+r_low.toFixed(2)).slice(-5) + 
                "," + (" "+r_high.toFixed(2)).slice(-5) + "]</span>" +
                "</tt> " +
                "<svg width=203 height=21 style=\\"vertical-align: middle\\">" +
                "<line x1="+tr(-1)+" y1=0 x2="+tr(-1)+" y2=21 stroke=#888 stroke-width=1 />" +
                "<line x1="+tr( 0)+" y1=0 x2="+tr( 0)+" y2=21 stroke=#888 stroke-width=1 />" +
                "<line x1="+tr( 1)+" y1=0 x2="+tr( 1)+" y2=21 stroke=#888 stroke-width=1 />" +
                "<line x1="+tr(r_low)+" y1=10 x2="+tr(r_high)+" y2=10 stroke=#000 stroke-width=2 />" +
                "<circle cx="+tr(r)+" cy=10 r=4 fill=#000 />" +
                "</svg>"
        }')
    )
}

n_coldef <- function(targets, maximum) {
    list(
        targets=targets,
        render=JS('
        function(data,type,row,meta) {
            var max = ',maximum,';
            var val = parseFloat(data);
            var size = 21*Math.sqrt(val/max);
            var pos = (21-size)/2;
            return "<tt>" + data + "</tt> " +
                "<svg width=21 height=21 style=\\"vertical-align: middle\\">" +
                "<rect x="+pos+" y="+pos+" width="+size+" height="+size+" style=\\"fill: #88f;\\" />"+
                "</svg>";
        }')
    )
}

rank_coldef <- function(target, fdr) {
    if (is.null(fdr))
        return(list(targets=target))

    list(
        targets=target,
        render=JS('function(data,type,row,meta) {
            var fdr = parseFloat(row[',fdr,']);
            var stars = "";
            var fdr_text, color;
            if (fdr <= 0.01)
                fdr_text = fdr.toExponential(0);
            else
                fdr_text = fdr.toPrecision(1);
            if (fdr <= 0.01)
                color = "#0a0"
            else if (fdr <= 0.05)
                color = "#cc0"
            else 
                color = "#888"
            return "<span style=\\"color:" + color + "\\">" + data + " " +
                "<div title=FDR style=\\"display: inline-block; font-size: 80%; width: 4em;\\">(" + 
                fdr_text + ")</div></span>";
        }')
    )
}



composable_shiny_panels_app <- function(panels, server, prefix="") {
    ui <- div(
        HTML('
<style>
table.dataTable.display tbody td { 
    border-top: 0; 
    border-bottom: 0; 
    border-style: none;
    padding-top: 0; 
    padding-bottom: 0; 
    line-height: 1;
    white-space: nowrap;
}
</style>
'),
        do.call(navlistPanel, c(list(id=paste0(prefix, "tabset"), widths=c(2,10),well=FALSE), panels)),
        div(style="height: 4em")
    )

    app <- composable_shiny_app(ui, server)
    app$component_panels <- panels
    app
}


#'
#' Show a report about an end shift analysis.
#'
#' @param result output from end_shift()
#'
#' @export
shiny_end_shift <- function(result) {
    result <- ensure_reactable(result)

    mr_plot <- shiny_plot(function(env) {
        df <- env$vars()$df
        if (is.null(df)) return()
        
        plot <- ggplot(df,aes(x=mean_reads, y=r)) +
            scale_x_log10() +
            scale_color_discrete("Significant by") +
            geom_segment(aes(xend=mean_reads, y=r_low, yend=r_high), color="#bbbbbb") +
            geom_point() +
            theme_bw() +
            labs(x = "Mean reads per sample (log scale)",
                 y = "r")
        
        size <- 3
        if (!is.null(df$fdr) && any(df$fdr <= 0.05)) {
            plot <- plot + geom_point(
                data=filter(df, fdr <= 0.05),aes(color="r"),
                shape=1,stroke=1.5,size=size)
            size <- size + 2
        }

        if (!is.null(df$edger_fdr) && any(df$edger_fdr <= 0.05)) {
            plot <- plot + geom_point(
                data=filter(df, edger_fdr <= 0.05),aes(color="edgeR"),
                shape=1,stroke=1.5,size=size)
            size <- size + 2
        }
        
        if (!is.null(df$limma_fdr) && any(df$limma_fdr <= 0.05)) {
            plot <- plot + geom_point(
                data=filter(df, limma_fdr <= 0.05),aes(color="limma"),
                shape=1,stroke=1.5,size=size)
            size <- size + 2
        }
        
        plot
    }, prefix="mr_plot_", width=800, brush=brushOpts(id="mr_plot_brush", delay=600000))

    panels <- list(
        tabPanel("Overview", 
            uiOutput("overview")),
        tabPanel("Results",
            DT::dataTableOutput("table"),
            downloadButton("table_download", "Download CSV file"),
            uiOutput("blurb"),
            DT::dataTableOutput("gene_table")
        )
    )
        
    server <- function(env) {
        output <- env$output
        input <- env$input
        session <- env$session
                
        env$vars <- reactive({
            result_val <- result(env)
            
            have_fdr <- "fdr" %in% colnames(result_val$table)
            
            get_cols <- c("rank","r","r_low","r_high","fdr",
                          "edger_rank","edger_fdr","limma_rank","limma_fdr",
                          "gene","biotype","product","parent","mean_reads")
            get_cols <- get_cols[ get_cols %in% colnames(result_val$table) ]
            df <- result_val$table[,get_cols]
            
            # column numbers as referred to by DT, counting from 0
            cols <- as.list(seq_len(ncol(df))-1)
            names(cols) <- colnames(df)

            list(result=result_val, df=df, have_fdr=have_fdr, cols=cols)
        })

        env$selected_rows <- reactive({
            df <- env$vars()$df
            
            if (!is.null(input$mr_plot_brush)) {
                df <- df %>%
                    filter(!is.na(r)) %>%
                    brushedPoints(input$mr_plot_brush, "mean_reads", "r")
            }
            
            df
        })
        
        observe({
            if (!is.null(input$mr_plot_brush))
                updateTabsetPanel(session, "tabset", selected="Results")
        })
    
        output$overview <- renderUI({
            vars <- env$vars()
            
            items <- list(
                tags$h2( vars$result$title ),
                mr_plot$component_ui,
                tags$p("Drag to select a subset of genes."),
                HTML("<p>Genes significant by various methods with FDR &le; 0.05 are circled.</p>"),
                tags$p( nrow(vars$df), "genes with multiple peaks and averaging at least", vars$result$min_reads, "reads per sample.")
            )
            
            if (!is.null(vars$result$edger))
                items <- c(items,list( 
                    tags$p( sprintf("%.1f",vars$result$edger$dge$prior.df), "edgeR prior degrees of freedom." )
                ))
            
            if (!is.null(vars$result$limma))
                items <- c(items,list(
                    tags$p( sprintf("%.1f",vars$result$limma$fit$df.prior), "limma prior degrees of freedom." )
                ))
            
            items <- c(items,list(
                tags$h2("Interpretation"),
                
                tags$p("r values lie between -1 and 1. Positive r values indicate a shift from proximal to distal peak usage, and negative vice versa."),
                
                tags$p( sprintf("%g%%",vars$result$ci*100), "confidence intervals on r account only for technical variance estimating r (ie genes with few reads will have a wide interval). Genes are ranked by the end of the confidence interval closest to 0."),
                
                tags$p("FDR values are shown in brackets next to rankings."),
                
                tags$p("The FDR values associated with r are based on a permuation test, and are quite conservative, espcially if the number of possible permutations of samples (respecting grouping) is limited."),
                
                tags$p("EdgeR and limma rankings and FDR values are based on differential exon usage tests applied to this data. The edgeR \"gene\" and limma \"F\" methods were used. These will pick up a shift in peak usage, but if there are more than two peaks no attention is paid to the order of peaks (in contrast to the r statistic, where the order is important)."),

                tags$p("More prior degrees of freedom in edgeR and limma indicates more consistent expression levels between genes, and generally leads to more significant results."),
                
                tags$h2("Details"), 

                tags$p("r definition: Consider each pairing of reads with one read from condition - and one read from condition +. Assign the value 1 if the + read is further 3' than the - read, -1 if the opposite, and 0 if they are from the same peak. r is then the average of all these values."),
                
                tags$p("r can be viewed as a rescaled Mann-Whitney-Wilcoxon U statistic. This has an associated normal approximation to the error distribution, so we can give an accuracy to which r has been estimated, at least in terms of technical variation."),
                
                tags$p("When there are multiple groups, r is estimated for each of the groups and then these r values are averaged weighted by their precision (1/variance).")
            ))
            
            do.call(tags$div, items)
        })
    
    
        output$table <- reactive({
            vars <- env$vars()
            cols <- vars$cols            
            df <- env$selected_rows()
            
            DT::renderDataTable(
                df, 
                rownames=F, 
                selection=list(mode="single", selection=c(1)),
                options=list(
                    pageLength=20,
                    columnDefs=c(
                        list(list(
                            targets=c(cols$r_low, cols$r_high, cols$fdr, cols$edger_fdr, cols$limma_fdr),
                            visible=FALSE
                        )),
                        list(rank_coldef(cols$rank, cols$fdr)),
                        list(r_coldef(cols$r, cols$r_low, cols$r_high)),
                        { if (is.null(cols$edger_rank)) list()
                          else list(rank_coldef(cols$edger_rank, cols$edger_fdr)) },
                        { if (is.null(cols$limma_rank)) list()
                          else list(rank_coldef(cols$limma_rank, cols$limma_fdr)) }
                    )
                )
            )(shinysession=session, name="table")
        })
        
        
        output$table_download <- downloadHandler(
            filename="alternative-utrs.csv",
            content=function(file) {
                write_csv(env$selected_rows(), file)
            }
        )
        
        
        output$blurb <- renderUI({
            vars <- env$vars()
            df <- env$selected_rows()
            
            row <- as.integer(input$table_rows_selected)
            if (length(row) != 1 || row > nrow(df)) {
                return(tags$div(style="height:2em"))
            }
            
            row <- as.integer(row)
            tags$h3(df$gene[row], df$parent[row])
        })
        
        output$gene_table <- reactive({
            vars <- env$vars()
            df <- env$selected_rows()
        
            row <- as.integer(input$table_rows_selected)
            if (length(row) != 1 || row > nrow(df)) {
                return(DT::renderDataTable(
                    data_frame("Select a gene to view details"=character(0)),
                    rownames=F,
                    options=list(dom="")
                )(shinysession=session, name="gene_table"))
            }
            
            parent <- df$parent[row]
            peaks <- vars$result$splitter[[parent]]
            n_peaks <- length(peaks)
            
            out <- data_frame(sample=colnames(vars$result$counts))
            if (!is.null(vars$result$group))
                out$group <- vars$result$group
            out$condition <- ifelse(vars$result$condition, "+","-")
            #escape <- ncol(out)
            
            mat <- t(vars$result$counts[peaks,,drop=F]) 
            colnames(mat) <- paste(vars$result$peak_info$id[peaks], vars$result$peak_info$relation[peaks])
            max_mat <- max(mat)
            mat <- mat %>% as.data.frame      
            
            df <- 
                bind_cols(out, mat) %>%
                arrange(condition, group)
    
            DT::renderDataTable(
                options=list(
                    pageLength=100,
                    dom="t",
                    columnDefs=list(
                        n_coldef( seq_len(n_peaks)+(ncol(df)-n_peaks-1), max_mat)
                    )
                ), 
                rownames=F,
                selection=list(mode="single",selection=integer(0)),
                df
            )(shinysession=session, name="gene_table") 
        })

        
        mr_plot$component_server(env)
    }

    app <- composable_shiny_panels_app(panels,server)
}


#'
#' Perform end shift analyses on pipeline output and report results.
#'
#' @param tests Named list of lists giving parameters to end_shift_pipeline. Each name should be unique.
#'
#' @export
shiny_end_shift_pipeline <- function(tests, cache_prefix="cache_") {
    get <- function(name, peak_set) {
        param <- tests[[name]]        
        param$antisense <- peak_set == "all"
        param$non_utr <- peak_set != "utr"
        
        filename <- paste0(cache_prefix,name,"_",peak_set,".rds")
        
        cached_value <- NULL
        if (file.exists(filename)) {
            cached_value <- readRDS(filename)
            if (identical(cached_value$param, param))
                return(cached_value$result)
        }
        
        result <- withProgress(
            message=paste0("Computing ",name," ",peak_set),
            do.call(end_shift_pipeline, param)
        )
        
        saveRDS(list(param=param,result=result), filename)
        
        result
    }
    
    assert_that(is.list(tests))
    assert_that(!is.null(names(tests)))
    assert_that(!any(duplicated(names(tests))))
    
    titles <- map2_chr(names(tests), tests, function(name,param) {
        if (!is.null(param$title)) param$title 
        else name
    })
    choices <- as.list(names(tests))
    names(choices) <- titles
    
    peak_choices = list(
        "All" = "all",
        "Sense strand only" = "sense",
        "3' UTR only" = "utr"
    )
        
    subapp <- shiny_end_shift(function(env) env$test_result())
    
    panels <- list(
        tabPanel("Select test",
            h2("Select test"),
            selectizeInput("test", "Test", choices=choices, selected=NULL, width="100%"),
            selectizeInput("peak_set", "Peaks to use", choices=peak_choices, selected="sense"),
            div(style="height: 4em"),
            actionButton("cache_all", "Ensure all tests are cached") 
        )
    )
    
    panels <- c(panels, subapp$component_panels)
        
    server <- function(env) {
        observe({
            if (env$input$cache_all == 0) return()
            
            for(name in names(tests))
                for(peak_set in unlist(peak_choices))
                    get(name, peak_set)
        })
    
        env$test_result <- reactive({
            get(env$input$test, env$input$peak_set)
        })    
    
        subapp$component_server(env)
    }
    
    composable_shiny_panels_app(panels, server)
}







