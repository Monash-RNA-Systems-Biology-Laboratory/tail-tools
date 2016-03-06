

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
                color = "#0c0"
            else if (fdr <= 0.05)
                color = "#cc0"
            else 
                color = "#888"
            return "<span style=\\"color:" + color + "\\">" + data + " " +
                "<div title=FDR style=\\"display: inline-block; font-size: 66%; width: 3em;\\">(" + 
                fdr_text + ")</div></span>";
        }')
    )
}


#'
#' Show a report about an end shift analysis.
#'
#' @param result output from end_shift()
#'
#' @export
shiny_end_shift <- function(result) {
    have_fdr <- "fdr" %in% colnames(result$table)
    
    get_cols <- c("rank","r","r_low","r_high","fdr",
                  "edger_rank","edger_fdr","limma_rank","limma_fdr",
                  "gene","biotype","product","parent")
    get_cols <- get_cols[ get_cols %in% colnames(result$table) ]
    df <- result$table[,get_cols]

    # column numbers as referred to by DT, counting from 0
    cols <- as.list(seq_len(ncol(df))-1)
    names(cols) <- colnames(df)
    
    ui <- fluidPage(
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
        navlistPanel(
            widths=c(2,10),
            well=FALSE,
            tabPanel("Overview", 
                uiOutput("overview")),
            tabPanel("Results",
                DT::dataTableOutput("table"),
                downloadButton("table_download", "Download CSV file"),
                uiOutput("blurb"),
                DT::dataTableOutput("gene_table")
            )
        )
    )
    
    server <- function(input,output,session) {
        output$overview <- renderUI({
            tags$div(
                tags$h2("Overview"),
                tags$p( nrow(df), "genes with multiple peaks." ),
                tags$p( sprintf("%.1f",result$edger$dge$prior.df), "edgeR prior degrees of freedom." ),
                tags$p( sprintf("%.1f",result$limma$fit$df.prior), "limma prior degrees of freedom." ),
                tags$p("(More prior degrees of freedom indicates more consistent expression levels between genes, and generally leads to more significant results.)"),
                
                tags$h2("Interpretation"),
                
                tags$p("r values lie between -1 and 1. Positive r values indicate a shift from proximal to distal peak usage, and negative vice versa."),
                
                tags$p( sprintf("%g%%",result$ci*100), "confidence intervals on r account only for technical variance estimating r (ie genes with few reads will have a wide interval). Genes are ranked by the end of the confidence interval closest to 0."),
                
                tags$p("The FDR values associated with r are based on a permuation test, and are quite conservative, espcially if the number of possible permutations of samples (respecting grouping) is limited."),
                
                tags$p("EdgeR and limma rankings and FDR values are based on differential exon usage tests applied to this data. The edgeR \"gene\" and limma \"F\" methods were used. These will pick up a shift in peak usage, but if there are more than two peaks no attention is paid to the order of peaks (in contrast to the r statistic, where the order is important).")
            )
        })
    
    
        output$table <- DT::renderDataTable(
            df, 
            rownames=F, 
            selection="single",
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
        )
        
        
        output$table_download <- downloadHandler(
            filename="alternative-utrs.csv",
            content=function(file) {
                write_csv(df, file)
            }
        )
        
        
        output$blurb <- renderUI({
            row <- input$table_row_last_clicked
            if (is.null(row)) {
                return("")
            }
            
            row <- as.integer(row)
            tags$h3(df$gene[row], df$parent[row])
        })
        
        output$gene_table <- reactive({
            row <- input$table_row_last_clicked
            if (is.null(row)) {
                return(NULL)
            }
            
            row <- as.integer(row)
            parent <- result$table$parent[row]
            peaks <- result$splitter[[parent]]
            n_peaks <- length(peaks)
            
            out <- data_frame(sample=colnames(result$counts))
            if (!is.null(result$group))
                out$group <- result$group
            out$condition <- ifelse(result$condition, "+","-")
            #escape <- ncol(out)
            
            mat <- t(result$counts[peaks,,drop=F]) 
            colnames(mat) <- paste(result$peak_info$id[peaks], result$peak_info$relation[peaks])
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
                df
            )(shinysession=session, name="gene_table") 
        })
    }

    shinyApp(ui,server)
}


#'
#' Perform an end shift analysis and report results.
#'
#' @param path Directory containing tail-tools pipeline output
#'
#' @export
shiny_end_shift_pipeline <- function(path) {

}


