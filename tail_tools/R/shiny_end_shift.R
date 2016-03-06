

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
                "<span style=\\"font-size: 75%;\\"> [" + (" "+r_low.toFixed(2)).slice(-5) + 
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

rankstar_coldef <- function(target, fdr) {
    list(
        targets=target,
        render=JS('function(data,type,row,meta) {
            var fdr = parseFloat(row[',fdr,']);
            var stars = "";
            if (fdr <= 0.0001) stars = "&lowast;&lowast;&lowast;&lowast;"
            else if (fdr <= 0.001) stars = "&lowast;&lowast;&lowast;"
            else if (fdr <= 0.01) stars = "&lowast;&lowast;"
            else if (fdr <= 0.05) stars = "&lowast;" 
            else if (fdr <= 0.1) stars = "+";
            return "<span title=\\"fdr=" + fdr + "\\">" + 
                data + " " +
                "<div style=\\"display: inline-block; width: 4em;\\">" + stars + "</div>" +
                "</span>";
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
    have_fdr <- "fdr" %in% colnames(result)
    if (have_fdr)
        df <- result$table %>%        
            select(rank,r,r_low,r_high,fdr,gene,biotype,product,parent)
    else
        df <- result$table %>%        
            select(rank,r,r_low,r_high,gene,biotype,product,parent)
    
    cols <- as.list(seq_len(ncol(df))-1)
    names(cols) <- colnames(df)
    
    ui <- tags$div(
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
        dataTableOutput("table"),
        downloadButton("table_download", "Download CSV file"),
        uiOutput("blurb"),
        dataTableOutput("gene_table")
    )
    
    server <- function(input,output,session) {
        output$table <- renderDataTable(
            df, 
            rownames=F, 
            selection="single",
            options=list(
                pageLength=20,
                columnDefs=list(
                    list(
                        targets=c(cols$r_low, cols$r_high),
                        visible=FALSE
                    ),
                    #rankstar_coldef(cols$rank, cols$fdr),
                    r_coldef(cols$r, cols$r_low, cols$r_high)
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
            colnames(mat) <- paste(colnames(mat), result$peak_info$relation[peaks])
            max_mat <- max(mat)
            mat <- mat %>% as.data.frame      
            
            df <- 
                bind_cols(out, mat) %>%
                arrange(condition, group)
    
            renderDataTable(
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


