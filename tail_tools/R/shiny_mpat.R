

meltdown <- function(mat, rows, columns, values) {
    result <- reshape2::melt(t(as.matrix(mat)))
    colnames(result) <- c(rows,columns,values)
    result
}

#'
#' Shiny report for mPAT results
#'
#' @export
shiny_mpat <- function(filename, normalizing_gene=NULL, title="mPAT results") {
    library(shiny)
    library(nesoni)
    library(reshape2)
    library(varistran)
    library(ggplot2)
    library(gridExtra)
    library(dplyr)

    tables <- read.grouped.table(filename)
    genes <- rownames(tables$Count)
    samples <- colnames(tables$Count)
    
    raw_data <-
        meltdown(tables$Count, "sample","gene","count") %>%
        left_join(meltdown(tables$Tail, "sample","gene","tail"), 
                  c("sample","gene")) %>%
        left_join(meltdown(tables$Tail_count, "sample", "gene", "tail_count"), 
                  c("sample","gene"))
    

    # Plots

    overview <- shiny_plot(prefix="overview", dlname="overview", width=800,height=800, function(env) {
        normed <- env$normed()
        
        print( 
            ggplot(normed$norm_data,aes(x=sample,y=log2_fold_norm_count)) + 
            geom_point() + 
            facet_wrap(~ gene, drop=F) + 
            ylab(ifelse(normed$is_differential,"log2 fold change","log2 normalized count")) +
            xlab("") +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )        
    })
    
    overview_tail <- shiny_plot(prefix="overview_tail", dlname="overview_tail", width=800,height=800, function(env) {
        normed <- env$normed()
        df <- filter_(normed$norm_data, ~tail_count > 0)
        df_bad <- filter_(df, ~is_low_tail_count)
        print(
            ggplot(df,aes(x=sample,y=relative_tail)) + 
            geom_point() + 
            geom_point(data=df_bad,color="red") +
            facet_wrap(~ gene, drop=F) +
            ylab(ifelse(normed$is_differential,"change in poly(A) tail length","poly(A) tail length")) +
            xlab("") + 
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )
    })
    
    individual <- shiny_plot(prefix="individual", dlname="individual", width=800, function(env) {
        normed <- env$normed()
        individual_gene <- env$input$individual_gene
        df <- filter_(normed$norm_data, ~gene == individual_gene)
        df_bad <- filter_(df, ~is_low_tail_count)
        
        plot.new()
        grid.arrange(
            ggplot(df, aes(x=sample, y=norm_count)) +
            geom_bar(stat="identity", color="#000000", fill="#cccccc") +
            ylab("Normlized count") +
            xlab("") +
            title(paste(env$input$individual_gene,"expression")) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
            
            ggplot(df, aes(x=sample, y=tail)) +
            geom_bar(stat="identity", color="#000000", fill="#ccffcc") +
            geom_bar(data=df_bad, 
                     stat="identity", color="#000000", fill="#ff0000") +
            ylab("poly(A) tail length") +
            xlab("") +
            title(paste(env$input$individual_gene,"poly(A) tail"))+
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
            
            ncol=2, newpage=F
        )    
    })


    # App

    ui <- div(
        titlePanel(title),
        navlistPanel(
            widths=c(2,10),
            well=FALSE,
            tabPanel("Select",
                selectInput("normalizing_gene","Normalizing gene",choices=sort(genes),selected=normalizing_gene),
                br(),
                
                selectInput("genes","Genes to display",choices=genes,selected=genes,multiple=TRUE,width="100%"),
                br(),
                
                selectInput("samples","Samples to display",choices=samples,selected=samples,multiple=TRUE,width="100%"),
                p(
                    actionButton("samples_all", "All")
                ),
                br(),
                
                selectInput("normalizing_samples","Samples to normalize to in overview",choices=samples,selected=c(),multiple=TRUE,width="100%"),
                p(
                    actionButton("normalizing_samples_all", "All"),
                    actionButton("normalizing_samples_none", "None")
                ),
                p("Tail lengths and log2 read counts are shown relative to the average of these samples, in the overview."), 
                p("Select none to show absolute log2 read counts and tail lengths in overview."),
                br(),
                
                numericInput("highlight_low", "Highlight in red tail lengths based on less than this many reads", 50),
                br(),
                
                h1("Normalization"),
                p("Read counts are divided by this to obtain normalized read counts."),
                htmlOutput("normalizer_table")
            ),
            tabPanel("Overview",
                h1("Expression levels"),
                p("Note these are shown on a log scale."),
                overview$component_ui,
                h1("Tail lengths"),
                p("Where there are few poly(A) reads, the length is shown in red."),
                overview_tail$component_ui
            ),
            tabPanel("Individual genes",
                selectInput("individual_gene", "Gene", choices=sort(genes)),
                p("Note expression levels are *not* log transformed in this plot."),
                individual$component_ui
            )
        )
    )
    
    server <- function(env) {
        input <- env$input
        output <- env$output
        
        observeEvent(input$samples_all, {
            updateSelectInput(env$session, "samples", selected=samples)
        })
    
        observeEvent(input$normalizing_samples_all, {
            updateSelectInput(env$session, "normalizing_samples", selected=samples)
        })
    
        observeEvent(input$normalizing_samples_none, {
            updateSelectInput(env$session, "normalizing_samples", selected=character(0))
        })
    
        env$normed <- reactive({
             normalizing_gene <- input$normalizing_gene
             selected_samples <- input$samples
             selected_genes <- input$genes
             highlight_low <- input$highlight_low
             
             normalizer <- raw_data %>%
                 filter_(~ gene == normalizing_gene) %>%
                 mutate_(normalizer =~ count / mean(count)) %>%
                 select_(~sample, ~normalizer)
             
             norm_data <- raw_data %>%
                 left_join(normalizer, "sample") %>%
                 mutate_(
                     norm_count =~ count / normalizer,
                     log2_norm_count =~ log2(norm_count + 0.5),
                     log2_fold_norm_count =~ log2_norm_count,
                     relative_tail =~ tail,
                     is_low_tail_count =~ tail_count < highlight_low
                 )
             
             is_differential <- length(env$input$normalizing_samples) > 0
             
             if (is_differential) {
                 baseline <- norm_counts %>%
                     select_(~sample %in% env$input$normalizing_samples) %>%
                     group_by_(~gene) %>%
                     summarize_(
                         baseline =~ mean(log2_fold_norm_count),
                         baseline_tail =~ mean(tail, na.rm=TRUE)
                     ) %>%
                     select_(~gene, ~baseline, ~baseline_tail)
             
                 norm_data <- norm_data %>%
                     left_join(baseline, "gene") %>%
                     mutate_(
                         log2_fold_norm_count =~ log2_fold_norm_count - baseline,
                         relative_tail =~ tail - baseline_tail
                     )
             }     
             
             norm_data <- norm_data %>%
                 filter_(
                     ~sample %in% selected_samples,
                     ~gene %in% selected_genes
                 ) %>%
                 mutate_(
                     sample =~ factor(sample, selected_samples),
                     gene =~ factor(gene, selected_genes)
                 )
                 
             
             list(
                 norm_data=norm_data,
                 normalizer=normalizer, 
                 is_differential=is_differential
             )
        })
    
        output$normalizer_table <- renderTable({
             env$normed()$normalizer
        }, digits=6)
        
        overview$component_server(env)
        overview_tail$component_server(env)
        individual$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}



#' @export
shiny_repat <- shiny_mpat

