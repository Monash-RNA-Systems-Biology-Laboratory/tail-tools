
#'
#' Shiny report for REPAT results
#'
#' @export
shiny_repat <- function(filename, normalizing_gene=NULL, title="3'TAP results") {
    library(shiny)
    library(nesoni)
    library(reshape2)
    library(varistran)
    library(ggplot2)
    library(gridExtra)

    tables <- read.grouped.table(filename)
    counts <- as.matrix(tables$Count)
    genes <- rownames(counts)
    samples <- colnames(counts)

    # Plots

    overview <- shiny_plot(prefix="overview", dlname="overview", width=800,height=800, function(env) {
        normed <- env$normed()
        m <- subset(normed$melted, sample %in% env$input$samples)
        m$sample <- factor(m$sample, levels = env$input$samples)
        m$gene <- factor(m$gene, levels = sort(unique(as.character(m$gene))))
        
        print( 
            ggplot(m,aes(x=sample,y=log2_fold_norm_count)) + 
            geom_point() + 
            facet_wrap(~ gene, drop=F) + 
            ylab(ifelse(normed$normalized_by_sample,"log2 fold change","log2 normalized count")) +
            xlab("") +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )        
    })
    
    overview_tail <- shiny_plot(prefix="overview_tail", dlname="overview_tail", width=800,height=800, function(env) {
        normed <- env$normed()
        m <- subset(normed$melted, sample %in% env$input$samples)
        m$sample <- factor(m$sample, levels = env$input$samples)
        m$gene <- factor(m$gene, levels = sort(unique(as.character(m$gene))))
        m <- subset(m, tail_count>0)
        print(
            ggplot(m,aes(x=sample,y=relative_tail_mean)) + 
            geom_point() + geom_point(data=subset(m,tail_count<env$input$highlight_low),color="red") +
            facet_wrap(~ gene, drop=F) +
            ylab(ifelse(normed$normalized_by_sample,"change in poly(A) tail length","poly(A) tail length")) +
            xlab("") + 
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
        )
    })
    
    individual <- shiny_plot(prefix="individual", dlname="individual", width=800, function(env) {
        m <- subset(env$normed()$melted, (sample %in% env$input$samples) & (gene == env$input$individual_gene))
        m$sample <- factor(m$sample, levels = env$input$samples)
        
        plot.new()
        grid.arrange(
            ggplot(m, aes(x=sample, y=norm_count)) +
            geom_bar(stat="identity", color="#000000", fill="#cccccc") +
            ylab("Normlized count") +
            xlab("") +
            title(paste(env$input$individual_gene,"expression")) +
            theme_bw() +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)),
            
            ggplot(m, aes(x=sample, y=tail_mean)) +
            geom_bar(stat="identity", color="#000000", fill="#ccffcc") +
            geom_bar(data=subset(m, tail_count<env$input$highlight_low), 
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

             normalizer <- counts[normalizing_gene,]
             normalizer <- normalizer / mean(normalizer)
             
             norm_counts <- t(t(counts)/normalizer)             
             log2_fold_norm_counts <- log2(norm_counts + 0.5)
             
             tail_mean <- as.matrix(tables$Tail)
             relative_tail_mean <- tail_mean
             
             normalized_by_sample <- length(env$input$normalizing_samples) > 0
             
             if (normalized_by_sample) {
                 sample_normalizer <- rowMeans(log2_fold_norm_counts[,env$input$normalizing_samples,drop=F])
                 log2_fold_norm_counts <- log2_fold_norm_counts - sample_normalizer
             
                 sample_tail_normalizer <- rowMeans(tail_mean[,env$input$normalizing_samples,drop=F])
                 relative_tail_mean <- relative_tail_mean - sample_tail_normalizer
             }
             
             melted <- melt(norm_counts)
             colnames(melted) <- c("gene", "sample", "norm_count")
             melted$log2_fold_norm_count <- melt(log2_fold_norm_counts)$value
             melted$tail_count <- melt(as.matrix(tables$Tail_count))$value
             melted$tail_mean <- melt(tail_mean)$value
             melted$relative_tail_mean <- melt(relative_tail_mean)$value
             melted$tail_median <- melt(as.matrix(tables$Tail_quantile_50))$value

             list(
                 normalizer=normalizer, 
                 counts=norm_counts, 
                 melted=melted,
                 normalized_by_sample=normalized_by_sample
             )
        })
    
        output$normalizer_table <- renderTable({
             v <- env$normed()$normalizer[input$samples]
             data.frame(sample=names(v), normalizer=v, row.names=NULL)
        }, digits=6)
        
        overview$component_server(env)
        overview_tail$component_server(env)
        individual$component_server(env)
    }
    
    composable_shiny_app(ui, server)
}

