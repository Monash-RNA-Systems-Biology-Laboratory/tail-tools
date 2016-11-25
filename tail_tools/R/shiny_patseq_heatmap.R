#' @title Produces detailed heatmap
#' @details 
#' Workhorse function for this package.
#'      
#' @param datfr List of dataframes 
#' Takes a read.grouped.table() as input or a list of four dataframes (more data frames are ok but it only uses these):
#' Counts - Genewise counts of expression.
#' Tail - Mean tail length 
#' Tail_counts - Number of poly-A tails counted.
#' Annotation - Information regarding the annotation information
#'      Gene name, chromosome, gene product, biotype etc...
#'  
#' @param sample_labels Sample labels
#' @param sample_labels2 Sample labels (second plot)
#' @param feature_labels Feature labels
#' @param prefix Prefix for plot
#' @param species Species of the data. Currently supports Human (Hs), Saccharomyces cerevisiae (Sc), Caenorhabditis elegans (Ce), Mus musculus (Mm)
#' Now optional. Not entering this disables GO term analysis
#' 
#' @return 
#' Returns a composable shiny app object
#' 
#' @export

shiny_patseq_heatmap <- function(datfr, sample_labels=NULL, sample_labels2=NULL, feature_labels=NULL, prefix="", species=NULL) {
    if(!("Annotation" %in% names(datfr)) || !("Count" %in% names(datfr)) || !("Tail_count" %in% names(datfr)) || !("Tail" %in% names(datfr))){
        cat("This a summary of your dataframe: \n")
        print(summary(datfr))
        stop("List is missing one or more of: Tail, Count, Tail_count, Annotation")
    }
    dmrow <- c(nrow(datfr$Count), nrow(datfr$Tail_count), nrow(datfr$Tail), nrow(datfr$Annotation))
    if(length(levels(as.factor(dmrow))) != 1){
        cat("Dimensions of your data: \n")
        for(i in 1:length(datfr)){
            cat(dim(datfr[[i]]),"\t", (names(datfr))[i],"\n")
        }
        stop("Number of rows in Tail, Count, Tail_count, Annotation not equal")
    }
    
    if(is.null(species) || sum(species == c("Hs", "Sc", "Ce", "Mm"))!=1){
        GOenable <- FALSE
        if (!is.null(species)) cat("Species is not supported. Not enabling GO Term Analysis\n")
    } else {
        GOenable <- TRUE
    }
    
    p <- NS(prefix)
    sample_labels <- ensure_reactable(sample_labels)
    sample_labels2 <- ensure_reactable(sample_labels2)
    feature_labels <- ensure_reactable(feature_labels)
    
    plot <- shiny_patseq_heatmap_inner(
        callback = function(env) {
            print(env[[p("grob")]]())
        },
        rorder = function(env) {
            env[[p("grob")]]()$info$row_order$order
        },
        width=1000,
        height=700,
        dlname="heatmap",
        prefix=p("plot"),
        goenabl=GOenable,
        selin = function(env){
            env[[p("seldat")]]()$ann
        },
        spp=species,
        hgc = function(env){
            env[[p("hgc")]]()
        },
        otype = function(env){
            env[[p("otype")]]()
        }
    )
    
    # Shiny's UI layout 
    colvec <- names(datfr$Tail)
    panels <- list(
        function(request) shiny::tabPanel("Select and filter",
            shiny::h2("Features"),
            shiny::radioButtons(p("selFeat"), 
                                label="Show features with greatest:", 
                                choices=list("Span of tail length"=1, "Span of expression"=2, "Average expression"=3), 
                                selected=1,
                                inline=TRUE),
            shiny::radioButtons(p("roword"), 
                                label="Features ordered by: ", 
                                choices=list("Tail Length"=1, "Expression"=2, "Chromosomal location"=3), 
                                selected=1,
                                inline=TRUE),
            shiny::numericInput(p("n"), "Number of features to show", 50, min=10,max=2000,step=10),
            shiny::br(),
            
            shiny::h2("Samples"),
            selectizeInput(p("choosecol"), "Choose samples to display",multiple=T,colvec, selected=colvec, width="100%"),
            shiny::radioButtons(p("clusterby"), 
                                label="Order samples by: ", 
                                choices=list("Tail length" = 2, "Expression" = 3, "Retain existing order" = 1), 
                                selected = 1,
                                inline=TRUE),
            shiny::br(),
            
            shiny::h2("Options"),
            shiny::checkboxInput(p("tail_relative"), "Tail lengths relative to mean for feature.", value=TRUE),
            shiny::numericInput(p("nmin"), "Trim Tail Counts below value to NA", 50, min=0,max=1000,step=1),
            shiny::numericInput(p("expmin"), "Exclude genes with all expression counts below: ", 0, min=0,max=1500,step=1)
        ),
        function(request) shiny::tabPanel("Tail heatmap", call_ui(plot$component_ui, request)),
        function(request) shiny::tabPanel("GO analysis of selection",
            shiny::fluidRow(
                shiny::column(3, shiny::numericInput(p("fdrcutoff"), "p value cutoff", 0.05, min=0.0001, max=1, step=0.01)),
                shiny::column(3, shiny::radioButtons(p("ontype"), label="Ontology search type",
                                                     choices=list("Biological Processes"=1, "Cellular Component"=2,"Molecular Function"=3),
                                                     selected=1,
                                                     inline=FALSE)),
                shiny::column(3, shiny::tags$label("Download analysis as .csv"), shiny::tags$br(),
                              shiny::downloadButton(p("dlanalysis"), ".csv"))
            ),
            shiny::textOutput(p("goerror")),
            DT::dataTableOutput(p("gotab")),
            
            shiny::p("Selected genes are compared to remaining gene in heatmap. Significance testing is by hypergeometric test (Fisher's Exact Test), no correction for multiple testing.")
            #shiny::tableOutput("tabout"))
        )
    )
    
    # Shiny's server
    server <- function(env) {
        # Nice idea, doesn't work well with DT or brushes
        #shinyURL::shinyURL.server(env$session)
        
        # Processes the input from datfr into the correct length and returns a list of 4 data frames
        wproc <- reactive({
            
            datfr2 <- list()
            
            #colvec <- which(names(datfr$Tail) %in% env$input[[p("choosecol")]])
            colvec <- match(env$input[[p("choosecol")]], names(datfr$Tail))
            
            #datfr$Tail[datfr$Tail_count < env$input[[p("nmin")]]] <- NA
            
            datfr2$Tail <- datfr$Tail[,colvec,drop=F]
            datfr2$Tail[datfr$Tail_count[,colvec,drop=F] < env$input[[p("nmin")]]] <- NA
            
            datfr2$Count <- datfr$Count[,colvec,drop=F]
            
            hldvec2 <- apply(datfr2$Tail,1,function(row) {
                any(sum(!is.na(row))>= 2)
            })
            
            datfr2$Tail <- datfr2$Tail[hldvec2,]
            datfr2$Count <- datfr2$Count[hldvec2,]
            
            #if(env$input[[p("clusterby")]] == 4){
            #    datfr$Tail <- datfr$Tail[env$input[[p("choosecol")]]]
            #    datfr2$Count <- datfr$Count[env$input[[p("choosecol")]]]
            #}
            
            #Append names for neatness in graphic
            #colnames(datfr2$Tail) <- paste(colnames(datfr2$Tail), "Tail")
            #colnames(datfr2$Count) <- paste(colnames(datfr2$Count), "Count")
            datfr3 <- list()
            
            hldvec3 <- apply(datfr2$Count,1,function(row) {
                any(row >= (env$input[[p("expmin")]]))
            })
            datfr2$Tail <- datfr2$Tail[hldvec3,]
            datfr2$Count <- datfr2$Count[hldvec3,]
            datfr3$Tail <- datfr2$Tail
            datfr3$Count <- varistran::vst(datfr2$Count)
            
            datfr3$annotate<- datfr$Annotation[rownames(datfr$Annotation) %in% rownames(datfr3$Tail),]
            
            datfr3$Tail_count <- datfr$Tail_count[rownames(datfr$Tail_count) %in% rownames(datfr3$Tail),]
            
            orderedvec <- order(datfr3$annotate$chromosome, datfr3$annotate$start)
            
            datfr3$annotate <- datfr3$annotate[orderedvec,]
            datfr3$Tail <- datfr3$Tail[orderedvec,]
            datfr3$Count <- datfr3$Count[orderedvec,]
            datfr3$Tail_count <- datfr3$Tail_count[orderedvec,]
            datfr3$Tail_count <- datfr3$Tail_count[,colvec,drop=F]
            
            #Truncate names for neatness
           # rownames(datfr3$Count) <- substr(rownames(datfr3$Count), 1, 17)
           # rownames(datfr3$Tail) <- substr(rownames(datfr3$Tail), 1, 17)
           # rownames(datfr3$Tail_count) <- substr(rownames(datfr3$Tail_count), 1, 17)
           # rownames(datfr3$annotate) <- substr(rownames(datfr3$annotate), 1, 17)
           # datfr3$annotate$product <- substr(datfr3$annotate$product , 1, 76)
           # datfr3$annotate$gene <- substr(datfr3$annotate$gene, 1, 11)
            
            datfr3
        })
        
        # Means to extract universe for selection
        env[[p("rows")]] <- reactive({ rownames(wproc()$annotate) })
        
        # Processes the selection of rows by calculating maximum span in either tail length of expression
        selproc <- reactive({
            a1 <- wproc()
            if(env$input[[p("selFeat")]] == 1){
                y <- a1$Tail
                y <- as.matrix(y)
                interest <- apply(y,1,max,na.rm=T) - apply(y,1,min,na.rm=T)
            } else if(env$input[[p("selFeat")]] == 2){
                y <- a1$Count
                y <- as.matrix(y)
                interest <- apply(y,1,max,na.rm=T) - apply(y,1,min,na.rm=T)
            } else {
                interest <- as.matrix(rowMeans(a1$Count, na.rm=TRUE), ncol=1)
            }
            
            n <- env$input[[p("n")]]
            if (n > 2000) stop("Drawing large heatmaps uses excessive system resources. Sorry.")
            
            selection <- rep(FALSE,length(interest))
            selection[ order(-interest)[ seq_len(n) ] ] <- TRUE
            
            if (sum(selection) < 1) stop("No features to show.")
            
            rtVal <- list()
            rtVal$sel <- selection
            rtVal$ann <- a1$annotate[selection,,drop=FALSE]
            rtVal$a1 <- a1
            return(rtVal)
        })        
        
        goDBret <- reactive({
            return(env$output[[p("GOsrch")]])
        })
        hgcret <- reactive({
            return(env$input[[p("fdrcutoff")]])
        })
        otype <- reactive({
            return(env$input[[p("ontype")]])
        })
        # Enables env to hold the data from sh_hmap_detailed
        env[[p("seldat")]] <- selproc
        #env[[p("goDB")]] <- goDBret
        env[[p("hgc")]] <- hgcret
        env[[p("otype")]] <- otype
        
        env$output[[p("debugOut")]] <- shiny::renderText({
            names(env)
        })
        
        env[[p("grob")]] <- reactive({
            
            selrt <- selproc()
            splitDF <- selrt$a1
            selection <- selrt$sel
            
            plot_patseq_heatmap(
                matf1=splitDF$Tail[selection,,drop=FALSE],
                matf2=splitDF$Count[selection,,drop=FALSE],
                gmatf=splitDF$annotate[selection,,drop=FALSE],
                sample_labels=sample_labels(env),
                sample_labels2=sample_labels(env),
                feature_labels=feature_labels(env)[selection],
                clusterby=env$input[[p("clusterby")]],
                row_ord=env$input[[p("roword")]],
                tail_relative=env$input[[p("tail_relative")]]   
            )
        })
        plot$component_server(env)
        if(GOenable == TRUE){
            env$output[[p("gotab")]] <- DT::renderDataTable(env[[p("gotab")]](), 
                server=F
                #extensions = 'TableTools',
                #options = list(searchHighlight = TRUE,
                #               dom = 'T<"clear">lfrtip',
                #               tableTools = list(
                #                   sSwfPath = DT::copySWF(),
                #                   aButtons = list('print'))
                #)
            )
        } else {
            env$output[[p("goerror")]] <- shiny::renderText({
                paste("GO Term analysis disabled")
            })
        }
        
        env$output[[p("dlanalysis")]] <- shiny::downloadHandler(
            paste0("Analysis.csv"),
            function(filename) {
                write.csv(env[[p("gotab")]](), filename)
            }
        )
    }

    composable_shiny_panels_app(panels, server, title="Tail heatmap")
}



