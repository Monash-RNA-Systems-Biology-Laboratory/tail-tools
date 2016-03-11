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
        cat("Species is missing or incorrect. Not enabling GO Term Analysis\n")
    } else {
        GOenable <- TRUE
    }
    
    p <- function(name) paste0(prefix,name)
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
        width=1250,
        height=900,
        dlname="heatmap",
        prefix=p("plot_"),
        goenabl=GOenable,
        selin = function(env){
            env$seldat()$ann
        },
        spp=species,
        hgc = function(env){
            env$hgc()
        },
        otype = function(env){
            env$otype()
        }
    )
    
    # Shiny's UI layout 
    ui <- shiny::tags$div(
        shiny::titlePanel("Tail heatmap"),
        shiny::tabsetPanel(
            shiny::tabPanel("Options",
                            shiny::fluidRow(
                                shiny::column(3,
                                              shiny::p("Features are selected based on span of:"),
                                              shiny::radioButtons(p("selFeat"), 
                                                                  label="Select features by:", 
                                                                  choices=list("Tail length"=1, "Expression"=2), 
                                                                  selected=1,
                                                                  inline=TRUE),
                                              shiny::uiOutput(p("chrs"))
                                ),
                                shiny::column(3,
                                              shiny::numericInput(p("n"), "Number of features to show", 50, min=10,max=2000,step=10),
                                              shiny::numericInput(p("nmin"), "Trim Tail Counts below value to NA", 50, min=0,max=1000,step=1),
                                              shiny::numericInput(p("expmin"), "Exclude genes with expression counts below: ", 0, min=0,max=1500,step=1),
                                              shiny::radioButtons(p("roword"), 
                                                                  label="Features ordered by: ", 
                                                                  choices=list("Tail Length"=1, "Expression"=2, "Group by location"=3), 
                                                                  selected=1,
                                                                  inline=TRUE),
                                              shiny::radioButtons(p("clusterby"), 
                                                                  label="Order samples by: ", 
                                                                  choices=list("Tail length" = 2, "Expression" = 3, "Manually order from select columns" = 1), 
                                                                  selected = 1,
                                                                  inline=TRUE)
                                ),
                                shiny::column(3,
                                              shiny::uiOutput(p("selCol"))
                                              #shinyURL::shinyURL.ui()
                                )
                                
                            )),
            shiny::tabPanel("Plot",plot$component_ui),
            shiny::tabPanel("GO analysis of selection",
                            shiny::fluidRow(
                                shiny::column(3, shiny::numericInput(p("fdrcutoff"), "Hyper geometric test cutoff", 0.05, min=0.0001, max=1, step=0.01)),
                                shiny::column(3, shiny::radioButtons(p("ontype"), label="Ontology search type",
                                                                     choices=list("Biological Processes"=1, "Cellular Component"=2,"Molecular Function"=3),
                                                                     selected=1,
                                                                     inline=TRUE)),
                                shiny::column(3, shiny::tags$label("Download analysis as .csv"), shiny::tags$br(),
                                              shiny::downloadButton(p("dlanalysis"), ".csv"))
                            ),
                            shiny::textOutput(p("goerror")),
                            DT::dataTableOutput(p("gotab")))
                            #shiny::tableOutput("tabout"))
        ),
        parenthetically("This plot is produced by tailtools::shiny_patseq_heatmap.")
        
    )
    
    # Shiny's server
    server <- function(env) {
        # Nice idea, doesn't work well with DT or brushes
        #shinyURL::shinyURL.server(env$session)
        
        # Processes the input from datfr into the correct length and returns a list of 4 data frames
        wproc <- reactive({
            
            datfr2 <- list()
            
            colvec <- -which(names(datfr$Tail) %in% env$input[[p("choosecol")]])
            
            datfr$Tail[datfr$Tail_count < env$input[[p("nmin")]]] = NA
            
            datfr2$Tail <- datfr$Tail[,-colvec]
            datfr2$Count <- datfr$Count[,-colvec]
            
            hldvec2 <- apply(datfr2$Tail,1,function(row) {
                any(sum(!is.na(row))>= 2)
            })
            
            datfr2$Tail <- datfr2$Tail[hldvec2,]
            datfr2$Count <- datfr2$Count[hldvec2,]
            
            if(env$input[[p("clusterby")]] == 4){
                datfr$Tail <- datfr$Tail[env$input[[p("choosecol")]]]
                datfr2$Count <- datfr$Count[env$input[[p("choosecol")]]]
            }
            
            #Append names for neatness in graphic
            colnames(datfr2$Tail) <- paste(colnames(datfr2$Tail), "Tail")
            colnames(datfr2$Count) <- paste(colnames(datfr2$Count), "Count")
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
            datfr3$Tail_count <- datfr3$Tail_count[,-colvec]
            
            #Sort into order
            cvec <- datfr3$annotate$chromosome
            cvec <- cvec %in% env$input[[p("choosechr")]]
            
            datfr3$annotate <- datfr3$annotate[cvec,]
            datfr3$Tail <- datfr3$Tail[cvec,]
            datfr3$Tail_count <- datfr3$Tail_count[cvec,]
            datfr3$Count <- datfr3$Count[cvec,]
            
            #Truncate names for neatness
           # rownames(datfr3$Count) <- substr(rownames(datfr3$Count), 1, 17)
           # rownames(datfr3$Tail) <- substr(rownames(datfr3$Tail), 1, 17)
           # rownames(datfr3$Tail_count) <- substr(rownames(datfr3$Tail_count), 1, 17)
           # rownames(datfr3$annotate) <- substr(rownames(datfr3$annotate), 1, 17)
           # datfr3$annotate$product <- substr(datfr3$annotate$product , 1, 76)
           # datfr3$annotate$gene <- substr(datfr3$annotate$gene, 1, 11)
            
            return(datfr3)
        })
        
        # Processes the selection of rows by calculating maximum span in either tail length of expression
        selproc <- reactive({
            a1 <- wproc()
            if(env$input[[p("selFeat")]] == 1){
                y <- a1$Tail
            } else if(env$input[[p("selFeat")]] == 2){
                y <- a1$Count
            }
            y <- ensure_reactable(y)
            
            n <- env$input[[p("n")]]
            if (n > 2000) stop("Drawing large heatmaps uses excessive system resources. Sorry.")
            
            y_val <- as.matrix(y(env))
            y_span <- apply(y_val,1,max,na.rm=T) - apply(y_val,1,min,na.rm=T)
            selection <- rep(FALSE,nrow(y_val))
            selection[ order(-y_span)[ seq_len(n) ] ] <- TRUE
            
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
        env$seldat <- selproc
        env$goDB <- goDBret
        env$hgc <- hgcret
        env$otype <- otype
        
        env$output[[p("debugOut")]] <- shiny::renderText({
            names(env)
        })
        # RenderUI output for shiny, dynamically generate two selctize elements ---
        env$output[[p("chrs")]] <- shiny::renderUI({
            tmpvec <- levels(datfr$Annotation$chromosome)
            selectizeInput("choosechr", "Choose chromosomes to display",multiple=T,tmpvec, selected=tmpvec)
        })
        env$output[[p("selCol")]] <- shiny::renderUI({
            colvec <- names(datfr$Tail)
            selectizeInput("choosecol", "Choose samples to display",multiple=T,colvec, selected=colvec)
        })
        #---
        
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
                row_ord=env$input[[p("roword")]]   
            )
        })
        plot$component_server(env)
        if(GOenable == TRUE){
            env$output[[p("gotab")]] <- DT::renderDataTable(env$gotab(), 
                                        server=F,
                                        extensions = 'TableTools',
                                        options = list(searchHighlight = TRUE,
                                                       dom = 'T<"clear">lfrtip',
                                                       tableTools = list(
                                                           sSwfPath = DT::copySWF(),
                                                           aButtons = list('print'))
                                        )
            )
        } else {
            env$output[[p("goerror")]] <- shiny::renderText({
                paste("GO Term analysis disabled")
            })
        }
        
        env$output[[p("dlanalysis")]] <- shiny::downloadHandler(
            paste0("Analysis.csv"),
            function(filename) {
                write.csv(env$gotab(), filename)
            }
        )
    }
    composable_shiny_app(ui, server)
}
