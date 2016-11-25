#' @title Integrates heatmap into shiny
#' 
#' @details 
#' Shiny wrapper for sh_hmap_detailed
#' Uses gridBase to produce brush for interactivity with plot
#' 
#' @param callback Heatmap Grob from plot_patseq_heatmap
#' @param width Default width of the heatmap grob when the shiny app loads
#' @param height Default height of the heatmap grob when the shiny app loads
#' @param dlname Default filename when downloading the heatmap image
#' @param prefix Prefix
#' @param selin Selection of genes from sh_hmap_detailed
#' @param rorder Row order from sh_hmap_detailed
#' @param spp Species of the data (Hs, Sc, Ce, Mm)
#' @param goenabl Whether GO Term analysis is enabled. TRUE by default. This might need changing in future
#' 
#' @export

shiny_patseq_heatmap_inner <- function(callback, width=500, height=500, dlname="plot", prefix="", selin, rorder,spp, hgc=0.05, otype=1, goenabl=TRUE) {
    selin <- ensure_reactable(selin)
    rorder <- ensure_reactable(rorder)
    hgc <- ensure_reactable(hgc)
    otype <- ensure_reactable(otype)
    #gosrch <- ensure_reactable(gosrch)
    
    p <- NS(prefix)
    
    # Shiny's UI layout
    ui <- function(request)
      shiny::tags$div(
        shiny::fluidRow(
            shiny::column(3, shiny::numericInput(p("width"), "Plot width", width, min=100, max=10000, step=50)),
            shiny::column(3, shiny::numericInput(p("height"), "Plot height", height, min=100, max=10000, step=50)),
            shiny::column(4, shiny::tags$label("Download Plot"), shiny::tags$br(),
                          shiny::downloadButton(p("pdf"), "PDF"),
                          shiny::downloadButton(p("eps"), "EPS")),
            shiny::column(3, shiny::tags$label("Download Selection as .csv"), shiny::tags$br(),
                          shiny::downloadButton(p("dlsel"), ".csv"))
        ),
        shiny::uiOutput(p("plotui"), width="auto", height="auto"),
        DT::dataTableOutput(p('datab'))
      )
    
    # Shiny's server 
    server <- function(env) { 
        output <- env$output
        
        
        i <- function(name) env$input
        
        # Renders plot and enables a brush object
        output[[p("plotui")]] <- shiny::renderUI({
            shiny::plotOutput(p("plot"),
                              brush = brushOpts(
                                  id =p("plot_brush"),  
                                  direction = 'y',
                                  resetOnNew = TRUE,
                                  clip = TRUE,
                                  delay = 600000
                              ),
                              width=env$input[[p("width")]],
                              height=env$input[[p("height")]]
            )
        })
        
        gosrch <- reactive({
            library(GOstats)
        
            selrows <- calcdt()
            allrows <- selin(env)
            if(nrow(selrows)==nrow(allrows) || nrow(selrows)==0)
                stop("Inappropriate selection for GO term analysis (One set is of size 0)")
            
            if(spp == "Sc"){
                library(org.Sc.sgd.db)
                annlabel <- "org.Sc.sgd.db"
                lib <- org.Sc.sgd.db
            }else if(spp == "Hs"){
                library(org.Hs.eg.db)
                annlabel <- "org.Hs.eg.db"
                lib <- org.Hs.eg.db
            }else if(spp == "Ce"){
                library(org.Ce.eg.db)
                annlabel <- "org.Ce.eg.db"
                lib <- org.Ce.eg.db
            }else if(spp == "Mm"){
                library(org.Mm.eg.db)
                annlabel <- "org.Mm.eg.db"
                lib <- org.Mm.eg.db
            }
            
            colkp <- c("GENENAME", "ENTREZID")
            remrows <- rownames(allrows[-which(rownames(selrows) %in% rownames(allrows)),])
            selrowsname <- rownames(selrows)
            unisel <- select(lib, keys=remrows, columns=colkp, keytype="ENSEMBL")
            intsel <- select(lib, keys=selrowsname, columns=colkp, keytype="ENSEMBL")
            #Deduplicate rows
            unisel <- unisel[!duplicated(unisel$ENTREZID),]
            intsel <- intsel[!duplicated(intsel$ENTREZID),]
            #Remove rows where ETREZID is NA
            unisel <- unisel[!is.na(unisel$ENTREZID),]
            intsel <- intsel[!is.na(intsel$ENTREZID),]
            
            if(otype(env) == 1){
                ot = "BP"
            } else if (otype(env)==2){
                ot = "CC"
            } else {
                ot = "MF"
            }
            
            #Setup hyperG or fischer's exact test params
            gc <- hgc(env)
            hgparam <- new("GOHyperGParams",annotation=annlabel,geneIds=intsel,universeGeneIds=unisel,ontology=ot, testDirection="over",pvalueCutoff=gc)
            hg <- hyperGTest(hgparam)
            hg.pv <- pvalues(hg)
            hgadjpv <- p.adjust(hg.pv,'fdr')
            sigGOID <- names(hgadjpv[hgadjpv < gc])
            df <- GOstats::summary(hg)   
            #df <- summary(hg)
            df[,3:4] <- format(signif(df[,3:4], digits=6), format="fg")
            df$Pvalue <- format(df$Pvalue, scientific=T)
            return(df)
        })
        if(goenabl){
            env[[p("gotab")]] <- gosrch
        } else {
            env[[p("gotab")]] <- data.frame()
        }
        
        # Works out if any rows are selected and returns selected rows
        # Calculates rows based on y-coordinates from brush output
        ## If no rows are selected, outputs all rows shown in the heatmap
        calcdt <- reactive({
            numrows <- isolate( nrow(selin(env)) )
            
            ytop <- round((env$input[[p("plot_brush")]]$ymax + 2)/(54/numrows))
            ybot <- round((env$input[[p("plot_brush")]]$ymin + 2)/(54/numrows)+1)
            if(length(ytop) == 0){
                return( isolate( selin(env) )[c(),])
                #return(selin(env)[rorder(env),])
            } else {
                sel <- isolate( selin(env)[rorder(env),] )
                sel <- sel[(ytop):ybot,]
                return(sel)
            }
        })
                
        # Means to extract selected rows
        env[[p("rows-selected")]] <- reactive({ rownames(calcdt()) })
        
        # Data table output with selected rows
        output[[p("datab")]] <- DT::renderDataTable(
            calcdt(), 
            server=F #,
            #extensions = 'TableTools',
            #options = list(searchHighlight = TRUE,
            #               dom = 'T<"clear">lfrtip',
            #                   tableTools = list(
            #                       sSwfPath = DT::copySWF(),
            #                       aButtons = list('print'))
            #                                                       )
            #                                        
        )

        # Produces plot output
        output[[p("plot")]] <- shiny::renderPlot(withProgress(message="Plotting", { 
            vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
            plot.new()
            callback(env)
            seekViewport("prodVP")
            pltvec <- gridPLT()
            pltvec <- c(0, 1, pltvec[3], pltvec[4])
            par(new=T, plt=pltvec)
            plot(1,type="n", axes=F, xlab ="",ylab="",xlim=c(0,50),ylim=c(0,50))
            
            popViewport()            
        }))
        
        # Download handlers---
        output[[p("pdf")]] <- shiny::downloadHandler(
            paste0(dlname,".pdf"),
            function(filename) {
                pdf(filename, width=env$input[[p("width")]]/72, height=env$input[[p("height")]]/72)
                callback(env)
                dev.off()
            }
        )
        output[[p("eps")]] <- shiny::downloadHandler(
            paste0(dlname,".eps"),
            function(filename) {
                postscript(filename, width=env$input[[p("width")]]/72, height=env$input[[p("height")]]/72,
                           paper="special", onefile=FALSE, horizontal=FALSE)
                callback(env)
                dev.off()
            }
        )
        output[[p("dlsel")]] <- shiny::downloadHandler(
            paste0("Selection.csv"),
            function(filename) {
                write.csv(calcdt(), filename)
            }
        )
        
        #---
    }
    
    composable_shiny_app(ui, server)
}