
#
# features must have seqinfo
#

#' @export
read_gff <- function(filename, unlist_parent=TRUE) {
    result <- import(
        filename, 
        sequenceRegionsAsSeqinfo=T,
        colnames=c("ID","Parent","Name"))
    if (unlist_parent) {
        unlisted <- BiocGenerics::unlist(result$Parent)
        assert_that(length(unlisted) == length(result))
        mcols(result)$Parent <- unlisted 
    }
    
    result
}
    

#' Load features needed by plot_genome().
#'
#' @export
load_genome_features <- function(reference_dir) {
    transcripts <- read_gff(file.path(reference_dir,"transcript.gff"))
    exons <- read_gff(file.path(reference_dir,"exon.gff"))
    cds <- read_gff(file.path(reference_dir,"cds.gff"))
    seqinfos <- seqinfo(transcripts)

    reorder <- order(
        seqnames(transcripts) %>% as.character, 
        strand(transcripts) %>% as.character, 
        (start(transcripts)+end(transcripts)) %>% as.integer)
    transcripts <- transcripts[ reorder ]
    
    y <- (seq_along(transcripts)*(2/(1+sqrt(5)))) %% 1
    mcols(transcripts)$y <- y
    mcols(exons)$y <- y[ match(as.character(exons$Parent), as.character(transcripts$ID)) ]
    mcols(cds)$y <- y[ match(as.character(cds$Parent), as.character(transcripts$ID)) ]
    
    transcripts$type <- "transcript"
    transcripts$thickness <- 0.5
    exons$type <- "exon"
    exons$thickness <- 2
    cds$type <- "CDS"
    cds$thickness <- 4

    c(transcripts, exons, cds)
}


#' Produce a depth of coverage plot of a set of samples.
#'
#' @export
plot_genome <- function(features, samples, pos, marks=c()) {
    pos_forward <- is_forward(pos)
    
    ggtemplate <- function(...) {
        template <- ggplot(...) + theme_minimal()
        if (pos_forward) 
            template <- template + 
                scale_x_continuous(limits=c(start(pos),end(pos)), expand=c(0,0), oob=function(a,b) a)
        else
            template <- template + 
                scale_x_reverse(limits=c(end(pos),start(pos)), expand=c(0,0), oob=function(a,b) a)
        template
    }
    
    if (pos_forward) {
        this_samples <- samples %>% mutate(
            cover_same = cover_fwd, cover_opp = cover_rev,
            span_same = span_fwd, span_opp = span_rev)
    } else {
        this_samples <- samples %>% mutate(
            cover_same = cover_rev, cover_opp = cover_fwd,
            span_same = span_rev, span_opp = span_fwd)
    }

    n <- min(1000, width(pos))
    
    lines <- purrr::map_df(seq_len(nrow(samples)), function(i) {
            cs <- rtracklayer::summary(BigWigFile(this_samples$cover_same[i]), pos, size=n, type="max")[[1]]
            co <- rtracklayer::summary(BigWigFile(this_samples$cover_opp[i]), pos, size=n, type="max")[[1]]
            ss <- rtracklayer::summary(BigWigFile(this_samples$span_same[i]), pos, size=n, type="max")[[1]]
            so <- rtracklayer::summary(BigWigFile(this_samples$span_opp[i]), pos, size=n, type="max")[[1]]
            scale <- samples$depth_normalizer[i] #max(1,cs$score,co$score)
            
            data_frame(
                name = this_samples$name[i],
                label = this_samples$label[i],
                cover_same = cs$score / scale,
                cover_opp = co$score / scale,
                span_same = ss$score / scale,
                span_opp = so$score / scale,
                xstart = start(cs)-0.5,
                xend = end(cs)+0.5,
                x = (start(cs)+end(co))*0.5
            )
        }) %>% 
        mutate(name = factor(name, levels=this_samples$name))
        
    alpha <- 0.333
    line_plot <- ggtemplate(lines,aes(x=x,group=name,color=label)) + 
        geom_density(aes(y=span_same), stat="identity", color=NA, fill="#dddddd", alpha=alpha) + 
        geom_density(aes(y=-span_opp), stat="identity", color=NA, fill="#dddddd", alpha=alpha) +
        geom_line(aes(y=cover_same), lwd=0.5) + 
        geom_line(aes(y=-cover_opp), lwd=0.5) +
        #geom_line(aes(y=cover_same), lwd=0.25) + 
        #geom_line(aes(y=-cover_opp), lwd=0.25) +
        geom_hline(yintercept=0) +
        #scale_color_hue(h=c(0,270), l=c(40,65,90,65)) +
        labs(y="", color="") +
        scale_y_continuous(
            limits=c(-max(lines$cover_opp), max(lines$cover_same)), 
            oob=function(a,b) a) +
        theme_minimal()

    if (length(marks))    
        line_plot <- line_plot + geom_vline(xintercept=marks)


    relevant <- subsetByOverlaps(features, pos, ignore.strand=TRUE)
    
    feature_plot <- function(forward, items) {
        p <- ggtemplate() + 
            ylab("") +
            geom_point(data=data_frame(x=c(start(pos),end(pos)),y=c(0,1)),aes(x=x,y=y),color=NA,fill=NA) +
            theme_minimal() +
            scale_y_continuous(limits=c(0,1),breaks=c())
        
        if (length(items)) {
            df_items <- GenomicRanges::as.data.frame(items)
            p <- p + 
                ggplot2::geom_segment(data=df_items, 
                    aes(x=start,xend=end,y=y,yend=y), 
                    lwd=df_items$thickness,lineend="butt")
                    
            df_items2 <- df_items[df_items$type == "transcript",,drop=FALSE]
            df_items2$x <- if (forward == pos_forward) df_items2$end else df_items2$start
            p <- p + geom_text(data=df_items2, 
                aes(x=x, y=y, label=paste0(" ",Name, " ")), 
                hjust=if (forward) 0 else 1,
                vjust=0.5)
        }
        
        p
    }

    
    ggbio::tracks(
         Same = feature_plot(T, relevant[ is_forward(relevant) == pos_forward ]),
         Coverage = line_plot,
         Opposite = feature_plot(F, relevant[ is_forward(relevant) != pos_forward ]),
         heights = c(1,2,1),
         xlim=c(start(pos),end(pos)),
         scale.height = unit(2,"lines")
    )
}




#' @export
shiny_genome_browser_ui <- function(id, location="") {
    ns <- NS(id)
    
    tagList(
        fluidRow(
            column(3, textInput(ns("location"), "Location", location)),
            column(9,
                tags$label("Go"),
                tags$br(),
                actionButton(ns("go_in"), "+"),
                actionButton(ns("go_out"), "-"),
                actionButton(ns("go_left"), HTML("&larr;")),
                actionButton(ns("go_right"), HTML("&rarr;"))
            )
        ),
        
        shiny_plot_ui(ns("plot"), width=1000)
    )
}

#' @export
shiny_genome_browser_server <- function(input, output, session, callback) {
    location <- reactive({
        text <- gsub(",", "", input$location)
        if (text == "") return(NULL)
        as(text, "GRanges")
    })
    
    observeEvent(input$go_in, {
        pos <- location()
        step <- width(pos) %/% 4
        start(pos) <- start(pos) + step
        end(pos) <- end(pos) - step
        updateTextInput(session, "location", value=as.character(pos))
    })

    observeEvent(input$go_out, {
        pos <- location()
        step <- width(pos) %/% 2
        start(pos) <- start(pos) - step
        end(pos) <- end(pos) + step
        updateTextInput(session, "location", value=as.character(pos))
    })
    
    observeEvent(input$go_left, {
        pos <- location()
        pos <- shift(pos, ifelse(is_forward(pos),-1,1) * max(1, width(pos) %/% 3))
        updateTextInput(session, "location", value=as.character(pos))        
    })
    
    observeEvent(input$go_right, {
        pos <- location()
        pos <- shift(pos, ifelse(is_forward(pos),1,-1) * max(1, width(pos) %/% 3))
        updateTextInput(session, "location", value=as.character(pos))        
    })
    
    plotter <- function() {
        pos <- location()
        if (is.null(pos)) return()
        callback(pos)
    }
    
    callModule(shiny_plot_server, "plot", plotter, "browser", session=session)
}

#' @export
shiny_genome_browser <- function(location, callback, prefix="") {
    composable_shiny_app(
        div( shiny_genome_browser_ui(prefix, location=location) ),
        function(env) {
            callModule(shiny_genome_browser_server, prefix, callback, session=env$session)
        }
    )
}




