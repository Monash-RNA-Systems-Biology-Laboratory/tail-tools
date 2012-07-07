
import nesoni
from nesoni import config, runr

@config.help(
'(Deprecated) Create a file of proportions of poly-A reads.'
)
@config.String_flag('norm_file', 'Normalization file.')
@config.Positional('all', 'Output from "nesoni count:" using all reads.')
@config.Positional('polya', 'Output from "nesoni count:" using only poly-A reads.')
class Proportions(runr.R_action, config.Action_with_prefix):
    norm_file = None
    all = None
    polya = None

    script = """
        library(nesoni)
        
        all.dgelist <- read.counts(all, norm.file=norm_file)
        polya.dgelist <- read.counts(polya)
        
        all <- all.dgelist$counts
        polya <- polya.dgelist$counts
        proportions <- polya / all
        colnames(proportions) <- colnames(all)
        
        result <- data.frame(
            Feature = rownames(all.dgelist$genes),
            check.names = FALSE
        )
        
        for(name in colnames(all.dgelist$genes))
            if (!all(is.na(all.dgelist$genes[,name])))
                result[,name] <- all.dgelist$genes[,name]
        
        for(name in colnames(proportions))
            result[,paste(name,'prop')] <- proportions[,name]

        for(name in colnames(all))
            result[,paste(name,'norm')] <- all[,name] * all.dgelist$samples[name,'normalizing.multiplier']
        
        sink(sprintf("%s.csv", prefix))
        
        cat('# prop columns give the proportion of reads with poly-A tails\n')
        cat('# norm columns give the normalized number of reads\n') 
        cat('#\n')        
        write.csv(result, row.names=FALSE)
        
        sink()        
    """
 

@config.help(
'(Deprecated) Plot heatmap of proportions of reads with poly(A) tails.',
)
@config.Int_flag('min_min', 'For a gene to be included, all samples must have this many reads aligning.')
@config.Int_flag('min_max', 'For a gene to be included, at least one sample must have this many reads aligning.')
@config.Float_flag('min_diff', 
    'For a gene to be included the proportion of poly(A) reads in one sample '
    'must be this many times larger than the ratio in some other sample.')
@config.String_flag('norm_file', 'Normalization file.')
@config.Positional('all', 'Output from "nesoni count:" using all reads.')
@config.Positional('polya', 'Output from "nesoni count:" using only poly(A) reads.')
class Proportions_heatmap(runr.R_action, config.Action_with_prefix):
    norm_file = None
    all = None
    polya = None
    min_min = 50
    min_max = 50
    min_diff = 0.25

    script = r"""
    library(nesoni)
    
    
    draw.prop.heatmap <- function(prefix, dend, expr, prop, labels) {
        n.rows <- nrow(expr)
        n.cols <- ncol(expr)
        
        res <- 150
        row.labels <- n.rows <= 300        
        height <- if(row.labels) (25*n.rows+600)*res/150 else 2500*res/150    
        png(sprintf('%s.png',prefix), width=2000*res/150, height=height, res=res)
        
        if (n.rows < 2 || n.cols < 2) {
            plot.new()
            title(sprintf('Can\'t plot %d x %d heatmap', n.rows, n.cols))
            dev.off()
            return()
        }
        
        y1 <- 4 / par()$fin[2]
        y2 <- 0.99
        
        dend.x1 <- 0
        dend.x2 <- 1/9
        
        prop.x1 <- 1/9
        prop.x2 <- 2/9 -1/100
        
        expr.x1 <- 2/9
        expr.x2 <- 3/9 -1/100
        
        legend.y2 <- y1*0.75
        legend.y1 <- legend.y2 - 0.03*aspect.ratio()
        
        prop.legend.x1 <- 11/30
        prop.legend.x2 <- 14/30
        expr.legend.x1 <- 15/30
        expr.legend.x2 <- 18/30
        
        expr.col <- signed.col
        expr.extreme <- max(0.0,abs(expr),na.rm=TRUE)
        expr.breaks <- seq(-expr.extreme,expr.extreme, length=length(expr.col)+1)
        
        prop.col <- unsigned.col
        prop.breaks <- seq(0.0, 1.0, length=length(prop.col)+1)
        
        plot.new()
        
        if (!all(is.na(dend$dendrogram))) {
            put.plot(dend.x1,dend.x2, y1,y2)
            plot(dend$dendrogram, horiz=TRUE, axes=FALSE, yaxs="i", leaflab="none")
        }
        
        put.plot(prop.x1,prop.x2, y1,y2)
        basic.image(1:n.cols,1:n.rows, t(prop[dend$order,]), col=prop.col, breaks=prop.breaks)        
        axis(1, at=1:n.cols, labels=colnames(prop), las=2, tick=FALSE)
        
        put.plot(expr.x1,expr.x2, y1,y2)
        basic.image(1:n.cols,1:n.rows, t(expr[dend$order,]), col=expr.col, breaks=expr.breaks)
        axis(1, at=1:n.cols, labels=colnames(expr), las=2, tick=FALSE)
        
        if (row.labels) {
            line <- 0
            for(i in basic.seq(length(labels))) {
                if (i < length(labels))
                    l <- trim.labels(labels[[i]])
                else
                    l <- as.character(labels[[i]])
                axis(4, at=1:n.rows, labels=l[dend$order], las=2, tick=FALSE, line=line)
                line <- line + 0.5 * max(0,nchar(l))
            }
        }
        
        put.plot(expr.legend.x1,expr.legend.x2, legend.y1,legend.y2)
        color.legend(expr.col, expr.breaks, 'log2 expression\ndifference from row average\n')
        
        put.plot(prop.legend.x1,prop.legend.x2, legend.y1,legend.y2)
        color.legend(prop.col, prop.breaks, 'proportion of reads\nwith poly(A) tail\n')
        
        dev.off()
    }
    
    all.dgelist <- read.counts(all, norm.file=norm_file)
    polya.dgelist <- read.counts(polya)
    
    all <- all.dgelist$counts
    polya <- polya.dgelist$counts
    proportions <- polya / all
    colnames(proportions) <- colnames(all)
    
    good <- (
        row.apply(all,min) >= min_min & 
        row.apply(all,max) >= min_max &
        row.apply(proportions,max) - row.apply(proportions,min) >= min_diff 
    )  
    all.good <- basic.subset(all,good)
    prop.good <- basic.subset(proportions,good)
    
    expression <- voom(all.dgelist)$E
    expression.good <- basic.subset(expression,good)
    expression.good.centered <- t(scale(t(expression.good), center=TRUE,scale=FALSE))
    
    dend <- do.dendrogram( t(scale(t(prop.good))) )
             
    label.list <- list( rownames(all.good) )
    annotation <- c('gene', 'product')
    for(colname in annotation)
        if (!all(is.na(all.dgelist$gene[,colname])))
            label.list[[ length(label.list)+1 ]] <- all.dgelist$gene[rownames(all.good),colname]
    
    draw.prop.heatmap(prefix, dend, expression.good.centered, prop.good, label.list)
    
    #res <- 150
    #row.labels <- nrow(heatmap.matrix) <= 500        
    #height <- if(row.labels) (16*nrow(heatmap.matrix)+800)*res/150 else 2500*res/150    
    #png(sprintf('%s.png',prefix), width=1500*res/150, height=height, res=res)
    #
    #heatmap <- nesoni.heatmap(heatmap.matrix.renamed, 
    #    col=unsigned.col,
    #    margins=c(20,20),
    #    labRow=(if(row.labels) NULL else NA),
    #    main='Proportion of reads with poly-A tail'
    #)
    #
    #dev.off()
    
    
    table.filename <- sprintf('%s.csv', prefix)
    
    sink(table.filename)
    cat('# prop columns give the proportion of reads with poly-A tails\n')
    cat('# norm columns give the normalized number of reads\n') 
    cat('#\n')
    cat(sprintf('# %d genes shown\n', nrow(data)))
    cat('#\n')
    
    frame <- data.frame(name=rownames(all.good), row.names=rownames(all.good), check.names=FALSE) 
    for(colname in annotation)
        if (!all(is.na(all.dgelist$gene[,colname])))
            frame[,colname] <- all.dgelist$gene[rownames(frame),colname]
    
    frame[dend$order,'cluster hierarchy'] <- dend$paths

    for(name in colnames(prop.good))
        frame[,paste(name,'prop')] <- prop.good[,name]

    for(name in colnames(all.good))
        frame[,paste(name,'norm')] <- all.good[,name] * all.dgelist$samples[name,'normalizing.multiplier']
    
    write.csv(frame[rev(dend$order),], row.names=FALSE)
    sink()
    """


