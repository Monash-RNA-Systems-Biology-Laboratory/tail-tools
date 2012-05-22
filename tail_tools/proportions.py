
import nesoni
from nesoni import config, runr

@config.help(
'Create a file of proportions of poly-A reads.'
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
'Plot heatmap of proportions of reads with poly-A tails.',
)
@config.Int_flag('min_min', 'For a gene to be included, all samples must have this many reads aligning.')
@config.Int_flag('min_max', 'For a gene to be included, at least one sample must have this many reads aligning.')
@config.Float_flag('min_fold', 
    'For a gene to be included the ratio of poly-A:non-poly-A reads in one sample '
    'must be this many times larger than the ratio in some other sample.')
@config.Positional('all', 'Output from "nesoni count:" using all reads.')
@config.Positional('polya', 'Output from "nesoni count:" using only poly-A reads.')
class Proportions_heatmap(runr.R_action, config.Action_with_prefix):
    all = None
    polya = None
    min_min = 10
    min_max = 50
    min_fold = 2.0

    script = """
        library(nesoni)
        
        all.dgelist <- read.counts(all)
        polya.dgelist <- read.counts(polya)
        
        all <- all.dgelist$counts
        polya <- polya.dgelist$counts
        
        good <- row.apply(all, function (x) (min(x) >= min_min && max(x) >= min_max))
        all.good <- all[good,]
        polya.good <- polya[good,]
        prop.good <- polya.good / all.good     
        colnames(prop.good) <- colnames(all.good)
        
        nonpolya.good <- all.good - polya.good
        ratio.good <- polya.good / nonpolya.good        
        interesting <- row.apply(ratio.good, function(x) 
            (min(x) * min_fold <= max(x) )
        )
        
        heatmap.matrix <- prop.good[interesting,]
        
        heatmap.matrix.renamed <- heatmap.matrix
        
        annotation <- c('gene', 'product')
        for(colname in annotation)
            if (!all(is.na(all.dgelist$gene[,colname])))
                rownames(heatmap.matrix.renamed) <- paste(
                    rownames(heatmap.matrix.renamed), 
                    all.dgelist$gene[rownames(heatmap.matrix),colname]
                )

        res <- 150
        row.labels <- nrow(heatmap.matrix) <= 500        
        height <- if(row.labels) (16*nrow(heatmap.matrix)+800)*res/150 else 2500*res/150    
        png(sprintf('%s.png',prefix), width=1500*res/150, height=height, res=res)
        
        heatmap <- nesoni.heatmap(heatmap.matrix.renamed, 
            col=unsigned.col,
            margins=c(20,20),
            labRow=(if(row.labels) NULL else NA),
            main='Proportion of reads with poly-A tail'
        )
        
        dev.off()


        shuffled.matrix <- heatmap.matrix[rev(heatmap$rowInd),]
        
        table.filename <- sprintf('%s.csv', prefix)
        
        sink(table.filename)
        cat('# poly-A proportion data\n')
        cat('#\n')
        cat('# Values given are proportion of reads with poly-A tail\n')
        cat('#\n')
        cat(sprintf('# %d genes shown\n', nrow(data)))
        cat('#\n')
        
        frame <- data.frame(name=rownames(shuffled.matrix), row.names=rownames(shuffled.matrix), check.names=FALSE) 
        for(colname in annotation)
            if (!all(is.na(all.dgelist$gene[,colname])))
                frame[,colname] <- all.dgelist$gene[rownames(frame),colname]
             
        frame[,'cluster hierarchy'] <- rev(dendrogram.paths(heatmap$rowDendrogram))
        frame <- data.frame(frame, shuffled.matrix, check.names=FALSE)
        
        write.csv(frame, row.names=FALSE)
        sink()


        
    """


