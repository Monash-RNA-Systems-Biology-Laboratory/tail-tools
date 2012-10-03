
import itertools, collections, math

import nesoni
from nesoni import annotation, sam, span_index, config, working_directory, io, runr

@config.help(
'Create file to be used by "aggregate-tail-lengths:".'
"""
""")
@config.String_flag('types', 'Comma separated list of feature types to use.')
class Tail_lengths(config.Action_with_working_dir):
     types = 'gene'

     _workspace_class = working_directory.Working

     def run(self):
         workspace = self.get_workspace()
         reference = workspace.get_reference()
         
         types = [ item.lower() for item in self.types.split(',') ]
     
         annotations = [ 
             item 
             for item in annotation.read_annotations(reference.annotations_filename())
             if item.type.lower() in types
         ]
         
         self.log.log('%d annotations\n' % len(annotations))
         
         index = { }
         
         for item in annotations:
             if item.seqid not in index:
                 index[item.seqid] = span_index.Span_index()
             index[item.seqid].insert(item)
             
             item.hits = [] # [ (rel_start, rel_end, tail_length) ]
         
         for item in index.itervalues(): 
             item.prepare()
         
         for read_name, fragment_alignments, unmapped in sam.bam_iter_fragments(workspace/'alignments_filtered.bam'):
             for fragment in fragment_alignments:
                 start = min(item.pos-1 for item in fragment)
                 end = max(item.pos+item.length-1 for item in fragment)
                 alignment_length = end-start
                 strand = -1 if fragment[0].flag&sam.FLAG_REVERSE else 1
                 
                 tail_length = 0
                 for item in fragment[0].extra:
                     if item.startswith('AN:i:'):
                         tail_length = int(item[5:])
                 
                 if fragment[0].rname in index:
                     for gene in index[fragment[0].rname].get(start,end):
                         if gene.strand != strand: continue
                         
                         if strand > 0:
                             rel_start = start - gene.start
                             rel_end = end - gene.start
                         else:
                             rel_start = gene.end - end
                             rel_end = gene.end - start
                         
                         gene.hits.append( (rel_start,rel_end,tail_length) )
                                  
         workspace.set_object(annotations, 'tail-lengths.pickle.gz')


@config.help(
'Aggregate data collected by "tail-lengths:" and produce various CSV tables.',
"""\
"""
)
@config.Int_flag('saturation', 
     'Duplicate start position saturation level. '
     'Reads that start at the same position will only '
     'count as up to this many. Zero for no staturation.'
)
@config.Main_section('working_dirs')
class Aggregate_tail_lengths(config.Action_with_prefix):         
    saturation = 1
    working_dirs = [ ]
         
    def run(self):
        workspaces = [ 
            working_directory.Working(item) 
            for item in self.working_dirs 
        ]
        
        data = [ 
            item.get_object('tail-lengths.pickle.gz') 
            for item in workspaces
        ]        
        
        names = [ item.name for item in workspaces ]
        
        annotations = data[0]
        
        max_length = max( 
            tail_length
            for sample in data
            for feature in sample
            for rel_start,rel_end,tail_length in feature.hits
        )+1
        
        for i, sample in enumerate(data):
            n_alignments = 0
            n_duplicates = 0
            n_good = 0
            for feature in sample:
                feature.tail_counts = [ 0.0 ] * max_length
                
                buckets = collections.defaultdict(list)
                for rel_start,rel_end,tail_length in feature.hits:
                    buckets[ (rel_start,rel_end) ].append(tail_length)
                for item in buckets.values():
                    n_alignments += len(item)
                    n_good += 1
                    if self.saturation < 1 or len(item) <= self.saturation:
                        weight = 1.0
                    else:
                        weight = float(self.saturation) / len(item)
                        n_duplicates += len(item)
                    for item2 in item:
                        feature.tail_counts[item2] += weight

            self.log.datum(workspaces[i].name, 'Alignments to genes', n_alignments)
            if self.saturation >= 1:
                self.log.datum(workspaces[i].name, 'Proportion of alignments with duplicate start and end position', float(n_duplicates)/max(1,n_alignments))
                self.log.datum(workspaces[i].name, 'Alignments to genes after deduplication', n_good)
                
        
        counts = [ ]
        
        for item in data: 
            assert len(item) == len(data[0])
        for row in itertools.izip(*data):
            this_counts = [ item.tail_counts for item in row ]
            counts.append(this_counts)
        
        #max_length = max(max(len(item) for item in row) for row in counts)
        #
        #for row in counts:
        #    for item in row:
        #        while len(item) < max_length:
        #            item.append(0)
        
        
        

        def counts_iter():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    row[('Count',names[j])] = '%d' % sum(counts[i][j]) 

                row[('Annotation','Length')] = annotations[i].end - annotations[i].start
                row[('Annotation','gene')] = annotations[i].attr.get('gene','')
                row[('Annotation','product')] = annotations[i].attr.get('product','')
                #row[('Annotation','Strand')] = str(annotations[i].strand)
                yield row
        io.write_csv(self.prefix + '-counts.csv', counts_iter())

        def raw_columns():
            for i in xrange(len(names)):
                row = collections.OrderedDict()
                row['Sample'] = names[i]
                for j in xrange(max_length):
                    row['length-%d' % j] = str(i*max_length+j+1) #For R+, so 1 based
                yield row
        io.write_csv(self.prefix + '-raw-columns.csv', raw_columns())

        #Somewhat inefficient        
        def raw():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    for k in xrange(max_length):
                        row['%d %s' % (k,names[j])] = str( counts[i][j][k] )
                yield row
        io.write_csv(self.prefix + '-raw.csv', raw())
        
        def pooled():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(max_length):
                    row[str(j)] = str( sum( counts[i][k][j] for k in xrange(len(names)) ) )
                yield row
        io.write_csv(self.prefix + '-pooled.csv', pooled())

        #Note: zero if zero count
        def mean_length():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    a = 0
                    b = 0
                    for pos,count in enumerate(counts[i][j]):
                        a += pos*count
                        b += count
                    row[names[j]] = str( float(a)/max(1,b) )
                yield row
        io.write_csv(self.prefix + '-mean-length.csv', mean_length())

        def mean4_length():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    a = 0
                    b = 0
                    for pos,count in enumerate(counts[i][j]):
                        if pos < 4: continue
                        a += pos*count
                        b += count
                    row[names[j]] = str( float(a)/max(1,b) )
                yield row
        io.write_csv(self.prefix + '-mean-length-beyond4.csv', mean4_length())


@config.help(
    'Create a spreadsheet containing various statistics, using the output of "aggregate-tail-lengths:"',
    'The prefix should be the same as that given to "aggregate-tail-lengths:".'
)
class Tail_stats(config.Action_with_prefix, runr.R_action):
    def log_filename(self):
        if self.prefix is None: return None
        return self.prefix + '_tail_stats_log.txt'
        
    script = r"""
    
    library(nesoni)

    # === Load data ==
    
    data <- read.grouped.table(sprintf('%s-counts.csv',prefix),require=c('Count','Annotation'))
    annotation <- data$Annotation
    expression <- voom( as.matrix(data$Count), normalize.method='quantile' )$E
    
    raw <- as.matrix(read.grouped.table(sprintf('%s-raw.csv',prefix))$All)
    
    cols <- as.matrix(read.grouped.table(sprintf('%s-raw-columns.csv',prefix))$All)   
    
    row.names <- rownames(annotation)
    n.samples <- nrow(cols)
    n.genes <- nrow(raw)
    n.lengths <- ncol(cols)
    lengths <- basic.seq(n.lengths)-1


    # === Calculate mean tail lengths for each sample ===
    
    # Omit columns for tail length 0,1,2,3    
    skip <- 4
    good.cols <- cols[,(skip+1):n.lengths,drop=FALSE]
    good.lengths <- lengths[(skip+1):n.lengths]
    
    totals <- matrix(nrow=n.genes,ncol=n.samples) 
    x.totals <- matrix(nrow=n.genes,ncol=n.samples) 
    xx.totals <- matrix(nrow=n.genes,ncol=n.samples) 
    means <- matrix(nrow=n.genes,ncol=n.samples)
    sds <- matrix(nrow=n.genes,ncol=n.samples)
    colnames(totals) <- rownames(cols)
    rownames(totals) <- row.names
    colnames(means) <- rownames(cols)
    rownames(means) <- row.names
    colnames(sds) <- rownames(cols)
    rownames(sds) <- row.names
    for(c in basic.seq(n.samples)) {
        totals[,c] <- rowSums(raw[,good.cols[c,],drop=FALSE])
        x.totals[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * good.lengths) )
        xx.totals[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * (good.lengths*good.lengths) ) )

        means[,c] <- x.totals[,c] / totals[,c]
        sds[,c] <- sqrt( totals[,c]/(totals[,c]-1) * (xx.totals[,c]/totals[,c] - means[,c]*means[,c]) )
        
        means[totals[,c] < 1,c] <- NA
        sds[totals[,c] < 2,c] <- NA
    }


    grand.totals <- rowSums(totals)
    grand.x.totals <- rowSums(x.totals)
    grand.xx.totals <- rowSums(xx.totals)
    grand.means <- grand.x.totals / grand.totals
    grand.sds <- sqrt( grand.totals/(grand.totals-1) * (grand.xx.totals/grand.totals - grand.means*grand.means) )
    grand.means[ grand.totals < 1 ] <- NA
    grand.sds[ grand.totals < 2 ] <- NA
    pool.frame <- data.frame(
        "Reads" = rowSums(raw),
        "Tailed reads" = grand.totals,
        "Mean tail length" = grand.means,
        "StdDev tail length" = grand.sds,
        row.names = row.names,
        check.names = FALSE,
        check.rows = FALSE
    )

    # === Output ===
    
    write.grouped.table(
        list(
            "Annotation" = annotation,
            "Count" = as.data.frame( data$Count ),
            "Expression" = as.data.frame( expression ),
            "Pooled" = pool.frame,
            "Tail count" = as.data.frame( totals ),
            "Mean tail length" = as.data.frame( means ),
            "StdDev tail length" = as.data.frame( sds )
        ),
        sprintf("%s-statistics.csv", prefix),
        comments = c(
            '',
            'Expression columns are quantile normalized log2 Reads Per Million as produced by the limma voom function',
            ''
        )
    )
    
    """



@config.help(
"""\
Show a heatmap of the number of reads with different tail lengths in different genes.
""",
"""\
- Reads with an estimated tail length of 3 or less are ignored.

- Reads are pooled over all samples.
"""
)
@config.Positional('aggregate', 'Prefix of output from "aggregate-tail-lengths:"')
@config.Int_flag('min_tails', 'Minimum number of reads with tails in order to include in heatmap.')
@config.Float_flag('min_svd', 'Attempt to pick a sample of genes representative of all the different types of variation.')
class Plot_pooled(config.Action_with_prefix, runr.R_action):
    aggregate = None
    min_tails = 1000
    min_svd = 0

    script = r"""
    library(nesoni)
    
    # === Load data ===
        
    annotation <- read.grouped.table(sprintf('%s-counts.csv',aggregate),require='Annotation')$Annotation
    pooled <- as.matrix( read.grouped.table(sprintf('%s-pooled.csv',aggregate))$All )
    
    maxtail <- ncol(pooled)
    
    
    # === Normalize by row maximum ===

    # Omit columns for tail length 0,1,2,3    
    skip <- 4

    norm.pooled <- pooled[,(skip+1):maxtail,drop=FALSE]
    for(i in basic.seq(ncol(norm.pooled)))
        if (i%%10 != 1)
            colnames(norm.pooled)[i] <- ''
    
    col.length <- basic.seq(ncol(norm.pooled)) + (skip-1)

    pool.total <- rep(0,nrow(pooled))
    pool.mean <- rep(0,nrow(pooled))
    pool.sd <- rep(0,nrow(pooled))        

    for(i in basic.seq(nrow(pooled))) {
        pool.total[i] <- sum(norm.pooled[i,])

        pool.mean[i] <- sum(norm.pooled[i,] * col.length) / pool.total[i]
        
        if (pool.total[i] > 1) {
            deviation <- col.length - pool.mean[i]
            pool.sd[i] <- sqrt( sum( deviation*deviation*norm.pooled[i,] ) / (pool.total[i]-1) )
        }
        
        norm.pooled[i,] <- norm.pooled[i,] / max(norm.pooled[i,],1)
    }


    # === Pick interesting genes ===
        
    keep <- (pool.total > min_tails)
    
    if (min_svd > 0) {
        keep[keep] <- svd.gene.picker( scale(norm.pooled[keep,,drop=FALSE]), min.svd=min_svd )
    }
    
    kept.annotation <- annotation[keep,,drop=FALSE]    
    kept.norm.pooled <- norm.pooled[keep,,drop=FALSE]
    kept.pool.total <- pool.total[keep]
    kept.pool.mean <- pool.mean[keep]
    kept.pool.sd <- pool.sd[keep]

    # === Output ===

    n.rows <- nrow(kept.norm.pooled)
    row.labels <- n.rows <= 300        
    height <- if(row.labels) (25*n.rows+600) else 2500        
    png(sprintf("%s.png",prefix), width=2000, height=height, res=150)
    
    if (row.labels) {
       labels <- list(
           rownames(kept.annotation), 
           kept.annotation[,'gene'],
           kept.annotation[,'product']
       )
    } else {
       labels <- NULL
    }
    
    result <- nesoni.heatmap( 
        kept.norm.pooled, labels, 
        sort.mat=kept.norm.pooled, signed=FALSE, legend='number of reads\nscaled by row maximum')

    dev.off()


    # === Output CSV ===

    shuffle <- rev(result$dend.row$order)

    sink(sprintf('%s.csv', prefix))
    cat('# Pooled tail-length heatmap data\n')
    cat('#\n')
    cat('# Reads are counted as having a tail if there is a poly-A sequence of length 4 or greater\n')
    cat('#\n')

    frame <- data.frame( 
        Name = rownames(kept.annotation)[shuffle],
        Gene = kept.annotation$gene[shuffle],
        Product = kept.annotation$product[shuffle],
        'Cluster hierarchy' = rev(result$dend.row$paths),
        'Total reads with tails' = kept.pool.total[shuffle],
        'Mean tail length' = kept.pool.mean[shuffle],
        'Standard deviation tail length' = kept.pool.sd[shuffle],
        check.names=FALSE
    )

    write.csv(frame, row.names=FALSE)
    sink()

    """

@config.help(
'Produce a heatmap comparing average tail lengths between samples for different genes.'
)
@config.Positional('aggregate', 'Prefix of output from "aggregate-tail-lengths:"')
@config.Int_flag('min_tails', 'Minimum number of reads with tails *in each sample* required in order to include a gene.')
@config.Float_flag('min_span', 'Minimum difference in average tail lengths required in order to include a gene.')
class Plot_comparison(config.Action_with_prefix, runr.R_action):
    aggregate = None
    min_tails = 1000
    min_span = 4

    script = r"""
    library(nesoni)
    
    # === Load data ==
    
    data <- read.grouped.table(sprintf('%s-counts.csv',aggregate),require=c('Count','Annotation'))
    annotation <- data$Annotation
    expression <- voom( as.matrix(data$Count), normalize.method='quantile' )$E
    
    raw <- as.matrix(read.grouped.table(sprintf('%s-raw.csv',aggregate))$All)
    
    cols <- as.matrix(read.grouped.table(sprintf('%s-raw-columns.csv',aggregate))$All)   
    
    n.samples <- nrow(cols)
    n.genes <- nrow(raw)
    n.lengths <- ncol(cols)
    lengths <- basic.seq(n.lengths)-1


    # === Calculate mean tail lengths for each sample ===
    
    # Omit columns for tail length 0,1,2,3    
    skip <- 4
    good.cols <- cols[,(skip+1):n.lengths,drop=FALSE]
    good.lengths <- lengths[(skip+1):n.lengths]
    
    totals <- matrix(nrow=n.genes,ncol=n.samples) 
    means <- matrix(nrow=n.genes,ncol=n.samples)
    colnames(means) <- rownames(cols)
    for(c in basic.seq(n.samples)) {
        totals[,c] <- rowSums(raw[,good.cols[c,],drop=FALSE])
        means[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * good.lengths) ) / totals[,c]
        means[totals[,c] < min_tails,c] <- NA
    }


    # === Choose interesting genes to keep ===
    
    keep <- (
        row.apply(totals,min) >= min_tails &
        row.apply(means,min) <= row.apply(means,max) - min_span
    )
    
    kept.annotation <- annotation[keep,,drop=FALSE]
    kept.raw <- raw[keep,,drop=FALSE]
    kept.totals <- totals[keep,,drop=FALSE]    
    kept.means <- means[keep,,drop=FALSE]
    kept.expression <- expression[keep,,drop=FALSE]

    # === Cluster and travelling-salesman by tail length ===
        
    dend.row <- do.dendrogram( t(scale(t(kept.means))) )

       
    # === Output heatmap ===
    
    n.rows <- sum(keep)
    row.labels <- n.rows <= 300        
    height <- if(row.labels) (25*n.rows+700) else 2500        
    png(sprintf("%s.png",prefix), width=2000, height=height, res=150)

    if (row.labels) {
       labels <- list(
           rownames(kept.annotation)[dend.row$order], 
           kept.annotation[dend.row$order,'gene'],
           kept.annotation[dend.row$order,'product']
       )
    } else {
       labels <- NULL
    }
    
    plots <- list(
        list(
            weight=1,
            type='dendrogram',
            data=dend.row
        ),
        list(
            weight=1,
            type='heatmap',
            data=t(scale(t(kept.means[dend.row$order,,drop=FALSE]),scale=FALSE,center=TRUE)),
            signed=TRUE,
            legend='Tail length\n vs row average',
            title='Tail length  -vs-'
        ),
        list(
            weight=0.3,
            type='scatter',
            data=rowMeans( kept.means[dend.row$order,,drop=FALSE] ),
            title='row\naverage'
        ),
        list(
            weight=1,
            type='heatmap',
            data=t(scale(t( kept.expression[dend.row$order,,drop=FALSE] ),scale=FALSE,center=TRUE)),
            signed=TRUE,
            legend='log2 Reads Per Million\n(quantile normalized)\nvs row average',
            title='log2 RPM  -vs-'
        ),
        list(
            weight=0.3,
            type='scatter',
            data=rowMeans( kept.expression[dend.row$order,,drop=FALSE] ),
            title='row\naverage'
        )
    )
    
    #for(i in basic.seq(n.samples)) {
    #    data <- (kept.raw[,good.cols[i,],drop=FALSE] / kept.totals[,i])[dend.row$order,,drop=FALSE]
    #    
    #    #if (nrow(data) > 0) {
    #    #    for(j in basic.seq(ncol(data)))
    #    #        data[,j] <- data[,j]*ncol(data)
    #    #}
    #    
    #    colnames(data) <- rep('', ncol(data))
    #    colnames(data)[1] <- rownames(cols)[i]
    #    plots[[length(plots)+1]] <- list(
    #        weight=1/n.samples,
    #        type='heatmap',
    #        data=data,
    #        signed=FALSE
    #    )
    #}
    
    multiplot(
        plots,
        labels
    )
    
    dev.off()


    # === Output CSV ===

    shuffle <- rev(dend.row$order)

    sink(sprintf('%s.csv', prefix))
    cat('# Tail-length heatmap data\n')
    cat('#\n')
    cat('# log2 counts are quantile normalized\n')
    cat('#\n')
    
    frame <- data.frame( 
        Name = rownames(kept.annotation)[shuffle],
        Gene = kept.annotation$gene[shuffle],
        Product = kept.annotation$product[shuffle],
        'Cluster hierarchy' = rev(dend.row$paths),
        check.names=FALSE
    )
    
    for(i in basic.seq(n.samples))
        frame[, paste(colnames(kept.means)[i],'mean tail')] <- kept.means[shuffle,i,drop=FALSE]

    for(i in basic.seq(n.samples))
        frame[, paste(colnames(kept.means)[i],'log2 count')] <- kept.expression[shuffle,i,drop=FALSE]

    write.csv(frame, row.names=FALSE)
    sink()
    
    """



@config.help(
'Analyse tail lengths of a set of samples.'
"""\
"""
)
@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.Int_flag('saturation',
     'Duplicate start position saturation level. '
     'Reads that start at the same position will only '
     'count as up to this many. Zero for no saturation.'
)
@config.Main_section('working_dirs')
class Analyse_tail_lengths(config.Action_with_prefix):         
    types = 'gene'
    saturation = 1
    working_dirs = [ ]

    def log_filename(self):
        if self.prefix is None: return None
        return self.prefix + '_analyse_tail_lengths_log.txt'
    
    def get_plot_pooleds(self):
        return [ 
        
            Plot_pooled(
                prefix = '%s-pooled-min-tails-%d' % (self.prefix, min_tails),
                aggregate = self.prefix,
                min_tails = min_tails,
            )
            
            for min_tails in (20,50,100,200,500,1000,2000)
        ]
    
    def get_plot_comparisons(self):
        return [ 
        
            Plot_comparison(
                prefix = '%s-comparison-min-tails-%d-min-span-%.1f' % (self.prefix, min_tails, min_span),
                aggregate = self.prefix,
                min_tails = min_tails,
                min_span = min_span,
            )
            
            for min_tails in (20,50,100,200)
            for min_span in (2,4,8)
        ]

    def get_heatmaps(self):
        return [
            
            nesoni.Heatmap(
                prefix = '%s-heatmap-min-fold-%.1f-min-max-%d' % (self.prefix, fold, min_count),
                counts = self.prefix + '-counts.csv',
                min_total = 0,
                min_max = min_count,
                min_span = math.log(fold)/math.log(2.0),
            )
            
            for fold in [ 1.5, 2.0, 4.0, 6.0, 8.0 ]
            for min_count in [ 10, 50, 250 ]
        ]

    def run(self):
         stage = nesoni.Stage()

         for dir in self.working_dirs:
             Tail_lengths(
                 working_dir=dir, 
                 types=self.types,
             ).process_make(stage)    

         stage.barrier() #=================================================
         
         Aggregate_tail_lengths(
             prefix=self.prefix, 
             working_dirs=self.working_dirs,
             saturation=self.saturation,
         ).make()        
         
         Tail_stats(
             prefix=self.prefix
         ).make()

         for action in self.get_plot_pooleds() + self.get_plot_comparisons() + self.get_heatmaps():
             action.process_make(stage)
         
         nesoni.Plot_counts(
             prefix=self.prefix,
             counts=self.prefix+'-counts.csv',
         ).process_make(stage)
         
         stage.barrier() #=================================================











