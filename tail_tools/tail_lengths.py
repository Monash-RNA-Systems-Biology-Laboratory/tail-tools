
import itertools, collections

import nesoni
from nesoni import annotation, sam, span_index, config, working_directory, io, runr

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


@config.Int_flag('saturation', 
     'Duplicate start position saturation level. '
     'Reads that start at the same position will only '
     'count as up to this many.'
)
@config.Main_section('working_dirs')
class Aggregate_tail_lengths(config.Action_with_prefix):         
    saturation = None
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
        
        for sample in data:
            for feature in sample:
                feature.tail_counts = [ 0.0 ] * max_length
                
                buckets = collections.defaultdict(list)
                for rel_start,rel_end,tail_length in feature.hits:
                    buckets[ (rel_start,rel_end) ].append(tail_length)
                for item in buckets.values():
                    if self.saturation is None or len(item) <= self.saturation:
                        weight = 1.0
                    else:
                        weight = float(self.saturation) / len(item)
                    for item2 in item:
                        feature.tail_counts[item2] += weight
                    
                
        
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
        
        
        

        def annotation_table():
            for item in annotations:
                row = collections.OrderedDict()
                row['Feature'] = item.get_id()
                row['Length'] = item.end - item.start
                row['Gene'] = item.attr.get('gene','')
                row['Product'] = item.attr.get('product','')
                row['Strand'] = str(item.strand)
                yield row
        io.write_csv(self.prefix + '-annotation.csv', annotation_table())

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
        
        def total():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    row[names[j]] = str( sum(counts[i][j]) )
                yield row
        io.write_csv(self.prefix + '-total.csv', total())

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


@config.Positional('aggregate', 'Prefix of output from "aggregate-tail-lengths:"')
@config.Int_flag('min_tails', 'Minimum number of reads with tails in order to include in heatmap.')
class Plot_pooled(runr.R_action, config.Action_with_prefix):
    aggregate = None
    min_tails = 1000

    script = r"""
    library(nesoni)
    
    ann <- nesoni.read.table(sprintf('%s-annotation.csv',aggregate))
    pooled <- as.matrix( nesoni.read.table(sprintf('%s-pooled.csv',aggregate)) )
    
    maxtail <- ncol(pooled)
    skip <- 4

    norm.pooled <- pooled[,(skip+1):maxtail,drop=FALSE]
    for(i in basic.seq(ncol(norm.pooled)))
        if (i%%10 != 1)
            colnames(norm.pooled)[i] <- ''

    pool.total <- rep(0,nrow(pooled))        

    for(i in basic.seq(nrow(pooled))) {
        pool.total[i] <- sum(norm.pooled[i,])
        norm.pooled[i,] <- norm.pooled[i,] / max(norm.pooled[i,],1)
    }
    
    keep <- pool.total > min_tails
    kept.ann <- ann[keep,,drop=FALSE]
    
    mat <- norm.pooled[keep,]

    n.rows <- nrow(mat)
    row.labels <- n.rows <= 300        
    height <- if(row.labels) (25*n.rows+600) else 2500        
    png(sprintf("%s.png",prefix), width=2000, height=height, res=150)
    
    if (row.labels) {
       labels <- list(rownames(kept.ann), kept.ann[,'Gene'],kept.ann[,'Product'])
    } else {
       labels <- NULL
    }
    
    nesoni.heatmap( 
        mat, labels, 
        sort.mat=mat, signed=FALSE, legend='number of reads\nscaled by row maximum\n')

    dev.off()
    """


@config.Positional('aggregate', 'Prefix of output from "aggregate-tail-lengths:"')
@config.Int_flag('min_tails', 'Minimum number of reads with tails in order to include in heatmap.')
class Plot_comparison(runr.R_action, config.Action_with_prefix):
    aggregate = None
    min_tails = 1000

    script = r"""
    library(nesoni)
    
    ann <- nesoni.read.table(sprintf('%s-annotation.csv',aggregate))
    raw <- as.matrix(nesoni.read.table(sprintf('%s-raw.csv',aggregate)))
    
    expression <- voom( as.matrix(nesoni.read.table(sprintf('%s-total.csv',aggregate))), normalize.method='quantile' )$E
    
    cols <- as.matrix(nesoni.read.table('test-raw-columns.txt'))   
    
    n.samples <- nrow(cols)
    n.genes <- nrow(raw)
    n.lengths <- ncol(cols)
    lengths <- basic.seq(n.lengths)-1
    
    skip <- 4
    good.cols <- cols[,(skip+1):n.lengths,drop=FALSE]
    good.lengths <- lengths[(skip+1):n.lengths]
    
    totals <- matrix(nrow=n.genes,ncol=n.samples) 
    means <- matrix(nrow=n.genes,ncol=n.samples)
    colnames(means) <- rownames(cols)
    for(c in basic.seq(n.samples)) {
        totals[,c] <- rowSums(raw[,good.cols[c,],drop=FALSE])
        means[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * good.lengths) ) / totals[,c]
    }
    
    keep <- (
        row.apply(totals,min) > min_tails &
        row.apply(means,min) <= row.apply(means,max) - 8
    )
    
    kept.ann <- ann[keep,,drop=FALSE]
    kept.raw <- raw[keep,,drop=FALSE]
    kept.totals <- totals[keep,,drop=FALSE]
    
    kept.means <- means[keep,,drop=FALSE]
    dend.row <- do.dendrogram( t(scale(t(kept.means))) )
       
    n.rows <- sum(keep)
    row.labels <- n.rows <= 300        
    height <- if(row.labels) (25*n.rows+700) else 2500        
    png(sprintf("%s.png",prefix), width=2000, height=height, res=150)

    if (row.labels) {
       labels <- list(rownames(kept.ann)[dend.row$order], kept.ann[dend.row$order,'Gene'],kept.ann[dend.row$order,'Product'])
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
            data=t(scale(t( expression[keep,,drop=FALSE][dend.row$order,,drop=FALSE] ),scale=FALSE,center=TRUE)),
            signed=TRUE,
            legend='log2 count\nvs row average',
            title='log2 count  -vs-'
        ),
        list(
            weight=0.3,
            type='scatter',
            data=rowMeans( expression[keep,,drop=FALSE][dend.row$order,,drop=FALSE] ),
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
    
    """




@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.Int_flag('saturation',
     'Duplicate start position saturation level. '
     'Reads that start at the same position will only '
     'count as up to this many.'
)
@config.Main_section('working_dirs')
class Analyse_tail_lengths(config.Action_with_prefix):         
    types = 'gene'
    saturation = None
    working_dirs = [ ]

    def run(self):
         stage = nesoni.Stage()

         for dir in self.working_dirs:
             Tail_lengths(
                 working_dir=dir, 
                 types=self.types,
             ).process_make(stage)    

         stage.barrier()
         
         Aggregate_tail_lengths(
             prefix=self.prefix, 
             working_dirs=self.working_dirs,
             saturation=self.saturation,
         ).make()






