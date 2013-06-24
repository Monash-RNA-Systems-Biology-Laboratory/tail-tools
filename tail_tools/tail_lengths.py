
import itertools, collections, math

import nesoni
from nesoni import annotation, sam, span_index, config, working_directory, workspace, io, runr, reporting

import cPickle as pickle

@config.help(
'Create file to be used by "aggregate-tail-lengths:".'
)
@config.String_flag('annotations', 'Filename containing annotations. Defaults to annotations in reference directory.')
@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.Positional('working_dir', 'Working directory to use as input.')
class Tail_lengths(config.Action_with_prefix):
     annotations = None
     types = 'gene'
     working_dir = None

     #_workspace_class = working_directory.Working

     def run(self):
         #workspace = self.get_workspace()
         workspace = working_directory.Working(self.working_dir, must_exist=True)
         if self.annotations == None:
             reference = workspace.get_reference()
             annotations_filename = reference.annotations_filename()
         else:
             annotations_filename = self.annotations
         
         types = [ item.lower() for item in self.types.split(',') ]
     
         annotations = [ 
             item 
             for item in annotation.read_annotations(annotations_filename)
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
                     hits = [ 
                         gene
                         for gene in index[fragment[0].rname].get(start,end)
                         if gene.strand == strand
                         ]                         
                     if hits:
                         hits.sort(key=lambda gene: gene.get_id()) #Ensure a deterministic order
                         gene = max(hits, key=lambda gene: 
                             (min(end,gene.end)-max(start,gene.start),
                              abs((start+end)-(gene.start+gene.end)))
                             ) 
                             #Maximum overlap hit
                             #(failing that, closest middles)
                             
                         if strand > 0:
                             rel_start = start - gene.start
                             rel_end = end - gene.start
                         else:
                             rel_start = gene.end - end
                             rel_end = gene.end - start
                         
                         gene.hits.append( (rel_start,rel_end,tail_length) )

         #workspace.set_object(annotations, 'tail-lengths.pickle.gz')
         f = io.open_possibly_compressed_writer(self.prefix + '.pickle.gz')
         pickle.dump((workspace.name, workspace.get_tags(), annotations), f, pickle.HIGHEST_PROTOCOL)
         f.close()
         


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
@config.Int_flag('tail',
     'Minimum tail length to count as having a tail.'
     )
@config.Int_flag('ntail',
     'Minimum number of observed tails lengths needed to call a mean tail length.'
     )
@config.Int_flag('nprop',
     'Minimum number of observed reads to call a proportion of reads with tails.'
     )
@config.Main_section('pickles')
class Aggregate_tail_lengths(config.Action_with_prefix):             
    saturation = 1
    tail = 4
    ntail = 10
    nprop = 10
    pickles = [ ]

    def log_filename(self):
        if self.prefix is None: return None
        return self.prefix + '_aggregate_log.txt'

    def state_filename(self):
        return self.prefix + '_aggregate.state'
         
    def run(self):
        #workspaces = [ 
        #    working_directory.Working(item) 
        #    for item in self.working_dirs 
        #]
        #
        #data = [ 
        #    item.get_object('tail-lengths.pickle.gz') 
        #    for item in workspaces
        #]        
        #
        #names = [ item.name for item in workspaces ]
        
        data = [ ]
        names = [ ]
        sample_tags = [ ]
        
        for item in self.pickles:
            f = io.open_possibly_compressed_file(item)
            name, tags, datum = pickle.load(f)
            f.close()
            data.append(datum)
            names.append(name)
            sample_tags.append(tags)
        
        annotations = data[0]
        
        all_lengths = [ 
            tail_length
            for sample in data
            for feature in sample
            for rel_start,rel_end,tail_length in feature.hits
            ]
        if all_lengths: 
            max_length = max(all_lengths)+1
        else:
            max_length = 1
        del all_lengths
        
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

            self.log.datum(names[i], 'Alignments to genes', n_alignments)
            if self.saturation >= 1:
                self.log.datum(names[i], 'Proportion of alignments with duplicate start and end position', float(n_duplicates)/max(1,n_alignments))
                self.log.datum(names[i], 'Alignments to genes after deduplication', n_good)
                
        
        counts = [ ]  # [feature][sample][taillength]
        
        for item in data: 
            assert len(item) == len(data[0])
        for row in itertools.izip(*data):
            this_counts = [ item.tail_counts for item in row ]
            counts.append(this_counts)
        
        totals = [ ]          # [feature][sample]  Total count
        sample_prop = [ ]     # [feature][sample]  Proportion of reads with tail
        sample_tail = [ ]     # [feature][sample]  Mean tail length in each sample
        overall_prop = [ ]    # [feature]          Overall proportion with tail
        overall_tail = [ ]    # [feature]          Overall mean tail length
        overall_n_tail = [ ]
        for row in counts:
            this_totals = [ ]
            this_prop = [ ]
            this_tail = [ ]
            this_n = 0.0
            this_n_tail = 0.0
            this_total_tail = 0.0
            for item in row:
                this_this_n = sum(item)
                this_totals.append( this_this_n )
                this_this_n_tail = sum(item[self.tail:])
                this_this_total_tail = sum( item[i]*i for i in xrange(self.tail,max_length) )
                if this_this_n < self.nprop:
                    this_prop.append(None)
                else:
                    this_prop.append(float(this_this_n_tail)/this_this_n)
                if this_this_n_tail < self.ntail:
                    this_tail.append(None)
                else:
                    this_tail.append(this_this_total_tail/this_this_n_tail)
                this_n += this_this_n
                this_n_tail += this_this_n_tail
                this_total_tail += this_this_total_tail

            totals.append(this_totals)
            sample_prop.append(this_prop)
            sample_tail.append(this_tail)
            overall_n_tail.append(this_n_tail)
            if this_n < self.nprop:
                overall_prop.append(None)
            else:
                overall_prop.append(float(this_n_tail)/this_n)
            if this_n_tail < self.ntail:
                overall_tail.append(None)
            else:
                overall_tail.append(this_total_tail/this_n_tail)
             
        for i, name in enumerate(names):
            this_total = 0.0            
            this_n = 0
            for row in sample_tail:
                if row[i] is not None:
                    this_total += row[i]
                    this_n += 1
            if this_n:
                self.log.datum(name, 'Average poly-A tail', this_total/this_n)
                
        for i, name in enumerate(names):
            this_total = 0.0
            this_n = 0
            for row in sample_prop:
                if row[i] is not None:
                    this_total += row[i]
                    this_n += 1
            if this_n:
                self.log.datum(name, 'Average proportion of reads with tail', this_total/this_n)
            
        
        #max_length = max(max(len(item) for item in row) for row in counts)
        #
        #for row in counts:
        #    for item in row:
        #        while len(item) < max_length:
        #            item.append(0)
                
        
        with open(self.prefix+'.gff','wb') as f:
            annotation.write_gff3_header(f)
            for i, item in enumerate(annotations):
                item.attr['mean_tail'] = '%.1f'%overall_tail[i] if overall_tail[i] else 'NA'
                item.attr['proportion_with_tail'] = '%.2f'%overall_prop[i] if overall_prop[i] else 'NA'
                
                if overall_tail[i] is None:
                    item.attr['color'] = '#444444'
                else:
                    a = (overall_tail[i]-self.tail)/max(1,max_length-self.tail)
                    item.attr['color'] = '#%02x%02x%02x' % (int(a*255),int((1-abs(a*2-1))*255),255-int(a*255))
                #item.attr['color'] = ...                
                print >> f, item.as_gff()
        
        
        comments = [ '#Counts' ] + [
            '#sampleTags='+','.join(tags)
            for tags in sample_tags
            ] + [
            'n-with-tail = Number of reads with a poly-A tail of length at least %d' % self.tail,
            'mean-tail = Mean poly-A tail length, for reads with a tail of length at least %d, if there were at least %d such reads' % (self.tail, self.ntail),
            'proportion-with-tail = Proportion of reads with a tail of length at least %d, if there were at least %d total reads' % (self.tail, self.nprop),
            '"Tail" group is mean tail per sample, as above',
            '"Proportion" group is proportion of reads with tail per sample, as above',
            ]

        def counts_iter():
            for i in xrange(len(counts)):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(len(names)):
                    row[('Count',names[j])] = '%d' % totals[i][j] #sum(counts[i][j]) 

                row[('Annotation','Length')] = annotations[i].end - annotations[i].start
                row[('Annotation','gene')] = annotations[i].attr.get('gene','')
                row[('Annotation','product')] = annotations[i].attr.get('product','')
                #row[('Annotation','Strand')] = str(annotations[i].strand)
                row[('Annotation','n-with-tail')] = str(overall_n_tail[i])
                row[('Annotation','mean-tail')] = str(overall_tail[i]) if overall_tail[i] is not None else 'NA'
                row[('Annotation','proportion-with-tail')] = str(overall_prop[i]) if overall_prop[i] is not None else 'NA'
                for j in xrange(len(names)):
                    row[('Tail',names[j])] = str(sample_tail[i][j]) if sample_tail[i][j] is not None else 'NA'
                for j in xrange(len(names)):
                    row[('Proportion',names[j])] = str(sample_prop[i][j]) if sample_prop[i][j] is not None else 'NA'
                yield row
        io.write_csv(self.prefix + '-counts.csv', counts_iter(), comments=comments)

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
        
        ##Note: zero if zero count
        #def mean_length():
        #    for i in xrange(len(counts)):
        #        row = collections.OrderedDict()
        #        row['Feature'] = annotations[i].get_id()
        #        for j in xrange(len(names)):
        #            a = 0
        #            b = 0
        #            for pos,count in enumerate(counts[i][j]):
        #                a += pos*count
        #                b += count
        #            row[names[j]] = str( float(a)/max(1,b) )
        #        yield row
        #io.write_csv(self.prefix + '-mean-length.csv', mean_length())
        #
        #def mean4_length():
        #    for i in xrange(len(counts)):
        #        row = collections.OrderedDict()
        #        row['Feature'] = annotations[i].get_id()
        #        for j in xrange(len(names)):
        #            a = 0
        #            b = 0
        #            for pos,count in enumerate(counts[i][j]):
        #                if pos < 4: continue
        #                a += pos*count
        #                b += count
        #            row[names[j]] = str( float(a)/max(1,b) )
        #        yield row
        #io.write_csv(self.prefix + '-mean-length-beyond4.csv', mean4_length())


#@config.help(
#    'Create a spreadsheet containing various statistics, using the output of "aggregate-tail-lengths:"',
#    'The prefix should be the same as that given to "aggregate-tail-lengths:".'
#)
#class Tail_stats(config.Action_with_prefix, runr.R_action):
#    def state_filename(self):
#        return self.prefix + '_tail_stats.state'
#
#    def log_filename(self):
#        if self.prefix is None: return None
#        return self.prefix + '_tail_stats_log.txt'
#        
#    script = r"""
#    
#    library(nesoni)
#
#    # === Load data ==
#    
#    data <- read.grouped.table(sprintf('%s-counts.csv',prefix),require=c('Count','Annotation'))
#    annotation <- data$Annotation
#    expression <- voom( as.matrix(data$Count), normalize.method='quantile' )$E
#    
#    raw <- as.matrix(read.grouped.table(sprintf('%s-raw.csv',prefix))$All)
#    
#    cols <- as.matrix(read.grouped.table(sprintf('%s-raw-columns.csv',prefix))$All)   
#    
#    row.names <- rownames(annotation)
#    n.samples <- nrow(cols)
#    n.genes <- nrow(raw)
#    n.lengths <- ncol(cols)
#    lengths <- basic.seq(n.lengths)-1
#
#
#    # === Calculate mean tail lengths for each sample ===
#    
#    # Omit columns for tail length 0,1,2,3    
#    skip <- 4
#    good.cols <- cols[,(skip+1):n.lengths,drop=FALSE]
#    good.lengths <- lengths[(skip+1):n.lengths]
#    
#    totals <- matrix(nrow=n.genes,ncol=n.samples) 
#    x.totals <- matrix(nrow=n.genes,ncol=n.samples) 
#    xx.totals <- matrix(nrow=n.genes,ncol=n.samples) 
#    means <- matrix(nrow=n.genes,ncol=n.samples)
#    sds <- matrix(nrow=n.genes,ncol=n.samples)
#    colnames(totals) <- rownames(cols)
#    rownames(totals) <- row.names
#    colnames(means) <- rownames(cols)
#    rownames(means) <- row.names
#    colnames(sds) <- rownames(cols)
#    rownames(sds) <- row.names
#    for(c in basic.seq(n.samples)) {
#        totals[,c] <- rowSums(raw[,good.cols[c,],drop=FALSE])
#        x.totals[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * good.lengths) )
#        xx.totals[,c] <- rowSums( t(t(raw[,good.cols[c,],drop=FALSE]) * (good.lengths*good.lengths) ) )
#
#        means[,c] <- x.totals[,c] / totals[,c]
#        sds[,c] <- sqrt( totals[,c]/(totals[,c]-1) * (xx.totals[,c]/totals[,c] - means[,c]*means[,c]) )
#        
#        means[totals[,c] < 1,c] <- NA
#        sds[totals[,c] < 2,c] <- NA
#    }
#
#
#    grand.totals <- rowSums(totals)
#    grand.x.totals <- rowSums(x.totals)
#    grand.xx.totals <- rowSums(xx.totals)
#    grand.means <- grand.x.totals / grand.totals
#    grand.sds <- sqrt( grand.totals/(grand.totals-1) * (grand.xx.totals/grand.totals - grand.means*grand.means) )
#    grand.means[ grand.totals < 1 ] <- NA
#    grand.sds[ grand.totals < 2 ] <- NA
#    pool.frame <- data.frame(
#        "Reads" = rowSums(raw),
#        "Tailed reads" = grand.totals,
#        "Mean tail length" = grand.means,
#        "StdDev tail length" = grand.sds,
#        row.names = row.names,
#        check.names = FALSE,
#        check.rows = FALSE
#    )
#
#    # === Output ===
#    
#    write.grouped.table(
#        list(
#            "Annotation" = annotation,
#            "Count" = as.data.frame( data$Count ),
#            "Expression" = as.data.frame( expression ),
#            "Pooled" = pool.frame,
#            "Tail count" = as.data.frame( totals ),
#            "Mean tail length" = as.data.frame( means ),
#            "StdDev tail length" = as.data.frame( sds )
#        ),
#        sprintf("%s-statistics.csv", prefix),
#        comments = c(
#            '',
#            'Expression columns are quantile normalized log2 Reads Per Million as produced by the limma voom function',
#            ''
#        )
#    )
#    
#    """



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
@config.Float_flag('glog_moderation', 'Amount of moderation used in log transformation of expression level.')
class Plot_comparison(config.Action_with_prefix, runr.R_action):
    aggregate = None
    min_tails = 100
    min_span = 4.0
    glog_moderation = 5.0

    script = r"""
    library(nesoni)
    
    # === Load data ==
    
    #data <- read.grouped.table(sprintf('%s-counts.csv',aggregate),require=c('Count','Annotation'))
    #annotation <- data$Annotation
    #expression <- voom( as.matrix(data$Count), normalize.method='quantile' )$E
    
    dgelist <- read.counts(sprintf('%s-counts.csv',aggregate))
    expression <- glog2.rpm.counts(dgelist, glog_moderation)$E
    annotation <- dgelist$genes
    
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
@config.String_flag('annotations', 'Filename containing annotations. Defaults to annotations in reference directory.')
@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.Int_flag('saturation',
     'Duplicate start position saturation level. '
     'Reads that start at the same position will only '
     'count as up to this many. Zero for no saturation.'
)
@config.Main_section('working_dirs')
class Analyse_tail_lengths(config.Action_with_prefix):         
    annotations = None
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
                #prefix = '%s-heatmap-min-fold-%.1f-min-max-%d' % (self.prefix, fold, min_count),
                prefix = '%s-heatmap-min-fold-%.1f' % (self.prefix, fold),
                counts = self.prefix + '-counts.csv',
                min_total = 0,
                #min_max = min_count,
                min_span = math.log(fold)/math.log(2.0),
            )
            
            for fold in [ 1.5, 2.0, 4.0, 6.0, 8.0, 10.0, 20.0 ]
            #for min_count in [ 10, 50, 250 ]
        ]

    def run(self):
         stage = nesoni.Stage()
         
         names = [ ]
         
         pickle_workspace = workspace.Workspace(self.prefix + '-pickles')
         pickle_filenames = [ ]

         for dir in self.working_dirs:
             working = working_directory.Working(dir, must_exist=True)
             pickle_filenames.append(pickle_workspace/working.name+'.pickle.gz')
             Tail_lengths(
                 pickle_workspace/working.name,
                 working_dir=dir, 
                 annotations=self.annotations,
                 types=self.types,
             ).process_make(stage)    

         stage.barrier() #=================================================
         
         Aggregate_tail_lengths(
             prefix=self.prefix, 
             #working_dirs=self.working_dirs,
             pickles=pickle_filenames,
             saturation=self.saturation,
         ).make()        
         
         #Tail_stats(
         #    prefix=self.prefix
         #).make()

         for action in self.get_plot_pooleds() + self.get_plot_comparisons() + self.get_heatmaps():
             action.process_make(stage)
         
         nesoni.Plot_counts(
             prefix=self.prefix,
             counts=self.prefix+'-counts.csv',
         ).process_make(stage)
         
         stage.barrier() #=================================================

         r = reporting.Reporter(self.prefix + '-report', 
             '3\'seq expression analysis ' + ('with read deduplication' if self.saturation else 'without read deduplication'), 
             '')         
         self.report_tail_lengths(r)
         r.close()



    def report_tail_lengths(self, r):        
        saturation = self.saturation
        if saturation:
            r.p(
                'Duplicate read removal: Sets of reads aligning with exactly the start and end position in the reference were counted as a single read.'
            )

        r.p(
            'This scatterplot show the number of reads aligning to each gene between each pair of samples. '
            'This can be used to discover poor samples.'
        )        
        
        r.p( r.get(self.prefix + '-count.png', image=True) )
        
        r.p('Counts are converted to '
            'log2 Reads Per Million using a variance stabilised transformation. '
            'Let the generalised logarithm with moderation m be')
        
        r.p('glog(x,m) = log2((x+sqrt(x*x+4*m*m))/2)')
        
        r.p('then the transformed values will be')
        
        r.p('glog( count/library_size*1e6, log_moderation/mean_library_size*1e6 )')
        
        r.p('where log_moderation is a parameter, here 5. '
            'The library sizes used are effective library sizes after TMM normalization.')
        
        
        r.subheading('Spreadsheet with statistics for all genes and all samples')
        
        if saturation:
            r.p( 'Note: Reads with the same start and end position sometimes don\'t have the same tail length. '
                 'After deduplication these can contribute fractionally to the number of reads with tails.' )
        
        #r.p( r.get(self.prefix + '-statistics.csv') )
        r.p( r.get(self.prefix + '-counts.csv') )
        
        r.subheading('Poly(A) tail length in reads')
        
        r.p(
            'Some reads contain a poly(A) sequence not found in the reference. '
            'It is hoped that the lengths of these poly(A) sequences gives an indication of the true length of polyadenylation.'
        )
        
        r.p(
            'Only reads with a poly(A) sequence of four or more bases are used.'
        )
        
        for heatmap in self.get_plot_pooleds():
            r.report_heatmap(heatmap)
            
        r.subheading('Average poly(A) tail length and its relation to expression levels')
        
        r.p(
            'Only reads with a poly(A) sequence of four or more bases was included in the averages.'
        )
        
        r.p(
            'Genes were selected based on there being at least some number of reads with poly(A) sequence in <i>each</i> sample (min-tails), '
            'and on there being at least some amount of difference in average tail length between samples (min-span).'
        )
        
        for heatmap in self.get_plot_comparisons():
            r.report_heatmap(heatmap)
        
        r.subheading('Heatmaps')
        
        r.p(
            'Genes were selected based on there being '
            #'at least some number of reads '
            #'in at least one of the samples (min-max), '
            #'and '
            'on there being at least some fold change difference between '
            'some pair of samples.'
        )
        
        #r.p(
        #    'The log2 counts have been quantile normalized.'
        #)
        
        for heatmap in self.get_heatmaps():
            r.report_heatmap(heatmap)

        r.heading('Raw data')
        
        r.p(r.get(self.prefix + '-counts.csv'))
        








