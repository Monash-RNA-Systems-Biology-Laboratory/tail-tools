"""

Test for shift of expression from one tail to another within a gene.

"""

import json, collections

import nesoni
from nesoni import config, io, bio, annotation, runr, reference_directory

def _float_or_none(text):
    if text == 'NA':
        return None
    return float(text)

def _text(value):
    if value is None:
        return 'None'
    return str(value)

def _annotation_sorter(item):
    if item.strand < 0:
        return -item.end
    else:
        return item.start

CHISQ = """

cat('Performing chi-square tests.\n\n')

n <- length(data)
result <- data.frame(
    name = names(data), 
    peaks = rep(NA, n),
    df = rep(NA, n),
    chisq = rep(NA, n),
    log.p = rep(NA, n),
    p.value = rep(NA, n),
    interestingness = rep(NA, n),
    gene = GENES,
    product = PRODUCTS,
    mean.tail = MEAN_TAILS,
    proportion.with.tail = PROP_TAILS,
    peak.names = PEAK_NAMES,
    chromosome = CHROMOSOME_NAMES,
    strand = STRANDS,
    transcription.stops = TRANSCRIPTION_STOPS,
    prepeak.seq = PREPEAK_SEQS,
    interpeak.seq = INTERPEAK_SEQS,
    row.names = names(data)
    )

for(gene in names(data)) {
    mat <- data[[gene]]    
    test <- chisq.test(mat)
    
    total <- sum(mat)
    expect <- outer(rowSums(mat),colSums(mat))/total
    interestingness <- sum(abs(mat-expect))/total
    
    # Subtract expected interestingness
    for(e in expect) {
        interestingness <- interestingness - sum(dbinom(0:total,total,e/total) * abs((0:total)-e))/total
    }
    
    result[gene,'peaks'] <- nrow(mat)
    result[gene,'df'] <- test$parameter
    result[gene,'chisq'] <- test$statistic    
    result[gene,'interestingness'] <- interestingness
    result[gene,'log.p'] <- pchisq(test$statistic, test$parameter, lower.tail=FALSE, log.p=TRUE)
    result[gene,'p.value'] <- test$p.value
}

result[,'FDR'] <- p.adjust(result$p.value, method='BH')

#reordering <- order(result$log.p)

reordering <- order(result$interestingness, decreasing=TRUE)

output <- result[
    reordering, 
    c('name','interestingness','peaks','df','chisq','FDR',
      'gene','product','mean.tail','proportion.with.tail',
      'peak.names', 'chromosome','strand','transcription.stops',
      'prepeak.seq','interpeak.seq'
      )
    ]

sink(OUTPUT_FILENAME)
cat('##Compare-peaks chi-square tests\n')
write.csv(output, na='', row.names=FALSE)
sink()

"""

@config.help(
    'Look for shifts of transcription end point within genes.',
    'Output is:\n'
    '\n'
    '1. "Interestingness" and chi-square tests on all genes, from a table of counts with samples as rows and peaks as columns.\n'
    '\n'
    '2. A "counts" file and accompanying normalization file '
    'with one row for each pair of peaks which share a common gene. '
    'This can be used as input to "nesoni test-counts:".\n'
    )
@config.String_flag('norm_file', 'Use normalization produced by "norm-from-counts:".')
@config.String_flag('parent_type', 'Type of parent genes in parents annotation file.')
@config.String_flag('utrs', '3\' UTR annotations.')
@config.Bool_flag('utr_only', 
    'Only use peaks in the 3\' UTR. '
    )
@config.Int_flag('top',
    'Only use top n most expressed peaks. 0 = use all.'
    )
@config.Int_flag('max_seq',
    'Maximum length of pre-peak and inter-peak sequence to include in output.'
    )
@config.Positional('reference', 'Reference sequences.')
@config.Positional('parents', 'Annotation file containing parent genes.')
@config.Positional('children', 'Annotation file containing child peaks. "Parent" property should be set.')
@config.Positional('counts', 'Counts file as produced by tail-stats.')
class Compare_peaks(config.Action_with_prefix):
    norm_file = None
    parent_type = 'gene'
    utrs = None
    utr_only = True
    top = 0
    max_seq = 50000
    reference = None
    parents = None
    children = None
    counts = None
    
    def run(self):
        #assert not self.utr_only or self.utrs, '--utrs-only yes but no --utrs given'
        
        # Reference genome
        
        #chromosome_lengths = reference_directory.Reference(self.reference, must_exist=True).get_lengths()
        chromosomes = collections.OrderedDict(io.read_sequences(self.reference))

        def get_interpeak_seq(peaks):
            start = min(item.transcription_stop for item in peaks)
            end = max(item.transcription_stop for item in peaks)
            if end-start > self.max_seq: return ''
            if peaks[0].strand >= 0:
                return chromosomes[peaks[0].seqid][start:end]
            else:
                return bio.reverse_complement(chromosomes[peaks[0].seqid][start:end])

        def get_prepeak_seq(gene,peaks):
            if gene.strand >= 0:
                start = gene.utr_pos
                end = min(item.transcription_stop for item in peaks)
                if end-start > self.max_seq: return ''
                return chromosomes[gene.seqid][start:end]
            else:
                start = max(item.transcription_stop for item in peaks)
                end = gene.utr_pos
                if end-start > self.max_seq: return ''
                return bio.reverse_complement(chromosomes[gene.seqid][start:end])
        
        # Normalization files
        
        if self.norm_file:
            norm_file = self.norm_file
        else:
            nesoni.Norm_from_counts(self.prefix+'-norm', self.counts).run()
            norm_file = self.prefix+'-norm.csv'

        norms = io.read_grouped_table(norm_file, [('All',str)])['All']
        pair_norm_names = [ ]
        pair_norms = [ ]
        for i in xrange(len(norms)):
            pair_norm_names.append(norms.keys()[i]+'-peak1')
            pair_norms.append(norms.values()[i])
        for i in xrange(len(norms)):
            pair_norm_names.append(norms.keys()[i]+'-peak2')
            pair_norms.append(norms.values()[i])
        io.write_grouped_csv(
            self.prefix+'-pairs-norm.csv',
            [('All',io.named_list_type(pair_norm_names)(pair_norms))],
            comments=['#Normalization'],
            )


        # Read data
        
        annotations = list(annotation.read_annotations(self.parents))
        if self.utrs:
            utrs = list(annotation.read_annotations(self.utrs))
        else:
            utrs = [ ]
        children = list(annotation.read_annotations(self.children))
        
        count_table = io.read_grouped_table(self.counts, [
            ('Count',int),
            ('Tail',_float_or_none),
            ('Proportion',_float_or_none),
            ('Annotation',str)
            ])
        counts = count_table['Count']
        proportions = count_table['Proportion']
        tails = count_table['Tail']
        
        samples = counts.value_type().keys()
        sample_tags = { }
        for line in count_table.comments:
            if line.startswith('#sampleTags='):
                parts = line[len('#sampleTags='):].split(',')
                assert parts[0] not in sample_tags
                sample_tags[parts[0]] = parts
        
        for item in children:
            item.weight = sum( counts[item.get_id()][name] * float(norms[name]['Normalizing.multiplier']) for name in samples )
        
        parents = [ ]
        id_to_parent = { }
        for item in annotations:
            if item.type != self.parent_type: continue
            assert item.get_id() not in id_to_parent, 'Duplicate id in parent file: '+item.get_id()
            parents.append(item)
            id_to_parent[item.get_id()] = item
            item.children = [ ]
            #item.cds = [ ]
            
            # Default utr
            if item.strand >= 0:
               item.utr_pos = item.end
            else:
               item.utr_pos = item.start
            
            if 'three_prime_UTR_start' in item.attr:
               if item.strand >= 0:
                  item.utr_pos = int(item.attr['three_prime_UTR_start'])-1
               else:
                  item.utr_pos = int(item.attr['three_prime_UTR_start'])
            
            
        
        for item in utrs:
            assert item.attr['Parent'] in id_to_parent, 'Unknown gene '+item.attr['Parent']    
            id_to_parent[item.attr['Parent']].utr_pos = (item.start if item.strand >= 0 else item.end)


        for item in children:
            item.transcription_stop = item.end if item.strand >= 0 else item.start #End of transcription, 0-based, ie between-positions based
            
            if 'Parent' in item.attr:
                for item_parent in item.attr['Parent'].split(','):
                    parent = id_to_parent[item_parent]
                    parent.children.append(item)
                    

        for item in parents:
            item.children.sort(key=_annotation_sorter)
            
            relevant = list(item.children)
            if self.utr_only:
                #if item.strand <= 0:
                #    relative_utr_start = item.end - int(item.attr['three_prime_UTR_start'])
                #else:
                #    relative_utr_start = int(item.attr['three_prime_UTR_start'])-1 - item.start
                #
                #def relative_start(peak):
                #    return item.end-peak.end if item.strand < 0 else peak.start-item.start
                #relevant = [ peak for peak in relevant if relative_start(peak) >= relative_utr_start ]
                relevant = [ 
                    peak for peak in relevant 
                    if (peak.end >= item.utr_pos if item.strand >= 0 else peak.start <= item.utr_pos)
                    ]
                    
            if self.top:
                relevant.sort(key=lambda peak:peak.weight, reverse=True)
                relevant = relevant[:self.top]
            relevant.sort(key=_annotation_sorter)
            item.relevant_children = relevant
        
        
        
        # JSON output
        
        j_data = { }
        j_genes = j_data['genes'] = { }
        
        j_genes['__comment__'] = 'start is 0-based'
        j_genes['name'] = [ ]
        j_genes['chromosome'] = [ ]
        j_genes['strand'] = [ ]
        j_genes['start'] = [ ]
        j_genes['utr'] = [ ]
        j_genes['end'] = [ ]
        j_genes['gene'] = [ ]
        j_genes['product'] = [ ]
        j_genes['peaks'] = [ ]
        j_genes['relevant_peaks'] = [ ]
        #j_genes['cds'] = [ ]
        #j_genes['cds_start'] = [ ]
        #j_genes['cds_end'] = [ ]
        for item in parents:
            j_genes['name'].append( item.get_id() )
            j_genes['chromosome'].append( item.seqid )
            j_genes['strand'].append( item.strand )
            j_genes['start'].append( item.start )
            j_genes['utr'].append( item.utr_pos )
            j_genes['end'].append( item.end )
            j_genes['gene'].append( item.attr.get('Name',item.attr.get('gene','')) )
            j_genes['product'].append( item.attr.get('Product',item.attr.get('product','')) )
            j_genes['peaks'].append( [ item2.get_id() for item2 in item.children ] )
            j_genes['relevant_peaks'].append( [ item2.get_id() for item2 in item.relevant_children ] )
            #j_genes['cds'].append( item.cds )
            #j_genes['cds_start'].append( item.cds_start )
            #j_genes['cds_end'].append( item.cds_end )
        
        j_peaks = j_data['peaks'] = { }
        j_peaks['__comment__'] = 'start is 0-based'
        j_peaks['name'] = [ ]
        j_peaks['chromosome'] = [ ]
        j_peaks['strand'] = [ ]
        j_peaks['start'] = [ ]
        j_peaks['end'] = [ ]
        j_peaks['parents'] = [ ]
        j_peaks['counts'] = [ ]
        j_peaks['tail_lengths'] = [ ]
        j_peaks['proportion_tailed'] = [ ]
        for item in children:
            j_peaks['name'].append( item.get_id() )
            j_peaks['chromosome'].append( item.seqid )
            j_peaks['strand'].append( item.strand )
            j_peaks['start'].append( item.start )
            j_peaks['end'].append( item.end )
            j_peaks['parents'].append( item.attr['Parent'].split(',') if 'Parent' in item.attr else [ ])
            j_peaks['counts'].append( counts[item.get_id()].values() )
            j_peaks['tail_lengths'].append( count_table['Tail'][item.get_id()].values() )
            j_peaks['proportion_tailed'].append( count_table['Proportion'][item.get_id()].values() )
        
        j_samples = j_data['samples'] = { }
        j_samples['name'] = [ ]
        j_samples['tags'] = [ ]
        j_samples['normalizing_multiplier'] = [ ]
        for name in samples:
            j_samples['name'].append(name)
            j_samples['tags'].append(sample_tags[name])
            j_samples['normalizing_multiplier'].append(float(norms[name]['Normalizing.multiplier']))
        
        j_chromosomes = j_data['chromosomes'] = { }
        j_chromosomes['name'] = [ ]
        j_chromosomes['length'] = [ ]
        for name, seq in chromosomes.iteritems():
            j_chromosomes['name'].append(name)
            j_chromosomes['length'].append(len(seq))        
        
        with open(self.prefix + '.json','wb') as f:
            json.dump(j_data, f)
        
        
        # Output paired peak file
        
        output_comments = [ '#Counts' ]
        output_samples = [ ]
        for item in samples:
            output_samples.append(item+'-peak1')
            output_comments.append('#sampleTags=' + ','.join([item+'-peak1','peak1']+sample_tags.get(item,[])))
        for item in samples:
            output_samples.append(item+'-peak2')
            output_comments.append('#sampleTags=' + ','.join([item+'-peak2','peak2']+sample_tags.get(item,[])))
        
        output_names = [ ]
        output_counts = [ ]
        output_proportions = [ ]
        output_tails = [ ]
        output_annotation_fields = [ 'gene', 'product', 'mean_tail_1', 'mean_tail_2', 'chromosome', 'strand', 
                                     'transcription_stops' ] #, 'interpeak_seq', ]
        output_annotations = [ ]
            
        for item in parents:
            peaks = item.relevant_children
            for i in xrange(len(peaks)-1):
                for j in xrange(i+1, len(peaks)):
                    id_i = peaks[i].get_id()
                    id_j = peaks[j].get_id()
                    id_pair = item.get_id() + '-'+id_i+'-'+id_j
                    output_names.append(id_pair)
                    
                    row = [ ]
                    row.extend(counts[id_i].values())
                    row.extend(counts[id_j].values())
                    output_counts.append(filter(_text,row))

                    row = [ ]
                    row.extend(proportions[id_i].values())
                    row.extend(proportions[id_j].values())
                    output_proportions.append(filter(_text,row))

                    row = [ ]
                    row.extend(tails[id_i].values())
                    row.extend(tails[id_j].values())
                    output_tails.append(filter(_text,row))
                    
                    output_annotations.append([
                        item.attr.get('Name',item.attr.get('gene','')),
                        item.attr.get('Product',item.attr.get('product','')),
                        count_table['Annotation'][id_i]['mean-tail'],
                        count_table['Annotation'][id_j]['mean-tail'],
                        
                        item.seqid,
                        str(item.strand),
                        '%d, %d' % (peaks[i].transcription_stop,peaks[j].transcription_stop),
                        #get_interpeak_seq([peaks[i],peaks[j]]),
                        ])
        
        output_count_table = io.named_matrix_type(output_names,output_samples)(output_counts)
        io.write_grouped_csv(
            self.prefix + '-pairs.csv',
            [ 
                ('Count',io.named_matrix_type(output_names,output_samples)(output_counts)),
                ('Proportion',io.named_matrix_type(output_names,output_samples)(output_proportions)),
                ('Tail',io.named_matrix_type(output_names,output_samples)(output_tails)),
                ('Annotation',io.named_matrix_type(output_names,output_annotation_fields)(output_annotations)),
                ],
            comments=output_comments,
            )
                        
        # Chi Sq tests
        
        #for id in relation:
        #    peaks = relation[id]
        #    if len(peaks) < 2: continue     
        
        mats = [ ]   
        genes = [ ]
        products = [ ]
        mean_tails = [ ]
        prop_tails = [ ]

        peak_names = [ ]
        chromosome_names = [ ]
        strands = [ ]
        transcription_stops = [ ]
        interpeak_seqs = [ ]
        prepeak_seqs = [ ]

        for parent in parents:
            id = parent.get_id()
            peaks = parent.relevant_children
            if len(peaks) < 2: continue
            
            matrix = [ ]
            for item in peaks:
                matrix.append(counts[item.get_id()].values())
            
            mats.append(
                runr.R_literal(id) + ' = ' + 
                runr.R_literal(matrix)
                )
            
            genes.append(parent.attr.get('Name',parent.attr.get('gene','')))
            products.append(parent.attr.get('Product',parent.attr.get('product','')))
            
            def format_mean(s):
                if s == 'NA': return 'NA'
                return '%.1f' % float(s)
            mean_tails.append(', '.join( format_mean(count_table['Annotation'][item.get_id()]['mean-tail']) for item in peaks ))
            
            def format_prop(s):
                if s == 'NA': return 'NA'
                return '%.2f' % float(s)
            prop_tails.append(', '.join( format_prop(count_table['Annotation'][item.get_id()]['proportion-with-tail']) for item in peaks ))
            
            peak_names.append(', '.join(item.get_id() for item in peaks))
            chromosome_names.append(parent.seqid)
            strands.append(parent.strand)
            transcription_stops.append(', '.join(str(item.transcription_stop) for item in peaks))
            interpeak_seqs.append(get_interpeak_seq(peaks))
            prepeak_seqs.append(get_prepeak_seq(parent,peaks))
            
            #if len(mats) >= 10: break
        
        text = 'cat("Loading data into R+\n")\n'
        text += 'data <- list(\n' + ',\n'.join(mats) + ')\n'        
        text += CHISQ
        
        runr.run_script(text,
            OUTPUT_FILENAME=self.prefix+'.csv',
            GENES = genes,
            PRODUCTS = products,
            MEAN_TAILS = mean_tails,
            PROP_TAILS = prop_tails,
            PEAK_NAMES = peak_names,
            CHROMOSOME_NAMES = chromosome_names,
            STRANDS = strands,
            TRANSCRIPTION_STOPS = transcription_stops,
            INTERPEAK_SEQS = interpeak_seqs,
            PREPEAK_SEQS = prepeak_seqs,
            )
        
            
        
        
        
        
        
        
        
        