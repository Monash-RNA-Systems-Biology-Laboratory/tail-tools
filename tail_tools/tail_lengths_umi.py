
# Variants for UMI counting

import os.path, itertools, collections, math

from nesoni import annotation, sam, span_index, config, grace, working_directory, workspace, io, legion

import cPickle as pickle


def str_na(value):
    if value is None:
        return "NA"
    else:
        return str(value)


@config.help("""\
Create file to be used by "aggregate-tail-lengths-umi:".\
""","""\
Reads are aligned to "parts" features. A parent is then sought of type "types", possibly several levels up, or possibly zero levels up.

If part features have a "max_extension" attribute, this is respected when extending them. Typically this is used to avoid extending into a following CDS.
""")
@config.String_flag('annotations', 'Filename containing annotations. Defaults to annotations in reference directory.')
@config.String_flag('types', 'Comma separated list of feature types to use.')
@config.String_flag('parts', 'Comma separated list of feature types that make up features. Defaults to types if blank.')
@config.Int_flag('extension', 'How far downstrand of the given annotations a read or peak belonging to a gene might be.')
@config.Positional('working_dir', 'Working directory to use as input.')
class Tail_count_umi(config.Action_with_prefix):
     annotations = None
     types = 'gene'
     parts = 'exon'
     working_dir = None
     
     extension = None
     
     ##Memory intensive, hack to run with reduced parallelism
     ## as nesoni doesn't have any memory usage management, only core usage
     #def cores_required(self):
     #    return min(4, legion.coordinator().get_cores())

     def run(self):
         assert self.extension is not None, '--extension must be specified'
     
         #workspace = self.get_workspace()
         workspace = working_directory.Working(self.working_dir, must_exist=True)
         if self.annotations == None:
             reference = workspace.get_reference()
             annotations_filename = reference.annotations_filename()
         else:
             annotations_filename = self.annotations
         
         types = [ item.lower() for item in self.types.split(',') ]
         
         parts = self.parts or self.types 
         parts = [ item.lower() for item in parts.split(',') ]
         
         
         all_annotations = list(annotation.read_annotations(annotations_filename))
         annotation.link_up_annotations(all_annotations)
         for item in all_annotations: 
             item.primary = None
     
         annotations = [ 
             item 
             for item in all_annotations
             if item.type.lower() in types
         ]
         
         part_annotations = [ ]
         seen = set()
         queue = [ (item,item) for item in annotations ]
         while queue:
             primary, item = queue.pop()
             if item.type.lower() in parts:
                 assert item.primary is None, "Feature with multiple parents"
                 item.primary = primary
                 key = (id(primary),item.start,item.end,item.seqid,item.strand)
                 # Ignore duplicate exons (many isoforms will have the same exons)
                 if key not in seen:
                     seen.add(key)
                     part_annotations.append(item)
             queue.extend( (primary, item2) for item2 in item.children )
         
         del seen
         del all_annotations
         
         self.log.log('%d annotations\n' % len(annotations))
         self.log.log('%d part annotations\n' % len(part_annotations))
         
         #assert annotations, 'No annotations of specified types in file'
         
         for item in part_annotations:
             this_extension = self.extension
             if "max_extension" in item.attr:
                 this_extension = min(this_extension,int(item.attr["max_extension"]))
                 
             if item.strand >= 0:
                 item.tail_pos = item.end
                 item.end += this_extension
             else:
                 item.tail_pos = item.start
                 item.start -= this_extension
         
         for item in annotations:    
             item.hits = [] # [ (tail_length, adaptor_bases) ]
         
         index = span_index.index_annotations(part_annotations)
         
         for alignment in sam.Bam_reader(workspace/'alignments_filtered_sorted.bam'):
             if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                 continue
             
             parts = alignment.qname.rsplit("_",1)
             assert len(parts) == 2, "Read name missing UMI: "+alignment.qname
             umi = parts[1]
             
             start = alignment.reference_start
             end = alignment.reference_end
             alignment_length = end-start
             strand = -1 if alignment.flag&sam.FLAG_REVERSE else 1
             fragment_feature = annotation.Annotation(
                 seqid=alignment.reference_name,
                 start=start,
                 end=end,
                 strand=strand
                 )
             
             if strand >= 0:
                 tail_pos = end
             else:
                 tail_pos = start
             
             tail_length = 0
             adaptor_bases = 0
             for item in alignment.extra:
                 if item.startswith('AN:i:'):
                     tail_length = int(item[5:])
                 elif item.startswith('AD:i:'):
                     adaptor_bases = int(item[5:])
             
             hits = index.get(fragment_feature, same_strand=True)
             if hits:
                 gene = min(hits, key=lambda gene: 
                     (abs(tail_pos - gene.tail_pos), gene.primary.get_id()))
                     # Nearest by tail_pos
                     # failing that, by id to ensure a deterministic choice
                 
                 gene.primary.hits.append( (tail_length,adaptor_bases,umi) )
         
         for item in annotations:
             del item.parents
             del item.children
             del item.primary
         
         f = io.open_possibly_compressed_writer(self.prefix + '.pickle.gz')
         pickle.dump((workspace.name, workspace.get_tags(), annotations), f, pickle.HIGHEST_PROTOCOL)
         f.close()
         


@config.help(
'Aggregate data collected by "tail-lengths:" and produce various CSV tables, UMI version.',
"""\
"""
)
@config.Int_flag('tail',
    'Minimum tail length to count as having a tail.'
    )
@config.Int_flag('adaptor',
    'Minimum number of adaptor bases required, 0 for no filtering.'
    )
@config.Int_flag('clip_tail',
    'Tails longer than this will be reduced to this length. 0 for no clipping.'
    )
@config.Main_section('pickles')
class Aggregate_tail_counts_umi(config.Action_with_output_dir):             
    tail = 4
    adaptor = 0
    clip_tail = 0
    pickles = [ ]
     
    #Memory intensive, don't run in parallel
    def cores_required(self):
        return legion.coordinator().get_cores()
         
    def run(self):
        assert len(self.pickles) > 0, "No samples to count."
        
        work = self.get_workspace()
        
        data = [ ]
        names = [ ]
        sample_tags = [ ]
        
        old = grace.status("Loading pickles")
        
        max_length = 1
        for i, item in enumerate(self.pickles):
            grace.status("Loading "+os.path.basename(item))
            f = io.open_possibly_compressed_file(item)
            name, tags, datum = pickle.load(f)
            f.close()
            data.append(datum)
            names.append(name)
            sample_tags.append(tags)
            
            if self.clip_tail:
                for feature in datum:
                    feature.hits = [ 
                        (min(self.clip_tail,item2[0]),item2[1],item2[2]) 
                        for item2 in feature.hits 
                        ]
            
            try:
                max_length = max(max_length, max( 
                    item2[0] #tail_length
                    for feature in datum
                    for item2 in feature.hits
                    ) + 1)
            except ValueError:
                pass
            
            if i == 0:
               annotations = datum
        
        grace.status(old)
        
        self.log.log("Maximum tail length %d\n" % max_length)
        
        for i in xrange(len(names)):        
            n_alignments = 0.0
            n_umis = 0.0
            for feature in data[i]:
                
                # TODO: error correction?
                umi_bins = { }
                for tail_length, adaptor_bases, umi in feature.hits:
                    if adaptor_bases < self.adaptor:
                        continue
                    if umi not in umi_bins:
                        umi_bins[umi] = [ ]
                    umi_bins[umi].append(tail_length)
                
                feature.tail_counts = [ 0.0 ] * max_length
                for bin in umi_bins.values():
                    for tail_length in bin:
                        feature.tail_counts[tail_length] += 1.0/len(bin)
                
                feature.total_count = len(umi_bins)
                feature.alignment_count = len(feature.hits)
                
                n_umis += feature.total_count
                n_alignments += feature.alignment_count
                
                del feature.hits
            
            self.log.datum(names[i], 'Alignments to features', int(n_alignments))
            self.log.datum(names[i], 'UMIs to features', int(n_umis))
        
        
        counts = [ ]  # [feature][sample](total_count, [taillength])
        
        for item in data: 
            assert len(item) == len(data[0])
        for row in itertools.izip(*data):
            this_counts = [ (item.total_count, item.alignment_count, item.tail_counts) for item in row ]
            counts.append(this_counts)
        
        n_features = len(counts)
        n_samples = len(data)
        
        sample_n = [ [0.0]*n_samples for i in xrange(n_features) ]        # [feature][sample]  Total UMI count
        sample_alignments = [ [0.0]*n_samples for i in xrange(n_features) ]        # [feature][sample]  Total alignment count
        sample_n_tail = [ [0.0]*n_samples for i in xrange(n_features) ]   # [feature][sample]  Polya UMI count
        sample_prop = [ [None]*n_samples for i in xrange(n_features) ]    # [feature][sample]  Proportion of UMIs with tail (deprecated)
        sample_tail = [ [None]*n_samples for i in xrange(n_features) ]    # [feature][sample]  Mean tail length in each sample
        sample_sd_tail = [ [None]*n_samples for i in xrange(n_features) ] # [feature][sample]  Std dev tail length in each sample
        sample_total_tail = [ [0.0]*n_samples for i in xrange(n_features) ]
        
        sample_quantile_tail = collections.OrderedDict( 
            (item, [ [None]*n_samples for i in xrange(n_features) ]) 
            for item in [25,50,75,100]
            )
        
        overall_n = [ 0.0 ]*n_features       # [feature]          Overall UMI count
        overall_prop = [ None ]*n_features   # [feature]          Overall proportion with tail
        overall_tail = [ None ]*n_features   # [feature]          Overall mean tail length
        overall_n_tail = [ 0.0 ]*n_features  # [feature]          Overall polya count
        for i, row in enumerate(counts):
            for j, (this_this_n, this_this_alignments, item) in enumerate(row):
                sample_n[i][j] = this_this_n
                sample_alignments[i][j] = this_this_alignments
                sample_n_tail[i][j] = sum(item[self.tail:])
                sample_total_tail[i][j] = sum( item[k]*k for k in xrange(self.tail,max_length) )

                if sample_n[i][j] >= 1:
                    sample_prop[i][j] = sample_n_tail[i][j]/sample_n[i][j]
                
                if sample_n_tail[i][j] >= 1:
                    sample_tail[i][j] = sample_total_tail[i][j]/sample_n_tail[i][j]
                
                    for quantile in sample_quantile_tail:
                        counter = sample_n_tail[i][j] * quantile / 100.0
                        for k in xrange(self.tail, max_length):
                            counter -= item[k]
                            if counter <= 0.0: break
                        sample_quantile_tail[quantile][i][j] = k
                
                if sample_n_tail[i][j] >= 2:
                    sample_sd_tail[i][j] = math.sqrt(
                        sum( item[k]*((k-sample_tail[i][j])**2) for k in xrange(self.tail,max_length) )
                        / (sample_n_tail[i][j]-1)
                        )
                    
            overall_n[i] = sum(sample_n[i])
            overall_n_tail[i] = sum(sample_n_tail[i])
            if overall_n[i] >= 1:
                overall_prop[i] = sum(sample_n_tail[i])/overall_n[i]
            if overall_n_tail[i] >= 1:
                overall_tail[i] = sum(sample_total_tail[i])/overall_n_tail[i]
             
        for i, name in enumerate(names):
            this_total = sum( item[i] for item in sample_total_tail )
            this_n = sum( item[i] for item in sample_n_tail )
            if this_n:
                self.log.datum(name, 'Average poly-A tail', this_total/this_n)
                
        for i, name in enumerate(names):
            this_total = sum( item[i] for item in sample_n_tail )
            this_n = sum( item[i] for item in sample_n )
            if this_n:
                self.log.datum(name, 'Average proportion of UMIs with tail', this_total/this_n)
        
        with open(work/'features-with-data.gff','wb') as f:
            annotation.write_gff3_header(f)
            for i, item in enumerate(annotations):
                item.attr['reads'] = str(overall_n[i])
                item.attr['reads_with_tail'] = str(overall_n_tail[i])
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
            '"Tail_count" group is number of UMIs with tail',
            '"Tail" group is mean tail per sample',
            '"Proportion" group is proportion of UMIs with tail',
            ]
            
        have_biotype = any("Biotype" in item.attr for item in annotations)
        have_parent = any("Parent" in item.attr for item in annotations)
        have_relation = any("Relation" in item.attr for item in annotations)
        have_antisense = any("Antisense_parent" in item.attr for item in annotations)

        def counts_iter():
            for i in xrange(n_features):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(n_samples):
                    row[('Count',names[j])] = '%d' % sample_n[i][j]
                
                row[('Annotation','Length')] = annotations[i].end - annotations[i].start
                row[('Annotation','gene')] = annotations[i].attr.get('Name','')
                row[('Annotation','product')] = annotations[i].attr.get('Product','')
                if have_biotype:
                    row[('Annotation','biotype')] = annotations[i].attr.get('Biotype','')
                if have_parent:
                    row[('Annotation','parent')] = annotations[i].attr.get('Parent','')
                if have_relation:
                    row[('Annotation','relation')] = annotations[i].attr.get('Relation','')
                
                if have_antisense:
                    row[('Annotation','antisense_gene')] = annotations[i].attr.get('Antisense_name','')
                    row[('Annotation','antisense_product')] = annotations[i].attr.get('Antisense_product','')
                    row[('Annotation','antisense_biotype')] = annotations[i].attr.get('Antisense_biotype','')
                    row[('Annotation','antisense_parent')] = annotations[i].attr.get('Antisense_parent','')
                
                row[('Annotation','chromosome')] = str(annotations[i].seqid)
                row[('Annotation','strand')] = str(annotations[i].strand)
                row[('Annotation','start')] = str(annotations[i].start+1)
                row[('Annotation','end')] = str(annotations[i].end)
                
                row[('Annotation','reads')] = str(overall_n[i])
                row[('Annotation','reads-with-tail')] = str(overall_n_tail[i])
                row[('Annotation','mean-tail')] = str_na(overall_tail[i])
                row[('Annotation','proportion-with-tail')] = str_na(overall_prop[i])
                for j in xrange(n_samples):
                    row[('Tail_count',names[j])] = '%d' % sample_n_tail[i][j]
                for j in xrange(n_samples):
                    row[('Tail',names[j])] = str_na(sample_tail[i][j])
                for j in xrange(n_samples):
                    row[('Tail_sd',names[j])] = str_na(sample_sd_tail[i][j])
                
                for quantile in sample_quantile_tail:
                    for j in xrange(n_samples):
                        row[('Tail_quantile_%d'%quantile,names[j])] = str_na(sample_quantile_tail[quantile][i][j])                    
                
                for j in xrange(len(names)):
                    row[('Proportion',names[j])] = str_na(sample_prop[i][j])
                
                for j in xrange(n_samples):
                    row[('Alignments',names[j])] = '%d' % sample_alignments[i][j]
                
                yield row
        io.write_csv(work/'counts.csv', counts_iter(), comments=comments)
        
        
        def write_csv_matrix(filename, matrix):
            def emitter():
                for i in xrange(n_features):
                    row = collections.OrderedDict()
                    row["Feature"] = annotations[i].get_id()
                    for j in xrange(n_samples):
                        row[names[j]] = str_na(matrix[i][j])
                    yield row
            io.write_csv(filename, emitter())
        
        write_csv_matrix(work/'umi_count.csv', sample_n)
        write_csv_matrix(work/'alignment_count.csv', sample_alignments)
        write_csv_matrix(work/'tail_count.csv', sample_n_tail)
        write_csv_matrix(work/'tail.csv', sample_tail)
        write_csv_matrix(work/'tail_sd.csv', sample_sd_tail)
        for quantile in sample_quantile_tail:
            write_csv_matrix(work/('tail_quantile_%d.csv'%quantile), sample_quantile_tail[quantile])
        
        
        #def raw_columns():
        #    for i in xrange(n_samples):
        #        row = collections.OrderedDict()
        #        row['Sample'] = names[i]
        #        for j in xrange(max_length):
        #            row['length-%d' % j] = str(i*max_length+j+1) #For R+, so 1 based
        #        yield row
        #io.write_csv(work/'raw-columns.csv', raw_columns())
        #
        ##Somewhat inefficient        
        #def raw():
        #    for i in xrange(n_features):
        #        row = collections.OrderedDict()
        #        row['Feature'] = annotations[i].get_id()
        #        for j in xrange(n_samples):
        #            for k in xrange(max_length):
        #                row['%d %s' % (k,names[j])] = str( counts[i][j][1][k] )
        #        yield row
        #io.write_csv(work/'raw.csv', raw())
        
        def pooled():
            for i in xrange(n_features):
                row = collections.OrderedDict()
                row['Feature'] = annotations[i].get_id()
                for j in xrange(max_length):
                    row[str(j)] = str( sum( counts[i][k][2][j] for k in xrange(n_samples) ) )
                yield row
        io.write_csv(work/'pooled.csv', pooled())

