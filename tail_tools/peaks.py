
import collections, os

import nesoni
from nesoni import config, sam, workspace, legion, grace, annotation, span_index


@config.help(
    'Call peaks that are higher than everything within "radius", '
    'and higher than some minimum depth.\n'
    '\n'
    'Note: The --lap parameter can be used to apply some smoothing, '
    'if you are looking for modes that have a bit of width as well as pure height.\n'
    )
@config.Int_flag(
    'lap',
    'Add this many bases to the end of each fragment when calculating depth, '
    'then subtract this many bases from called peaks, back to original size. '
    'Positive values will tend to join up small gaps. '
    'Negative values may allow calling of slightly overlapping transcripts, '
    'and enhance calling of transcripts that are close together. '
    )
@config.String_flag(
    'type',
    'Type of feature to produce in GFF output.'
    )
@config.Int_flag(
    'min_depth',
    'Minimum depth.'
    )
@config.Int_flag(
    'radius',
    'Smaller mode suppression radius.'
    )
@config.Bool_flag(
    'polya', 
    'Only use poly(A) reads.'
    )
@config.Main_section(
    'filenames',
    'Working directories or BAM files (sorted by read name).'
    )
class Find_peaks(config.Action_with_prefix):
    lap = 0
    type = 'peak'
    polya = True
    filenames = [ ]

    min_depth = 5
    radius = 20
    
    def run(self):
        spans = collections.defaultdict(list)
        
        #for item in legion.parallel_imap(self._load_bam, self.filenames):
        #    for key,value in item.items():
        for filename in self.filenames:
            for key,value in self._load_bam(filename).items():
                spans[key].extend(value)

        grace.status('Calling peaks')

        f = open(self.prefix+'.gff', 'wb')
        annotation.write_gff3_header(f)
        
        n = 0

        for (rname, strand), span_list in spans.items():
            length = 1+max( item[1] for item in span_list )
            depth = [ 0.0 ] * length
            AN_total = [ 0.0 ] * length
            AG_total = [ 0.0 ] * length 
            for start, end, AN, AG in span_list:
                depth[start] += 1.0
                depth[end] -= 1.0
                AN_total[start] += AN
                AN_total[end] -= AN
                AG_total[start] += AG
                AG_total[end] -= AG
            
            for i in xrange(1,length):
                depth[i] += depth[i-1]
                AN_total[i] += AN_total[i-1]
                AG_total[i] += AG_total[i-1]

            for start, end in self._find_spans(depth):
                if end-self.lap-start <= 0: continue
                
                n += 1
                
                id = 'peak%d' % n
                
                ann = annotation.Annotation()
                ann.source = 'tailtools'
                ann.type = self.type
                ann.seqid = rname
                ann.start = start
                ann.end = end - self.lap
                
                if ann.end != ann.start+1:
                    self.log.log("%s odd: start %d end %d\n" % (id, ann.start, ann.end))

                ann.strand = strand
                ann.score = None
                ann.phase = None
                ann.attr = { 
                    'id' : id,
                    'n' : str(depth[start+self.lap//2]),
                    'mean_tail' : str(AN_total[start+self.lap//2]/depth[start+self.lap//2]),
                    'mean_genomic' : str(AG_total[start+self.lap//2]/depth[start+self.lap//2]),
                    'color' : '#00ff00' if strand > 0 else '#0000ff' if strand < 0 else '#008080',
                    }
                print >> f, ann.as_gff()
            f.flush()

        f.close()
        
        self.log.datum('-','called peaks',n)
        
        grace.status('')


    def _load_bam(self, filename):
        spans = { }

        if os.path.isdir(filename):
            filename = os.path.join(filename, "alignments_filtered_sorted.bam")
        
        for alignment in sam.Bam_reader(filename):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
        
            if self.polya and not any( item.startswith("AA:i:") for item in alignment.extra ):
                continue
                
            AN = 0.0
            AG = 0.0
            for item in alignment.extra:
                if item.startswith("AN:i:"): AN = float(item[5:])
                if item.startswith("AG:i:"): AG = float(item[5:])
        
            strand = -1 if alignment.flag&sam.FLAG_REVERSE else 1
        
            start = alignment.reference_start
            end = alignment.reference_end

            # 3' end            
            if strand >= 0:
                start = end-1
            else:
                end = start+1
            
            if end+self.lap-start <= 0: continue
            
            rname = alignment.reference_name
            if (rname,strand) not in spans: 
                spans[(rname,strand)] = [ ]          
            spans[(rname, strand)].append((start,end+self.lap, AN,AG))
                
        return spans


    def _find_spans(self, depth):
        result = [ ]
        
        i = 0
        while i < len(depth):
            j = i+1
            while j < len(depth) and depth[j] == depth[i]:
                j += 1
            
            lap = max(0,self.lap-(j-i)+1)
            lap_back = lap//2
            lap_forward = lap-lap_back

            if depth[i] >= self.min_depth:
                for k in xrange(max(0,i-self.radius),min(len(depth),j+self.radius)):
                    if depth[k] > depth[i] or (depth[k] == depth[i] and k < i): #Resolve ties arbitrarily
                        break
                else:
                    result.append((i-lap_back,j+lap_forward))                
            i = j
            
        return result



def join_descriptions(seq, joiner='/'):
    result = [ ]
    for item in seq:
        if item and item not in result: 
            result.append(item)
    return joiner.join(result)


def _extend(feature, extension):
    result = feature.shifted(0, min(extension,int(feature.attr.get("max_extension","10000"))))
    result.parents = feature.parents
    return result

def _three_prime(feature):
    result = feature.three_prime()
    result.parents = feature.parents
    return result


# Note: this now also does some filtering.
@config.String_flag('parent')
@config.String_flag('child')
@config.Int_flag('extension', 'How far downstrand of the gene can the peak be.')
@config.Float_flag('min_tail', 'Minimum tail length to retain peak.')
class Filter_and_relate_peaks_to_genes(config.Action_with_prefix):
    parent = None
    child = None
    extension = None
    min_tail = 0.0
    
    def run(self):
        items = list(annotation.read_annotations(self.parent))
                
        annotation.link_up_annotations(items)
        for item in items:
            assert len(item.parents) <= 1
        
        genes = [ item for item in items if item.type == "gene" ]
        downstrand_genes = [ _extend(_three_prime(item), self.extension) for item in genes ]
        exons = [ item for item in items if item.type == "exon" ]
        utrs = [ _extend(item, self.extension) for item in items if item.type == "three_prime_utr" ]
        
        gene_index = span_index.index_annotations(genes)
        downstrand_gene_index = span_index.index_annotations(downstrand_genes)
        exon_index = span_index.index_annotations(exons)
        utr_index = span_index.index_annotations(utrs)
        
        peaks = [ ]
        for peak in annotation.read_annotations(self.child):
            if float(peak.attr.get("mean_tail","0.0")) < self.min_tail:
                continue
            peaks.append(peak)
        
        for peak in peaks:
            # Query is final base in genome before poly(A) starts
            query = peak.three_prime().shifted(-1,0)
            
            hit_to = "3'UTR"
            hits = [ item.parents[0].parents[0] for item in
                     utr_index.get(query, True) ]

            if not hits:
                hit_to = "Exon"
                hits = [ item.parents[0].parents[0] for item in
                         exon_index.get(query, True) ]
            
            # For non-coding RNAs, which don't have a 3' UTR
            if not hits:
                hit_to = "Downstrand"
                hits = downstrand_gene_index.get(query, True)
            
            if not hits:
                hit_to = "Intron"
                hits = gene_index.get(query, True)

            antisense_hits = gene_index.get(query.reversed(), True)
            if not hits:
                hit_to = "Antisense"
                hits = antisense_hits
            
            if hits:
                peak.attr["Parent"] = join_descriptions([ item.get_id() for item in hits ], ",")
                peak.attr["Relation"] = hit_to
                peak.attr["Name"] = join_descriptions([ item.attr.get("Name","") for item in hits ])
                peak.attr["Product"] = hit_to + " " + join_descriptions([ item.attr.get("Product","") for item in hits ])
                peak.attr["Biotype"] = join_descriptions([ item.attr.get("Biotype","") for item in hits ])
            
            if antisense_hits:
                peak.attr["Antisense_parent"] = join_descriptions([ item.get_id() for item in antisense_hits ], ",")
                peak.attr["Antisense_name"] = join_descriptions([ item.attr.get("Name","") for item in antisense_hits ])
                peak.attr["Antisense_product"] = "Antisense " + join_descriptions([ item.attr.get("Product","") for item in antisense_hits ])
                peak.attr["Antisense_biotype"] = join_descriptions([ item.attr.get("Biotype","") for item in antisense_hits])
                
        
        annotation.write_gff3(self.prefix+"-parent.gff", genes) #Hmm
        annotation.write_gff3(self.prefix+"-child.gff", peaks)



@config.help(
    'Call peaks in depth of coverage indicating transcription end sites.',
    'Note that you must specify --shift-start, and --shift-end.'
    '\n\n'
    'Operates as follows:'
    '\n\n'
    '- Call end points using "nesoni modes: --what 3prime".\n'
    '- Extend called modes back by --peak-length.\n'
    '- Relate resultant peaks to the given annotations.\n'
    '\n'
    'Note: As this is a workflow, you may need to specify "--make-do all" to force everything to recompute if an input file is changed. '
    'Use "--make-do -modes" to recompute everything but the peak calling.'
    )
@config.Int_flag('lap', '--lap value for "nesoni modes:". How fuzzy the pileup of 3\' ends can be when calling a peak.')
@config.Int_flag('radius', '--radius value for "nesoni modes:". How close peaks can be to one another.')
@config.Int_flag('min_depth', '--min-depth value for "nesoni modes:".')
@config.Int_flag('peak_length', 'Number of bases to extend peak back from 3\' end point.')
@config.String_flag('annotations', 'Annotation file. A GFF file containing genes, which contain mRNAs, which contain exons and a three_prime_utr.')
@config.Int_flag('extension', 'How far downstrand of the gene can the peak be.')
@config.Bool_flag('polya', 'Only use poly(A) reads.')
@config.Float_flag('min_tail', 'Minimum average tail length to retain peak.')
@config.Main_section('samples', 'List of sample directories as produced by "analyse-polya:" or "analyse-polya-batch:".')
class Call_peaks(config.Action_with_output_dir):
    lap = 10
    radius = 50
    min_depth = 50
    peak_length = 100
    polya = True
    min_tail = 15.0
    
    annotations = None
    extension = None
    
    samples = [ ]
    
    def run(self):
        assert self.extension is not None, '--extension must be specified'
        assert self.annotations is not None, '--annotations must be specified'
    
        outspace = self.get_workspace()
        working = workspace.Workspace(outspace / 'working', must_exist=False)
        
        Find_peaks(
            working/'modes',
            filenames = self.samples,
            lap = self.lap,
            radius = self.radius,
            min_depth = self.min_depth,
            polya = self.polya,
            ).make()
        
        nesoni.Modify_features(
            working/'peaks',
            working/'modes.gff',
            shift_start = str(-self.peak_length),
            ).make()
        
        Filter_and_relate_peaks_to_genes(
            outspace/'relation',
            parent = self.annotations,
            child = working/'peaks.gff',
            extension = self.extension,
            min_tail = self.min_tail,
            ).make()
            




