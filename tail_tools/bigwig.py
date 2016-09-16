

import heapq, subprocess, os, itertools, json, pickle

import nesoni
from nesoni import io, bio, grace, config, working_directory, sam

from . import web

def pile(spanners, initial=0):
    n = len(spanners)
    heap = [ (0,i,0) for i in xrange(n) ]
    
    pos = 0
    value = initial
    result = [ ]
    while heap:
        new_pos, i, j = heapq.heappop(heap)
        
        if new_pos != pos:
            if result and result[-1][1] == value:
                result[-1] = (result[-1][0]+new_pos-pos, value)
            else:
                result.append( (new_pos-pos, value) )
        pos = new_pos

        spanner = spanners[i]
        if j:
            value -= spanner[j-1][1]
        if j < len(spanner):
            item = spanner[j]
            value += item[1]
            heapq.heappush(heap, (pos+item[0], i,j+1))
    
    return result


def map_spanner(a_func, a_spanner):
    result = [ ]
    for length, value in a_spanner:
        result.append((length, a_func(value)))
    return result


def scale_spanner(scale, a_spanner):
    if scale == 1: return a_spanner
    return map_spanner(lambda x: x*scale, a_spanner)


def clip_spanner(length, spanner):
    remaining = length
    result = [ ]
    for item_length, item_value in spanner:
        if remaining <= 0: break
        if item_length > remaining:
            item_length = remaining
        result.append((item_length, item_value))
        remaining = remaining - item_length
    return result


class Piler(object):
    def __init__(self, length):
        self.length = length
        self.spanners = [[(length,0)]]
    
    def squash_from(self, i):
        if i < len(self.spanners)-1:
            #print "squash", i, len(self.spanners)
            self.spanners[i:] = [ pile(self.spanners[i:]) ]
            #print "done"
    
    def add(self, spanner):
        if not spanner: return
        
        self.spanners.append(spanner)
        
        #if len(self.spanners) < 1000000: return
        
        i = len(self.spanners)-1
        n = len(self.spanners[i])
        while i > 0 and len(self.spanners[i-1]) <= n:
            i -= 1
            n += len(self.spanners[i])
        self.squash_from(i)
   
    def get(self):
        self.squash_from(0)
        return clip_spanner(self.length, self.spanners[0])


def bedgraph(filename, spanners):
    with open(filename,"wb") as f:
        print >> f, "track type=bedGraph"
        for name, spanner in spanners:
            pos = 0
            for length, value in spanner:
                #if value != 0:
                print >> f, "{}\t{}\t{}\t{}".format(name,pos,pos+length,value)
                pos = pos + length



def iter_fragments(alf):
    pool = { }
    for item in alf:
        if item.is_unmapped or item.is_secondary or item.is_supplementary:
            continue
    
        if not item.is_proper_pair or item.mate_is_unmapped:
            yield [item]
            continue
        
        my_key = (item.query_name, item.reference_name, item.reference_start, item.is_reverse)
        if my_key in pool:
            mate = pool[my_key].pop()
            if not pool[my_key]: del pool[my_key]
            yield [ mate, item ]
            continue
        
        mate_key = (item.query_name, item.next_reference_name, item.next_reference_start, item.mate_is_reverse)
        if mate_key not in pool: pool[mate_key] = [ ]
        pool[mate_key].append(item)
    
    if len(pool): print len(pool), "items left in pool"
    for items in pool.itervalues():
        for item in items:
            yield [item]



def make_bigwig(prefix, bam_filenames, make_spanner, fragments=False, stop_after=None, scale=1.0): 
    have_pysam = False
    try:
        import pysam
        have_pysam = True
    except ImportError:
        pass
    
    #alf = pysam.AlignmentFile(bam_filenames[0])
    #header = alf.header
    header = sam.parsed_bam_headers(bam_filenames[0])
    
    with open(prefix+"-chrom.sizes","wb") as f:
        for entry in header["SQ"]:
            f.write("{}\t{}\n".format(entry["SN"],entry["LN"]))
    
    chrom_names = [ entry["SN"] for entry in header["SQ"] ]
    chrom_sizes = [ int(entry["LN"]) for entry in header["SQ"] ]

    #alf.close()

    forward = dict([ (i,Piler(j)) for i,j in zip(chrom_names,chrom_sizes) ])
    reverse = dict([ (i,Piler(j)) for i,j in zip(chrom_names,chrom_sizes) ])
    
    old = grace.status("Bigwig")

    for filename in bam_filenames:
        if have_pysam:
            alf = pysam.AlignmentFile(filename)
        else:
            alf = sam.Bam_reader(filename)
        
        n = 0
        
        if not fragments:
            for item in alf:
                if item.is_unmapped or item.is_secondary or item.is_supplementary:
                    continue
            
                # Assume --> <-- oriented read pairs
                which = forward if bool(item.is_reverse) == bool(item.is_read2) else reverse
                which[item.reference_name].add( make_spanner(item) )
                    
                n += 1
                if stop_after is not None and n > stop_after: break
                if n % 1000000 == 0: grace.status(os.path.basename(prefix)+" "+filename+" "+grace.pretty_number(n))
        
        else:
            for item in iter_fragments(alf):
                # Assume --> <-- oriented read pairs
                which = forward if bool(item[0].is_reverse) == bool(item[0].is_read2) else reverse
                which[item[0].reference_name].add( make_spanner(item) )
                    
                n += 1
                if stop_after is not None and n > stop_after: break
                if n % 1000000 == 0: grace.status(os.path.basename(prefix)+" "+filename+" "+grace.pretty_number(n))
        
        if have_pysam:
            alf.close()

    bedgraph(prefix+"-fwd.bedgraph", zip(chrom_names, [ scale_spanner(scale, forward[item].get()) for item in chrom_names ]))
    subprocess.check_call([
        "wigToBigWig",prefix+"-fwd.bedgraph",prefix+"-chrom.sizes",prefix+"-fwd.bw"])
    os.unlink(prefix+"-fwd.bedgraph")
        
    bedgraph(prefix+"-rev.bedgraph", zip(chrom_names, [ scale_spanner(scale, reverse[item].get()) for item in chrom_names ]))    
    subprocess.check_call([
        "wigToBigWig",prefix+"-rev.bedgraph",prefix+"-chrom.sizes",prefix+"-rev.bw"])
    os.unlink(prefix+"-rev.bedgraph")
    
    os.unlink(prefix+"-chrom.sizes")
    
    grace.status(old)


def read2_starts(item):
    if not item.is_read2:
        return None

    pos = item.reference_end-1 if item.is_reverse else item.reference_start
    return ((pos,0),(1,1))

def read1_starts(item):
    if item.is_read2:
        return None

    pos = item.reference_end-1 if item.is_reverse else item.reference_start
    return ((pos,0),(1,1))

def read_starts(item):
    pos = item.reference_end-1 if item.is_reverse else item.reference_start
    return ((pos,0),(1,1))

def read_ends(item):
    pos = item.reference_start if item.is_reverse else item.reference_end-1
    return ((pos,0),(1,1))


def coverage(item):
    pos = 0
    result = [ ]
    blocks = item.get_blocks()
    for start, end in blocks:
        assert start >= pos
        result.append((start-pos,0))
        result.append((end-start,1))
        pos = end
    return result


def fragment_split_coverage(items):
    """ Don't count twice for overlapping reads in a pair. """
    pos = 0
    result = [ ]
    blocks = [ ]
    for item in items: 
        blocks.extend(item.get_blocks())
    blocks.sort()
    for start, end in blocks:
        if start > pos:
            result.append((start-pos,0))
            pos = start
        if end > pos:
            result.append((end-pos,1))
            pos = end
    return result

def fragment_coverage(items):
    start = min( item.reference_start for item in items )
    end = max( item.reference_end for item in items )
    return ((start,0),(end-start,1))



def make_ambiguity_bigwig(prefix, bam_filenames, stop_after=None, subsample=1): 
    import pysam
    
    #alf = pysam.AlignmentFile(bam_filenames[0])
    #header = alf.header
    header = sam.parsed_bam_headers(bam_filenames[0])
    
    with open(prefix+"-chrom.sizes","wb") as f:
        for entry in header["SQ"]:
            f.write("{}\t{}\n".format(entry["SN"],entry["LN"]))
    
    chrom_names = [ entry["SN"] for entry in header["SQ"] ]
    chrom_sizes = [ int(entry["LN"]) for entry in header["SQ"] ]

    #alf.close()

    unambiguous = [ Piler(i) for i in chrom_sizes ]
    total = [ Piler(i) for i in chrom_sizes ]

    for filename in bam_filenames:
        alf = pysam.AlignmentFile(filename)
        n = 0
        
        sub = subsample-1
        for item in alf:
            if item.is_unmapped or item.is_supplementary:
                continue
        
            sub = (sub + 1) % subsample
            if sub: continue
        
            spanner = fragment_split_coverage([item])
            total[item.reference_id].add(spanner)
            if item.get_tag("NH") == 1:
                unambiguous[item.reference_id].add(spanner)
                
            n += 1
            if stop_after is not None and n > stop_after: break
            if n % 1000000 == 0: print prefix, filename, grace.pretty_number(n)
        
        alf.close()


    ambiguities = [ ]
    for i in xrange(len(total)):
        u = unambiguous[i].get()
        t = map_spanner(lambda x: x*1j, total[i].get())
        c = pile([u,t],initial=0.0)
        c = map_spanner(lambda x: max(0.0,x.imag-x.real)/max(x.imag,1.0), c)
        ambiguities.append(c)

    bedgraph(prefix+".bedgraph", zip(chrom_names, [ item for item in ambiguities ]))
    subprocess.check_call([
        "wigToBigWig",prefix+".bedgraph",prefix+"-chrom.sizes",prefix+".bw"])
    os.unlink(prefix+".bedgraph")
    os.unlink(prefix+"-chrom.sizes")



def make_ambiguity_bigwig_by_readname(prefix, bam_filenames, stop_after=None, subsample=1): 
    #import pysam
    
    #alf = pysam.AlignmentFile(bam_filenames[0])
    #header = alf.header
    header = sam.parsed_bam_headers(bam_filenames[0])
    
    with open(prefix+"-chrom.sizes","wb") as f:
        for entry in header["SQ"]:
            f.write("{}\t{}\n".format(entry["SN"],entry["LN"]))
    
    chrom_names = [ entry["SN"] for entry in header["SQ"] ]
    chrom_sizes = [ int(entry["LN"]) for entry in header["SQ"] ]

    #alf.close()

    unambiguous = dict([ (i,Piler(j)) for i,j in zip(chrom_names,chrom_sizes) ])
    total = dict([ (i,Piler(j)) for i,j in zip(chrom_names,chrom_sizes) ])
    
    old = grace.status("Ambiguity bigwig")

    for filename in bam_filenames:
        #alf = pysam.AlignmentFile(filename)
        alf = sam.Bam_reader(filename)
        n = 0
        
        sub = subsample-1
        for (key,items) in itertools.groupby(alf, lambda item: item.query_name):
            sub = (sub + 1) % subsample
            if sub: continue
        
            items = [ item for item in items if not item.is_unmapped and not item.is_supplementary ]
            if not items:
                continue
            
            # Only use top scoring alignments
            AS = [ item.get_AS() for item in items ]
            best_AS = max(AS)
            items = [ item for item, this_AS in zip(items,AS) if this_AS >= best_AS ]
        
            for item in items:
                #spanner = fragment_split_coverage([item])
                spanner = fragment_coverage([item])        #TODO fixme when blocks available
                spanner = scale_spanner(1.0/len(items), spanner)
                total[item.reference_name].add(spanner)
                if len(items) == 1:
                    unambiguous[item.reference_name].add(spanner)
                
            n += 1
            if stop_after is not None and n > stop_after: break
            if n % 1000000 == 0: grace.status(os.path.basename(prefix)+" "+filename+" "+grace.pretty_number(n))
        
        alf.close()


    ambiguities = [ ]
    for i in xrange(len(total)):
        u = unambiguous[chrom_names[i]].get()
        t = map_spanner(lambda x: x*1j, total[chrom_names[i]].get())
        c = pile([u,t],initial=0.0)
        c = map_spanner(lambda x: max(0.0,x.imag-x.real)/max(x.imag,1.0), c)
        ambiguities.append(c)

    bedgraph(prefix+".bedgraph", zip(chrom_names, [ item for item in ambiguities ]))
    subprocess.check_call([
        "wigToBigWig",prefix+".bedgraph",prefix+"-chrom.sizes",prefix+".bw"])
    os.unlink(prefix+".bedgraph")
    os.unlink(prefix+"-chrom.sizes")
    
    grace.status(old)



@config.help("Produce various bigwig depth of coverage files from a BAM file.", """\
This tool uses PySam, so can not currently be used with PyPy.

Bigwig files produced:

cover - Depth of coverage by actually sequenced bases.
span - Depth of coverage of region spanned by reads or fragments.
start - Fragment start locations (5' end of read 1).
end - Fragment end locations (5' end of read 2).
5p - Read start locations.
3p - Read end locations.
ambiguity - What proportion of reads are multi-mappers at each base, using NH attribute in BAM file.

Note: "cover" alone requires the pysam library. 

""")
@config.Main_section("bam_files")
@config.String_flag("what", "What bigwig files to actually produce. Comma separated list.")
@config.Int_flag("subsample", "(currently for ambiguity plots only) Subsample alignments by this factor.")
@config.Float_flag("scale", "Scale output by this (eg a normalizing multiplier).")
class Bam_to_bigwig(config.Action_with_prefix):
    what = "cover,span,start,end,ambiguity"
    subsample = 1
    scale = 1.0
    bam_files = [ ]
    
    def run(self):
        with nesoni.Stage() as stage:
            for item in self.what.split(","):
                if item == "cover":
                    stage.process(make_bigwig,
                        self.prefix + "-cover", self.bam_files, fragment_split_coverage, True, scale=self.scale)
                elif item == "span":
                    stage.process(make_bigwig,
                        self.prefix + "-span", self.bam_files, fragment_coverage, True, scale=self.scale)
                elif item == "start":
                    stage.process(make_bigwig,
                        self.prefix + "-start", self.bam_files, read1_starts, False, scale=self.scale)
                elif item == "end":
                    stage.process(make_bigwig,
                        self.prefix + "-end", self.bam_files, read2_starts, False, scale=self.scale)
                elif item == "5p":
                    stage.process(make_bigwig,
                        self.prefix + "-5p", self.bam_files, read_starts, False, scale=self.scale)
                elif item == "3p":
                    stage.process(make_bigwig,
                        self.prefix + "-3p", self.bam_files, read_ends, False, scale=self.scale)
                elif item == "ambiguity":
                    stage.process(make_ambiguity_bigwig,
                        self.prefix + "-ambiguity", self.bam_files, subsample=self.subsample)
                else:
                    raise config.Error("Don't know how to make: "+item)


@config.help("Produce ambiguity bigwig from a BAM file.", """\
""")
@config.Main_section("bam_files")
@config.Int_flag("subsample", "Subsample alignments by this factor.")
@config.Bool_flag("by_readname", "Determine ambiguity by read name, in a BAM file sorted by readname. If no, the NH attribute is used.")
class Bam_ambiguity(config.Action_with_prefix):
    what = "cover,span,start,end,ambiguity"
    subsample = 1
    by_readname = True
    bam_files = [ ]
    
    def run(self):
        (make_ambiguity_bigwig_by_readname if self.by_readname else make_ambiguity_bigwig) \
            (self.prefix, self.bam_files, subsample=self.subsample)


@config.help("Create a set of bigwig files based on pipeline output.")
@config.String_flag("norm_file", "File of normalizations produced by \"nesoni norm-from-counts:\" or \"nesoni norm-from-samples:\".")
@config.String_flag("peaks_file", "For convenience the generated page can also load a peaks GFF file.")
@config.String_flag("title", "Title for HTML page.")
@config.Main_section("working_dirs", "Working directories or pipeline output directory (from \"analyse-polya-batch:\").")
class Polya_bigwigs(config.Action_with_output_dir):
    pipeline_dir = None
    norm_file = ""
    peaks_file = None
    title = "IGV tracks"
    
    def run(self):
        working_dirs = [ ] 
        peaks_file = self.peaks_file       
        for item in self.working_dirs:
            state_filename = os.path.join(item,'analyse-polya-batch.state')
            if not os.path.exists(state_filename):
                working_dirs.append(item)
            else:
                with open(state_filename,'rb') as f:
                    state = pickle.load(f)

                for sample in state.samples:
                    working_dirs.append(os.path.join(item,'samples',sample.output_dir))
                
                if not peaks_file:
                    peaks_file = os.path.join(self.pipeline_dir, "peaks", "relation-child.gff")

        #state_filename = os.path.join(self.pipeline_dir,'analyse-polya-batch.state')
        #with open(state_filename,'rb') as f:
        #    state = pickle.load(f)

        #sample_names = [ ]
        #working_dirs = [ ]
        #for sample in state.samples:
        #    sample_names.append(sample.output_dir)
        #    working_dirs.append(os.path.join(self.pipeline_dir,'samples',sample.output_dir))
        
        sample_names = [ os.path.split(dirname)[1] for dirname in working_dirs ]
        workspaces = [ working_directory.Working(dirname, must_exist=True) for dirname in working_dirs ]
        workspaces_polya = [ working_directory.Working(dirname+"-polyA", must_exist=True) for dirname in working_dirs ]
        
        workspace = self.get_workspace()
        
        with open(workspace/"index.html","wb") as f:
            web.emit(f, "igv.html", dict(
                SAMPLES = json.dumps(sample_names),
                HAVE_NORM = json.dumps(bool(self.norm_file)),
                TITLE = self.title,
            ))
        
        bams = [ item/"alignments_filtered_sorted.bam" for item in workspaces ]
        bams_raw = [ item/"alignments.bam" for item in workspaces ]
        bams_polya = [ item/"alignments_filtered_sorted.bam" for item in workspaces_polya ]

        
        for i in xrange(len(sample_names)):
            io.symbolic_link(bams[i], workspace/(sample_names[i]+".bam"))
            io.symbolic_link(bams[i]+".bai", workspace/(sample_names[i]+".bam.bai"))
        
        io.symbolic_link(peaks_file, workspace/"peaks.gff")

        
        #print "SKIP!  TODO: bams"
        #return
                
        if self.norm_file:
            mults = io.read_grouped_table(self.norm_file)['All']
            norm_mult = [ float(mults[name]['Normalizing.multiplier']) for name in sample_names ]
        
        with nesoni.Stage() as stage:
            Bam_ambiguity(workspace/"ambiguity", bam_files=bams_raw).process_make(stage)
            
            Bam_to_bigwig(workspace/"total-all", bam_files=bams, what="span,3p",
                ).process_make(stage)
            Bam_to_bigwig(workspace/"total-polya", bam_files=bams_polya, what="span,3p",
                ).process_make(stage)
            
            for i in xrange(len(sample_names)):
                for scale_desc, scale in [("raw",1.0)]+([("norm",norm_mult[i])] if self.norm_file else []):
                    for bam_desc, bam in [("all", bams[i]), ("polya", bams_polya[i])]:
                        Bam_to_bigwig(
                            workspace/(sample_names[i]+"-"+bam_desc+"-"+scale_desc), 
                            bam_files=[bam], what="span,3p", scale=scale
                            ).process_make(stage)









