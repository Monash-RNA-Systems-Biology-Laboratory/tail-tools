

import heapq, subprocess, os

import nesoni
from nesoni import bio, config

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

class Piler(object):
    def __init__(self, length):
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
        return self.spanners[0]


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
        if not item.is_proper_pair or item.mate_is_unmapped:
            yield [item]
        
        my_key = (item.query_name, item.reference_id, item.reference_start, item.is_reverse)
        if my_key in pool:
            mate = pool[my_key].pop()
            if not pool[my_key]: del pool[my_key]
            yield [ mate, item ]
        else:
            mate_key = (item.query_name, item.next_reference_id, item.next_reference_start, item.mate_is_reverse)
            if mate_key not in pool: pool[mate_key] = [ ]
            pool[mate_key].append(item)
    
    print len(pool), "items left in pool"
    for items in pool.itervalues():
        for item in items:
            yield [item]



def make_bigwig(prefix, bam_filenames, make_spanner, fragments=False, stop_after=None): 
    import pysam
    
    alf = pysam.AlignmentFile(bam_filenames[0])
    
    with open(prefix+"-chrom.sizes","wb") as f:
        for entry in alf.header["SQ"]:
            f.write("{}\t{}\n".format(entry["SN"],entry["LN"]))
    
    chrom_names = [ entry["SN"] for entry in alf.header["SQ"] ]
    chrom_sizes = [ entry["LN"] for entry in alf.header["SQ"] ]

    alf.close()

    forward = [ Piler(i) for i in chrom_sizes ]
    reverse = [ Piler(i) for i in chrom_sizes ]

    for filename in bam_filenames:
        alf = pysam.AlignmentFile(filename)
        n = 0
        
        if not fragments:
            for item in alf:
                if item.is_secondary or item.is_unmapped:
                    continue
            
                # Assume --> <-- oriented read pairs
                which = forward if bool(item.is_reverse) == bool(item.is_read2) else reverse
                which[item.reference_id].add( make_spanner(item) )
                    
                n += 1
                if stop_after is not None and n > stop_after: break
                if n % 1000000 == 0: print prefix, n
        
        else:
            for item in iter_fragments(alf):
                if item[0].is_secondary or item[0].is_unmapped:
                    continue
            
                # Assume --> <-- oriented read pairs
                which = forward if bool(item[0].is_reverse) == bool(item[0].is_read2) else reverse
                which[item[0].reference_id].add( make_spanner(item) )
                    
                n += 1
                if stop_after is not None and n > stop_after: break
                if n % 1000000 == 0: print prefix, n
        
        alf.close()

    bedgraph(prefix+"-fwd.bedgraph", zip(chrom_names, [ item.get() for item in forward ]))
    subprocess.check_call([
        "wigToBigWig",prefix+"-fwd.bedgraph",prefix+"-chrom.sizes",prefix+"-fwd.bw"])
    os.unlink(prefix+"-fwd.bedgraph")
        
    bedgraph(prefix+"-rev.bedgraph", zip(chrom_names, [ item.get() for item in reverse ]))    
    subprocess.check_call([
        "wigToBigWig",prefix+"-rev.bedgraph",prefix+"-chrom.sizes",prefix+"-rev.bw"])
    os.unlink(prefix+"-rev.bedgraph")
    
    os.unlink(prefix+"-chrom.sizes")


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



@config.help("Produce various bigwig depth of coverage files from a BAM file.", """\
This tool uses PySam, so can not currently be used with PyPy.

Bigwig files produced:

cover - Depth of coverage by actually sequenced bases.
span - Depth of coverage of region spanned by reads or fragments.
start - Fragment start locations.
end - Fragment end locations.

""")
@config.Main_section("bam_files")
class Bam_to_bigwig(config.Action_with_prefix):
    bam_files = [ ]
    
    def run(self):
        with nesoni.Stage() as stage:
            stage.process(make_bigwig,
                self.prefix + "-cover", self.bam_files, fragment_split_coverage, True)
            stage.process(make_bigwig,
                self.prefix + "-span", self.bam_files, fragment_coverage, True)
            stage.process(make_bigwig,
                self.prefix + "-start", self.bam_files, read1_starts, False)
            stage.process(make_bigwig,
                self.prefix + "-end", self.bam_files, read2_starts, False)






