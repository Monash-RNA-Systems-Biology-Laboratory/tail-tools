"""

Functions to help play with output from tail-tools.

Added to as necessity and whim dictates, don't expect any kind of sense.

"""


import mmap, collections
from os.path import join
from nesoni import io, annotation, reference_directory, span_index

import numpy

COLORS = {
    'A':(0,1,0),
    'T':(1,0,0),
    'C':(0.33,0.33,1),
    'G':(1,1,0),
    }
    
def index(filename, type=None, modify = lambda item: item, name = lambda item: item.get_id()):
    result = { }
    for item in annotation.read_annotations(filename):
        if type is not None and item.type != type: continue
        item = modify(item)
        assert name(item) not in result
        result[name(item)] = item
    return result    


class Reference(object): pass

def load_ref(dirname):
    ref = reference_directory.Reference(dirname, must_exist=True)

    result = Reference()
    
    #result.seqs = dict(io.read_sequences(join(dirname,'reference.fa')))
    result.seqs = { }
    for name, length in ref.get_lengths():
        with open(join(dirname,ref.name,name+'.txt'),'rb') as f:
            #result.seqs[name] = f.read()
            result.seqs[name] = mmap.mmap(f.fileno(), 0, access=mmap.PROT_READ)
    
    result.genes = index(join(dirname,'reference.gff'), 'gene')
    result.utrs = index(join(dirname,'utr.gff'))
    
    return result


class Analysis(object): pass

def load_analysis(dirname):
    result = Analysis()
    result.peaks = index(join(dirname,'peaks','relation-child.gff'),
        modify=lambda item: item.three_prime())
    result.peak_index = span_index.index_annotations(result.peaks.itervalues())
    return result



class Piler(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.x = numpy.arange(start,end)
        self.pile = numpy.zeros((end-start), 'int64')
    
    def add(self, offset, amount=1):
        if offset < self.start or offset >= self.end: return
        self.pile[offset-self.start] += amount

def peak_pile(ref,ana, start, end):
    fwd = Piler(start,end)
    rev = Piler(start,end)
    for gene in ref.genes.itervalues():
        three_prime = gene.three_prime()
        for peak in ana.peak_index.get(three_prime.shifted(start,end), False):
            rel = peak.relative_to(three_prime)
            if rel.strand > 0:
                fwd.add(rel.start)
            else:
                rev.add(rel.start)
    return fwd, rev


class Kmer_pile(object): pass

def kmers(k):
    if not k: return [ '' ]
    return [
        item+base
        for item in kmers(k-1)
        for base in ['A','C','G','T']
        ]
          
def kmer_pile(ref,feat,k,start,end):
    result = Kmer_pile()
    result.k = k
    result.n = len(feat)
    result.start = start
    result.end = end
    result.x = numpy.arange(start,end)
    
    result.kmers = kmers(k)
    result.kmer_index = dict( (item,i) for i,item in enumerate(result.kmers) )
    result.piles = numpy.zeros((len(result.kmers),end-start), 'int32')
    for item in feat:
        seq = item.three_prime().shifted(start,end+k).get_seq(ref.seqs).upper()
        for i in xrange(end-start):
            kmer = seq[i:i+k]
            if kmer in result.kmer_index:
               result.piles[result.kmer_index[kmer]][i] += 1
    return result

#def show_piles(p):
#    import pylab
#    for name, pile in p.items():
#        pylab.plot(pile)

def stack_kmer_pile(p):
    import pylab
    from matplotlib.patches import Rectangle
    pylab.figure(figsize=(20.0,6.0))
    
    k = p.k
    
    items = [ ]
    for i in xrange(len(p.kmers)):
        for j in xrange(p.end-p.start):
            if p.piles[i][j] >= 0.01*p.n:
                items.append((p.piles[i][j], p.kmers[i], j))
    
    heights = [ 0 ] * (p.end-p.start)
    
    while items:
        i = min(xrange(len(items)),
            key=lambda i: max(heights[items[i][2]:items[i][2]+k])  -items[i][0])
        height, kmer, i = items.pop(i)
    
        offset = max(heights[i:i+k])
        heights[i:i+k] = [ offset+height ] * k
        for j in xrange(k):
            pylab.gca().add_patch(Rectangle((p.start-0.5+i+j,offset),1,height,
                linewidth=0.1,
                facecolor = COLORS[kmer[j]]))
            #if height > 0.02*p.n:
            #    pylab.annotate(kmer[j],(p.start-0.5+i+j+0.5,offset+height*0.5), size=8, ha='center', va='center')
        #if height > 0.05*p.n:
        #    pylab.annotate('%.0f%%' % (height*100.0/p.n),(p.start-0.5+i+k,offset), size=8,ha='right',va='bottom')
        pylab.gca().add_patch(Rectangle((p.start-0.5+i,offset),k,height, fill=False,linewidth=1.5))

    for i,c in enumerate('ACGT'):
        pylab.figtext(0.95,0.9-i*0.05, c, size=15.0, color=COLORS[c], ha='center')
    
    pylab.xlim(p.start-1,p.end+1)
    pylab.ylim(0,max(heights))   
    
    pylab.yticks([ 0.1*p.n ], [ '10%' ])
    #pylab.yticks([])
    
    pylab.gca().set_position([0.05,0.05, 0.85,0.93])
    #pylab.tight_layout() 
    
    
    #pylab.xlim(0,n)    
    #pylab.ylim(0,1)
    #for i in xrange(n):
    #    vec = numpy.array([ item[i] for item in p.itervalues() ], 'float64')
    #    norm = vec/sum(vec)
    #    y = 0.0
    #    for j in sorted(xrange(len(norm)), key=lambda j: -norm[j]):
    #        if norm[j] > 0.01:
    #            pylab.gca().add_patch(Rectangle((i,y),1,norm[j],fill=False))
    #            pylab.annotate(kmers[j],(i,y),size=10.0,rotation=90.0,va='bottom')
    #            y += norm[j]
    
