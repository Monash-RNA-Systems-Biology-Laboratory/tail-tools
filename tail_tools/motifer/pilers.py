
from __future__ import division

from tail_tools import env
import numpy, re
from numpy import random



class Piler(object):
    """
    
    n - n piles
    ticks - [ (x, label) ]
    
    """
    
    def __init__(self, n):
        self.n = n
        self.x = numpy.arange(n) + 0.5
        self.ticks = [ ]

    
    def pile(self, seqs, recognizer, weights=None):
        if weights is None:
            weights = [ 1.0/len(self.fetchers) ]*len(self.fetchers)
    
        pile = numpy.zeros(self.n)
        
        for weight, (fetcher, bins) in zip(weights,self.fetchers):
            if not weight: continue
            
            seq = fetcher.shifted(0,recognizer.length).get_seq(seqs).upper()
            
            for i in xrange(len(bins)-1):
                a = bins[i]
                b = max(a+1,bins[i+1])
                for j in xrange(a,b):
                    pile[i] += float(recognizer(seq[j:j+recognizer.length])) * (weight/(b-a))
        
        return pile


    def pile_features(self, index, weights=None):
        if weights is None:
            weights = [ 1.0/len(self.fetchers) ]*len(self.fetchers)
    
        pile = numpy.zeros(self.n)
        
        for weight, (fetcher, bins) in zip(weights, self.fetchers):
            hits = index.get(fetcher, True)
            hits = [ item.relative_to(fetcher) for item in hits ]
            
            for i in xrange(len(bins)-1):
                a = bins[i]
                b = max(a+1,bins[i+1])
                for item in hits:
                    n = max(a,min(b,item.end)) - max(a,min(b,item.start))
                    pile[i] += n * weight / (b-a)
        
        return pile


    def setup_figure(self):
        import pylab
        pylab.xlim(0,self.n-1)
        pylab.xticks([item[0] for item in self.ticks],[item[1] for item in self.ticks])
        for item in self.ticks:
            if item[2]:
                pylab.axvline(item[0], color="black")


class Anchored_piler(Piler):
    def __init__(self, before, after, anchor_name, locations, stride=1):
        for item in locations:
            assert item.start == item.end
        
        Piler.__init__(self, ((before+after) // stride) ) 
        
        self.offset = -before
        self.ticks = [(before // stride, anchor_name, True)]
        
        tick_stride = 1
        i = 0
        while self.n > tick_stride * 20:
            if i%3 == 1:
               tick_stride = tick_stride * 5 // 2
            else:
               tick_stride *= 2
            i += 1
        for i in xrange(tick_stride,before//stride+1,tick_stride):
            self.ticks.append((before//stride-i,str(-i*stride), False))
        for i in xrange(tick_stride,after//stride+1,tick_stride):
            self.ticks.append((before//stride+i,str(i*stride), False))
        self.ticks.sort()
        
        fetch_bins = range(0,before+after+1,stride)
        assert len(fetch_bins) == self.n+1
        self.fetchers = [ 
            (item.shifted(-before,after), fetch_bins)
            for item in locations
            ]


#class Stretched_piler(Piler):
#    def __init__(self, n, locations):
#        Piler.__init__(self, n)
#        
#        self.ticks = [ ]
#        for i in xrange(11):
#            self.ticks.append((self.n*i/10.0, "%.0f%%" % (i*10.0), False))
#        
#        self.fetchers = [ ]
#        for item in locations:
#            length = item.end-item.start
#            bins = [ int(i*length/n) for i in xrange(n+1) ]
#            self.fetchers.append((item, bins))
        


class Stretched_piler(Piler):
    def __init__(self, first_label, sections):
        Piler.__init__(self, sum([ item[0] for item in sections ]))
        
        self.fetchers = [ ]
        for a,b in zip(sections[0][1],sections[-1][1]):
            self.fetchers.append((a.span_with(b),[0]))
        
        self.ticks = [ ]
        self.ticks.append((0,first_label,True))
        i = 0
        for n,features,show_percent,label in sections:
            if show_percent:
                for j in xrange(1,5):
                    self.ticks.append((i+n*j/5.0, "%.0f%%" % (j*20.0), False))
            self.ticks.append((i+n,label,True))

            for j, feature in enumerate(features):
                assert feature.strand == self.fetchers[j][0].strand
                feature = feature.relative_to(self.fetchers[j][0])
                assert feature.strand == 1
                assert feature.start >= 0
                assert feature.end <= self.fetchers[j][0].end-self.fetchers[j][0].start
                start = feature.start
                length = feature.end-feature.start
                for k in xrange(1,n+1):
                    self.fetchers[j][1].append(int(start+k*length/n))

            i += n
        

