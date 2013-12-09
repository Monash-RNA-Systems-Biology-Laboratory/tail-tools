
import collections

from nesoni import config, io, grace

import sys


State = collections.namedtuple('State','state score a_start a_end') 


@config.help(
'Clip low quality sequence and poly-A runs from the end of colorspace reads.',
"""\
The run of 0s can contain up to one fifth other colors, to allow for \
sequencing errors.

Reads should be in CSFASTQ format.
""")
@config.String_flag('sample', 'Sample name (for logging of statistics).')
@config.Int_flag('quality', 'Minimum quality.')
#@config.String_flag('adaptor', 'Adaptor sequence expected after poly-A tail (basespace only).')
@config.Int_flag('length', 'Minimum length.')
@config.Bool_flag('debug', 'Show detected poly-A region and adaptor location in each read.')
@config.Main_section('filenames', 'Input FASTQ files.')
class Clip_runs_colorspace(config.Action_with_prefix):
    sample = 'sample'
    quality = 20
    length = 25
    debug = False
    filenames = [ ]

    def run(self):
        min_quality = chr(33+self.quality)

        with io.open_possibly_compressed_writer(self.prefix+'.csfastq.gz') as out_file:        
            n = 0
            n_discarded = 0
            n_clipped = 0
            total_before = 0
            total_clipped = 0
            for filename in self.filenames:
                for name, seq, qual in io.read_sequences(filename, qualities='required'):
                    
                    score = 0
                    start = 0
                    for i in xrange(len(seq)-1):
                        if qual[i] >= min_quality:
                            if seq[i+1] == '0':
                                score += 1
                            else:
                                score = max(0, score-4)
                                if not score: start = i+2
                    
                    n += 1
                    total_before += len(seq)
                
                    if start > self.length+1:
                        if start < len(seq):
                            n_clipped += 1
                            total_clipped += len(seq)-start
                        
                        print >> out_file, '@'+name
                        print >> out_file, seq[:start]
                        print >> out_file, '+'
                        print >> out_file, qual[:start-1]
                    else:
                        n_discarded += 1
        
        self.log.datum(self.sample,'reads',n)
        if n:
            self.log.datum(self.sample,'mean length before poly-A clipping',float(total_before)/n)
        self.log.datum(self.sample,'reads discarded as too short after poly-A clipping',n_discarded)
        self.log.datum(self.sample,'reads poly-A clipped and kept',n_clipped)
        if n_clipped:
            self.log.datum(self.sample,'mean length clipped',float(total_clipped)/n_clipped)





@config.help(
'Clip low quality sequence and poly-A runs from the end of basespace reads.',
"""\
Looks for a run of As, followed by an adaptor sequence, followed by anything. \
The run of As and the adaptor sequence can contain up to one fifth errors, to allow \
for sequencing errors.

Reads should be in FASTQ format.
""")
@config.String_flag('sample', 'Sample name (for logging of statistics).')
@config.Int_flag('quality', 'Minimum quality, for base-space.')
@config.String_flag('adaptor', 'Adaptor sequence expected after poly-A tail (basespace only).')
@config.Int_flag('length', 'Minimum length.')
@config.Bool_flag('debug', 'Show detected poly-A region and adaptor location in each read.')
@config.Main_section('filenames', 'Input FASTQ files.')
class Clip_runs_basespace(config.Action_with_prefix):
    sample = 'sample'
    adaptor = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    quality = 10
    length = 20
    debug = False
    filenames = [ ]

    def run(self):
        """
        
        <sequence> <poly-A> <adaptor> <anything>
        
        """
        min_quality = chr(33+self.quality)
        
        with io.open_possibly_compressed_writer(self.prefix+'.fastq.gz') as out_file, \
             io.open_possibly_compressed_writer(self.prefix+'.clips.gz') as out_clips_file:
            print >> out_clips_file, '#Read\tread length\tpoly-A start\tpoly-A end\tpoly-A start, ignoring adaptor\rpoly-A end, ignoring adaptor'
             
            n = 0
            n_discarded = 0
            n_clipped = 0
            total_before = 0
            total_clipped = 0

            for filename in self.filenames:
                for name, seq, qual in io.read_sequences(filename, qualities='required'):
                    best_score = 0
                    best_a_start = len(seq)
                    best_a_end = len(seq)
                    best_aonly_score = 0
                    best_aonly_start = len(seq)
                    best_aonly_end = len(seq)
                    for a_start in xrange(len(seq)):
                        if a_start and seq[a_start-1] == 'A': continue
                        
                        a_end = a_start
                        aonly_score = 0
                        while True:
                            if aonly_score > best_aonly_score:
                                best_aonly_score = aonly_score
                                best_aonly_start = a_start
                                best_aonly_end = a_end
                                            
                            score = aonly_score
                            for i in xrange(a_end,min(a_end+len(self.adaptor),len(seq))):
                                if qual[i] >= min_quality:
                                    if seq[i] == self.adaptor[i-a_end]:
                                        score += 1
                                    else:
                                        score -= 4
                            if score > best_score:
                                best_score = score
                                best_a_start = a_start
                                best_a_end = a_end
                        
                            if a_end >= len(seq): break
                            if qual[a_end] >= min_quality:
                                if seq[a_end] == 'A':
                                    aonly_score += 1
                                else:
                                    aonly_score -= 4
                                    if aonly_score <= 0: break
                            a_end += 1
                        
                    a_start = best_a_start
                    a_end = best_a_end
                    aonly_start = best_aonly_start
                    aonly_end = best_aonly_end                    
                        
                    if self.debug: # and a_end == a_start and a_end < len(seq)-10:        
                        print name
                        print ''.join( 'X' if item<min_quality else ' ' for item in qual )
                        print seq
                        print ' '*a_start + 'A'*(a_end-a_start) + self.adaptor
                        print ' '*aonly_start + 'A'*(aonly_end-aonly_start)
                        print 
                        sys.stdout.flush()

                    n += 1
                    total_before += len(seq)

                    print >> out_clips_file, '%s\t%d\t%d\t%d\t%d\t%d' % (name, len(seq) , a_start, a_end, aonly_start, aonly_end)
                    
                    if a_start > self.length:
                        if a_start < len(seq):
                            n_clipped += 1
                            total_clipped += a_start
                    
                        print >> out_file, '@'+name
                        print >> out_file, seq[:a_start]
                        print >> out_file, '+'
                        print >> out_file, qual[:a_start]
                    else:
                        n_discarded += 1
                    
                    if n%10000 == 0: 
                        grace.status('Clip-runs ' + self.sample + ' ' + grace.pretty_number(n)) # + ' (' + grace.pretty_number(len(dstates)) + ' dstates)')
        
        grace.status('')
        
        self.log.datum(self.sample,'reads',n)
        if n:
            self.log.datum(self.sample,'mean length before poly-A/adaptor clipping',float(total_before)/n)
        self.log.datum(self.sample,'reads discarded as too short after poly-A/adaptor clipping',n_discarded)
        self.log.datum(self.sample,'reads poly-A/adaptor clipped and kept',n_clipped)
        if n_clipped:
            self.log.datum(self.sample,'mean length clipped',float(total_clipped)/n_clipped)




if __name__ == '__main__':
    config.shell_run(Clip_runs(), sys.argv[1:], sys.executable + ' ' + __file__)
