
from nesoni import config, io

import sys

import util

@config.help(
'Clip low quality sequence and runs of 0s from the end of SOLiD reads.',
"""\
The run of 0s can contain up to one fifth other colors, to allow for \
sequencing errors.

Reads should be in FASTQ format.
""")
@config.String_flag('sample', 'Sample name (for logging of statistics).')
@config.Int_flag('quality', 'Minimum quality.')
@config.Int_flag('length', 'Minimum length.')
@config.Main_section('filenames', 'Input FASTQ files.')
class Clip_runs(config.Action_with_prefix):
    sample = 'sample'
    quality = 20
    length = 25
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
                for name, seq, qual in util.read_fastq(open(filename,'rb')):
                    
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


if __name__ == '__main__':
    config.shell_run(Clip_runs(), sys.argv[1:], sys.executable + ' ' + __file__)