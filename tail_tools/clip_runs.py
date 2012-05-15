
from nesoni import config

import sys

import util

@config.help(
'Clip low quality sequence and runs of 0s from the end of SOLiD reads.',
"""\
The run of 0s can contain up to one fifth other colors, to allow for \
sequencing errors.

Reads should be in FASTQ format.
""")
@config.Int_flag('quality', 'Minimum quality.')
@config.Int_flag('length', 'Minimum length.')
@config.Positional('filename', 'Input FASTQ file.')
class Clip_runs(config.Action_with_optional_output):
    quality = 20
    length = 25
    filename = None

    def run(self):
        out_file = self.begin_output()
        
        min_quality = chr(33+self.quality)
        
        for name, seq, qual in util.read_fastq(open(self.filename,'rb')):
            score = 0
            start = 0
            for i in xrange(len(seq)-1):
                if qual[i] >= min_quality:
                    if seq[i+1] == '0':
                        score += 1
                    else:
                        score = max(0, score-4)
                        if not score: start = i+2
        
            if start > self.length+1:
                print >> out_file, '@'+name
                print >> out_file, seq[:start]
                print >> out_file, '+'
                print >> out_file, qual[:start-1]

        self.end_output(out_file)


if __name__ == '__main__':
    config.shell_run(Clip_runs(), sys.argv[1:], sys.executable + ' ' + __file__)