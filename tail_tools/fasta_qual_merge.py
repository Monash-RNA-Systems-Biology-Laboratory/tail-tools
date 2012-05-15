
from nesoni import config

import sys, itertools

@config.help("""\
Merge FASTA and qual files into a single FASTQ file (streamed to stdout).
""")
@config.Positional('fasta_file', 'FASTA file.')
@config.Positional('qual_file', 'Qual file.')
class Fasta_qual_merge(config.Action_with_optional_output):
    fasta_file = None
    qual_file = None

    def run(self):
        fa = open(self.fasta_file, 'rb')
        fq = open(self.qual_file, 'rb')

        out_file = self.begin_output()
        
        while True:
            a1 = fa.readline()
            if not a1: break
            a1 = a1.strip()
            a2 = fa.readline().strip()
            q1 = fq.readline().strip()
            q2 = fq.readline().strip()
            
            assert a1.startswith('>')
            assert a1 == q1
            
            print >> out_file, '@' + a1[1:]
            print >> out_file, a2
            print >> out_file, '+'
            print >> out_file, ''.join( chr(33+max(0,int(item))) for item in q2.split() )

        self.end_output(out_file)
        
        fa.close()
        fq.close()
        
                   
if __name__ == '__main__':
    config.shell_run(Fasta_qual_merge(), sys.argv[1:], sys.executable + ' ' + __file__)
