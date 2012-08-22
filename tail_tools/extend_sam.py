
import nesoni
from nesoni import config

import util

import sys, os, array, itertools

FLAG_PAIRED = 1
FLAG_PROPER = 2
FLAG_UNMAPPED = 4
FLAG_MATE_UNMAPPED = 8
FLAG_REVERSE = 16
FLAG_MATE_REVERSE = 32
FLAG_FIRST = 64
FLAG_SECOND = 128
FLAG_NONPRIMARY = 256
FLAG_FAIL = 512 #Read failed platform/vendor quality checks
FLAG_DUP = 1024 #PCR or optical duplicate

class Alignment(object):
    def __init__(self, line):
        parts = line.rstrip('\n').split('\t')
        self.extra = parts[11:]
        (self.qname, 
         flag, 
         self.rname, 
         pos, 
         mapq, 
         self.cigar, 
         self.mrnm, 
         mpos, 
         isize, 
         self.seq, 
         self.qual) = parts[:11]
        
        self.flag = int(flag)
        self.pos = int(pos) # 1-based
        self.mapq = int(mapq)
        self.mpos = int(mpos)
        self.isize = int(isize)
                
        self.length = 0
        n = 0
        for value in array.array('B', self.cigar):
            if 48 <= value <= 57:
                n = n*10+(value-48)
            else:
                #if char in 'MDNP=X':
                if value == 77 or value == 68 or value == 78 or value == 80 or value == 61 or value == 88:
                    self.length += n
                n = 0

    def __repr__(self):
        #return '%s @ %d..%d %s %s' % (self.original_name(),self.pos,self.pos+self.length-1, '-' if self.flag&FLAG_REVERSE else '+',self.rname,)
        return '\t'.join([
            self.qname,
            str(self.flag),
            self.rname,
            str(self.pos),
            str(self.mapq),
            self.cigar,
            self.mrnm,
            str(self.mpos),
            str(self.isize),
            self.seq,
            self.qual,
        ] + self.extra) 
    
    def original_name(self):
        #Assuming it was Illumina
        if self.flag&FLAG_PAIRED:
            if self.flag&FLAG_FIRST:
                return self.qname+'/1'
            elif self.flag&FLAG_SECOND:
                return self.qname+'/2'
            else:
                raise grace.Error('Confused by SAM file')                
        else:
            return self.qname

    def get_qual(self):
        if self.qual == '*':
            return '^' * len(self.seq)
        else:
            return self.qual

    def set_flag(self, mask, value):
        if value:
            self.flag |= mask
        else:
            self.flag &= ~mask

    def get_mrnm(self):
        if self.mrnm == '=':
            return self.rname
        else:
            return self.mrnm
    
    def get_AS(self):
        """ Get alignment score.
            Prefer AS if present. Meaning of mapq is complicated. """
        for item in self.extra:
            if item.startswith('AS:i:'):
                return int(item[5:])
        
        return self.mapq        
        #raise grace.Error('SAM line lacks AS')


solid_encoding = {
   'AA':'0', 'CA':'1', 'GA':'2', 'TA':'3',
   'AC':'1', 'CC':'0', 'GC':'3', 'TC':'2',
   'AG':'2', 'CG':'3', 'GG':'0', 'TG':'1',
   'AT':'3', 'CT':'2', 'GT':'1', 'TT':'0',  
}

solid_decoding = { }
for dinuc in solid_encoding:
    solid_decoding[ dinuc[0] + solid_encoding[dinuc] ] = dinuc[1]

def solid_encode(seq):
    seq = seq.upper()
    enc = [ ]
    for i in xrange(len(seq)-1):
        enc.append(solid_encoding.get(seq[i:i+2], 'N'))
    
    return ''.join(enc)

def solid_decode(seed, seq):
    dec = [ seed ]
    for c in seq:
        dec.append( solid_decoding.get(dec[-1]+c, 'N') )
    return ''.join(dec[1:])


complement = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' }
def rev_comp(seq):
    return ''.join( complement.get(c,'N') for c in seq.upper()[::-1] )


def alignment_score(qual, seq1, seq2, qual_cutoff):
    """ Score an alignment starting from the start, 
        ending at an arbitrary point,
        no indels,
        low-quality scores 0,
        match scores 1,
        mismatch scores -4 
    """
    min_quality = chr(33+  qual_cutoff  )

    best_score = 0
    best_end = 0
    score = 0
    for i in xrange(len(qual)):
        if qual[i] < min_quality: continue
        if seq1[i] == seq2[i]:
            score += 1
        else:
            score -= 4
        if score > best_score:
            best_score = score
            best_end = i+1
    return best_score, best_end


def cigar_decode(cigar):
    dec = [ ]
    n = 0
    for char in cigar:
        if '0' <= char <= '9':
            n = n*10+ord(char)-ord('0')
        else:
            dec.append(char*n)
            n = 0
    return ''.join(dec)

def cigar_encode(cigin):
    enc = [ ]
    for key, subiterator in itertools.groupby(cigin):
        n = len(list(subiterator))
        enc.append('%d%s' % (n,key))
    return ''.join(enc)
      


@config.help(
'Having aligned clipped reads the reference, if the clipping was too enthusiastic and \
the reference truly does contain a poly-A run, extend the alignment.', 
"""\
Stream SAM records from stdin to stdout.

SAM records of reads that have a tail of at least three bases are tagged with an 'AA' attribute.
""")
@config.Int_flag('quality', 'Minimum quality.')
@config.Positional('reads_filename', 'Original reads in FASTQ format.')
@config.Main_section('reference_filenames', 'Reference sequences in FASTA format.')
class Extend_sam(config.Action_filter):
    quality = 20
    reads_filename = None
    reference_filenames = [ ]

    def cores_required(self):
        # This is memory intensive, so require it to be exclusive
        # (ideally there would be some sort of memory resource management...)
        return nesoni.coordinator().get_cores()

    def run(self):
        references = { }
        for filename in self.reference_filenames:
            print >> sys.stderr, 'Load', filename
            for name, seq in util.read_fasta(open(filename,'rb')):
                references[name] = seq
        
        print >> sys.stderr, 'Load', self.reads_filename
        reads = { }
        for name, seq, qual in util.read_fastq(open(self.reads_filename,'rb')):
            reads[name] = (seq, qual)
        
        print >> sys.stderr, 'Begin'
        
        in_file = self.begin_input()
        out_file = self.begin_output()
        
        for line in in_file:
            line = line.rstrip()
            if line.startswith('@'):
                print >> out_file, line
                continue
            
            al = Alignment(line)
            
            reverse = al.flag & FLAG_REVERSE
            if reverse:
                read_bases = rev_comp(al.seq)
                read_qual = al.qual[::-1]
                cigar = cigar_decode(al.cigar)[::-1]
            else:
                read_bases = al.seq
                read_qual = al.qual
                cigar = cigar_decode(al.cigar)
            
            al.extra = [ item for item in al.extra
                         if not item.startswith('CQ:Z:') and 
                            not item.startswith('CS:Z:') ] + [
                         'CQ:Z:'+reads[al.qname][1],
                         'CS:Z:'+reads[al.qname][0],
                       ]
            
            ref = references[al.rname]
            
            seq_tail = reads[al.qname][0][ len(al.seq)+1: ]
            qual_tail = reads[al.qname][1][ len(al.seq): ]
            n_tail = len(seq_tail)
            
            if reverse:
                if al.pos-1-n_tail < 0: continue #TODO: handle tail extending beyond end of reference
                bases_ref = rev_comp(ref[al.pos-1-n_tail:al.pos-1+1])    
            else:
                if al.pos-1+al.length+n_tail > len(ref): continue #TODO: handle tail extending beyond end of reference
                bases_ref = ref[al.pos-1+al.length-1:al.pos-1+al.length+n_tail]
                    
            seq_ref = solid_encode( bases_ref )
            
            basic_score = alignment_score(qual_tail, seq_tail, seq_ref, self.quality)
            
            if n_tail:
                tail_score, tail_pos = max(
                    (alignment_score(
                        qual_tail, seq_tail,
                        solid_encode(bases_ref[:1+i] + 'A'*(n_tail-i)), 
                        self.quality)[0],
                     i)
                    for i in xrange(n_tail+1)
                )
                
                baseline = max(0, alignment_score(
                    qual_tail[:tail_pos], seq_tail[:tail_pos], seq_ref[:tail_pos], self.quality)[0])
        
                if tail_score >= baseline + 3:
                    #Record position of end of transcript in 'AA' (1-based position)
                    if reverse:
                        tail_refpos = al.pos - tail_pos
                    else:
                        tail_refpos = al.pos+al.length + tail_pos - 1 
                    al.extra.append('AA:i:%d' % tail_refpos)
                
                #Record sequence's poly(A) tail length in AN: from the end of correspondence with the reference sequence to the end of good quality As
                estimated_tail_length = alignment_score(qual_tail,seq_tail,solid_encode(bases_ref[:1+tail_pos]+'A'*(n_tail-tail_pos)),self.quality)[1] - tail_pos
                if estimated_tail_length > 0:
                    al.extra.append('AN:i:%d' % estimated_tail_length)
                
                if tail_pos:    
                    read_bases += solid_decode(read_bases[-1], seq_tail[:tail_pos])
                    read_qual += qual_tail[:tail_pos]
                    cigar += 'M'*tail_pos
                    al.length += tail_pos
                    if reverse:
                        al.pos -= tail_pos
                        al.seq = rev_comp(read_bases)
                        al.qual = read_qual[::-1]
                        al.cigar = cigar_encode(cigar[::-1])
                    else: 
                        al.seq = read_bases
                        al.qual = read_qual
                        al.cigar = cigar_encode(cigar)
            
            print >> out_file, al
        
        self.end_output(out_file)
        self.end_input(in_file)
                   
if __name__ == '__main__':
    config.shell_run(Extend_sam(), sys.argv[1:], sys.executable + ' ' + __file__)


