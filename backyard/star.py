"""
The genome can contain runs of "A"s and sometimes reads are incorrectly clipped at these. A read will be clipped at a run of 10 or more "A"s (essentially).

Playing around with the idea of not clipping reads before aligning with STAR.

However, long runs of "A"s can attract reads. STAR does not realize this has lower information content than other sequences.

A possible ideal solution would be to align both clipped and non-clipped reads, and then choose the better alignment between them.
"""


import os
from nesoni import config, sam, bio

def cigar_parts(cigar):
    result = [ ]
    n = 0
    for char in cigar:
        if '0' <= char <= '9':
            n = n*10+ord(char)-ord('0')
        else:
            result.append((n,char))
            n = 0
    return result


def a_count(seq):
    score = 0
    n = 0
    this_score = 0
    this_n = 0
    while this_n < len(seq):
        if seq[this_n] == "A":
            this_score += 1
        else:
            this_score -= 4
        this_n += 1
        if this_score >= score:
            score = this_score
            n = this_n
    return n
        

def a_adaptor_count(seq, adaptor):
    score = 0
    n = 0
    nadapt = 0
    
    this_score = 0
    this_n = 0
    while this_n < len(seq):
        if seq[this_n] == "A":
            this_score += 1
        else:
            this_this_score = this_score
            this_this_nadapt = 0
            while this_this_nadapt < len(adaptor) and \
                  this_n+this_this_nadapt < len(seq):
                if seq[this_n+this_this_nadapt] == adaptor[this_this_nadapt]:
                    this_this_score += 1
                else:
                    this_this_score -= 4
                this_this_nadapt += 1
                if this_this_score >= score:
                    score = this_this_score
                    n = this_n
                    nadapt = this_this_nadapt
            
            this_score -= 4
                    
        this_n += 1
        if this_score >= score:
            score = this_score
            n = this_n
            nadapt = 0
    
    return n, nadapt
  

@config.help("Annotate alignments with poly(A) tail length. Mark alignments with too little non-poly(A) as unmapped. Designed for use with STAR.")
@config.Positional("input", "Input BAM filename.")
@config.String_flag('adaptor', 'Adaptor sequence expected after poly(A) tail.')
@config.String_flag('name', 'Name, for logging.')
class Annotate_bam(config.Action_with_prefix):
    input = None
    
    adaptor = 'GATCGGAAGAGCACACGTCTGAACTCCAGTCAC'    
    min_genomic = 20
    
    name = ''
    
    def run(self):
        adaptor = self.adaptor.upper()
        name = self.name or os.path.basename(self.prefix)
        
        headers = sam.bam_headers(self.input)
        
        writer = sam.Bam_writer(self.prefix+"_temp.bam", headers)
        
        n_kept = 0
        n_unaligned = 0
        n_discarded = 0
        n_multi = 0

        for i, al in enumerate(sam.Bam_reader(self.input)):
            if al.flag & sam.FLAG_UNMAPPED:
                writer.write(al)
                n_unaligned += 1
                continue
            
            reverse = al.flag & sam.FLAG_REVERSE
            if reverse:
                read_bases = bio.reverse_complement(al.seq)
                cigar = cigar_parts(al.cigar)[::-1]
            else:
                read_bases = al.seq.upper()
                cigar = cigar_parts(al.cigar)
            
            n_unaligned = 0
            if cigar and cigar[-1][1] == "S": 
                n_unaligned = cigar[-1][0]
            
            n_aligned = len(read_bases) - n_unaligned
            seq_unaligned = read_bases[n_aligned:]
            seq_aligned = read_bases[:n_aligned]
            AN, AD = a_adaptor_count(seq_unaligned, adaptor)
            AG = a_count(seq_aligned[::-1])
            
            if AN: al.extra.append("AN:i:%d" % AN)
            if AD: al.extra.append("AD:i:%d" % AN)
            if AG: al.extra.append("AG:i:%d" % AN)
            if AN >= 4: al.extra.append("AA:i:1")
            
            if n_aligned - AG < self.min_genomic:
                al.flag = al.flag | sam.FLAG_UNMAPPED
                n_discarded += 1
            else:
                n_kept += 1

                NH = 1
                for item in al.extra:
                    if item.startswith("NH:i:"):
                        NH = int(item[5:])
                if NH > 1: n_multi += 1
            
            if i % 10000 == 0:
                print al.rname, al.pos
                print cigar
                print " "*(len(seq_aligned)-AG)+"="*AG
                print seq_aligned#[::-1]
                print seq_unaligned
                print "="*AN+"D"*AD
            
            writer.write(al)

        self.log.datum(name, "reads", n_unaligned+n_kept+n_discarded)
        self.log.datum(name, "did not align", n_unaligned)
        self.log.datum(name, "short non-A alignments discarded", n_discarded)
        self.log.datum(name, "alignments kept", n_kept)
        self.log.datum(name, "multimappers", n_multi)

        writer.close()        
        sam.sort_and_index_bam(
            self.prefix+"_temp.bam",
            self.prefix)
        os.unlink(self.prefix+"_temp.bam")

    

