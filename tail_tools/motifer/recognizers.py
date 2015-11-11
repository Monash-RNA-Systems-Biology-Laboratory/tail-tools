
from __future__ import division


def kmers(letters, n):
    if not n: return [""]
    shorter = kmers(letters, n-1)
    result = [ ]
    for letter in letters:
        for item in shorter:
            result.append(letter+item)
    return result 



class Recognizer(object):
    def count(self, seq):
        n = 0
        for i in xrange(0,len(seq)-self.length+1):
            if self(seq[i:i+self.length]):
                n += 1
        return n


class Recognize_string(Recognizer):
    def __init__(self, string):
        self.string = string
        self.length = len(string)
    
    def __call__(self, string):
        return string == self.string


class Recognize_regex(Recognizer):
    def __init__(self, length, pattern):
        import re
        self.length = length
        self.regex = re.compile(pattern)
    
    def __call__(self, string):
        return self.regex.match(string) is not None



base_number = { "A":0,"C":1,"G":2,"T":3 }

class Bad_pwm_exception(Exception): pass

class Recognize_pwm(Recognizer):
    def __init__(self, filename, within=2.0):
        with open(filename,"rU") as f: 
            lines = f.readlines()
        if not lines:
           raise Bad_pwm_exception()
        assert lines[0].startswith("Pos")
        self.length = len(lines)-1
        matrix = numpy.zeros((self.length,4))
        for i in xrange(self.length):
            parts = lines[i+1].strip().split()
            for j in xrange(4):
                matrix[i,j] = float(parts[j+1])
        self.score_matrix = numpy.log(matrix)/numpy.log(2.0) + 2.0 #TODO: background base frequency
        
        max_score = numpy.sum(numpy.maximum.reduce(self.score_matrix,axis=1))
        self.cutoff = max_score - within
        
    def __call__(self, string):
        score = 0.0
        for i in xrange(self.length):
            score += self.score_matrix[i, base_number[string[i]]]
        
        #return 2.0 ** score
        return score >= self.cutoff



def kmer_recognizers(upto=3):
    recognizers = [ ]
    for n in xrange(1,upto+1):
        recognizers.append( "%d-mers" % n )
        for kmer in kmers("ACGT", n):
            recognizers.append((
                kmer, Recognize_string(kmer)
                ))

    return recognizers





