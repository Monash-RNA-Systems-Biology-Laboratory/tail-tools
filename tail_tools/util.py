
def filesystem_friendly_name(name):
    """ Remove special characters from a name """

    for char in '\'"<>&|/\\_ .':
        name = name.replace(char,'_')
    return name

def read_fasta(f):
    line = f.readline()
    while line:
        line = line.rstrip()
        assert line.startswith('>'), 'Not a FASTA file?'
        name = line[1:].split()[0]
        assert name, 'FASTA file contains record with no name'
        
        line = f.readline()
        parts = [ ]
        while line and not line.startswith('>'):
            parts.append(line.rstrip())
            line = f.readline()
        
        yield name, ''.join(parts)

def read_fastq(reads_file):
    while True:
        line1 = reads_file.readline()
        if not line1: break
        line2 = reads_file.readline()
        line3 = reads_file.readline()
        line4 = reads_file.readline()
            
        assert line1.startswith('@'), 'Not a FASTQ file?'
        assert line3.startswith('+'), 'Not a FASTQ file?'
            
        read_name = line1.rstrip('\n')[1:]
        read_seq = line2.rstrip('\n')
        read_qual = line4.rstrip('\n')
        yield read_name, read_seq, read_qual
