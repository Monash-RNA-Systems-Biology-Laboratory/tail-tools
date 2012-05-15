
import util

import sys, os

out_dir = sys.argv[1]
filenames = sys.argv[2:]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for filename in filenames:
    for name, seq in util.read_fasta(open(filename,'rU')):
        path = os.path.join(out_dir, util.filesystem_friendly_name(name))
        f = open(path,'wb')
        f.write(seq)
        f.close()
