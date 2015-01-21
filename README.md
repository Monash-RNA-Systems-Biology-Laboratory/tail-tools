
Tail Tools
==========

This is a Python 2 based suite of tools for analysing SOLiD or Illumina
sequencing reads with poly(A) tails.

Releases can be found in the Python Package Index:

https://pypi.python.org/pypi/tail-tools/


License:
========

This software is distributed under the terms of the GPL, version 2 or later,
excepting that:

- The third party javascript libraries included for convenience
  in directory tail_tools/web/third_party are covered by the terms of
  their respective licenses (also in that directory).

- The remaining files in the directory tail_tools/web are placed in the
  public domain.


Requirements:
=============

Use of PyPy is recommened for speed.

- "nesoni", available from https://github.com/Victorian-Bioinformatics-Consortium/nesoni or using

    pip install nesoni

  You don't need to install all of nesoni's dependencies, just Python 2.7 
  or later or PyPy.

- bowtie2
  (SHRiMP for legacy color-space data)

- samtools

- SplitsTree from http://www.splitstree.org/, note: v4.13.1 seems to be broken, v4.11.3 works

- The "convert" tool from ImageMagick.

- rsync (for downloads from UCSC browser)

- R+, with package "seriation" and BioConductor packages "limma" and "edgeR"

- degust.py from 
  https://victorian-bioinformatics-consortium.github.io/degust/dist/latest/degust.py


Installation:
=============

Easy way:

    pip install tail-tools

From source:

    python setup.py install

For PyPy it seems to be currently easiest to set up in a virtualenv:

    virtualenv -p pypy myenv
    myenv/bin/pip install tail-tools


Usage:
======

This package contains a number of tools, which can be listed by typing:

    tail-tools
  

The package can be used directly from the source directory with:

    python -m tail_tools


These tools may also be used as part of a nesoni-style workflow python script.

Typical usage of the pipeline is described below.


Reference format:
=================

Before processing any reads, you need to create a "tail-tools reference directory".

References are most easily downloaed from the UCSC browser using:

    tail-tools make-ucsc-reference: \
        <output_dir> \
        <ucsc_reference_name>

If creating your own reference, it needs to consist of:

- sequences, eg in FASTA format
- annotations in GFF3 format

The reference directory is then created with the command:

    tail-tools make-tt-reference: \
        <output_dir> \
        <sequence_file> \
        <annotations_file>

Annotations shall include the following feature types and attributes:

gene
* ID - unique identifier
* Name (optional) - nomenclature name
* Product (optional) - short description

mRNA
* ID      - unique identifier
* Parent  - gene ID

CDS
* Parent  - mRNA ID

exon
* Parent  - mRNA ID



Pipeline:
=========

Having created a reference directory, the next step is to run the pipeline,
"analyse-polya-batch". This can be done from the command line, but is more
usefully done from a python script. We suggest adapting the following example
to your data:

```python

import tail_tools, nesoni, glob

tags = [
    ('logRep1',         ['BY',   'rep1']),
    ('logRep2',         ['BY',   'rep2']),
    ('deltaccr4logRep1',['ccr4', 'rep1']),
    ('deltaccr4logRep2',['ccr4', 'rep2']),
    ('deltaccr4logRep3',['ccr4', 'rep3']),
    ('YPEGRep1',        ['ypeg', 'rep1']),
    ('YPEGRep2',        ['ypeg', 'rep2']),
    ('GALRep1',         ['gal',  'rep1']),
    ('GALRep2',         ['gal',  'rep2']),
    ('GLU10Rep1',       ['glu10','rep1']),
    ('GLU10Rep2',       ['glu10','rep2']),
    ('GLU20Rep1',       ['glu20','rep1']),
    ('GLU20Rep2',       ['glu20','rep2']),
]

filename_pattern = 'mydata/Sample_scBY4741%s/*.fastq.gz'

# For each sample we create a tail_tools.Analyse_polya instance
# Each sample is given a set of tags
samples = [ ]
for name, tags in tags:
    reads = sorted(glob.glob(filename_pattern % name))
    samples.append(tail_tools.Analyse_polya(
        name,
        reads = reads,
        tags = tags,            
        ))

action = tail_tools.Analyse_polya_batch(
        # Output directory
        'yeast-june-2013',
        
        # Title for report
        title = 'Yeast June 2013',
        
        # Files in report will have this prefix
        file_prefix = 'yeast-june-2013',
        
        # Reference directory you created earlier
        reference = 'sacCer3',
        
        # Allow reads/peaks this far downstrand of 
        # the annotated transcript end point
        # For sparser genomes than yeast, perhaps use 1000
        extension = 200,
        
        # Whether to include .genome file for IGV in plots tarball
        # Not necessary for model organisms where IGV 
        # already provides the genome.
        include_genome = False,
                
        # List of instances of tail_tools.Analyse_polya
        samples = samples,
        
        # List of sample groups
        # A sample group is specified as 
        # '<nesoni-selection-expression>=<name>'
        # (=<name> may be omitted)
        # See nesoni help for description of selection expressions,
        # this uses the tags given to each sample to concisely 
        # specify sets of samples.
        groups = [ 'BY', 'ccr4', 'ypeg', 'gal', 'glu10' ],
        
        # (Advanced)
        # Perform differential tests
        tests = [
            tail_tools.Test(
                'BY-ccr4',
                title='BY vs ccr4',
                null=['BY/ccr4'],
                alt=['ccr4'],
                ),
            #etc
            ],
        )



# A little boilerplate so that
# - multiprocessing works
# - you can control making
#   (see nesoni help on --make-* flags)

def main():
    action.make()

if __name__ == '__main__': 
    nesoni.run_script(main)

# If run again with adjusted parameters,
# only the parts that need to be run again will run.
#
# To force a complete re-run:
#     python myscript.py --make-do all
#
# To re-run everything but the alignment to reference
# (eg if there is a new version of tail-tools)
#     python myscript.py --make-do all --make-done analyse-polya
#

```

BAM file annotations
====================

AA:i:...  
- Will be present if the read is considered poly(A) (has at lest four non-templated As)

AN:i:...  
- Gives the observed non-templated poly(A) tail length

AD:i:...  
- Gives the number of adaptor sequence bases observed after the end of 
  the poly(A) sequence. If the adaptor sequence is observed, then we have
  sequenced the entirety of the poly(A) tail of this fragment.
  If absent assume zero. Not supported for colorspace reads.









