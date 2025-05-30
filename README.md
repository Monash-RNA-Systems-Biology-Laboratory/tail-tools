
# Tail Tools

This is a Python 2 based suite of tools for analysing Illumina or SOLiD sequencing reads with poly(A) tails, as produced using the PAT-Seq technique. The PAT-Seq technique was developed by Dr. Traude Beilharz, who heads the RNA Systems Biology Laboratory at Monash University.

Tail Tools is developed by Dr. Paul Harrison (paul.harrison@monash.edu) at Monash University. Development was initially under the auspices of the Victorian Bioinformatics Consortium and now continues with the Monash Genomics and Bioinformatics Platform. Michael See contributed R code to visualize output as an interactive heatmap.

Please feel free to email Paul any questions you have about getting Tail Tools up and running.

Links

* [RNA Systems Biology Laboratory](http://rnasystems.erc.monash.edu)
* [Monash Genomics and Bioinformatics Platform](https://www.monash.edu/researchinfrastructure/mgbp)


## License

This software is distributed under the terms of the GPL, version 2 or later,
excepting that:

- The third party javascript libraries included for convenience
  in directory tail_tools/web/third_party are covered by the terms of
  their respective licenses (also in that directory).

- The remaining files in the directory tail_tools/web are placed in the
  public domain.


## Installation

Conda is the easiest way to set up an environment to run Tail Tools. Otherwise install the requirements listed below, and then the Python and R parts of Tail Tools.

### Conda

This has been tested in May 2025.

```
conda create --name tail-tools
conda activate tail-tools

conda install -c bioconda -c conda-forge \
    'pypy2.7>=5.10.0' 'star>=2.7.11b' 'r-base>=4.4.3' 'r-pak>=0.8.0.2' \
    'samtools>=1.21' 'igvtools>=2.17.3' 'ucsc-wigtobigwig>=472' \
    'imagemagick>=7.1.1' 'pigz>=2.8'

pypy -m ensurepip --upgrade

pypy -m pip install --upgrade \
    'git+https://github.com/Victorian-Bioinformatics-Consortium/nesoni.git#egg=nesoni' \
    'git+https://github.com/Monash-RNA-Systems-Biology-Laboratory/tail-tools.git#egg=tail-tools'

Rscript -e 'pak::pkg_install(c(
    "MonashBioinformaticsPlatform/varistran",
    "Victorian-Bioinformatics-Consortium/nesoni/nesoni/nesoni-r",
    "Monash-RNA-Systems-Biology-Laboratory/tail-tools/tail_tools"))'
```


### Requirements

Tail Tools uses Python version 2. Use of PyPy is recommened for speed.

- [nesoni](https://github.com/Victorian-Bioinformatics-Consortium/nesoni), most easy installed with pip in Python and BiocManager in R:

```
# Install into standard python2:
pip install --user 'git+https://github.com/Victorian-Bioinformatics-Consortium/nesoni.git#egg=nesoni'

# Install into pypy:
pypy -m pip install --user 'git+https://github.com/Victorian-Bioinformatics-Consortium/nesoni.git#egg=nesoni'

R
install.packages("BiocManager")
BiocManager::install("Victorian-Bioinformatics-Consortium/nesoni", subdir="nesoni/nesoni-r")
```

  You don't need to install all of nesoni's dependencies, just Python 2.7 or later or PyPy. Do be sure to install the R component of nesoni.

- [STAR aligner](https://github.com/alexdobin/STAR)

- [samtools](http://www.htslib.org/)

- The "convert" tool from [ImageMagick](https://imagemagick.org/). Ubuntu users may need to further install `libmagickcore*-extra`, and also `inkscape` to ensure the RSVG renderer is available.

- The "wigToBigWig" tool from the [UCSC Genome Browser utilities](http://hgdownload.soe.ucsc.edu/admin/exe/).

- "pigz" parallel gzip utility.

- "igvtools" command line utility from IGV.

- [R](https://www.r-project.org/) 3.6 or higher.

- My R package [varistran](https://github.com/MonashBioinformaticsPlatform/varistran).

```
BiocManager::install("MonashBioinformaticsPlatform/varistran")
```

- My R package [topconfects](https://github.com/pfh/topconfects). Topconfects is in Bioconductor, but Tail Tools will often depend on "devel" features, so may need the github version.

```
BiocManager::install("pfh/topconfects", dependencies=TRUE)
```

### Optional requirement to support deprecated features:

- Optionally can also use bowtie2 instead of STAR for Illumina reads, or SHRiMP for SOLiD reads.

- [SplitsTree](http://www.splitstree.org/)
  Note: v4.13.1 seems to be broken, v4.11.3 works

- rsync (for downloads from UCSC browser)

- [Fitnoise](https://github.com/pfh/fitnoise) for old-style differential testing.

- [degust.py](https://victorian-bioinformatics-consortium.github.io/degust/dist/latest/degust.py) for old-style differential testing.

- topconfects branches to support alternative versions of differential tests (no longer recommended):

```
BiocManager::install("pfh/topconfects@wald", dependencies=TRUE)
BiocManager::install("pfh/topconfects@ql", dependencies=TRUE)
```

### Python part installation

Easy way:

```
pip install --user --upgrade 'git+https://github.com/Monash-RNA-Systems-Biology-Laboratory/tail-tools.git#egg=tail-tools'
```

From source:

```
python setup.py install
```

For PyPy, you can use pip by invoking the module:

```
pypy -m pip install --user --upgrade 'git+https://github.com/Monash-RNA-Systems-Biology-Laboratory/tail-tools.git#egg=tail-tools'
```

### R part installation

Easy way:

```R
BiocManager::install("Monash-RNA-Systems-Biology-Laboratory/tail-tools", subdir="tail_tools", dependencies=TRUE)
```

From source:

```
R
devtools::install("tail_tools")
```


## Usage

This package contains a number of tools, which can be listed by typing:

```
tail-tools
```

The package can be used directly from the source directory with:

```
python -m tail_tools
```

These tools may also be used from a python script (using the same system as my older genomics python package "nesoni"). A typical example of invoking the pipeline from python can be found below.


### R library usage

The tailtools R library can be loaded in R with:

```R
library(tailtools)
```

The two main uses of this are two Shiny apps (interactive web-based applications):

* The pipeline creates a Shiny app as part of its output (in subdirectory "shiny"). This can be served with ShinyServer or viewed from within R with `shiny::runApp("pipelineoutputdir/shiny")`.

* You can create a Shiny app for differential tests. A template is given below.


## Reference format

Before processing any reads, you need to create a "tail-tools reference directory".

References are most easily downloaded from Ensembl. It is also possible to use references from the UCSC browser, but not recommended.

References can be downloaded from Ensembl by downloading the "primary_assembly" version of the genome and a gene annotationn gff3 file from ftp://ftp.ensembl.org/pub/ and running:

```
tail-tools make-ensembl-reference: \
    <output_dir> \
    <assembly_file.fa.gz> \
    <gff3_file.gff3.gz>
```

If creating your own reference, it needs to consist of:

- sequences, eg in FASTA format
- annotations in GFF3 format

The reference directory is then created with the command:

```
tail-tools make-tt-reference: \
    <output_dir> \
    <sequence_file> \
    <annotations_file>
```

Annotations shall include the following feature types and attributes:

gene
* ID - unique identifier
* Name (optional) - nomenclature name
* Product (optional) - short description
* Biotype (optional) - what type of gene it is (protein_coding, rRNA, etc)

mRNA
* ID      - unique identifier
* Parent  - gene ID

CDS
* Parent  - mRNA ID

exon
* Parent  - mRNA ID


The pipeline assumes that two different genes do not have overlapping exons on the same strand. Is this too much to ask for in a reasonable genome annotation? Apparently the answer is yes. The UCSC and ENSEMBL genome downloaders merge genes with overlapping exons -- ids and names are concatenated with "/" as a separator. This can complicate downstream analysis, and Paul apologises for the pain this causes. In the case of ENSEMBL, higher confidence transcripts are given priority in order to avoid merging genes, and gene merging is rare.



## Pipeline

Having created a reference directory, the next step is to run the pipeline,
"analyse-polya-batch". This can be done from the command line, but is more
usefully done from a python script. We suggest adapting the following example
to your data:

```python

import tail_tools, nesoni, glob

tags = [
    ('wt1',   ['wt','rep1']),
    ('wt2',   ['wt','rep2']),
    ('wt3',   ['wt','rep3']),
    ('mutA1', ['mutA','rep1']),
    ('mutA2', ['mutA','rep2']),
    ('mutA3', ['mutA','rep3']),
    ('mutB1', ['mutB','rep1']),
    ('mutB2', ['mutB','rep2']),
    ('mutB3', ['mutB','rep3']),
]

# Where to find FASTQ files. %s becomes the name of the sample. Wildcards * and ? may be used.
filename_pattern = 'raw_data/%s.fastq.gz'

# For each sample we create a tail_tools.Analyse_polya instance
# Each sample is given a set of tags
samples = [ ]
for name, tags in tags:
    reads = sorted(glob.glob(filename_pattern % name))    
    assert reads, 'No reads for '+name
    
    samples.append(tail_tools.Analyse_polya(
        name,
        reads = reads,
        tags = tags,
        
        #To adjust clipping prior to alignment, modify these defaults:
        clip_runs_basespace = tail_tools.Clip_runs_basespace(
            # How many As to make up for one mismatch in the poly(A) sequence?
            # - Default is 4, which has been good with older four-color Illumina sequencing.
            #   (eg MiSeq)
            # - For two-color sequencing, completely disable mismatches
            #   with a value longer than the read length, eg 1000.
            #   (eg MiniSeq, NextSeq 550, NovaSeq 6000)
            a_mismatch_penalty = ...,
            
            # Clip to high quality region 
            # - Default is 0, no clipping, has been good for four-color sequencing.
            # - For two-color sequencing, tested with a NovaSeq dataset:
            #     clip_quality=20 gave good results, but we noticed it over-estimated 
            #                     tail lengths on 60-A Sequins.
            #     clip_quality=30 gave improved estimation of Sequin tail lengths, 
            #                     and is my current recommendation.
            # (note: high quality Gs are ignored)
            clip_quality = ...,
            
            # Adaptor used to clip reads. The default is correct for non-UMI PAT-Seq.
            # If using PAT-Seq with UMIs, ensure reads are named like: READNAME_BARCODE_UMI
            # Next check the orientation of barcodes and umis. Are they forward or reverse-complement?
            # Use 
            #   tail-tools detect-adaptor: my-reads.fastq.gz 
            # to determine correct choice.
            #adaptor = "umi_barcode" or "rc_umi_rc_barcode",
        ),

        #To discard multi-mapping reads. Default is to choose one location at random.
        #discard_multimappers=True,
        
        #To use bowtie2 rather than STAR
        #aligner="bowtie2",
        
        # Only use alignments with this many matching bases (STAR aligner only) 
        #min_match=30,
        
        #To allow for looser mispriming, lower this.
        #To only allow for stricter mispriming, raise this (maximum 1).
        #extension_prop_a = 0.6,
        ))


action = tail_tools.Analyse_polya_batch(
        # Output directory
        'pipeline',
                
        # List of instances of tail_tools.Analyse_polya
        samples = samples,
        
        # Title for report
        title = 'Pipeline output',
        
        # Reference directory you created earlier
        reference = '/data/reference/tail-tools/...',

        # Where will the shiny part of the html report be served from?
        # (In the output this part is pipeline/report/shiny)
        #shiny_report_url = "http://myserver:3838/...",

        # Allow reads/peaks this far downstrand of 
        # the annotated transcript end point
        # (however extension will not continue into CDS on the same strand)
        # For yeast use 400, for sparser genomes than yeast use 2000
        # (Left blank since it's easy to forget to change.)
        extension = ... ,
        
        # To be recognized as a poly(A) read there must be this many non-templated "A"s.
        # Default choice for older PAT-Seq is 4.
        # For more recent variants of PAT-Seq 12 As are always present when there is any 
        #   poly(A) tail, and you should specify 13.
        min_tail = ... ,
        
        # Minimum average tail length required to call a peak.
        # Set higher then 0.0 if there is mispriming.
        # 15.0 may be reasonable.
        peak_min_tail = 15.0,
        
        # Size of peak features generated
        # ie how far back from a site a read can end and still be counted towards it
        # Should be read length or a little shorter
        peak_length = 300,
        
        # Should UMI rather than read counting be used?
        # UMI should be part of the read name.
        # Ensure reads are named like: READNAME_BARCODE_UMI or READNAME_UMI
        umis = False,
        
        # Optional: Species to use in GO term analysis, choices are: Sc Ce Mm Hs
        #species="Sc",
        
        # The pipeline can crash on some steps if few or weird peaks are called.
        # Setting this to False will disable some of the fragile stages.
        do_fragile = True,
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

## Differential testing

Differential testing is most conveniently performed using a Shiny app in R.

Create a directory for a Shiny app containing an `app.R` file with something like:

```R
library(tailtools)
library(tidyverse)

pipeline_dir <- "...full path to.../pipeline"

# Construct a data frame of samples
# The tailtools package provides convenience functions for this, based on existing tags:
samples <- 
    pipeline_samples(pipeline_dir) %>%
    samples_group_tags("strain", c("wt","mutA","mutB")) %>%
    samples_group_tags("rep", c("rep1","rep2","rep3"))


# Construct a named list which defines the desired tests
tests <- list()

# Define desired tests in terms of samples to use, design matrix, and a coefficient.
tests[["wt_mutA"]] <- list(
    "test",
    title="wt to mutA",
    pipeline_dir=pipeline_dir,
    samples=samples$name,
    design=model.matrix(~ strain + rep, data=samples),
    coef="strainmutA")

# ... or in terms of a contrast.
tests[["mutA_to_mutB"]] <- list(
    "test",
    title="mutA to mutB",
    pipeline_dir=pipeline_dir,
    samples=samples$name,
    design=model.matrix(~ strain + rep, data=samples),
    contrast=c(0,-1,1,0,0))

# Multiple coefficients, or multiple contrasts (as columns in a matrix),
# may be specified to perfom an omnibus test.
tests[["strain_any"]] <- list(
    "test",
    title="Any differences between strains",
    pipeline_dir=pipeline_dir,
    samples=samples$name,
    design=model.matrix(~ strain + rep, data=samples),
    coef=c("strainmutA", "strainmutB"))


# A subset of samples may be used. 
# For example we could drop groups irrelevant to the test.
# Statistically this is safer but less powerful.
# Do this if samples are different enough that 
# noise levels won't be uniform across them.
keep <- 
    filter(samples, strain %in% c("wt","mutA")) %>% 
    droplevels()

tests[["wt_to_mutA"]] <- list(
    "test",
    title="wt to mutA",
    pipeline_dir=pipeline_dir,
    samples=keep$name,
    design=model.matrix(~ strain + rep, data=keep),
    coef="strainmutA")


# If neither coef or contrast is specified, test will look for 
# excess variability not explained by the design (which can be simply ~1).
#
# You can optionally specify a "calibration design" to be used in the
# weitrix calibration step. This sets a lower noise level that genes will need
# to exceed to be considered to have excess variability.
tests[["any_variability"]] <- list(
    "test",
    title="Any variability",
    pipeline_dir=pipeline_dir,
    samples=samples$name,
    design=model.matrix(~1, data=samples),
    calibration_design=model.matrix(~ strain + rep, data=samples))


app <- shiny_tests(tests, title="My shiny test app")
app
```

The app can be run from within R with `shiny::runApp("app-dir")` or using Shiny Server.


## BAM-file alignment attributes

* [BAM-file alignment attributes](doc/bam-files.md)


## Pipeline output

* [Directories and files produced by the pipeline](doc/output.md)


## Statistics

* [Statistics produced by `tail-tools analyse-tail-counts:`.](doc/statistics.md)
* [Differential expression and tail length, old method.](doc/differential.md)
* To understand new-style differential tests, see the vignettes in [topconfects](https://bioconductor.org/packages/devel/bioc/html/topconfects.html).





