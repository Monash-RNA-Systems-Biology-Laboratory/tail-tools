
Directories and files produced by the tail-tools pipeline
===

This output has grown in an organic fashion, for which I apologize.



### `report/`

An HTML report. 

Note that this contains symlinks to other files in the directory, so when copying it make sure to dereference symlinks (eg `cp -L`).



### `plotter_config.json`

JSON file containing various useful pieces of information, intended for visualization software that uses the pipeline output directory.



### `peaks/relation-child.gff`

Peaks called from poly(A) reads. Really only the 3' end of these features is meaningful, the 3' end is the polyadenylation site.



### `samples/<samplename>/`

Per-sample files.

Of particular importance:

* `alignments_filtered_sorted.bam` is the BAM file to use for viewing. These are also included as a .zip file in the html report.

The `samples/<samplename>-polyA` directories contain BAM files limited to reads with a poly(A) tail.



### `raw/`

Read count and tail length statistics in CSV format. These are also given as a .zip file in the html report.

These are produced for several types of feature:

* `genewise` - reads anywhere in annotated genes (or some number of bases downstrand, as controlled by the `extensions` parameter).

* `peakwise` - peaks called based on poly(A) reads.

* `pairwise` - the two most prominent peaks in the 3' UTR of each gene.

The files contain:

* `-count.csv` - read counts for each feature.

* `-tail-count.csv` - number of reads with a poly(A) tail for each feature.

* `-tail.csv` - average poly(A) tail length.

