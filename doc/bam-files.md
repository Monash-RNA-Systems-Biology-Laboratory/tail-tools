
BAM-file alignment attributes
===

Alignments in BAM files (output_dir/samples/sample_name/alignments_filtered_sorted.bam) contain the following attributes:


AA:i:...  
- Will be present if the read is considered poly(A) (has at least four non-templated As)

AN:i:...  
- Gives the observed non-templated poly(A) tail length

AD:i:...  
- Gives the number of adaptor sequence bases observed after the end of 
  the poly(A) sequence. If the adaptor sequence is observed, then we have
  sequenced the entirety of the poly(A) tail of this fragment.
  If absent assume zero. Not supported for colorspace reads.
