
BAM-file alignment attributes
===

Alignments in BAM files (output_dir/samples/sample_name/alignments_filtered_sorted.bam) contain the following attributes:


AA:i:...  
- Will be present if the read is considered poly(A) (has at least four non-templated As)

AN:i:...  
- Gives the observed non-templated poly(A) tail length

AG:i:...
- Gives the number of templated "A"s. Note this is the length in the genome, but may be as low as 60% actual "A"s in the genome. Non-"A"s are still counted in this length.

AD:i:...  
- Gives the number of adaptor sequence bases observed after the end of 
  the poly(A) sequence. If the adaptor sequence is observed, then we have
  sequenced the entirety of the poly(A) tail of this fragment.
  If absent assume zero. Not supported for colorspace reads.
