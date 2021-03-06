% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tail_weight.R
\name{weighted_log2_tails}
\alias{weighted_log2_tails}
\title{Weight tail lengths}
\usage{
weighted_log2_tails(
  tails,
  tail_counts,
  design,
  genes = NULL,
  min_reads = 10,
  biovar = TRUE
)
}
\description{
Tail lengths are log2 transformed, and given weights based on read counts.
}
\details{
Weights are calculated based on a variance model including a technical variation component inversely proportional to the number of reads plus a constant biological variance component.

Weighting is fine-tuned based on genes with at least min_reads in each sample.

Features are filtered so that only features that can be fitted using samples with min_reads reads are used. (Samples with less than min_reads reads will still be used in the calculation, so this is stricter than absolutely necessary.)

Returns an EList of weighted log2 tail lengths.
}
