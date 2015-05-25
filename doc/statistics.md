
Statistics produced by analyse-tail-counts
===

`tail-tools analyse-tail-counts:` produces a variety of CSV files as output.

`read_count.csv` contains read counts.

`tail_count.csv` contains counts of reads with tails. This is gives an indication of how accurate mean or median tail lengths are. Means and medians from low tail counts will have a high standard error.

`tail.csv` contains average tail lengths.

`tail_sd.csv` contains the standard deviation of tail lengths.

`tail_quantile_50.csv` contains the median tail lengths. 

`tail_quantile_25.csv` and `tail_quantile_75.csv` provide the other quartiles, and can be used as a more robust measure of tail length variability than the standard deviation.



counts.csv
---

`counts.csv` contains various outputs, in a special Nesoni format. This was the original way `analyse-tail-counts` produced its results. If you've installed the Nesoni R component, it can be loaded in R with:

```
library(nesoni)
x <- read.grouped.table("counts.csv")
```

This turned out to be annoying way to do things. Hence the CSV files described above are now also created.