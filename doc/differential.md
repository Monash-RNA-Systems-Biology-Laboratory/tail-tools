
Differential expression and tail length testing
===

Tail Tools provides differential expression, as in RNA-seq, and also differential tail length testing.

Tail Tools uses a package called [Fitnoise](https://github.com/pfh/fitnoise) to perform differential testing. This is very similar to the popular Limma package. Like Limma, Fitnoise2 is based on the idea of linear models, a general idea that encompasses many popular statistical tests, such as t-tests, ANOVA-style tests, linear regression, and testing in the presence of a batch effect.

In a differential test based on linear models, there are two models, a null hypothesis model and an alternative hypothesis model. The alternative hypothesis contains all the terms of the null hypothesis plus some additional terms. That is, the null hypothesis is nested within the alternative hypothesis. The alternative hypothesis will therefore always fit better than the null hypothesis. The significance test is whether how much better it fits is more than would have been expected just by chance (an F-test).

Differential testing with Tail Tools may be performed using the `tail-tools test:` tool from the command line. Instances of `tail_tools.Test( )` may also be passed to the `tests` parameter of the pipeline.

The basic form this takes is:

    tail-tools test: test-output-dir pipeline-dir \
        null: <terms in null linear model> \
        alt: <additional terms in alternative linear model>

Terms are expressions involving sample names and sample tags, which specify columns of a model matrix. Matching samples are given a value of 1 in this column, and non-matching samples given a value of 0. Type `nesoni` for documentation on these, but briefly:

* `tag1` - an expression selecting tag1
* `all` - an expression selecting select all samples
* `expr1/expr2` - an expression selecting either expr1 or expr2
* `expr1:expr2` - an expression selecting sampes matching both expr1 and expr2
* `-expr1` - not having expr1
* `[expression]` - square brackets can be used to group more complex expressions

Sample names may be used as well as tags.

If the model matrix needs to contain values other than 0 and 1, there is an alternative form for the terms. For example if sample1, sample2, and sample3 needed to be given values -0.5, 2 and 0.5 in a column, this could be specified with:

    {-0.5}sample1,{2}sample2,{0.5}sample3




### Example 1

For example, if you had samples which you have tagged as being either "experimental" or "control" group when you ran the pipeline, with two replicates within each group, a possible test would be:

    tail-tools test: test-output-dir pipeline-dir null: control/experimental alt: experimental
  
From this, our null hypothesis model matrix would be a 4x1 matrix:

    1    control-replicate-1
    1    control-replicate-2
    1    experimental-replicate-1
    1    experimental-replicate-2

and our alternative hypothesis model matrix would be a 4x2 matrix:

    1 0  control-replicate-1
    1 0  control-replicate-2
    1 1  experimental-replicate-1
    1 1  experimental-replicate-2

`tail-tools test` will print out the model matrix when it is run, so you can check that is what you wanted.

Here the `experimental` term is in addition to the baseline level given by `control`, and therefore is the difference between the experimental and control groups. Terms coefficients fitted for terms in `alt:` are reported in the output of `tail-tools test:`.

Hint: In R, you can check how the coefficients of a linear model will be calculated using the `ginv` function from the library `MASS`, which computes the pseudo-inverse.


### Example 2

Say you had samples from three experimental groups, "group1", "group2", and "group3" and you wanted to know if there were any differences between the groups. A suitable ANOVA-style test would be:

    tail-tools test: test-output-dir pipeline-dir null: group1/group2/group3 alt: group2 group3

From this our null hypothesis model would be:

    1      group1-rep1
    1      group1-rep2
    1      group2-rep1
    1      group2-rep2
    1      group3-rep1
    1      group3-rep2

and our alternative hypothesis model would be:

    1 0 0  group1-rep1
    1 0 0  group1-rep2
    1 1 0  group2-rep1
    1 1 0  group2-rep2
    1 0 1  group3-rep1
    1 0 1  group3-rep2

The two additional coefficients will be group2 minus group1 and group3 minus group1,
group1 acting as a "baseline". Again, you can check this in R with `ginv`.

### See also

See the [Fitnoise](https://github.com/pfh/fitnoise) documentation for further details of the mathematics behind this.


