
<div style="margin-top: 5px;">

<img src="man/figures/hex.png" align="right" width="150"/>

</div>

## <span style="font-size:60px;">`ondisc`</span>

Single-cell datasets are growing in size, posing challenges as well as
opportunities for genomics researchers. `ondisc` is an R package for
computing on large-scale single-cell data, enabling users to analyze
data out-of-core on a laptop or distributed across multiple nodes on a
computing cluster. In both of these settings `ondisc` requires only a
few gigabytes of memory, even if the input data are tens of gigabytes in
size. `ondisc` is oriented toward single-cell differential expression,
single-cell gene co-expresssion, and single-cell CRISPR screen analyses.
`ondisc` is powered by several novel, efficient algorithms for large,
sparse matrices.

`ondisc` is a companion package to
[`sceptre`](https://katsevich-lab.github.io/sceptre/), an R package for
statistically rigorous and user-friendly single-cell CRISPR screen
analysis. Although `ondisc` and `sceptre` work best in conjunction,
`ondisc` can be used independently of `sceptre` (and conversely,
`sceptre` can be used independently of `ondisc`).

## Get started

We recommend getting started by reading the `ondisc` vignette. Please
submit issues on the `ondisc` [Github
page](https://github.com/timothy-barry/ondisc/issues).
