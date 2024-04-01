
<div style="margin-top: 5px;">

<img src="man/figures/hex.png" align="right" width="150"/>

</div>

## <span style="font-size:60px;">`ondisc`</span>

Single-cell datasets are growing in size, posing challenges as well as
opportunities for genomics researchers. `ondisc` is an R package that
facilitates analysis of large-scale single-cell data out-of-core on a
laptop or distributed across tens to hundreds processors on a cluster or
cloud. In both of these settings, `ondisc` requires only a few gigabytes
of memory, even if the input data are tens of gigabytes in size.
`ondisc` mainly is oriented toward single-cell CRISPR screen analysis,
but ondisc also can be used for single-cell differential expression and
single-cell co-expression analyses. `ondisc` is powered by several new,
efficient algorithms for manipulating and querying large, sparse
expression matrices.

`ondisc` is a companion package to
[`sceptre`](https://katsevich-lab.github.io/sceptre/), an R package for
statistically rigorous and user-friendly single-cell CRISPR screen
analysis. Although `ondisc` and `sceptre` work best in conjunction,
`ondisc` can be used independently of `sceptre` (and conversely,
`sceptre` can be used independently of `ondisc`). Users can submit
issues on the `ondisc` [Github
page](https://github.com/timothy-barry/ondisc/issues).
