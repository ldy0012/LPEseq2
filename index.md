# LPEseq2

LPEseq2 is an R package for local pooled error-based ANOVA of RNA-seq
count data.

The package provides functions for preprocessing RNA-seq count data,
estimating intensity-dependent variance trends, and performing
multi-group differential expression analysis using a local pooled
error-based ANOVA framework.

## Installation

You can install the development version of LPEseq2 from GitHub:

``` r

install.packages("devtools")
devtools::install_github("ldy0012/LPEseq2")
```

## Optional dependencies

LPEseq2 supports several normalization methods.

For basic use, `library_size` normalization can be used without
additional Bioconductor packages.

For TMM and DESeq2 normalization, install the following packages:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("edgeR", "DESeq2"))
```

## Main functions

| Function | Description |
|----|----|
| [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md) | Checks count matrix and sample metadata, filters low-count genes, and performs normalization |
| [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md) | Performs local pooled error-based ANOVA for multi-group differential expression analysis |
| [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md) | Estimates an intensity-dependent variance trend and stores trimming information |

## Input format

The count matrix should contain genes as rows and samples as columns.

``` r

# Example count matrix
#        sample1 sample2 sample3 sample4
# gene1      100     120      80      95
# gene2       50      60      55      70
```

The sample metadata should contain samples as rows and variables as
columns.

``` r

# Example colData
#          group
# sample1 Control
# sample2 Control
# sample3 Treatment
# sample4 Treatment
```

The column names of the count matrix must match the row names of
`colData`.

## Example

``` r

library(LPEseq2)

set.seed(123)

counts <- matrix(
  rnbinom(1000, mu = 50, size = 10),
  nrow = 100,
  ncol = 10
)

rownames(counts) <- paste0("gene", 1:100)
colnames(counts) <- paste0("sample", 1:10)

colData <- data.frame(
  group = rep(c("Control", "Treatment"), each = 5),
  row.names = colnames(counts)
)

prep <- LPE_preprocess(
  counts = counts,
  colData = colData,
  design = ~ group,
  normalize.method = "library_size",
  verbose = FALSE
)

res <- LPE_ANOVA(
  object = prep,
  n.bin = 100,
  df = 10,
  use_weighted_between = FALSE,
  p.method = "chisq",
  verbose = FALSE
)

head(res)
```

### Advanced example with local fixed trimming

``` r

res_local <- LPE_ANOVA(
  object = prep,
  n.bin = 100,
  df = 10,
  trim.method = "local_fixed",
  d = 1.2,
  local.k = 3,
  min.local.bin.size = 10,
  use_weighted_between = TRUE,
  p.method = "chisq",
  verbose = FALSE
)

head(res_local)
attr(res_local, "trim.info")
```

## Trimming options

LPEseq2 supports three between-group trimming methods:

| Method | Description |
|----|----|
| `fixed` | Removes between-group differences whose absolute raw log2 difference is greater than or equal to the fixed threshold `d` |
| `local_fixed` | Uses an A-bin-specific local threshold estimated from within-group differences |
| `none` | Uses between-group-derived differences without outlier trimming |

Trimming is applied only to between-group-derived differences.  
Within-group pairwise differences are retained because they represent
replicate-based residual variation.

## Weighted between-group differences

The `use_weighted_between` option controls whether between-group-derived
differences are included in variance trend estimation.

When `use_weighted_between = FALSE`, LPEseq2 primarily uses within-group
pairwise differences for variance trend estimation.

When `use_weighted_between = TRUE`, LPEseq2 additionally uses weighted
between-group differences. This can be useful when replicate information
is limited, but between-group differences may contain true differential
expression signals.

Therefore, when `use_weighted_between = TRUE`, trimming methods such as
`fixed` or `local_fixed` can be used to reduce the influence of large
between-group differences.

### fixed

The `fixed` method applies one global threshold to raw between-group
log2 differences.

``` r

abs(D_between) < d
```

Here, `D_between` is the raw log2-scale difference between group means,
and `d` is the fixed trimming threshold.

### local_fixed

The `local_fixed` method estimates a local threshold for each
expression-intensity bin.

``` r

d_local = max(d, local.k * MAD(D_within_bin))
```

where:

- `D_within_bin` is the set of within-group raw log2 differences in the
  corresponding A-bin.
- `local.k` is a multiplier for the local MAD-based threshold.
- `d` is used as the minimum threshold.

If there are not enough within-group differences in an A-bin, the method
falls back to the fixed threshold `d`.

### none

The `none` method does not remove between-group-derived differences.

This option can be useful for diagnostic comparison, but true
differential expression signals may influence variance trend estimation
when many genes are differentially expressed.

## Output

[`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
returns a data frame with the following columns:

| Column       | Description                     |
|--------------|---------------------------------|
| `gene`       | Gene identifier                 |
| `mean`       | Mean expression value           |
| `var`        | Estimated local pooled variance |
| `MS_between` | Between-group mean square       |
| `F`          | Test statistic                  |
| `p.value`    | Raw p-value                     |
| `q.value`    | BH-adjusted p-value             |

Trimming information is stored as an attribute of the result object.

``` r

attr(res, "trim.info")
```

This includes the trimming method, threshold values, the number of
between-group-derived values before and after trimming, and the number
of removed values.

For `trim.method = "local_fixed"`, the threshold table contains
A-bin-specific local thresholds.

## Website

Package website: <https://ldy0012.github.io/LPEseq2/>

You can run LPEseq2 directly in the browser using the Shiny web tool:

[Launch LPEseq2 Shiny Web
Tool](https://ldy0012.shinyapps.io/lpeseq2app/)

## License

MIT
