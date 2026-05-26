# LPEseq2

LPEseq2 is an R package for differential expression analysis of RNA-seq
count data using local pooled error-based ANOVA and conventional one-way
ANOVA.

The package provides functions for preprocessing RNA-seq count data,
estimating intensity-dependent variance trends, and performing
multi-group differential expression analysis. LPEseq2 supports LPE-ANOVA
for small-sample RNA-seq data, standard gene-wise one-way ANOVA for
larger sample sizes, and an automatic mode that selects the analysis
method based on group sample size.

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

For upper-quantile regression-based variance trend estimation, install:

```r
install.packages("quantreg")
```


    ## Main functions

    | Function | Description |
    |---|---|
    | `LPE_preprocess()` | Checks count matrix and sample metadata, filters low-count genes, and performs normalization |
    | `LPE_ANOVA()` | Performs LPE-ANOVA, standard one-way ANOVA, or automatic method selection for multi-group differential expression analysis |
    | `LPE_ANOVA_var()` | Estimates an intensity-dependent variance trend and stores trimming information |

    ## Input format

    The count matrix should contain genes as rows and samples as columns.

    ```r
    # Example count matrix
    #        sample1 sample2 sample3 sample4
    # gene1      100     120      80      95
    # gene2       50      60      55      70

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
  analysis.method = "auto",
  standard.min.group.n = 5,
  n.bin = 100,
  df = 10,
  use_weighted_between = FALSE,
  p.method = "chisq",
  verbose = FALSE
)

head(res)
attr(res, "analysis.method")
```

### Advanced example with IQR-based trimming

The `iqr` trimming method applies the conventional 1.5 × IQR boxplot
rule within expression-intensity A-bins. Within-group and
between-group-derived pairwise values are pooled before trimming, and
outliers are detected on the M-value scale used for variance trend
estimation.

``` r

res_iqr <- LPE_ANOVA(
  object = prep,
  analysis.method = "LPE",
  n.bin = 100,
  df = 10,
  trim.method = "iqr",
  trend.method = "mean_spline",
  use_weighted_between = TRUE,
  p.method = "chisq",
  verbose = FALSE
)

head(res_iqr)
attr(res_iqr, "trim.info")
attr(res_iqr, "trend.info")
attr(res_iqr, "base.var")
```

## Analysis methods

LPEseq2 supports three analysis modes through the `analysis.method`
argument:

| Method | Description |
|----|----|
| `LPE` | Uses local pooled error-based ANOVA. This mode estimates an intensity-dependent pooled variance trend and is useful for small-sample RNA-seq data. |
| `standard_anova` | Uses conventional gene-wise one-way ANOVA based on within-group residual variance. This mode is more appropriate when each group has enough samples to estimate gene-wise variance. |
| `auto` | Automatically selects the analysis method based on the minimum group sample size. |

### LPE-ANOVA mode

``` r

res_lpe <- LPE_ANOVA(
  object = prep,
  analysis.method = "LPE",
  n.bin = 100,
  df = 10,
  use_weighted_between = FALSE,
  p.method = "chisq",
  verbose = FALSE
)
```

### Standard one-way ANOVA mode

``` r

res_standard <- LPE_ANOVA(
  object = prep,
  analysis.method = "standard_anova",
  verbose = FALSE
)
```

### Auto mode

In `auto` mode, LPEseq2 selects standard one-way ANOVA when every group
has at least `standard.min.group.n` samples. Otherwise, LPE-ANOVA is
used.

``` r

res_auto <- LPE_ANOVA(
  object = prep,
  analysis.method = "auto",
  standard.min.group.n = 5,
  n.bin = 100,
  df = 10,
  use_weighted_between = FALSE,
  p.method = "chisq",
  verbose = FALSE
)

attr(res_auto, "analysis.method")
```

The default threshold is:

``` r

standard.min.group.n = 5
```

This threshold is a practical heuristic. It can be adjusted depending on
the study design and the expected reliability of gene-wise variance
estimation.

## Trimming options

LPEseq2 supports boxplot/IQR-based outlier trimming for pairwise values
used in variance trend estimation.

| Method | Description |
|----|----|
| `iqr` | Applies the conventional 1.5 × IQR boxplot rule within expression-intensity A-bins after pooling within-group and between-group-derived pairwise values |
| `none` | Uses pairwise values without outlier trimming |

When `trim.method = "iqr"`, LPEseq2 first pools within-group pairwise
values and between-group-derived values. The pooled values are divided
into expression-intensity A-bins, and the conventional boxplot rule is
applied within each bin:

``` r

lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR
```

Outlier detection is performed on the M-value scale because M is the
scale used for local pooled variance estimation.

The `none` method does not remove pairwise values before variance trend
estimation.

## Variance trend methods

LPEseq2 supports two methods for estimating the intensity-dependent
variance trend.

| Method | Description |
|----|----|
| `mean_spline` | Fits a smoothing spline to bin-level local pooled variance estimates |
| `quantile_regression` | Fits an upper-quantile regression trend to reduce potential variance underestimation |

### Mean smoothing spline

The default method is:

``` r

trend.method = "mean_spline"
```

This method estimates the average local pooled variance trend across
expression-intensity bins.

### Quantile regression

The quantile regression method can be used when a more conservative
variance trend is desired.

``` r

trend.method = "quantile_regression"
tau = 0.75
```

This method estimates an upper-quantile variance trend. It may reduce
potential variance underestimation for high-variability genes, but it
can also reduce statistical power.

The `quantreg` package is required only when
`trend.method = "quantile_regression"` is used.

## Output

[`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
returns a data frame containing gene-level test statistics.

The output format depends on the selected analysis method, but the main
columns are shared across methods.

### Common output columns

| Column | Description |
|----|----|
| `gene` | Gene identifier from the row names of the input expression matrix |
| `mean` | Mean expression value of the gene across samples |
| `var` | Estimated variance used in the test. For LPE-ANOVA, this is the local pooled variance predicted from the intensity-dependent variance trend. For standard ANOVA, this corresponds to the within-group residual variance |
| `MS_between` | Between-group mean square |
| `F` | Test statistic |
| `p.value` | Raw p-value |
| `q.value` | Benjamini-Hochberg adjusted p-value |
| `method` | Analysis method used for the result |

### Additional columns for standard one-way ANOVA

When `analysis.method = "standard_anova"` is used, the result may also
include the following columns:

| Column      | Description                    |
|-------------|--------------------------------|
| `MS_within` | Within-group mean square       |
| `df1`       | Numerator degrees of freedom   |
| `df2`       | Denominator degrees of freedom |

### Analysis method information

The selected analysis method is stored as an attribute of the result
object.

``` r

attr(res, "analysis.method")
```

If `analysis.method = "auto"` is used, this attribute shows which method
was actually selected after checking the group sample sizes.

The originally requested analysis mode can also be stored as:

``` r

attr(res, "requested.analysis.method")
```

The sample-size threshold used for auto mode can be checked with:

``` r

attr(res, "standard.min.group.n")
```

### Trimming information

When LPE-ANOVA is used, trimming information is stored as an attribute
of the result object.

``` r

attr(res, "trim.info")
```

This information includes the trimming method, the IQR rule information,
the number of pairwise values before and after trimming, and the number
of removed values.

For `trim.method = "iqr"`, the threshold table may include:

| Column | Description |
|----|----|
| `bin` | Expression-intensity A-bin index |
| `A_low` | Lower boundary of the A-bin |
| `A_high` | Upper boundary of the A-bin |
| `Q1` | First quartile of M values within the bin |
| `Q3` | Third quartile of M values within the bin |
| `IQR` | Interquartile range, calculated as `Q3 - Q1` |
| `lower_bound` | Lower outlier boundary, calculated as `Q1 - 1.5 * IQR` |
| `upper_bound` | Upper outlier boundary, calculated as `Q3 + 1.5 * IQR` |
| `n_values` | Number of pairwise values in the bin |
| `n_within` | Number of within-group pairwise values in the bin |
| `n_between` | Number of between-group-derived values in the bin |
| `n_removed` | Total number of removed values in the bin |
| `n_within_removed` | Number of removed within-group values in the bin |
| `n_between_removed` | Number of removed between-group-derived values in the bin |

When standard one-way ANOVA is selected, trimming information is not
available because the LPE variance trend is not estimated.

``` r

attr(res, "trim.info")
```

In this case, the trimming information is expected to be `NULL`.

### Variance trend information

When LPE-ANOVA is used, variance trend information is stored as an
attribute of the result object.

``` r

attr(res, "trend.info")
```

This includes the variance trend fitting method, the quantile level used
for quantile regression, and the degrees of freedom used for spline
fitting when available.

The bin-level variance points used for trend fitting can be checked
with:

``` r

attr(res, "base.var")
```

This object contains the representative expression-intensity value and
the estimated local pooled variance for each bin. \## Website

Package website: <https://ldy0012.github.io/LPEseq2/>

You can run LPEseq2 directly in the browser using the Shiny web tool:

[Launch LPEseq2 Shiny Web
Tool](https://ldy0012.shinyapps.io/lpeseq2app/)

## License

MIT
