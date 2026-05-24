
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMetaTraj

`scMetaTraj` is an R package for metabolism-centered trajectory analysis
in single-cell data. The package scores curated metabolic modules,
embeds cells in a metabolic feature space, identifies metabolic
subclusters, and infers metabolic pseudotime (mPT) from graph distances
in that space.

## Installation

The development version can be installed from GitHub:

``` r
# install.packages("pak")
pak::pak("Dr-xyGreg/scMetaTraj")
```

After CRAN release, install with:

``` r
install.packages("scMetaTraj")
```

## Core workflow

The package is organized around five main steps:

1.  score pathway-level metabolic activity with `scMetaTraj_score()`
2.  embed cells with `scMetaTraj_embed()`
3.  identify metabolic subclusters with `scMetaTraj_cluster()`
4.  infer metabolic pseudotime with `scMetaTraj_infer()`
5.  quantify dynamic module changes with `scMetaTraj_trend()` and
    related plotting helpers

## Minimal example

The example below uses a small simulated expression matrix and toy
metabolic gene sets so that it can run without private data files.

``` r
library(scMetaTraj)

set.seed(123)

expr <- matrix(
  rexp(12 * 80, rate = 1),
  nrow = 12,
  ncol = 80,
  dimnames = list(
    c(
      "HK1", "PFKP", "LDHA", "GPI", "CS", "ACO2",
      "NDUFA1", "COX4I1", "G6PD", "PGD", "ACLY", "FASN"
    ),
    paste0("Cell", seq_len(80))
  )
)

gene_sets <- list(
  Glycolysis = c("HK1", "PFKP", "LDHA", "GPI"),
  TCA = c("CS", "ACO2"),
  OXPHOS = c("NDUFA1", "COX4I1"),
  PPP = c("G6PD", "PGD"),
  Lipid = c("ACLY", "FASN")
)

scores <- scMetaTraj_score(expr, gene_sets, min_genes = 2, scale = FALSE)
emb_pca <- scMetaTraj_embed(scores, method = "PCA", n_pcs = 4)
clusters <- scMetaTraj_cluster(emb_pca, k = 10, method = "louvain")
traj <- scMetaTraj_infer(emb_pca, k = 10, root_mode = "pc1_min")

trend_df <- scMetaTraj_trend(
  scores = scores[, "Glycolysis"],
  mPT = traj$mPT,
  n_bins = 20
)
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : pseudoinverse used at 0.425
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : neighborhood radius 0.1
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : reciprocal condition number 0
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : There are other near singularities as well. 0.01
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : Chernobyl! trL>n 14
#> Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
#> : Chernobyl! trL>n 14
#> Warning in sqrt(sum.squares/one.delta): NaNs produced

head(trend_df)
#>    mPT_bin    score score_smooth
#> 1    0.025 2.017809     2.017809
#> 7    0.325 1.388076     1.388076
#> 8    0.375 1.107979     1.107979
#> 9    0.425 1.875511     1.875511
#> 11   0.525 1.026035     1.026035
#> 12   0.575 1.029004     1.029004
table(clusters)
#> clusters
#>  1  2  3  4  5  6 
#> 22  6 14 11 15 12
```

## Package scope

`scMetaTraj` is intended to complement transcriptome-centered trajectory
tools. It is most useful when the primary biological question concerns
continuous metabolic remodeling rather than transcriptomic lineage
structure alone.

## Documentation

- The main workflow vignette is in `vignettes/scMetaTraj-workflow.Rmd`.
- Function-level help is available with `?scMetaTraj_score`,
  `?scMetaTraj_infer`, and related topics.
