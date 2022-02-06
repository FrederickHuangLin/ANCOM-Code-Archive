
# User Manual for [ANCOM](https://www.tandfonline.com/doi/full/10.3402/mehd.v26.27663) v2.1

[![R-CMD-check](https://github.com/xec-cm/ANCOM/workflows/R-CMD-check/badge.svg)](https://github.com/xec-cm/ANCOM/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3577802.svg)](https://doi.org/10.5281/zenodo.3577802)

The current code implements ANCOM in cross-sectional and longitudinal
datasets while allowing the use of covariates.

## Installation

To get a bug fix or to use a feature from the development version, you
can install the development version of ANCOM from GitHub.

``` r
if(!requireNamespace("devtools")){
  install.packages("devtools")
}

devtools::install_github("xec-cm/ANCOM")
```

## Usage

``` r
library(ANCOM)

# Step 0: Load or define input data

feature_table <- ANCOM:::moving_pics_table_taxa
meta_data <- ANCOM:::moving_pics_sample_metadata

# Step 1: Data Preprocessing

prepro <- feature_table_pre_process(
  feature_table = feature_table,
  meta_data = meta_data,
  sample_var = "Sample.ID",
  group_var = NULL,
  out_cut = 0.05,
  zero_cut = 0.90,
  lib_cut = 1000,
  neg_lb = FALSE
)

# Step 2: ANCOM

res <- ANCOM(
  feature_table = prepro$feature_table, 
  meta_data = prepro$meta_data,
  struc_zero = prepro$structure_zeros, 
  main_var = "Subject", 
  p_adj_method = "BH", 
  alpha = 0.05, 
  adj_formula = NULL, 
  rand_formula = NULL
)

res$fig
```

<img src="README_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

## Getting help

If you encounter a clear bug, please file an issue with a minimal
reproducible example on
[GitHub](https://github.com/xec-cm/ANCOM/issues).

## Code of Conduct

Please note that the ANCOM project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
