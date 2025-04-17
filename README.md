
<!-- README.md is generated from README.Rmd. Please edit that file -->

# foodSeq

<!-- badges: start -->
<!-- badges: end -->

The goal of foodSeq is to …

## Installation

You can install the development version of foodSeq like so:

``` r
# 1. Install bioconductor and devtools if needed
# This package depends on Bioconductor packages (phyloseq, microbiome) 
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)
library(BiocManager) 

# 2. Install the package
devtools:install_github("cr403/foodSeq") 

# 3. Load library
library(foodSeq)
```
