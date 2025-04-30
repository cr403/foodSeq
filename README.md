
<!-- README.md is generated from README.Rmd. Please edit that file -->

# foodSeq

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of foodSeq like so:

This chunk will be for future use when I can actually push to my github
rip :(

``` r
# # 1. Install bioconductor and devtools if needed
# # This package depends on Bioconductor packages (phyloseq, microbiome) 
# if (!requireNamespace("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
# }
# 
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# 
# library(devtools)
# library(BiocManager) 
# 
# # 2. Install the package
# devtools:install_github("cr403/foodSeq") 
# 
# # 3. Load library
# library(foodSeq)
```

For now, this is the workaround:

Within the foodSeq project, run the following chunk. Now, you should be
able to load this repo like a regular package in another project.

``` r
devtools::install_local() 
#> Skipping install of 'foodSeq' from a local remote, the SHA1 (0.0.0.90) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Loading Updates

To check for updates, you can run this chunk within your desired
project–**NOT within foodSeq project/repo**. This will uninstall and
reinstall the package without the need for restarting your R session.
It’s commented out here so that the Rmd file will knit.

``` r
# detach("package:foodSeq", unload = TRUE)
# devtools::install("~/Library/CloudStorage/Box-Box/project_davidlab/LAD_LAB_Personnel/Caroline_R/10_Repo/foodSeq")
# library(foodSeq)
```

ahh im having issues with loading Rmd and md files bc they got out of
sync
