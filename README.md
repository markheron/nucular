<!-- README.md is generated from README.Rmd. Please edit that file -->
nucular
-------

The **nucular** package provides convienience functions for handling genomic data as lists of ff vectors (or vectors). The focus lies on comparing the genomic data by computing correlations and other comparison scores between the vector lists.

### Installation

``` r
devtools::install_bitbucket(repo="markheron/nucular", auth_user="user_name", password="your_password", keep_source=TRUE)
```

### Usage

``` r
library(nucular)

a <- list("chr1"=ff::ff(1:10), "chr2"=ff::ff(1:10))
b <- list("chr1"=ff::ff(1:10), "chr2"=ff::ff(10:1))

cor_ff_list(a, b)
#> [1] 0
```

### License

The **nucular** package is licensed under the GPLv3 (<http://www.gnu.org/licenses/gpl.html>).
