
# Total empirical (TOTEM) statistics

## Overview

`totem` is a R package implementing the total empiricism (TOTEM)
framework. As a show case we have implemented the two-way $I$-test.

## Installation

`totem` is available at github

``` r
if (! requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
  
devtools::install_github(
  repo = "horyunchung/totem", 
  dependencies = TRUE, 
  build_vignettes = TRUE
)
```

The vignette shows how to use the package

``` r
vignette("I-test", package = "totem")
```
