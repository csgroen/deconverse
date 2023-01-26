
<!-- README.md is generated from README.Rmd. Please edit that file -->

## `deconverse`: bulk RNA-seq deconvolution benchmark using single-cell reference profiles <img src="man/figures/logo.png" align="right" width="120"/>

[![](https://img.shields.io/badge/devel%20version-0.2.0-blue.svg)](https://github.com/csgroen/deconverse)
[![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/csgroen/deconverse.svg)](https://github.com/csgroen/deconverse/commits/master)

**Note**: This package is still immature, APIs are not yet well
documented and may change

### Installation

To install, run:

``` r
remotes::install_github("csgroen/ggheatmapper")
remotes::install_github("jamesotto852/ggdensity")
remotes::install_github("Danko-Lab/BayesPrism/BayesPrism")
remotes::install_github("csgroen/deconverse")
```

`docker` or `singularity` must be available to run some deconvolution
methods. To install, see: <https://docs.docker.com/get-docker/>

To install CIBERSORTx docker, run:

``` r
library(deconverse)
install_cibersortx()
```

## Usage

[See the PBMC tutorial to get
started](http://csgroen.github.io/deconverse/articles/intro_pbmc.html).
