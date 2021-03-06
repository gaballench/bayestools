
<!-- README.md is generated from README.Rmd. Please edit that file -->
bayestools
==========

Overview
--------

An R package for pre- and post-processing in phylogenetic and divergence-times bayesian inference. You can use it for processing output of phylogenetic bayesian programs (such as [Beast](https://www.beast2.org/) or [MrBayes](http://nbisweden.github.io/MrBayes/index.html)) or to prepare input files and decision making. It also supports sensitivity calculations using intersection (Ballen, in prep.) for comparing two given densities (e.g. prior vs. posterior).

Installation
------------

For now only the github development version is available. You can install it using `devtools`:

``` r
install.packages("devtools")
devtools::install_github("gaballench/bayestools")
```

Once the package makes its way into CRAN the standard `install.packages` may be used for installing stable versions.

Problems?
---------

If you find a bug or unexpected behavior, please [file an issue](https://github.com/gaballench/bayestools/issues).
