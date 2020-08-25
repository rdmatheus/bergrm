# `bergrm`: BerG Regression Model for Count Data

[![Build Status](https://travis-ci.com/travis-ci/travis-web.svg?branch=master)](https://travis-ci.com/travis-ci/travis-web)

> Rodrigo M. R. Medeiros
> <rodrigo.matheus@live.com>, IME-USP

The `bergrm` package provide a set of functions for a complete regression analysis of count data, in which it is assumed that the dependent variable follows a BerG distribution. The BerG regression was proposed by Bourguignon and Medeiros (2020), and can be used to fit count data with overdispersion, equidispersion (but are not Poisson-distributed), and underdispersion.

## Installation

You can install the current development version of `bergrm` from [GitHub](https://github.com/rdmatheus/sdlrm) with:

``` r
devtools::install_github("rdmatheus/bergrm")
```
To run the above command, it is necessary that the `devtools` package is previously installed on R. If not, install it using the following command:

``` r
install.packages("devtools")
```
After installing the devtools package, if you are using Windows, install the most current [RTools](https://cran.r-project.org/bin/windows/Rtools/) program. Finally, run the command `devtools::install_github("rdmatheus/bergrm")`, and then the package will be installed on your computer.
