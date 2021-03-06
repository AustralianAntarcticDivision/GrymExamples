
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Grym Examples

This package provides a range of examples for the Grym package located
at <https://github.com/AustralianAntarcticDivision/Grym>

Unless otherwise stated all examples are works in progress.

## Installing

The package is easily installed from GitHub, using the remotes package.

``` r
remotes::install_github("AustralianAntarcticDivision/GrymExamples", build_vignettes=TRUE)
```

If you don’t have `remotes` installed already, install it first.

``` r
install.packages("remotes")
```

GrymExamples does not otherwise require `remotes` for normal use.

## TODO

  - Make examples

  - Check examples for consistency of implementation

# Grym

<!-- badges: start -->

<!-- badges: end -->

The goal of Grym is to have provide an R based implementation of the
Generalised Yield Model by Constable and de la Mare (1996). This package
provides specific functions which users can use in combination to build
their own projection functions and packages.

Like most stock assessment packages Grym wont warn you or stop you doing
something stupid. The purpose of GrymExamples is to provide a variety of
examples of basic use to minimize that likelihood.

## Example

To access a list of the examples run the following code:

``` r
vignette(package="GrymExamples")
```

To access a specific example run the following code with the example
name such as:

``` r
vignette(topic="Krill96_1",package="GrymExamples")
```

## References

Constable, AJ, and WK de la Mare. 1996. “A Generalised Model for
Evaluating Yield and the Long-Term Status of Fish Stocks Under
Conditions of Uncertainty.” CCAMLR Science 3: 31–54.
