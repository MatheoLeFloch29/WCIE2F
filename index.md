# <img src="man/figures/stickers_mixWCE.png" align="right" width=130 style="margin-right: 0px;vertical-align:middle"/> <span style="font-size:38px"> Estimation of the Weighted Cumulative Effect with Two-Level </span>


&nbsp;

<p align="justify">
The mixWCE package implements the WCIE approach to estimate the effect of a time-varying exposure on different types of outcomes such as binary, survival, or conditional logistic models. The WCIE allows modeling the influence of past exposure using a weighted cumulative framework.
All models are estimated using the default optimization algorithms implemented in standard functions such as glm.
The package also includes functions to simulate data, as well as utility tools for a complete statistical analysis.
</p>

A detailed companion paper isn't available in Journal of Statistical Software :

<p align="justify">
</p>


&nbsp;

## Install the package

The mixWCE package needs version >= 3.5.0 or newer of the R software.

To get the most recent update, install it from github (to update i think) :

``` r
remotes::install_github("MatheoLeFloch29/mixWCE")
```

The mixWCE package depends on other R package, namely :

- survival (>=2.37-2) for dealing with the survival outcomes
- MASS
- lcmm
- splines
- ggplot2
- dplyr


To run the examples proposed in this website, the following package are also needed :

- ggplot2
- dplyr
- splines

&nbsp;

## Documentation

<p align="justify">
This website is intended to help the mixWCE users in their statistical analyses. It provides an [overview](articles/mixWCE.html) of the package, several [vignettes](articles/index.html), and the [help pages](reference/index.html) of all functions included in the mixWCE package.
</p>

<p align="justify">
Further issues and questions about the use of the mixWCE package are reported on the github issue page <https://github.com/MatheoLeFloch29/mixWCE/issues>. (to update)
Please check both opened and closed issues to make sure that the topic has not already been treated before creating a new issue. To report a bug, please provide a reproducible example.
</p>
