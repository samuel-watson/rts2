# rts2
There are a wide variety of sources of case data that are finely spatially-resolved and time stamped to monitor disease epidemiology in real-time. 
These data may include daily hospital admissions, positive tests, or calls to public health telephone services. Geospatial statistical models provide 
a principled basis for generating predictions of disease epidemiology and its spatial and temporal distribution from these data. rts2 provides a set 
of functions to conduct disease surveillance using a real-time feed of spatio-temporal data of cases.

A detailed vignette with examples can be found at: <https://www.sam-watson.xyz/vignette.html>

## Installation and dependencies
The two main dependencies are `sf` and `rstan` (or `cmdstanr`). 

### sf
`sf` can be installed as normal:
```r
install.packages("sf")
```

### rstan
The `rst2` package uses Stan, which is a probabilistic programming language
used to draw posterior samples from Bayesian models. `rstan` is an R package
to interface with Stan. The trickiest part of the installation is getting a C++
compiler working with R. The easiest way of doing this on Windows is to install Rtools - see <https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started> for a detailed description. 


### cmdstanr
The package can be used with `cmdstanr` instead of `rstan`. `cmdstanr` is not 
necessary to run the package so if you've gotten `rstan` working then you
could skip this section. `cmdstanr` is lightweight, it can compile models much faster than `rstan`, and gets updates sooner, so it is worth looking into. 
 Installing cmdstan and cmdstanr is straightforward:
```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```
and then to install cmdstan:
```r
cmdstanr::install_cmdstan()
```
See <https://mc-stan.org/cmdstanr/> for more information.

### rts2
When all these packages are installed, rts2 can be installed:
```r
devtools::install_github("samuel-watson/rts2")
```
