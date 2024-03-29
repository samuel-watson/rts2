  <!-- badges: start -->
  [![R-CMD-check](https://github.com/samuel-watson/rts2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/samuel-watson/rts2/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# rts2

There are a wide variety of sources of case data that are finely spatially-resolved and time stamped to monitor disease epidemiology in real-time.
These data may include daily hospital admissions, positive tests, or calls to public health telephone services. Geospatial statistical models provide a
principled basis for generating predictions of disease epidemiology and its spatial and temporal distribution from these data.
`rts2` provides a set of functions to conduct disease surveillance using a real-time feed of spatial or spatio-temporal data of cases.

This readme is now out of date. It will be updated along with a manual describing use. The package supports LGCP with full Gaussian Process, Nearest Neigbour Gaussian Process, and Hibert Space Gaussian Process models, Bayesian (via Stan) and MCMC Maximum likelihood model fitting, "regional" data models for data aggregated to an irregular lattice, and more!
