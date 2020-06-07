wiqid
=====

[![CRAN status](https://www.r-pkg.org/badges/version/wiqid)](https://cran.r-project.org/web/packages/wiqid/index.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/last-month/wiqid)](https://www.r-pkg.org/services)


## Quick and Dirty functions for Wildlife simulations.

This package began life as a collection of little functions to estimate occupancy or abundance or survival directly in `R`, without having to export the data to `PRESENCE` or `MARK` or `EstimateS`, put together for a simulations workshop.

Recently we have started to use it in basic study design and data analysis workshops, where we need Quick and Dirty ways to do both MLE and Bayesian analysis for example data sets, without having to teach people to use the industry-standard software.

So we want to have data in `R` in ordinary data frames, and to use `R`'s formula interface and "twiddle" notation to define models, as in `lm`, `glm` and the `secr` package. Bayesian versions should produce output with the same look and feel as the `BEST` package. We also want it to work on all platforms: Windows, Linux, and Mac.

### Recent changes
The `Bwiqid` class for MCMC output has been scrapped and many of the functions for handling MCMC output have moved to the new package `mcmcOutput` on which `wiqid` now Depends.
