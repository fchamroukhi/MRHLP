
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

## Overview

**MRHLP** is an R package for flexible and user-friendly
**probabilistic** joint **segmentation** of multivariate time series (or
multivariate structured longitudinal data) with smooth and/or abrupt
regime changes by a mixture model-based multiple regression approach
with a hidden logistic process (Multiple Regression model with a Hidden
Logistic Process (MRHLP)). The model is fitted by the EM algorithm.

## Installation

You can install the MRHLP package from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/MRHLP")
```

To build *vignettes* for examples of usage, type the command below
instead:

``` r
# install.packages("devtools")
devtools::install_github("fchamroukhi/MRHLP", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
```

Use the following command to display vignettes:

``` r
browseVignettes("MRHLP")
```

## Usage

``` r
library(MRHLP)

data("toydataset")

K <- 5 # Number of regimes (mixture components)
p <- 3 # Dimension of beta (order of the polynomial regressors)
q <- 1 # Dimension of w (order of the logistic regression: to be set to 1 for segmentation)
variance_type <- "heteroskedastic" # "heteroskedastic" or "homoskedastic" model

n_tries <- 1
max_iter <- 1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE

mrhlp <- emMRHLP(toydataset$x, toydataset[,c("y1", "y2", "y3")], K, p, q, 
                 variance_type, n_tries, max_iter, threshold, verbose, verbose_IRLS)
#> EM: Iteration : 1 || log-likelihood : -4975.54177550763
#> EM: Iteration : 2 || log-likelihood : -3108.34368262058
#> EM: Iteration : 3 || log-likelihood : -3083.17524290617
#> EM: Iteration : 4 || log-likelihood : -3052.50226046505
#> EM: Iteration : 5 || log-likelihood : -3020.60866761548
#> EM: Iteration : 6 || log-likelihood : -2967.37662637476
#> EM: Iteration : 7 || log-likelihood : -2948.61300516787
#> EM: Iteration : 8 || log-likelihood : -2945.45995948196
#> EM: Iteration : 9 || log-likelihood : -2937.99296980136
#> EM: Iteration : 10 || log-likelihood : -2924.28973590932
#> EM: Iteration : 11 || log-likelihood : -2901.25080505023
#> EM: Iteration : 12 || log-likelihood : -2859.88249265728
#> EM: Iteration : 13 || log-likelihood : -2858.05147227319
#> EM: Iteration : 14 || log-likelihood : -2856.38015373797
#> EM: Iteration : 15 || log-likelihood : -2854.68196733762
#> EM: Iteration : 16 || log-likelihood : -2852.69581368828
#> EM: Iteration : 17 || log-likelihood : -2849.93140687413
#> EM: Iteration : 18 || log-likelihood : -2846.34467342533
#> EM: Iteration : 19 || log-likelihood : -2843.82658697638
#> EM: Iteration : 20 || log-likelihood : -2842.75921489778
#> EM: Iteration : 21 || log-likelihood : -2842.2361309076
#> EM: Iteration : 22 || log-likelihood : -2841.91343876731
#> EM: Iteration : 23 || log-likelihood : -2841.66202744546
#> EM: Iteration : 24 || log-likelihood : -2841.41784741157
#> EM: Iteration : 25 || log-likelihood : -2841.14668922972
#> EM: Iteration : 26 || log-likelihood : -2840.82033081985
#> EM: Iteration : 27 || log-likelihood : -2840.39141033072
#> EM: Iteration : 28 || log-likelihood : -2839.74532802897
#> EM: Iteration : 29 || log-likelihood : -2838.62532237046
#> EM: Iteration : 30 || log-likelihood : -2836.64319641069
#> EM: Iteration : 31 || log-likelihood : -2833.87378876047
#> EM: Iteration : 32 || log-likelihood : -2831.75584262499
#> EM: Iteration : 33 || log-likelihood : -2831.16293539695
#> EM: Iteration : 34 || log-likelihood : -2831.0646784204
#> EM: Iteration : 35 || log-likelihood : -2831.06467491195

mrhlp$summary()
#> ----------------------
#> Fitted MRHLP model
#> ----------------------
#> 
#> MRHLP model with K = 5 regimes
#> 
#>  log-likelihood nu       AIC       BIC       ICL
#>       -2831.065 98 -2929.065 -3149.921 -3149.146
#> 
#> Clustering table:
#>   1   2   3   4   5 
#> 100 120 200 100 150 
#> 
#> 
#> ------------------
#> Regime 1 (K = 1):
#> 
#> Regression coefficients:
#> 
#>       Beta(d = 1)   Beta(d = 2) Beta(d = 3)
#> 1       0.4466558     0.8104534    -2.36719
#> X^1   -25.5100013   -20.5995360    32.75195
#> X^2   413.8717640   498.0085618  -541.38904
#> X^3 -1811.4612012 -2477.5546420  2523.64723
#> 
#> Covariance matrix:
#>                                   
#>  1.17712613  0.1114059  0.07303969
#>  0.11140591  0.8394152 -0.02442220
#>  0.07303969 -0.0244222  0.85240361
#> ------------------
#> Regime 2 (K = 2):
#> 
#> Regression coefficients:
#> 
#>     Beta(d = 1) Beta(d = 2) Beta(d = 3)
#> 1      21.30187   -4.108239    1.838238
#> X^1  -199.86512  112.953325  112.257782
#> X^2   905.60445 -449.623857 -493.914613
#> X^3 -1316.42937  581.197948  694.872075
#> 
#> Covariance matrix:
#>                                     
#>   1.0409982 -0.180821350 0.137568024
#>  -0.1808214  1.042169409 0.009699162
#>   0.1375680  0.009699162 0.754147599
#> ------------------
#> Regime 3 (K = 3):
#> 
#> Regression coefficients:
#> 
#>     Beta(d = 1) Beta(d = 2) Beta(d = 3)
#> 1     4.4721830    9.349642    6.349724
#> X^1   0.7467282  -33.315977   17.837763
#> X^2 -11.9302818   96.730621  -51.086769
#> X^3  16.1571109  -85.951201   42.760070
#> 
#> Covariance matrix:
#>                                     
#>   1.02026230 -0.04094457 -0.02544812
#>  -0.04094457  1.15656511  0.02852275
#>  -0.02544812  0.02852275  0.99750511
#> ------------------
#> Regime 4 (K = 4):
#> 
#> Regression coefficients:
#> 
#>     Beta(d = 1) Beta(d = 2) Beta(d = 3)
#> 1      1267.288   -840.5119   -10.37768
#> X^1   -5458.816   3613.7273    19.40201
#> X^2    7813.122  -5184.1100    14.37103
#> X^3   -3718.619   2475.7168   -29.55020
#> 
#> Covariance matrix:
#>                                       
#>   0.822157811  0.006792726 -0.03667011
#>   0.006792726  1.093351047 -0.07477892
#>  -0.036670114 -0.074778924  0.85425249
#> ------------------
#> Regime 5 (K = 5):
#> 
#> Regression coefficients:
#> 
#>     Beta(d = 1) Beta(d = 2) Beta(d = 3)
#> 1      194.7894    12.88268    483.8383
#> X^1   -658.4685   -45.73544  -1634.9482
#> X^2    753.1086    61.92925   1858.1529
#> X^3   -286.1078   -27.37495   -702.9064
#> 
#> Covariance matrix:
#>                                 
#>  1.1282728 0.25684915 0.02034990
#>  0.2568491 1.21055927 0.04414336
#>  0.0203499 0.04414336 0.77644297

mrhlp$plot()
```

<img src="man/figures/README-unnamed-chunk-5-1.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-2.png" style="display: block; margin: auto;" /><img src="man/figures/README-unnamed-chunk-5-3.png" style="display: block; margin: auto;" />
