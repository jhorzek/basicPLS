
<!-- README.md is generated from README.Rmd. Please edit that file -->

# basicPLS

> basicPLS should not be used for any relevant data analysis; We
> recommend cSEM instead.

The objective of basicPLS is to provide a simple implementation of
PLS-SEM to better understand the algorithm underlying the modeling
approach. The implementation follows that found in
[plspm](https://github.com/gastonstat/plspm).

The main function is `basicPLS::PLS` and fits PLS-SEM with formative
measurement models.

## Example:

``` r
library(cSEM)
library(basicPLS)
data_set <- 10*cSEM::threecommonfactors + 2
PLS_result <- PLS(measurement = alist(eta1 ~ y11 + y12 + y13,
                                      eta2 ~ y21 + y22 + y23,
                                      eta3 ~ y31 + y32 + y33),
                  structure = alist(eta2 ~ eta1,
                                    eta3 ~ eta1 + eta2),
                  data = data_set)
PLS_result
#> 
#> #### PLS SEM Results ####
#> Component Weights:
#> eta1 = 0.344*y11 + 0.343*y12 + 0.540*y13 
#> eta2 = 0.171*y21 + 0.444*y22 + 0.574*y23 
#> eta3 = 0.554*y31 + 0.226*y32 + 0.393*y33 
#> 
#> Effects:
#> eta2 ~ 0.510*eta1 
#> eta3 ~ 0.353*eta1 + 0.309*eta2

# same thing with cSEM
model <- "
# Structural model
eta2 ~ eta1
eta3 ~ eta1 + eta2

# measurement model
eta1 <~ y11 + y12 + y13
eta2 <~ y21 + y22 + y23
eta3 <~ y31 + y32 + y33
"
fit_csem <- csem(data_set,
                 model,
                 .PLS_modes = "modeB",
                 .PLS_weight_scheme_inner = "centroid")

# Let's compare the results:
fit_csem$Estimates$Weight_estimates
#>            y11       y12       y13       y21       y22       y23       y31
#> eta1 0.3438924 0.3427876 0.5401396 0.0000000 0.0000000 0.0000000 0.0000000
#> eta2 0.0000000 0.0000000 0.0000000 0.1705062 0.4438504 0.5743282 0.0000000
#> eta3 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.5543436
#>            y32      y33
#> eta1 0.0000000 0.000000
#> eta2 0.0000000 0.000000
#> eta3 0.2261267 0.393099
PLS_result$weights
#> $eta1
#>       y11       y12       y13 
#> 0.3438924 0.3427876 0.5401396 
#> 
#> $eta2
#>       y21       y22       y23 
#> 0.1705062 0.4438504 0.5743282 
#> 
#> $eta3
#>       y31       y32       y33 
#> 0.5543436 0.2261267 0.3930990

fit_csem$Estimates$Path_estimates
#>           eta1      eta2 eta3
#> eta1 0.0000000 0.0000000    0
#> eta2 0.5096528 0.0000000    0
#> eta3 0.3534210 0.3092769    0
PLS_result$effects
#> $eta2
#>      eta1 
#> 0.5096528 
#> 
#> $eta3
#>      eta1      eta2 
#> 0.3534210 0.3092769

head(fit_csem$Estimates$Construct_scores -
       PLS_result$components)
#>               eta1          eta2          eta3
#> [1,] -6.938894e-18 -1.942890e-16 -2.012279e-16
#> [2,]  0.000000e+00  2.220446e-16  2.220446e-16
#> [3,]  0.000000e+00  0.000000e+00  3.885781e-16
#> [4,]  0.000000e+00  1.110223e-16 -2.220446e-16
#> [5,]  0.000000e+00 -1.110223e-15  0.000000e+00
#> [6,]  1.110223e-16 -1.554312e-15  4.440892e-16

assess(fit_csem, "r2")
#> ________________________________________________________________________________
#> 
#>  Construct        R2      
#>  eta2           0.2597    
#>  eta3           0.3320    
#> ________________________________________________________________________________
get_r2(PLS_result)
#> $eta2
#> [1] 0.2597459
#> 
#> $eta3
#> [1] 0.3319737
```

Using regression weights instead of the default centroid weights:

``` r
PLS_result <- PLS(measurement = alist(eta1 ~ y11 + y12 + y13,
                                      eta2 ~ y21 + y22 + y23,
                                      eta3 ~ y31 + y32 + y33),
                  structure = alist(eta2 ~ eta1,
                                    eta3 ~ eta1 + eta2),
                  data = data_set, 
                  path_estimation = "regression")
#> The algorithm took 5 iterations to converge.
PLS_result
#> 
#> #### PLS SEM Results ####
#> Component Weights:
#> eta1 = 0.344*y11 + 0.343*y12 + 0.540*y13 
#> eta2 = 0.172*y21 + 0.443*y22 + 0.574*y23 
#> eta3 = 0.555*y31 + 0.227*y32 + 0.392*y33 
#> 
#> Effects:
#> eta2 ~ 0.510*eta1 
#> eta3 ~ 0.353*eta1 + 0.309*eta2

fit_csem <- csem(data_set,
                 model,
                 .PLS_modes = "modeB",
                 .PLS_weight_scheme_inner = "path")

# Let's compare the results:
fit_csem$Estimates$Weight_estimates
#>           y11       y12       y13       y21       y22       y23      y31
#> eta1 0.343982 0.3427578 0.5400872 0.0000000 0.0000000 0.0000000 0.000000
#> eta2 0.000000 0.0000000 0.0000000 0.1715504 0.4431438 0.5743126 0.000000
#> eta3 0.000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.554802
#>            y32       y33
#> eta1 0.0000000 0.0000000
#> eta2 0.0000000 0.0000000
#> eta3 0.2270798 0.3917486
PLS_result$weights
#> $eta1
#>       y11       y12       y13 
#> 0.3439750 0.3427653 0.5400870 
#> 
#> $eta2
#>       y21       y22       y23 
#> 0.1715545 0.4431479 0.5743062 
#> 
#> $eta3
#>       y31       y32       y33 
#> 0.5548014 0.2270858 0.3917439

fit_csem$Estimates$Path_estimates
#>           eta1      eta2 eta3
#> eta1 0.0000000 0.0000000    0
#> eta2 0.5096772 0.0000000    0
#> eta3 0.3534568 0.3092132    0
PLS_result$effects
#> $eta2
#>      eta1 
#> 0.5096772 
#> 
#> $eta3
#>      eta1      eta2 
#> 0.3534571 0.3092129

head(fit_csem$Estimates$Construct_scores -
       PLS_result$components)
#>               eta1          eta2          eta3
#> [1,]  9.311063e-06  9.020625e-06 -4.743024e-06
#> [2,]  6.818080e-06  3.571668e-07  2.819106e-06
#> [3,] -2.689275e-07 -4.565194e-06  1.248594e-05
#> [4,] -1.167041e-05  3.850598e-06 -8.770794e-06
#> [5,] -2.524142e-06 -1.605349e-07 -7.439025e-07
#> [6,] -1.210029e-05  8.019773e-06  8.083453e-06

assess(fit_csem, "r2")
#> ________________________________________________________________________________
#> 
#>  Construct        R2      
#>  eta2           0.2598    
#>  eta3           0.3320    
#> ________________________________________________________________________________
get_r2(PLS_result)
#> $eta2
#> [1] 0.2597708
#> 
#> $eta3
#> [1] 0.3319534
```
