
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
                 .PLS_modes = "modeB")

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
#> 0.3438901 0.3427898 0.5401398 
#> 
#> $eta2
#>       y21       y22       y23 
#> 0.1705072 0.4438518 0.5743262 
#> 
#> $eta3
#>       y31       y32       y33 
#> 0.5543440 0.2261283 0.3930971

fit_csem$Estimates$Path_estimates
#>           eta1      eta2 eta3
#> eta1 0.0000000 0.0000000    0
#> eta2 0.5096772 0.0000000    0
#> eta3 0.3534568 0.3092132    0
PLS_result$effects
#> $eta2
#>      eta1 
#> 0.5096527 
#> 
#> $eta3
#>      eta1      eta2 
#> 0.3534210 0.3092768

head(fit_csem$Estimates$Construct_scores -
       PLS_result$components)
#>               eta1          eta2          eta3
#> [1,]  0.0001228359  0.0001530722  8.159066e-04
#> [2,]  0.0006256111  0.0021622667  1.204252e-03
#> [3,] -0.0007905130 -0.0002880882 -3.213786e-03
#> [4,] -0.0004487291  0.0008937740  1.579693e-03
#> [5,] -0.0011977423 -0.0034357471 -1.724097e-03
#> [6,]  0.0005104761 -0.0021602940 -2.691448e-05

assess(fit_csem, "r2")
#> ________________________________________________________________________________
#> 
#>  Construct        R2      
#>  eta2           0.2598    
#>  eta3           0.3320    
#> ________________________________________________________________________________
get_r2(PLS_result)
#> $eta2
#> [1] 0.2597459
#> 
#> $eta3
#> [1] 0.3319737
```
