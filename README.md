
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
#> 0.3438897 0.3427902 0.5401398 
#> 
#> $eta2
#>       y21       y22       y23 
#> 0.1705073 0.4438522 0.5743259 
#> 
#> $eta3
#>       y31       y32       y33 
#> 0.5543442 0.2261288 0.3930964

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
#> 0.3534211 0.3092768

head(fit_csem$Estimates$Construct_scores -
       PLS_result$components)
#>               eta1          eta2          eta3
#> [1,]  0.0001233242  0.0001536421  8.154323e-04
#> [2,]  0.0006259693  0.0021623806  1.204322e-03
#> [3,] -0.0007905236 -0.0002884053 -3.212220e-03
#> [4,] -0.0004493231  0.0008940692  1.578552e-03
#> [5,] -0.0011978506 -0.0034360563 -1.723869e-03
#> [6,]  0.0005098771 -0.0021601718 -2.639222e-05

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
