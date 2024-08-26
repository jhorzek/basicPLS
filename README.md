
<!-- README.md is generated from README.Rmd. Please edit that file -->

# basicPLS

The objective of basicPLS is to provide a reduced implementation of
PLS-SEM. The package is currently intentionally incomplete, as we are
only implementing the features we need. The implementation follows that
found in [plspm](https://github.com/gastonstat/plspm).

The main function is `basicPLS::PLS` and fits PLS-SEM with formative and
pseudo-reflective measurement models (pseudo-reflective because PLS-SEM
always fits a formative model).

## Example:

``` r
library(basicPLS)
data_set <- basicPLS::satisfaction

# Both, measurement and structural model are specified using R's formulas:
PLS_result <- PLS(
  measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
                      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
                      LOY ~ loy1 + loy2 + loy3 + loy4,
                      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
                      SAT ~ sat1 + sat2 + sat3 + sat4,
                      VAL ~ val1 + val2 + val3 + val4),
  structure = alist(QUAL ~ EXPE,
                    EXPE ~ IMAG,
                    SAT ~ IMAG + EXPE + QUAL + VAL,
                    LOY ~ IMAG + SAT,
                    VAL ~ EXPE + QUAL),
  data = data_set)
PLS_result
#> 
#> #### PLS SEM Results ####
#> Composite Weights:
#> EXPE = 0.095*expe1 + 0.433*expe2 + 0.164*expe3 + 0.315*expe4 + 0.238*expe5 
#> IMAG = -0.032*imag1 +  0.221*imag2 +  0.493*imag3 +  0.038*imag4 +  0.472*imag5 
#> LOY =  0.560*loy1 +  0.072*loy2 +  0.515*loy3 + -0.077*loy4 
#> QUAL = 0.216*qual1 + 0.381*qual2 + 0.163*qual3 + 0.211*qual4 + 0.242*qual5 
#> SAT = 0.502*sat1 + 0.353*sat2 + 0.026*sat3 + 0.221*sat4 
#> VAL = 0.494*val1 + 0.190*val2 + 0.116*val3 + 0.384*val4 
#> 
#> Effects:
#> QUAL ~ 0.851*EXPE 
#> EXPE ~ 0.609*IMAG 
#> SAT ~ 0.205*IMAG + 0.001*EXPE + 0.086*QUAL + 0.626*VAL 
#> LOY ~ 0.227*IMAG + 0.560*SAT 
#> VAL ~ 0.152*EXPE + 0.648*QUAL

# Use confidence_intervals to bootstrap confidence intervals for all parameters:
ci <- confidence_intervals(PLS_result,
                           # increase for actual use:
                           R = 50)
ci$confidence_intervals$effects
#>       Parameter     Estimate    lower_ci   upper_ci
#> 1  QUAL <- EXPE 0.8511298309  0.81659304 0.87109779
#> 2  EXPE <- IMAG 0.6085981452  0.54769258 0.68113772
#> 3   SAT <- IMAG 0.2046077647  0.12973748 0.30249271
#> 4   SAT <- EXPE 0.0005350609 -0.09381798 0.08434397
#> 5   SAT <- QUAL 0.0856953262 -0.04374396 0.28566479
#> 6    SAT <- VAL 0.6256122447  0.46587410 0.75732287
#> 7   LOY <- IMAG 0.2274903742  0.11825588 0.35169751
#> 8    LOY <- SAT 0.5601144411  0.44916450 0.68132001
#> 9   VAL <- EXPE 0.1524291435  0.04463562 0.30011672
#> 10  VAL <- QUAL 0.6478498478  0.50508385 0.77578808
```

To switch to mode_A (“reflective”) for a composite, use `as_reflective`:

``` r
PLS_result <- PLS(
  measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
                      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
                      LOY ~ loy1 + loy2 + loy3 + loy4,
                      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
                      SAT ~ sat1 + sat2 + sat3 + sat4,
                      VAL ~ val1 + val2 + val3 + val4),
  structure = alist(QUAL ~ EXPE,
                    EXPE ~ IMAG,
                    SAT ~ IMAG + EXPE + QUAL + VAL,
                    LOY ~ IMAG + SAT,
                    VAL ~ EXPE + QUAL),
  as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
  data = data_set)
#> The algorithm took 5 iterations to converge.
PLS_result
#> 
#> #### PLS SEM Results ####
#> Composite Weights:
#> EXPE = 0.237*expe1 + 0.282*expe2 + 0.224*expe3 + 0.259*expe4 + 0.265*expe5 
#> IMAG = 0.206*imag1 + 0.297*imag2 + 0.308*imag3 + 0.181*imag4 + 0.285*imag5 
#> LOY = 0.378*loy1 + 0.248*loy2 + 0.375*loy3 + 0.218*loy4 
#> QUAL = 0.242*qual1 + 0.269*qual2 + 0.225*qual3 + 0.245*qual4 + 0.247*qual5 
#> SAT = 0.322*sat1 + 0.308*sat2 + 0.246*sat3 + 0.268*sat4 
#> VAL = 0.349*val1 + 0.294*val2 + 0.250*val3 + 0.325*val4 
#> 
#> Effects:
#> QUAL ~ 0.846*EXPE 
#> EXPE ~ 0.560*IMAG 
#> SAT ~ 0.186*IMAG + 0.008*EXPE + 0.138*QUAL + 0.580*VAL 
#> LOY ~ 0.289*IMAG + 0.474*SAT 
#> VAL ~ 0.118*EXPE + 0.660*QUAL
```

## Weighted PLS

basicPLS can compute weighted PLS estimates:

``` r
# data set with missings:
data_set <- basicPLS::satisfaction_NA

PLS_result <- PLS(
  measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
                      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
                      LOY ~ loy1 + loy2 + loy3 + loy4,
                      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
                      SAT ~ sat1 + sat2 + sat3 + sat4,
                      VAL ~ val1 + val2 + val3 + val4),
  structure = alist(QUAL ~ EXPE,
                    EXPE ~ IMAG,
                    SAT ~ IMAG + EXPE + QUAL + VAL,
                    LOY ~ IMAG + SAT,
                    VAL ~ EXPE + QUAL),
  data = data_set,
  # specify missing value treatment:
  imputation_function = mean_impute,
  # specify sample weights:
  sample_weights = data_set$sample_weights
  )
#> The algorithm took 10 iterations to converge.
PLS_result
#> 
#> #### PLS SEM Results ####
#> Composite Weights:
#> EXPE = 0.077*expe1 + 0.405*expe2 + 0.246*expe3 + 0.395*expe4 + 0.161*expe5 
#> IMAG = -0.101*imag1 +  0.459*imag2 +  0.416*imag3 +  0.054*imag4 +  0.369*imag5 
#> LOY = 0.394*loy1 + 0.051*loy2 + 0.669*loy3 + 0.026*loy4 
#> QUAL = 0.202*qual1 + 0.208*qual2 + 0.283*qual3 + 0.309*qual4 + 0.245*qual5 
#> SAT = 0.364*sat1 + 0.395*sat2 + 0.182*sat3 + 0.298*sat4 
#> VAL = 0.425*val1 + 0.193*val2 + 0.172*val3 + 0.436*val4 
#> 
#> Effects:
#> QUAL ~ 0.848*EXPE 
#> EXPE ~ 0.598*IMAG 
#> SAT ~  0.137*IMAG + -0.023*EXPE +  0.083*QUAL +  0.672*VAL 
#> LOY ~ 0.285*IMAG + 0.550*SAT 
#> VAL ~ 0.081*EXPE + 0.726*QUAL

get_r2(PLS_result)
#> $QUAL
#> [1] 0.7196226
#> 
#> $EXPE
#> [1] 0.3577394
#> 
#> $SAT
#> [1] 0.6832522
#> 
#> $LOY
#> [1] 0.5892923
#> 
#> $VAL
#> [1] 0.6326801
```
