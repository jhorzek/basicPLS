test_that("Test consistent PLS", {
  rm(list = ls())
  library(plsR)
  data_set <- plsR::satisfaction

  # Both, measurement and structural model are specified using R's formulas:
  basic_sat <- PLS(
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
    data = data_set,
    sample_weights = data_set$sample_weights)

  # comparison to results from smartPLS
  load(system.file("testdata",
                   "expected_satisfaction_consistent_weighted.RData",
                   package="plsR"))
  expected <- expected_satisfaction_consistent_weighted
  rm(expected_satisfaction_consistent_weighted)
  for(i in expected$effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$effects$path[expected$effects$Parameter == i] -
                              basic_sat$effects[[dep]][pred]), .01)
  }

  for(i in expected$total_effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$total_effect$total_effect[expected$total_effect$Parameter == i] -
                              basic_sat$total_effects[[dep]][pred]), .01)
  }

  for(i in expected$weights$Parameter){
    # weights are incorrectly encoded in the PLS SEM output (the arrow shows the
    # wrong way).
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$weights$weight[expected$weights$Parameter == i] -
                              basic_sat$weights[[dep]][pred]), .01)
  }

  for(i in expected$loadings$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$loadings$loading[expected$loadings$Parameter == i] -
                              basic_sat$loadings[[dep]][pred]), .01)
  }
})

test_that("Test consistent PLS - Mean Impute", {
  rm(list = ls())
  library(plsR)
  data_set <- plsR::satisfaction_NA

  # smartPLS does not take sampling weights into account when imputing data.
  # Therefore, we also don't do this here
  data_set <- mean_impute(data_set, weights = rep(1, nrow(data_set)))

  # Both, measurement and structural model are specified using R's formulas:
  basic_sat <- PLS(
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
    data = data_set,
    sample_weights = data_set$sample_weights,
    imputation_function = mean_impute)

  # comparison to results from smartPLS
  load(system.file("testdata",
                   "expected_satisfaction_consistent_NA_weighted.RData",
                   package="plsR"))
  expected <- expected_satisfaction_consistent_NA_weighted
  rm(expected_satisfaction_consistent_NA_weighted)
  for(i in expected$effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$effects$path[expected$effects$Parameter == i] -
                              basic_sat$effects[[dep]][pred]), .01)
  }

  for(i in expected$total_effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$total_effect$total_effect[expected$total_effect$Parameter == i] -
                              basic_sat$total_effects[[dep]][pred]), .01)
  }

  for(i in expected$weights$Parameter){
    # weights are incorrectly encoded in the PLS SEM output (the arrow shows the
    # wrong way).
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$weights$weight[expected$weights$Parameter == i] -
                              basic_sat$weights[[dep]][pred]), .01)
  }

  for(i in expected$loadings$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$loadings$loading[expected$loadings$Parameter == i] -
                              basic_sat$loadings[[dep]][pred]), .01)
  }

})
