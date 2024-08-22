test_that("Test Satisfaction - Formative Regression", {
  library(cSEM)
  library(basicPLS)
  model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT ~ IMAG + EXPE + QUAL + VAL
LOY ~ IMAG + SAT
VAL ~ EXPE + QUAL
# Measurement model
EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY <~ loy1 + loy2 + loy3 + loy4
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT <~ sat1 + sat2 + sat3 + sat4
VAL <~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
                   .PLS_modes = "modeB",
                   .disattenuate = FALSE)

  ## Estimate model with basicPLS
  basic_sat <- PLS(
    measurement = alist(
      EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
      LOY ~ loy1 + loy2 + loy3 + loy4,
      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
      SAT ~ sat1 + sat2 + sat3 + sat4,
      VAL ~ val1 + val2 + val3 + val4),
    structure = alist(
      QUAL ~ EXPE,
      EXPE ~ IMAG,
      SAT ~ IMAG + EXPE + QUAL + VAL,
      LOY ~ IMAG + SAT,
      VAL ~ EXPE + QUAL),
    data = satisfaction,
    path_estimation = "regression")

  for(i in rownames(csem_sat$Estimates$Path_estimates)){
    for(j in colnames(csem_sat$Estimates$Path_estimates)){
      if(csem_sat$Estimates$Path_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$effects[[i]][j] - csem_sat$Estimates$Path_estimates[i, j]), .01)
    }
  }

  for(i in rownames(csem_sat$Estimates$Weight_estimates)){
    for(j in colnames(csem_sat$Estimates$Weight_estimates)){
      if(csem_sat$Estimates$Weight_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$weights[[i]][j] - csem_sat$Estimates$Weight_estimates[i, j]), .01)
    }
  }

  for(i in rownames(csem_sat$Estimates$Loading_estimates)){
    for(j in colnames(csem_sat$Estimates$Loading_estimates)){
      if(csem_sat$Estimates$Loading_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$loadings[[i]][j] - csem_sat$Estimates$Loading_estimates[i, j]), .01)
    }
  }
})

test_that("Test Satisfaction - Formative Centroid", {
  library(cSEM)
  library(basicPLS)
  model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT ~ IMAG + EXPE + QUAL + VAL
LOY ~ IMAG + SAT
VAL ~ EXPE + QUAL
# Measurement model
EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY <~ loy1 + loy2 + loy3 + loy4
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT <~ sat1 + sat2 + sat3 + sat4
VAL <~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
                   .PLS_modes = "modeB",
                   .PLS_weight_scheme_inner = "centroid",
                   .disattenuate = FALSE)

  ## Estimate model with basicPLS
  basic_sat <- PLS(
    measurement = alist(
      EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
      LOY ~ loy1 + loy2 + loy3 + loy4,
      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
      SAT ~ sat1 + sat2 + sat3 + sat4,
      VAL ~ val1 + val2 + val3 + val4),
    structure = alist(
      QUAL ~ EXPE,
      EXPE ~ IMAG,
      SAT ~ IMAG + EXPE + QUAL + VAL,
      LOY ~ IMAG + SAT,
      VAL ~ EXPE + QUAL),
    data = satisfaction,
    path_estimation = "centroid")

  for(i in rownames(csem_sat$Estimates$Path_estimates)){
    for(j in colnames(csem_sat$Estimates$Path_estimates)){
      if(csem_sat$Estimates$Path_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$effects[[i]][j] - csem_sat$Estimates$Path_estimates[i, j]), .01)
    }
  }

  for(i in rownames(csem_sat$Estimates$Weight_estimates)){
    for(j in colnames(csem_sat$Estimates$Weight_estimates)){
      if(csem_sat$Estimates$Weight_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$weights[[i]][j] - csem_sat$Estimates$Weight_estimates[i, j]), .01)
    }
  }
  for(i in rownames(csem_sat$Estimates$Loading_estimates)){
    for(j in colnames(csem_sat$Estimates$Loading_estimates)){
      if(csem_sat$Estimates$Loading_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(basic_sat$loadings[[i]][j] - csem_sat$Estimates$Loading_estimates[i, j]), .01)
    }
  }
})

test_that("Test Satisfaction - Bootstrap", {
  library(cSEM)
  library(basicPLS)
  set.seed(123)
  model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT ~ IMAG + EXPE + QUAL + VAL
LOY ~ IMAG + SAT
VAL ~ EXPE + QUAL
# Measurement model
EXPE <~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY <~ loy1 + loy2 + loy3 + loy4
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT <~ sat1 + sat2 + sat3 + sat4
VAL <~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
                   .PLS_modes = "modeB",
                   .resample_method = "bootstrap",
                   .R = 100,
                   .disattenuate = FALSE)
  summarized <- summarize(csem_sat)

  ## Estimate model with basicPLS
  basic_sat <- PLS(
    measurement = alist(
      EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
      LOY ~ loy1 + loy2 + loy3 + loy4,
      QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
      SAT ~ sat1 + sat2 + sat3 + sat4,
      VAL ~ val1 + val2 + val3 + val4),
    structure = alist(
      QUAL ~ EXPE,
      EXPE ~ IMAG,
      SAT ~ IMAG + EXPE + QUAL + VAL,
      LOY ~ IMAG + SAT,
      VAL ~ EXPE + QUAL),
    data = satisfaction,
    path_estimation = "regression")

  cis <- basicPLS::confidence_intervals(PLS_result = basic_sat,
                                        include_loadings = TRUE,
                                        R = 100)

  for(i in summarized$Estimates$Path_estimates$Name){
    if(!i %in% cis$confidence_intervals$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%L`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$lower_ci[cis$confidence_intervals$Parameter == i]), .1)
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%U`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$upper_ci[cis$confidence_intervals$Parameter == i]), .1)
  }
  for(i in summarized$Estimates$Weight_estimates$Name){
    if(!i %in% cis$confidence_intervals$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%L`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$lower_ci[cis$confidence_intervals$Parameter == i]), .1)
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%U`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$upper_ci[cis$confidence_intervals$Parameter == i]), .1)
  }

  for(i in summarized$Estimates$Loading_estimates$Name){
    if(!i %in% cis$confidence_intervals$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%L`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$lower_ci[cis$confidence_intervals$Parameter == i]), .1)
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%U`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$upper_ci[cis$confidence_intervals$Parameter == i]), .1)
  }
})
