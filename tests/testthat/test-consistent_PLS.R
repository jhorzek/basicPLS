test_that("Test reliability and debiased correlation", {
  rm(list = ls())
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
    as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
    data = data_set)

  reliability_basic <- PLS_result$reliability

  library(cSEM)
  model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT ~ IMAG + EXPE + QUAL + VAL
LOY ~ IMAG + SAT
VAL ~ EXPE + QUAL
# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY =~ loy1 + loy2 + loy3 + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT =~ sat1 + sat2 + sat3 + sat4
VAL =~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
                   .disattenuate = TRUE)

  for(i in names(csem_sat$Estimates$Reliabilities))
    testthat::expect_lt(abs(csem_sat$Estimates$Reliabilities[i] - reliability_basic[i]), 1e-4)

  testthat::expect_true(all(abs(PLS_result$debiased_correlation[rownames(csem_sat$Estimates$Construct_VCV), colnames(csem_sat$Estimates$Construct_VCV)] -
                                  csem_sat$Estimates$Construct_VCV) < 1e-3))

  # Estimates
  for(i in rownames(csem_sat$Estimates$Path_estimates)){
    for(j in colnames(csem_sat$Estimates$Path_estimates)){
      if(csem_sat$Estimates$Path_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(PLS_result$effects[[i]][j] - csem_sat$Estimates$Path_estimates[i, j]), .01)
    }
  }

  for(i in rownames(csem_sat$Estimates$Weight_estimates)){
    for(j in colnames(csem_sat$Estimates$Weight_estimates)){
      if(csem_sat$Estimates$Weight_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(PLS_result$weights[[i]][j] - csem_sat$Estimates$Weight_estimates[i, j]), .01)
    }
  }

  for(i in rownames(csem_sat$Estimates$Loading_estimates)){
    for(j in colnames(csem_sat$Estimates$Loading_estimates)){
      if(csem_sat$Estimates$Loading_estimates[i,j] == 0)
        next
      testthat::expect_lt(abs(PLS_result$loadings[[i]][j] - csem_sat$Estimates$Loading_estimates[i, j]), .01)
    }
  }

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(PLS_result$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(PLS_result$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test reliability and debiased correlation - Mixed", {
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
    as_reflective = c("EXPE", "IMAG", "LOY"),
    data = data_set)

  reliability_basic <- get_reliability(PLS_result)

  library(cSEM)
  model <- "
# Structural model
QUAL ~ EXPE
EXPE ~ IMAG
SAT ~ IMAG + EXPE + QUAL + VAL
LOY ~ IMAG + SAT
VAL ~ EXPE + QUAL
# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY =~ loy1 + loy2 + loy3 + loy4
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT <~ sat1 + sat2 + sat3 + sat4
VAL <~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
                   .disattenuate = TRUE)

  for(i in names(csem_sat$Estimates$Reliabilities))
    testthat::expect_lt(abs(csem_sat$Estimates$Reliabilities[i] - reliability_basic[i]), 1e-4)

  testthat::expect_true(all(abs(debias_reflective_correlation(PLS_result)[rownames(csem_sat$Estimates$Construct_VCV), colnames(csem_sat$Estimates$Construct_VCV)] -
                                  csem_sat$Estimates$Construct_VCV) < 1e-3))
})
