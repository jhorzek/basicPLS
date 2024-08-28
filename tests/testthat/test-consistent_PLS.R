test_that("Test consistent PLS", {
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
  csem_sat <- csem(data_set,
                   model,
                   .disattenuate = TRUE)

  for(i in names(csem_sat$Estimates$Reliabilities))
    testthat::expect_lt(abs(csem_sat$Estimates$Reliabilities[i] - reliability_basic[i]), 1e-4)

  testthat::expect_true(all(abs(PLS_result$debiased_correlation[rownames(csem_sat$Estimates$Construct_VCV), colnames(csem_sat$Estimates$Construct_VCV)] -
                                  csem_sat$Estimates$Construct_VCV) < 1e-3))

  testthat::expect_true(all(abs(PLS_result$R2[names(csem_sat$Estimates$R2)] -
                                  csem_sat$Estimates$R2) < 1e-3))

  testthat::expect_true(all(abs(PLS_result$R2adj[names(csem_sat$Estimates$R2adj)] -
                                  csem_sat$Estimates$R2adj) < 1e-3))

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

test_that("Test consistent PLS - Mixed", {
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
    as_reflective = c("EXPE", "IMAG", "LOY"),
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
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT <~ sat1 + sat2 + sat3 + sat4
VAL <~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(data_set,
                   model,
                   .disattenuate = TRUE)

  testthat::expect_true(all(abs(PLS_result$R2[names(csem_sat$Estimates$R2)] -
                                  csem_sat$Estimates$R2) < 1e-3))

  testthat::expect_true(all(abs(PLS_result$R2adj[names(csem_sat$Estimates$R2adj)] -
                                  csem_sat$Estimates$R2adj) < 1e-3))

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


test_that("Test consistent bootstrap", {
  # NOTE: We are using a reduced model here, because the full model has a very
  #       serious issue with multicollinearity, which makes the bootstrap results
  #       extremely variable across runs.

  rm(list = ls())
  set.seed(3453)
  library(cSEM)
  library(basicPLS)
  data_set <- basicPLS::satisfaction

  model <- "
# Structural model
EXPE ~ IMAG
SAT ~ IMAG + EXPE
LOY ~ IMAG + SAT
# Measurement model
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY =~ loy1 + loy2 + loy3 + loy4
SAT =~ sat1 + sat2 + sat3 + sat4
"
  ## Estimate model with csem
  csem_sat <- csem(data_set,
                   model,
                   .resample_method = "bootstrap",
                   .R = 100,
                   .disattenuate = TRUE)
  summarized <- summarize(csem_sat)

  ## Estimate model with basicPLS
  basic_sat <- PLS(
    measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
                        IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
                        LOY ~ loy1 + loy2 + loy3 + loy4,
                        SAT ~ sat1 + sat2 + sat3 + sat4),
    structure = alist(EXPE ~ IMAG,
                      SAT ~ IMAG + EXPE,
                      LOY ~ IMAG + SAT),
    as_reflective = c("EXPE", "IMAG", "LOY", "SAT"),
    data = data_set)

  testthat::expect_true(all(abs(basic_sat$R2[names(csem_sat$Estimates$R2)] -
                                  csem_sat$Estimates$R2) < 1e-3))

  testthat::expect_true(all(abs(basic_sat$R2adj[names(csem_sat$Estimates$R2adj)] -
                                  csem_sat$Estimates$R2adj) < 1e-3))

  cis <- basicPLS::confidence_intervals(PLS_result = basic_sat,
                                        R = 100)

  for(i in summarized$Estimates$Path_estimates$Name){
    # the naming differs between basicPLS and cSEM
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%L`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$effects$lower_ci[cis$confidence_intervals$effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%U`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$effects$upper_ci[cis$confidence_intervals$effects$Parameter == name_basic]), .2)
  }
  for(i in summarized$Estimates$Weight_estimates$Name){
    name_basic <- gsub("<~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$weights$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%L`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$weights$lower_ci[cis$confidence_intervals$weights$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%U`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$weights$upper_ci[cis$confidence_intervals$weights$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Loading_estimates$Name){
    name_basic <- gsub("=~", "->", i)
    if(!name_basic %in% cis$confidence_intervals$loadings$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%L`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$loadings$lower_ci[cis$confidence_intervals$loadings$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%U`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$loadings$upper_ci[cis$confidence_intervals$loadings$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Effect_estimates$Total_effect$Name){
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$total_effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Total_effect$`CI_percentile.95%L`[summarized$Estimates$Effect_estimates$Total_effect$Name == i] -
                              cis$confidence_intervals$total_effects$lower_ci[cis$confidence_intervals$total_effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Total_effect$`CI_percentile.95%U`[summarized$Estimates$Effect_estimates$Total_effect$Name == i] -
                              cis$confidence_intervals$total_effects$upper_ci[cis$confidence_intervals$total_effects$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Effect_estimates$Indirect_effect$Name){
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$indirect_effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Indirect_effect$`CI_percentile.95%L`[summarized$Estimates$Effect_estimates$Indirect_effect$Name == i] -
                              cis$confidence_intervals$indirect_effects$lower_ci[cis$confidence_intervals$indirect_effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Indirect_effect$`CI_percentile.95%U`[summarized$Estimates$Effect_estimates$Indirect_effect$Name == i] -
                              cis$confidence_intervals$indirect_effects$upper_ci[cis$confidence_intervals$indirect_effects$Parameter == name_basic]), .2)
  }

  # Testing multicore bootstrapping
  cis <- basicPLS::confidence_intervals(PLS_result = basic_sat,
                                        R = 100,
                                        parallel = "snow")

  for(i in summarized$Estimates$Path_estimates$Name){
    # the naming differs between basicPLS and cSEM
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%L`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$effects$lower_ci[cis$confidence_intervals$effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Path_estimates$`CI_percentile.95%U`[summarized$Estimates$Path_estimates$Name == i] -
                              cis$confidence_intervals$effects$upper_ci[cis$confidence_intervals$effects$Parameter == name_basic]), .2)
  }
  for(i in summarized$Estimates$Weight_estimates$Name){
    name_basic <- gsub("<~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$weights$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%L`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$weights$lower_ci[cis$confidence_intervals$weights$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Weight_estimates$`CI_percentile.95%U`[summarized$Estimates$Weight_estimates$Name == i] -
                              cis$confidence_intervals$weights$upper_ci[cis$confidence_intervals$weights$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Loading_estimates$Name){
    name_basic <- gsub("=~", "->", i)
    if(!name_basic %in% cis$confidence_intervals$loadings$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%L`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$loadings$lower_ci[cis$confidence_intervals$loadings$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Loading_estimates$`CI_percentile.95%U`[summarized$Estimates$Loading_estimates$Name == i] -
                              cis$confidence_intervals$loadings$upper_ci[cis$confidence_intervals$loadings$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Effect_estimates$Total_effect$Name){
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$total_effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Total_effect$`CI_percentile.95%L`[summarized$Estimates$Effect_estimates$Total_effect$Name == i] -
                              cis$confidence_intervals$total_effects$lower_ci[cis$confidence_intervals$total_effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Total_effect$`CI_percentile.95%U`[summarized$Estimates$Effect_estimates$Total_effect$Name == i] -
                              cis$confidence_intervals$total_effects$upper_ci[cis$confidence_intervals$total_effects$Parameter == name_basic]), .2)
  }

  for(i in summarized$Estimates$Effect_estimates$Indirect_effect$Name){
    name_basic <- gsub("~", "<-", i)
    if(!name_basic %in% cis$confidence_intervals$indirect_effects$Parameter)
      stop("Could not find parameter!")
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Indirect_effect$`CI_percentile.95%L`[summarized$Estimates$Effect_estimates$Indirect_effect$Name == i] -
                              cis$confidence_intervals$indirect_effects$lower_ci[cis$confidence_intervals$indirect_effects$Parameter == name_basic]), .2)
    testthat::expect_lt(abs(summarized$Estimates$Effect_estimates$Indirect_effect$`CI_percentile.95%U`[summarized$Estimates$Effect_estimates$Indirect_effect$Name == i] -
                              cis$confidence_intervals$indirect_effects$upper_ci[cis$confidence_intervals$indirect_effects$Parameter == name_basic]), .2)
  }
})

