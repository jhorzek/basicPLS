test_that("Test Satisfaction - Formative Regression", {
  rm(list = ls())
  set.seed(3453)
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

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test Satisfaction - Formative Regression Single Items", {
  rm(list = ls())
  set.seed(3453)
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
EXPE <~ expe1
IMAG <~ imag1 + imag2 + imag3 + imag4 + imag5
LOY <~ loy1
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
      EXPE ~ expe1,
      IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
      LOY ~ loy1,
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

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test Satisfaction - Formative Centroid", {
  rm(list = ls())
  set.seed(3453)
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

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test Satisfaction - Bootstrap", {
  rm(list = ls())
  set.seed(3453)
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


test_that("Test Satisfaction - Reflective Regression", {
  rm(list = ls())
  set.seed(3453)
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
    as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
    consistent = FALSE,
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

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test Satisfaction - Reflective Weighted", {
  rm(list = ls())
  set.seed(3453)
  library(basicPLS)
  satisfaction_data <- basicPLS::satisfaction
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
    data = satisfaction_data,
    as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
    consistent = FALSE,
    path_estimation = "regression",
    sample_weights = satisfaction_data$sample_weights)

  # comparison to results from smartPLS
  load(system.file("testdata",
                   "expected_weighted_satisfaction_results.RData",
                   package="basicPLS"))
  expected <- expected_weighted_satisfaction_results
  rm(expected_weighted_satisfaction_results)
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

test_that("Test Satisfaction - Reflective Weighted Mean Imputation", {
  rm(list = ls())
  set.seed(3453)
  library(basicPLS)
  data("satisfaction_NA")
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
    data = satisfaction_NA,
    as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
    consistent = FALSE,
    path_estimation = "regression",
    imputation_function = mean_impute,
    sample_weights = satisfaction_NA$sample_weights)

  # comparison to results from smartPLS
  load(system.file("testdata",
                   "expected_weighted_satisfaction_NA_results.RData",
                   package="basicPLS"))
  expected <- expected_weighted_satisfaction_NA_results
  rm(expected_weighted_satisfaction_NA_results)
  for(i in expected$effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$effects$path[expected$effects$Parameter == i] -
                              basic_sat$effects[[dep]][pred]), .01)
  }

  for(i in expected$total_effects$Parameter){
    pred <- stringr::str_extract(i, "^[a-zA-Z0-9]+")
    dep  <- stringr::str_extract(i, "[a-zA-Z0-9]+$")
    testthat::expect_lt(abs(expected$total_effects$total_effect[expected$total_effects$Parameter == i] -
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

test_that("Test Satisfaction - Mixed Regression", {
  rm(list = ls())
  set.seed(3453)
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
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY =~ loy1 + loy2 + loy3 + loy4
QUAL <~ qual1 + qual2 + qual3 + qual4 + qual5
SAT =~ sat1 + sat2 + sat3 + sat4
VAL =~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  csem_sat <- csem(satisfaction,
                   model,
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
    as_reflective = c("IMAG", "LOY", "SAT", "VAL"),
    consistent = FALSE,
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

  summarized_csem <- summarize(csem_sat)
  for(i in summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$indirect_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Indirect_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Indirect_effect$Name == i]), .01)
  }

  for(i in summarized_csem$Estimates$Effect_estimates$Total_effect$Name){
    # split in dependent and predictor
    dep <- stringr::str_extract(i, "^[a-zA-Z]+")
    pred <- stringr::str_extract(i, "[a-zA-Z]+$")

    testthat::expect_lt(abs(basic_sat$total_effects[[dep]][pred] -
                              summarized_csem$Estimates$Effect_estimates$Total_effect$Estimate[summarized_csem$Estimates$Effect_estimates$Total_effect$Name == i]), .01)
  }
})

test_that("Test Satisfaction - Imputation", {
  rm(list = ls())
  set.seed(3453)
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
EXPE =~ expe1 + expe2 + expe3 + expe4 + expe5
IMAG =~ imag1 + imag2 + imag3 + imag4 + imag5
LOY =~ loy1 + loy2 + loy3 + loy4
QUAL =~ qual1 + qual2 + qual3 + qual4 + qual5
SAT =~ sat1 + sat2 + sat3 + sat4
VAL =~ val1 + val2 + val3 + val4
"
  ## Estimate model with csem
  satisfaction_with_missings <- satisfaction
  missings <- matrix(sample(c(TRUE, FALSE),
                            nrow(satisfaction_with_missings)*ncol(satisfaction_with_missings),
                            prob = c(.1, .9),
                            replace = TRUE),
                     nrow = nrow(satisfaction_with_missings),
                     ncol = ncol(satisfaction_with_missings))
  satisfaction_with_missings[missings] <- NA

  satisfaction_mean_imputed <- satisfaction_with_missings
  for(i in 1:ncol(satisfaction_mean_imputed)){
    satisfaction_mean_imputed[[i]] <- ifelse(is.na(satisfaction_mean_imputed[[i]]),
                                             mean(satisfaction_with_missings[[i]], na.rm = TRUE),
                                             satisfaction_mean_imputed[[i]])
  }

  csem_sat <- csem(satisfaction_mean_imputed,
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
    data = satisfaction_with_missings,
    imputation_function = basicPLS::mean_impute,
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

  # test bootstrapping
  ci <- confidence_intervals(basic_sat, R = 10)
})
