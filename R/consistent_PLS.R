#' get_reliability
#'
#' Compute the reliabilities of formative constructs in a PLS model
#' @param PLS_result model estimated with PLS()
#' @returns vector with reliabilities. For composites, the reliability is 1.0
#' @importFrom stats cov.wt
#' @export
#' @examples
#' library(basicPLS)
#' data_set <- basicPLS::satisfaction
#'
#' # Both, measurement and structural model are specified using R's formulas:
#' PLS_result <- PLS(
#'   measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
#'                       IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
#'                       LOY ~ loy1 + loy2 + loy3 + loy4,
#'                       QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
#'                       SAT ~ sat1 + sat2 + sat3 + sat4,
#'                       VAL ~ val1 + val2 + val3 + val4),
#'   structure = alist(QUAL ~ EXPE,
#'                     EXPE ~ IMAG,
#'                     SAT ~ IMAG + EXPE + QUAL + VAL,
#'                     LOY ~ IMAG + SAT,
#'                     VAL ~ EXPE + QUAL),
#'   as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
#'   data = data_set)
#' PLS_result
#'
#' get_reliability(PLS_result)
get_reliability <- function(PLS_result){
  weights <- PLS_result$weights
  as_reflective <- PLS_result$input$as_reflective

  reliability <- rep(1.0, length(weights))
  names(reliability) <- names(weights)

  for(i in as_reflective){
    weights_vector <- matrix(weights[[i]], ncol = 1)

    # we use the correlation because we are using standardized data in PLS
    cov_mat <- stats::cov.wt(PLS_result$input$data[, names(weights[[i]])],
                             wt = PLS_result$input$sample_weights,
                             cor = TRUE,
                             method = "ML")$cor

    reliability[i] <- ((t(weights_vector)%*%weights_vector)^2) *
      (
        (t(weights_vector) %*% (cov_mat - diag(diag(cov_mat))) %*% weights_vector) /
          (t(weights_vector) %*% (weights_vector %*% t(weights_vector) - diag(diag(weights_vector %*% t(weights_vector))))%*% weights_vector)
      )
  }

  return(reliability)
}

#' debias_reflective_correlation
#'
#' Debias the correlations of reflective constructs
#' @param PLS_result model estimated with PLS()
#' @returns matrix with debiased correlations for reflective constructs
#' @importFrom stats cov.wt
#' @export
#' @examples
#' library(basicPLS)
#' data_set <- basicPLS::satisfaction
#'
#' # Both, measurement and structural model are specified using R's formulas:
#' PLS_result <- PLS(
#'   measurement = alist(EXPE ~ expe1 + expe2 + expe3 + expe4 + expe5,
#'                       IMAG ~ imag1 + imag2 + imag3 + imag4 + imag5,
#'                       LOY ~ loy1 + loy2 + loy3 + loy4,
#'                       QUAL ~ qual1 + qual2 + qual3 + qual4 + qual5,
#'                       SAT ~ sat1 + sat2 + sat3 + sat4,
#'                       VAL ~ val1 + val2 + val3 + val4),
#'   structure = alist(QUAL ~ EXPE,
#'                     EXPE ~ IMAG,
#'                     SAT ~ IMAG + EXPE + QUAL + VAL,
#'                     LOY ~ IMAG + SAT,
#'                     VAL ~ EXPE + QUAL),
#'   as_reflective = c("EXPE", "IMAG", "LOY", "QUAL", "SAT", "VAL"),
#'   data = data_set)
#' PLS_result
#'
#' debias_reflective_correlation(PLS_result)
debias_correlation <- function(PLS_result){

  indicators <- unlist(sapply(PLS_result$weights, names))

  reliability <- rep(1, length(indicators))
  names(reliability) <- indicators
  reliability <- c(reliability,
                   get_reliability(PLS_result))

  imputed_data <- PLS_result$input$imputation_function(data = PLS_result$input$data[,indicators],
                                                       weights = PLS_result$input$sample_weights)

  biased_correlation <- stats::cov.wt(cbind(PLS_result$composites, imputed_data),
                                      wt = PLS_result$input$sample_weights,
                                      cor = TRUE,
                                      method = "ML")$cor
  debiased_correlation <- biased_correlation
  for(i in rownames(biased_correlation)){
    for(j in colnames(biased_correlation)){
      if(i == j)
        next
      debiased_correlation[i, j] <- biased_correlation[i, j] / sqrt(reliability[i]*reliability[j])
    }
  }

  # make symmetric
  debiased_correlation <- .5*(debiased_correlation + t(debiased_correlation))

  return(debiased_correlation)
}


#' debias_PLS
#'
#' Debiases the regressions and loadings
#' in a PLS model with reflective constructs.
#' @param PLS_result a model fitted with PLS()
#' @returns the model with updated loadings and regressions
#' @keywords internal
debias_PLS <- function(PLS_result){
  PLS_result$reliability <- get_reliability(PLS_result)
  PLS_result$debiased_correlation <- debias_correlation(PLS_result)

  # debias effects
  for(i in names(PLS_result$effects)){
    PLS_result$effects[[i]] <- regression_from_correlation(R = PLS_result$debiased_correlation,
                                                           dependent = i,
                                                           predictors = names(PLS_result$effects[[i]]))
  }

  indirect_total <- compute_summarized_effects(PLS_result)
  # we already have the direct effects; we only need the indirect and total effects
  PLS_result$total_effects <- indirect_total$total_effects
  PLS_result$indirect_effects <- indirect_total$indirect_effects

  # debias loadings
  for(i in names(PLS_result$loadings)){
    loading_names <- names(PLS_result$loadings[[i]])
    PLS_result$loadings[[i]] <- sapply(loading_names,
                                       function(x)
                                         regression_from_correlation(R = PLS_result$debiased_correlation,
                                                                     dependent = x,
                                                                     predictors = i)
    )
    names(PLS_result$loadings[[i]]) <- loading_names
  }

  return(PLS_result)
}
