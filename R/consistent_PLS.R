#' get_reliability
#'
#' Compute the reliabilities of formative constructs in a PLS model.
#' See Dijkstra, T. K., & Henseler, J. (2015). Consistent partial least squares path modeling. MIS quarterly, 39(2), 297-316.
#' @param PLS_result model estimated with PLS()
#' @returns vector with reliabilities. For composites, the reliability is 1.0
#' @importFrom stats cov.wt
#' @export
#' @examples
#' library(plsR)
#' data_set <- plsR::satisfaction
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

  reliability <- rep(1.0, length(weights))
  names(reliability) <- names(weights)

  imputed_data <- PLS_result$input$imputation_function(data = PLS_result$input$data[, unlist(sapply(weights, names)), drop = FALSE],
                                                       weights = PLS_result$input$sample_weights)

  for(i in PLS_result$input$as_reflective){
    weights_vector <- matrix(weights[[i]], ncol = 1)

    # we use the correlation because we are using standardized data in PLS
    cov_mat <- stats::cov.wt(imputed_data[, names(weights[[i]]), drop = FALSE],
                             wt = PLS_result$input$sample_weights,
                             cor = TRUE,
                             method = "ML")$cor
    # See Eq 3 in: Dijkstra, T. K., & Henseler, J. (2015). Consistent partial least squares path modeling. MIS quarterly, 39(2), 297-316.
    reliability[i] <- ((t(weights_vector)%*%weights_vector)^2) *
      (
        (t(weights_vector) %*% (cov_mat - diag(diag(cov_mat))) %*% weights_vector) /
          (t(weights_vector) %*% (weights_vector %*% t(weights_vector) - diag(diag(weights_vector %*% t(weights_vector))))%*% weights_vector)
      )
  }

  return(reliability)
}

#' debias_correlation
#'
#' Debias the correlations of reflective constructs
#' See Dijkstra, T. K., & Henseler, J. (2015). Consistent partial least squares path modeling. MIS quarterly, 39(2), 297-316.
#' @param PLS_result model estimated with PLS()
#' @returns matrix with debiased correlations for reflective constructs
#' @importFrom stats cov.wt
debias_correlation <- function(PLS_result){

  reliability <- get_reliability(PLS_result)

  biased_correlation <- stats::cov.wt(PLS_result$composites,
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
  diag(debiased_correlation) <- 1.0

  return(debiased_correlation)
}


#' debias_PLS
#'
#' Debiases the regressions and loadings
#' in a PLS model with reflective constructs.
#' See Dijkstra, T. K., & Henseler, J. (2015). Consistent partial least squares path modeling. MIS quarterly, 39(2), 297-316.
#' @param PLS_result a model fitted with PLS()
#' @returns the model with updated loadings and regressions
#' @keywords internal
debias_PLS <- function(PLS_result){
  PLS_result$reliability <- get_reliability(PLS_result)
  PLS_result$debiased_correlation <- debias_correlation(PLS_result)

  # Debias effects
  # See Eq 6 in Dijkstra, T. K., & Henseler, J. (2015). Consistent partial
  # least squares path modeling. MIS quarterly, 39(2), 297-316.
  R2    <- rep(NA, length(PLS_result$effects))
  R2adj <- rep(NA, length(PLS_result$effects))
  names(R2) <- names(PLS_result$effects)
  names(R2adj) <- names(PLS_result$effects)
  for(i in names(PLS_result$effects)){
    PLS_result$effects[[i]] <- regression_from_correlation(R = PLS_result$debiased_correlation,
                                                           dependent = i,
                                                           predictors = names(PLS_result$effects[[i]]),
                                                           N = nrow(PLS_result$composites))
    R2[i] <- attr(PLS_result$effects[[i]], "R2")
    R2adj[i] <- attr(PLS_result$effects[[i]], "R2adj")
  }

  PLS_result$R2 <- R2
  PLS_result$R2adj <- R2adj

  indirect_total <- compute_summarized_effects(PLS_result)
  # we already have the direct effects; we only need the indirect and total effects
  PLS_result$total_effects <- indirect_total$total_effects
  PLS_result$indirect_effects <- indirect_total$indirect_effects

  # Debias loadings
  # See Eq 5 in Dijkstra, T. K., & Henseler, J. (2015). Consistent partial
  # least squares path modeling. MIS quarterly, 39(2), 297-316.
  for(i in PLS_result$input$as_reflective){
    PLS_result$loadings[[i]] <- PLS_result$weights[[i]] *
      (sqrt(PLS_result$reliability[i])/(t(PLS_result$weights[[i]]) %*% PLS_result$weights[[i]])[1,1])
  }

  return(PLS_result)
}
