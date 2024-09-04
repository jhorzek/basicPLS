#' mean_impute
#'
#' Basic mean imputation for missing data
#' @param data data set with missings
#' @param weights vector with weights for each person in the data set
#' @returns data set with imputed missings
#' @export
#' @examples
#' library(plsR)
#' satisfaction_with_missings <- plsR::satisfaction
#' missings <- matrix(sample(c(TRUE, FALSE),
#'                    nrow(plsR::satisfaction)*ncol(plsR::satisfaction),
#'                    prob = c(.1, .9),
#'                    replace = TRUE),
#'                    nrow = nrow(plsR::satisfaction),
#'                    ncol = ncol(plsR::satisfaction))
#' satisfaction_with_missings[missings] <- NA
#' mean_impute(data = satisfaction_with_missings,
#'             weights = rep(1, nrow(satisfaction_with_missings)))
mean_impute <- function(data, weights){
  for(i in 1:ncol(data)){
    data[[i]] <- ifelse(is.na(data[[i]]),
                        sum(weights * data[[i]], na.rm = TRUE) / sum(weights[!is.na(data[[i]])]),
                        data[[i]])
  }
  return(data)
}

#' fail_on_NA
#'
#' Fails if there is any missing data
#' @param data data set with missings
#' @param weights vector with weights for each person in the data set
#' @returns data in case of no error
#' @export
#' @examples
#' library(plsR)
#' satisfaction_with_missings <- plsR::satisfaction
#' missings <- matrix(sample(c(TRUE, FALSE),
#'                    nrow(plsR::satisfaction)*ncol(plsR::satisfaction),
#'                    prob = c(.1, .9),
#'                    replace = TRUE),
#'                    nrow = nrow(plsR::satisfaction),
#'                    ncol = ncol(plsR::satisfaction))
#' satisfaction_with_missings[missings] <- NA
#' try(fail_on_NA(data = satisfaction_with_missings,
#'                weights = rep(1, nrow(satisfaction_with_missings))))
fail_on_NA <- function(data, weights){
  if(anyNA(data))
    stop("Data has missings. Specify an imputation method to address the missingness (e.g., imputation_function = mean_impute).")

  return(data)
}
