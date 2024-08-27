#' check_formulas
#'
#' R's formulas are much more flexible than the current framework of basicPLS
#' supports. We have to make sure that users are not trying to fit more complex
#' models than what basicPLS allows.
#' @param measurement alist with the measurement formulas
#' @param structure alist with the structural formulas
#' @param allowed pattern of supported formulas
#' @returns nothing. Throws error in case of unsupported formulas
#' @examples
#' basicPLS:::check_formulas(measurement = alist(C ~ x1 + x2,
#'                                               D ~ x3 + x4),
#'                           structure = alist(D ~ C))
#'
check_formulas <- function(measurement,
                           structure,
                           allowed = "^[a-zA-Z0-9_]+[ ]*~[a-zA-Z0-9_\\+ ]+$"){
  for(form in c(measurement, structure)){
    if(!grepl(pattern = allowed,
              x = paste0(format(form),
                         collapse = "")))
      stop(paste0("basicPLS currently only supports formulas with the pattern z ~ x + y. ",
                  "The following is not allowed: ",
                  paste0(format(form),
                         collapse = ""), "."))
  }
}

# We will have to run multiple regressions and extract the coefficients. The
# following makes this a bit easier to read:
#' regression_coef
#'
#' Computes regression coefficients and allows weighting. This function is necessary
#' to re-assign the environment of the formula and make sure that R finds the weights.
#' @param formula regression formula
#' @param data data set
#' @param wt sample weights vector
#' @returns regression coefficients
#' @keywords internal
regression_coef <- function(formula, data, wt = NULL){
  formula <- as.formula(formula)
  environment(formula) <- environment()
  return(coef(lm(formula = formula, data = data , weights = wt))[-1])
}

#' For debiasing, we need to compute regressions based on correlations. This
#' function provides that feature. See https://wviechtb.github.io/metafor/reference/matreg.html
#'
#' @param R correlation matrix
#' @param dependent name of the dependent variable
#' @param predictors vector with names of the predictors
#' @returns regression coefficients
#' @keywords internal
regression_from_correlation <- function(R,
                                        dependent,
                                        predictors){
  # see https://wviechtb.github.io/metafor/reference/matreg.html
  b <- solve(R[predictors, predictors]) %*% R[predictors, dependent]
  # make the output consistent with that of regression_coef
  reg <- c(b)
  names(reg) <- predictors
  return(reg)
}

#' compute_summarized_effects
#'
#' Computes direct, indirect, and total effects
#' for a PLS model
#' @param PLS_result PLS model fitted with PLS() (see ?PLS)
#' @returns list with direct, indirect, and total effects.
#' @keywords internal
compute_summarized_effects <- function(PLS_result){
  effects <- PLS_result$effects
  measurement <- PLS_result$input$measurement
  # The following is adapted from plspm. See
  # https://github.com/gastonstat/plspm/blob/master/R/get_effects.r
  composites <- sapply(measurement, function(x) all.vars(x)[1], simplify = TRUE)
  structure_matrix <- matrix(0,
                             nrow = length(composites),
                             ncol = length(composites),
                             dimnames = list(composites, composites))
  for(dep in composites){
    if(!dep %in% names(effects))
      next
    structure_matrix[dep, names(effects[[dep]])] <- effects[[dep]]
  }

  total_effects <- structure_matrix
  current_effects <- structure_matrix

  if(length(composites) > 2){

    for(i in 2:(length(composites)-1)){
      current_effects <- current_effects %*% structure_matrix
      total_effects <- total_effects + current_effects
    }
  }

  indirect_effects <- total_effects - structure_matrix

  # because everything else is organized as lists, we will do the same here:
  total_effects_list <- list()
  indirect_effects_list <- list()
  for(dep in rownames(total_effects)){
    if(!all(total_effects[dep,] == 0))
      total_effects_list[[dep]] <- total_effects[dep,][total_effects[dep,] != 0]
    if(!all(indirect_effects[dep,] == 0))
      indirect_effects_list[[dep]] <- indirect_effects[dep, ][indirect_effects[dep, ] != 0]
  }

  return(list(direct_effects = effects,
              total_effects = total_effects_list,
              indirect_effects = indirect_effects_list))
}

#' flatten_effects
#'
#' Given a list of direct or indirect effects, this function flattens the list
#' to a vector with relabeled parameters
#' @param effects summarized total or indirect effects
#' @param separator separator used when renaming the parameters
#' @returns vector with effects
#' @keywords internal
flatten_effects <- function(effects, separator = "<-"){
  flattened <- sapply(names(effects),
                      function(x) {
                        vals <- effects[[x]]
                        names(vals) <- paste0(x, " ", separator, " ", names(vals))
                        vals},
                      simplify = FALSE) |>
    unname() |>
    unlist()
  return(flattened)
}

#' get_r2
#'
#' Compute the R squared value for a PLS-SEM
#' @param PLS_result fitted PLS-SEM
#' @return list with R squared
#' @importFrom stats lm
#' @importFrom methods is
#' @export
get_r2 <- function(PLS_result){
  r2 <- sapply(PLS_result$input$structure,
               function(x) summary(lm(as.formula(x),
                                      data = as.data.frame(PLS_result$composites),
                                      weights = PLS_result$input$sample_weights))$r.squared,
               simplify = FALSE)
  names(r2) <- names(PLS_result$effects)
  return(r2)
}
