#' PLS
#'
#' Basic function estimating a PLS-SEM with formative measurement model.
#' The algorithm is adapted from https://github.com/gastonstat/plspm.
#'
#' @param measurement alist with formulas defining the measurements
#' @param structure alist with formulas defining structural model
#' @param data data set (data.frame)
#' @param as_reflective vector with names of composites that should be treated
#' as reflective. Internally, basicPLS will use mode A estimation for these
#' composites, and mode B estimation for all other composites.
#' @param imputation_function a function that imputes missing data. The function
#' will be given the raw data set and the weights.
#' @param sample_weights data weights for weighted PLS-SEM
#' @param path_estimation how should the inner paths be estimated? Available
#' are "centroid" and "regression" based estimation.
#' @param max_iterations maximal number of iterations for the estimation
#' algorithm
#' @param convergence convergence criterion. If the maximal change in the
#' weights falls below this value, the estimation is stopped.
#' @return list with weights, effects, and components
#' @importFrom stats runif
#' @importFrom stats cor
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats as.formula
#' @importFrom corpcor wt.scale
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
#'   data = data_set)
#' PLS_result
#'
#' # R squared:
#' get_r2(PLS_result)
#'
#' # Use confidence_intervals to bootstrap confidence intervals for all parameters:
#' ci <- confidence_intervals(PLS_result,
#'                            # Increase for actual use:
#'                            R = 10)
PLS <- function(measurement,
                structure,
                data,
                as_reflective = NULL,
                imputation_function = fail_on_NA,
                sample_weights = NULL,
                path_estimation = "regression",
                max_iterations = 1000,
                convergence = 1e-5){

  check_formulas(measurement = measurement,
                 structure = structure)

  # create weights
  if(!is.null(sample_weights)){
    if(length(sample_weights) != nrow(data))
      stop("There must be one weight for each row in the data set.")
  }else{
    sample_weights <- rep(1, nrow(data))
  }

  # save input; this will make computing the standard errors easier
  input <- list(
    measurement = measurement,
    structure = structure,
    data = data,
    imputation_function = imputation_function,
    path_estimation = path_estimation,
    sample_weights = sample_weights,
    max_iterations = max_iterations,
    convergence = convergence
  )

  #### Preparation ####
  # We will first extract all components from the measurement structure. To this
  # end, we leverage that each component must be defined as the first element in
  # the measurement models:
  components <- sapply(measurement, function(x) all.vars(x)[1], simplify = TRUE)
  names(measurement) <- components

  # Next, we extract all observed variables. Those are on the right hand side of the
  # measurement models. The result will be a list and we will assign the
  # names of the constructs to the list elements.
  observables <- sapply(measurement, function(x) all.vars(x)[-1], simplify = FALSE)
  names(observables) <- components

  # Check if all observables are in the data set:
  if(any(!unlist(observables) %in% colnames(data)))
    stop("The following observables were not found in the data set: ",
         paste0(unlist(observables)[!unlist(observables) %in% colnames(data)], collapse = ", "), ".")
  # check if requested latents are in the components
  for(as_refl in as_reflective){
    if(!as_refl %in% components)
      stop(paste0("Could not find ", as_refl, " in measurements."))
  }

  # impute and scale:
  data_std <- corpcor::wt.scale(imputation_function(data[, unique(unlist(observables))],
                                                    sample_weights),
                                w = sample_weights)

  # Now, we construct a matrix indicating which of the components has an effect on
  # which other component. To this end, we will have to check the structure.
  # Extract the names of all variables in the structural part:
  structures <- sapply(structure, all.vars, simplify = FALSE)
  if(any(!unlist(structures) %in% components))
    stop("Could not find a measurement model for ",
         unlist(structures)[!unlist(structures) %in% components], ".")
  # The structure matrix will have ones for each directed effect.
  structure_matrix <- matrix(0,
                             nrow = length(components),
                             ncol = length(components),
                             dimnames = list(components, components))
  for(component in components){
    for(str in structures){
      # check if the component is on the left hand side of a structure
      if(component == str[1]){
        # then all right hand side variables are predictors of the component
        structure_matrix[component, str[-1]] <- 1
      }
    }
  }
  # The adjacency matrix below is similar to the structure matrix, but it indicates
  # with a 1 if two variables are related or not, irrespective of the direction of the
  # effect. We will need this for the inner model.
  adjacency_matrix <- structure_matrix + t(structure_matrix)

  # Finally, we initialize the weights. We could be using a full weights
  # matrix, but here we will keep things simpler (and slower) by having
  # component-specific weights.
  weights <- sapply(observables,
                    function(x) {
                      w <- rep(1, length(x))
                      names(w) <- x
                      return(w)
                    },
                    simplify = FALSE)

  # Now that we have everything prepared, we can get started with the
  # estimation process.
  for(i in seq_len(max_iterations)){

    #### Outer Model ####
    # First, we predict the components using the weights and the data.
    # As we don't estimate any parameters here, no weighting is needed.
    composite_values <- sapply(weights,
                               function(x) data_std[, names(x), drop = FALSE] %*% matrix(x, ncol = 1))
    # Next, we scale the components. Here, we use the weights to compute the
    # means and covariances.
    composite_values <- corpcor::wt.scale(composite_values,
                                          w = sample_weights)

    #### Inner Model ####
    # We compute the correlations of the components.
    component_correlation <- stats::cov.wt(composite_values,
                                           wt = sample_weights,
                                           cor = TRUE,
                                           method = "ML")$cor
    # Now, we take the component structure into account by creating a structure
    # matrix.
    if(path_estimation == "centroid"){
      # For the centroid approach, we simply replace all correlations between
      # unrelated constructs with 0. If constructs are related, we put a 1 in the
      # as effect if the correlation is positive and a -1 if the correlation is negative.
      effects_matrix <- sign(component_correlation) *
        adjacency_matrix[colnames(composite_values), colnames(composite_values)]
    }else if(path_estimation == "regression"){
      # For the path estimation, we first set all correlations of unrelated
      # variables to zero
      effects_matrix <- component_correlation *
        adjacency_matrix[colnames(composite_values), colnames(composite_values)]
      # Now, we replace all correlations with regression weights for endogenous
      # variables:
      effects <- sapply(structure,
                        function(x) regression_coef(formula = x,
                                                    data = as.data.frame(composite_values),
                                                    wt = sample_weights)[-1],
                        simplify = FALSE)
      names(effects) <- sapply(structure, function(x) all.vars(x)[1])
      for(effect in names(effects)){
        # Note that we have to put the effects into the columns as we are using
        # the matrix multiplication composite_values %*% effects_matrix below.
        effects_matrix[names(effects[[effect]]), effect] <- effects[[effect]]
      }
    }else{
      stop("Unknown path_estimation selected. Only centroid and regression are supported.")
    }
    # With this effects matrix, we update the components:
    composite_values <- composite_values %*% effects_matrix

    # Finally, we update the weights by predicting the component values with the
    # observed data using linear regressions.
    weights_upd <- list()
    for(comp in names(measurement)){
      if(comp %in% as_reflective){
        # Mode A: We predict the items using the composites
        weights_upd[[comp]] <- mode_A(measurement = measurement[[comp]],
                                      composite = comp,
                                      data_std = data_std,
                                      composite_values = composite_values,
                                      sample_weights = sample_weights)
      }else{
        # Mode B: We predict the composite using the items
        weights_upd[[comp]] <- mode_B(measurement = measurement[[comp]],
                                      composite = comp,
                                      data_std = data_std,
                                      composite_values = composite_values,
                                      sample_weights = sample_weights)
      }
    }

    #### Check Convergence ####
    # The algorithm converged, if the maximal change in the weights falls below
    # the pre-defined threshold:
    max_diff <- max(sapply(names(weights), function(x) max(abs(weights_upd[[x]] - weights[[x]]))))
    if(max_diff < convergence)
      break
    weights <- weights_upd
  }
  if(i == max_iterations){
    warning("Reached maximal number of iterations. Consider increasing max_iterations.")
  }else{
    message("The algorithm took ", i, " iterations to converge.")
  }

  # Finally, we compute the actual structural effects. To this end, we go one last
  # time through the steps above.
  # Predict components:
  composite_values <- sapply(weights,
                             function(x) data_std[, names(x), drop = FALSE] %*% matrix(x, ncol = 1))
  # Scale components:
  composite_values <- corpcor::wt.scale(composite_values,
                                        w = sample_weights)

  # Compute the final weights. In this round, we do not distinguish between
  # mode A and mode B. This ensures that we can always compute the scores as
  # data * weights
  weights <- list()
  for(comp in names(measurement)){
    # Mode B: We predict the composite using the items
    weights[[comp]] <- mode_B(measurement = measurement[[comp]],
                              composite = comp,
                              data_std = data_std,
                              composite_values = composite_values,
                              sample_weights = sample_weights)
  }

  # unstandardized weights
  weights_unstandardized <-  sapply(measurement,
                                    function(x) regression_coef(x,
                                                                data = as.data.frame(cbind(data, composite_values)),
                                                                wt = sample_weights)[-1],
                                    simplify = FALSE)
  names(weights_unstandardized) <- names(weights)

  # Compute effects as linear regressions between the components:
  effects <- sapply(structure,
                    function(x) regression_coef(x,
                                                data = as.data.frame(composite_values),
                                                wt = sample_weights)[-1],
                    simplify = FALSE)
  names(effects) <- sapply(structure, function(x) all.vars(x)[1])

  # Finally, we can also compute the pseudo-loadings of PLS-SEM.
  loadings <- list()
  for(m in measurement){
    # we have to reverse the direction of the measurements
    component <- all.vars(m)[1]
    measurements <-  all.vars(m)[-1]
    if(length(measurements) == 1){
      loadings[[component]] <- 1
      names(loadings[[component]]) <- measurements
      next
    }
    loadings[[component]] <- sapply(measurements,
                                    function(x) coef(lm(as.formula(paste0(x, " ~ ", component)),
                                                        data = as.data.frame(cbind(data_std, composite_values)),
                                                        weights = sample_weights))[2],
                                    simplify = TRUE)
    names(loadings[[component]]) <- measurements
  }

  # Combine results:
  result <- list(components = composite_values,
                 weights = weights,
                 loadings = loadings,
                 weights_unstandardized = weights_unstandardized,
                 effects = effects,
                 input = input)

  # compute total and indirect effects
  indirect_total <- compute_summarized_effects(result)
  # we already have the direct effects; we only need the indirect and total effects
  result$total_effects <- indirect_total$total_effects
  result$indirect_effects <- indirect_total$indirect_effects

  result$input <- input
  class(result) <- "PLS_SEM"
  return(result)
}

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

#' mode_A
#'
#' Mode A computes the weights as pseudo loadings. Say, composite C is a combination
#' of c_1, c_2, c_3 (e.g., C = .4*c_1 + .6*c_2 + .2*c_3). Then, mode A turns the
#' regression around and predicts each of the items c_1, c_2, c_3 with composite C:
#' c_1 = b_01 + b_11 * C + e, c_2 = b_02 + b_12 * C + e, c_3 = b_03 + b_13 * C + e.
#' The new weights are given by b_11, b_12, b_13.
#' @param measurement formula specifying the measurement (e.g., C ~ c_1 + c_2 + c_3)
#' @param composite the name of the composite which should be computed as (e.g., "C")
#' @param data_std standardized data set
#' @param composite_values current composite values
#' @param sample_weights weight for each observation
#' @returns a vector with updated weights
mode_A <- function(measurement,
                   composite,
                   data_std,
                   composite_values,
                   sample_weights){
  # Mode A: We predict the items using the composites
  current_items <- all.vars(measurement)[-1]
  wt <- sapply(current_items,
               function(x) coef(lm(as.formula(paste0(x, " ~ ", composite)),
                                   data = as.data.frame(cbind(data_std, composite_values)),
                                   weights = sample_weights))[-1],
               simplify = TRUE)
  names(wt) <- current_items

  return(wt)
}

#' mode_B
#'
#' Mode B computes the weights using a regression of the composite on the items.
#' Say, composite C is a combination of c_1, c_2, c_3 (e.g., C = .4*c_1 + .6*c_2 + .2*c_3).
#' Then, mode B estimates the new weights by estimating the weights of the regression
#' C = b_0 + b_1*c_1 + b_2*c_2 + b_3*c_3 + e.
#' The new weights are given by b_1, b_2, b_3.
#' @param measurement formula specifying the measurement (e.g., C ~ c_1 + c_2 + c_3)
#' @param composite the name of the composite which should be computed as (e.g., "C")
#' @param data_std standardized data set
#' @param composite_values current composite values
#' @param sample_weights weight for each observation
#' @returns a vector with updated weights
mode_B <- function(measurement,
                   composite,
                   data_std,
                   composite_values,
                   sample_weights){
  # Mode B: We predict the composite using the items
  wt <- regression_coef(formula = measurement,
                        data = as.data.frame(cbind(data_std, composite_values)),
                        wt = sample_weights)[-1]
  return(wt)
}

#' print.PLS_SEM
#'
#' Print results of PLS SEM
#' @param x PLS SEM model
#' @param ... not used
#' @returns nothing
#' @export
print.PLS_SEM <- function(x, ...){
  cat("\n#### PLS SEM Results ####\n")
  cat("Component Weights:\n")
  for(comp in names(x$weights)){
    cat(paste0(comp, " = ",
               paste0(paste0(format(round(x$weights[[comp]], 3), nsmall = 3), "*", names(x$weights[[comp]])),
                      collapse = " + ")),
        "\n")
  }
  cat("\nEffects:\n")
  for(comp in names(x$effects)){
    cat(paste0(comp, " ~ ",
               paste0(paste0(format(round(x$effects[[comp]], 3), nsmall = 3), "*", names(x$effects[[comp]])),
                      collapse = " + ")),
        "\n")
  }
}

#' coef.PLS_SEM
#'
#' Coefficient of PLS SEM
#' @param object PLS_SEM object
#' @param include_loadings if set to TRUE, loadings will also be shown. Note that
#' loadings are not proper parameters of the model as they are not optimized for.
#' @param ... not used
#' @returns vector with parameters
#' @importFrom stats coef
#' @export
coef.PLS_SEM <- function(object, include_loadings = FALSE, ...){
  par_values <- c()
  par_names  <- c()
  for(comp in names(object$weights)){
    par_values <- c(par_values, object$weights[[comp]])
    par_names <- c(par_names, paste0(comp, " <- ", names(object$weights[[comp]])))
    if(include_loadings){
      par_values <- c(par_values, object$loadings[[comp]])
      par_names <- c(par_names, paste0(comp, " -> ", names(object$loadings[[comp]])))
    }
  }
  for(comp in names(object$effects)){
    par_values <- c(par_values, object$effects[[comp]])
    par_names <- c(par_names, paste0(comp, " <- ", names(object$effects[[comp]])))
  }
  names(par_values) <- par_names
  return(par_values)
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
  return(coef(lm(formula = formula, data = data , weights = wt)))
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
                                      data = as.data.frame(PLS_result$components),
                                      weights = PLS_result$input$sample_weights))$r.squared,
               simplify = FALSE)
  names(r2) <- names(PLS_result$effects)
  return(r2)
}

#' confidence_intervals
#'
#' Compute confidence intervals for the parameters using bootstrap samples
#' @param PLS_result results from PLS
#' @param alpha significance level
#' @param R number of repetitions
#' @param ... additional arguments passed to boot. See ?boot::boot
#' @returns bootstrap results
#' @importFrom boot boot
#' @importFrom stats sd
#' @importFrom stats quantile
#' @import stringr
#' @export
confidence_intervals <- function(PLS_result,
                                 alpha = .05,
                                 R = 1000,
                                 ...){

  bootstrap_pls <- function(dat,
                            indices,
                            measurement,
                            structure,
                            imputation_function,
                            sample_weights,
                            path_estimation,
                            max_iterations,
                            convergence){
    fit_results <- suppressMessages(basicPLS::PLS(measurement = measurement,
                                                  structure = structure,
                                                  data = dat[indices,, drop = FALSE],
                                                  imputation_function = imputation_function,
                                                  sample_weights = sample_weights[indices],
                                                  path_estimation = path_estimation,
                                                  max_iterations = max_iterations,
                                                  convergence = convergence))
    # We have to do some trickery because boot expects that we only return
    # a single vector.
    return_vector <- c()
    for(p_type in c("effects", "weights", "loadings", "total_effects", "indirect_effects")){
      add_to_return <- flatten_effects(fit_results[[p_type]],
                                       separator = ifelse(p_type == "loadings", "->", "<-"))
      names(add_to_return) <- paste0(p_type, ":", names(add_to_return))
      return_vector <- c(return_vector, add_to_return)
    }

    return(return_vector)
  }

  # bootstrapping
  results <- boot::boot(PLS_result$input$data,
                        statistic = bootstrap_pls,
                        R = R,
                        measurement = PLS_result$input$measurement,
                        structure = PLS_result$input$structure,
                        imputation_function = PLS_result$input$imputation_function,
                        sample_weights = PLS_result$input$sample_weights,
                        path_estimation = PLS_result$input$path_estimation,
                        max_iterations = PLS_result$input$max_iterations,
                        convergence = PLS_result$input$convergence,
                        ...)
  colnames(results$t) <- names(results$t0)
  lower_ci <- apply(results$t, 2, quantile, alpha)
  upper_ci <- apply(results$t, 2, quantile, 1-alpha)

  confidence_intervals <- data.frame(Parameter = names(results$t0),
                                     Estimate = unname(results$t0),
                                     lower_ci = unname(lower_ci),
                                     upper_ci = unname(upper_ci))

  # the results combine weights, loadings, direct, indirect, and total effects.
  # We want to split those up again.
  parameter_class <- stringr::str_extract(string = confidence_intervals$Parameter,
                                          pattern = "^.*:") |>
    stringr::str_remove(pattern = ":")
  confidence_intervals$Parameter <- stringr::str_extract(string = confidence_intervals$Parameter,
                                                         pattern = ":.*$") |>
    stringr::str_remove(pattern = ":")

  confidence_intervals <- split(confidence_intervals, f = parameter_class)

  return(list("confidence_intervals" = confidence_intervals,
              "full_results" = results))
}

#' mean_impute
#'
#' Basic mean imputation for missing data
#' @param data data set with missings
#' @param weights vector with weights for each person in the data set
#' @returns data set with imputed missings
#' @export
#' @examples
#' library(basicPLS)
#' satisfaction_with_missings <- basicPLS::satisfaction
#' missings <- matrix(sample(c(TRUE, FALSE),
#'                    nrow(basicPLS::satisfaction)*ncol(basicPLS::satisfaction),
#'                    prob = c(.1, .9),
#'                    replace = TRUE),
#'                    nrow = nrow(basicPLS::satisfaction),
#'                    ncol = ncol(basicPLS::satisfaction))
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
#' library(basicPLS)
#' satisfaction_with_missings <- basicPLS::satisfaction
#' missings <- matrix(sample(c(TRUE, FALSE),
#'                    nrow(basicPLS::satisfaction)*ncol(basicPLS::satisfaction),
#'                    prob = c(.1, .9),
#'                    replace = TRUE),
#'                    nrow = nrow(basicPLS::satisfaction),
#'                    ncol = ncol(basicPLS::satisfaction))
#' satisfaction_with_missings[missings] <- NA
#' try(fail_on_NA(data = satisfaction_with_missings,
#'                weights = rep(1, nrow(satisfaction_with_missings))))
fail_on_NA <- function(data, weights){
  if(anyNA(data))
    stop("Data has missings. Specify an imputation method to address the missingness (e.g., imputation_function = mean_impute).")

  return(data)
}
