#' PLS
#'
#' Basic function estimating a PLS-SEM with formative measurement model.
#' The algorithm is adapted from https://github.com/gastonstat/plspm.
#'
#' @param measurement alist with formulas defining the measurements
#' @param structure alist with formulas defining structural model
#' @param data data set (data.frame)
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
#' library(cSEM)
#' library(basicPLS)
#' data_set <- 10*cSEM::threecommonfactors + 2
#' PLS_result <- PLS(measurement = alist(eta1 ~ y11 + y12 + y13,
#'                                       eta2 ~ y21 + y22 + y23,
#'                                       eta3 ~ y31 + y32 + y33),
#'                   structure = alist(eta2 ~ eta1,
#'                                     eta3 ~ eta1 + eta2),
#'                   data = data_set,
#'                   path_estimation = "centroid")
#'
#' # same thing with cSEM
#' model <- "
#' # Structural model
#' eta2 ~ eta1
#' eta3 ~ eta1 + eta2
#'
#' # measurement model
#' eta1 <~ y11 + y12 + y13
#' eta2 <~ y21 + y22 + y23
#' eta3 <~ y31 + y32 + y33
#' "
#' fit_csem <- csem(data_set,
#'                  model,
#'                  .PLS_modes = "modeB",
#'                  .PLS_weight_scheme_inner = "centroid")
#' # Let's compare the results:
#' fit_csem$Estimates$Weight_estimates
#' PLS_result$weights
#'
#' fit_csem$Estimates$Path_estimates
#' PLS_result$effects
#'
#' head(fit_csem$Estimates$Construct_scores -
#'        PLS_result$components)
#'
#' assess(fit_csem, "r2")
#' get_r2(PLS_result)
#'
#' # We can fit the same model with the regression scheme:
#' PLS_result <- PLS(measurement = alist(eta1 ~ y11 + y12 + y13,
#'                                       eta2 ~ y21 + y22 + y23,
#'                                       eta3 ~ y31 + y32 + y33),
#'                   structure = alist(eta2 ~ eta1,
#'                                     eta3 ~ eta1 + eta2),
#'                   data = data_set,
#'                   path_estimation = "regression")
#'
#' # same thing with cSEM
#' fit_csem <- csem(data_set,
#'                  model,
#'                  .PLS_modes = "modeB",
#'                  .PLS_weight_scheme_inner = "path")
#' # Let's compare the results:
#' fit_csem$Estimates$Weight_estimates
#' PLS_result$weights
#'
#' fit_csem$Estimates$Path_estimates
#' PLS_result$effects
#'
#' head(fit_csem$Estimates$Construct_scores -
#'        PLS_result$components)
#'
#' assess(fit_csem, "r2")
#' get_r2(PLS_result)
PLS <- function(measurement,
                structure,
                data,
                sample_weights = NULL,
                path_estimation = "regression",
                max_iterations = 1000,
                convergence = 1e-5){
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

  # Next, we extract all observed variables. Those are on the right hand side of the
  # measurement models. The result will be a list and we will assign the
  # names of the constructs to the list elements.
  observables <- sapply(measurement, function(x) all.vars(x)[-1], simplify = FALSE)
  names(observables) <- components
  # Check if all observables are in the data set:
  if(any(!unlist(observables) %in% colnames(data)))
    stop("The following observables were not found in the data set: ",
         paste0(unlist(observables)[!unlist(observables) %in% colnames(data)], collapse = ", "), ".")


  data_std <- corpcor::wt.scale(data[, unique(unlist(observables))],
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
    component_values <- sapply(weights,
                               function(x) data_std[, names(x), drop = FALSE] %*% matrix(x, ncol = 1))
    # Next, we scale the components. Here, we use the weights to compute the
    # means and covariances.
    component_values <- corpcor::wt.scale(component_values,
                                          w = sample_weights)

    #### Inner Model ####
    # We compute the correlations of the components.
    component_correlation <- stats::cov.wt(component_values,
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
        adjacency_matrix[colnames(component_values), colnames(component_values)]
    }else if(path_estimation == "regression"){
      # For the path estimation, we first set all correlations of unrelated
      # variables to zero
      effects_matrix <- component_correlation *
        adjacency_matrix[colnames(component_values), colnames(component_values)]
      # Now, we replace all correlations with regression weights for endogenous
      # variables:
      effects <- sapply(structure,
                        function(x) regression_coef(formula = x,
                                                    data = as.data.frame(component_values),
                                                    wt = sample_weights)[-1],
                        simplify = FALSE)
      names(effects) <- sapply(structure, function(x) all.vars(x)[1])
      for(effect in names(effects)){
        # Note that we have to put the effects into the columns as we are using
        # the matrix multiplication component_values %*% effects_matrix below.
        effects_matrix[names(effects[[effect]]), effect] <- effects[[effect]]
      }
    }else{
      stop("Unknown path_estimation selected. Only centroid and regression are supported.")
    }
    # With this effects matrix, we update the components:
    component_values <- component_values %*% effects_matrix

    # Finally, we update the weights by predicting the component values with the
    # observed data using linear regressions.
    weights_upd <- sapply(measurement,
                          function(x) regression_coef(formula = x,
                                                      data = as.data.frame(cbind(data_std, component_values)),
                                                      wt = sample_weights)[-1],
                          simplify = FALSE)
    names(weights_upd) <- names(weights)

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
  component_values <- sapply(weights,
                             function(x) data_std[, names(x), drop = FALSE] %*% matrix(x, ncol = 1))
  # Scale components:
  component_values <- corpcor::wt.scale(component_values,
                                        w = sample_weights)

  # Compute the final weights
  weights_upd <- sapply(measurement,
                        function(x) regression_coef(x,
                                                    data = as.data.frame(cbind(data_std, component_values)),
                                                    wt = sample_weights)[-1],
                        simplify = FALSE)
  names(weights_upd) <- names(weights)
  weights <- weights_upd

  # unstandardized weights
  weights_unstandardized <-  sapply(measurement,
                                    function(x) regression_coef(x,
                                                                data = as.data.frame(cbind(data, component_values)),
                                                                wt = sample_weights)[-1],
                                    simplify = FALSE)
  names(weights_unstandardized) <- names(weights)

  # Compute effects as linear regressions between the components:
  effects <- sapply(structure,
                    function(x) regression_coef(x,
                                                data = as.data.frame(component_values),
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
      next
    }
    loadings[[component]] <- sapply(measurements,
                                    function(x) coef(lm(as.formula(paste0(x, " ~ ", component)),
                                                        data = as.data.frame(cbind(data_std, component_values)),
                                                        weights = sample_weights))[2],
                                    simplify = TRUE)
    names(loadings[[component]]) <- measurements
  }

  # Return results:
  result <- list(components = component_values,
                 weights = weights,
                 loadings = loadings,
                 weights_unstandardized = weights_unstandardized,
                 effects = effects,
                 input = input)
  class(result) <- "PLS_SEM"
  return(result)
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
    par_names <- c(par_names, paste0(comp, " <~ ", names(object$weights[[comp]])))
    if(include_loadings){
      par_values <- c(par_values, object$loadings[[comp]])
      par_names <- c(par_names, paste0(comp, " =~ ", names(object$loadings[[comp]])))
    }
  }
  for(comp in names(object$effects)){
    par_values <- c(par_values, object$effects[[comp]])
    par_names <- c(par_names, paste0(comp, " ~ ", names(object$effects[[comp]])))
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

#' confidence_intervals
#'
#' Compute confidence intervals for the parameters using bootstrap samples
#' @param PLS_result results from PLS
#' @param include_loadings should loadings also be returned?
#' @param alpha significance level
#' @param R number of repetitions
#' @returns boostrap results
#' @importFrom boot boot
#' @importFrom stats sd
#' @importFrom stats quantile
#' @export
confidence_intervals <- function(PLS_result,
                                 include_loadings = FALSE,
                                 alpha = .05,
                                 R = 1000){

  bootstrap_pls <- function(dat,
                            indices,
                            measurement,
                            structure,
                            sample_weights,
                            path_estimation,
                            max_iterations,
                            convergence,
                            include_loadings){
    fit_results <- suppressMessages(basicPLS::PLS(measurement = measurement,
                                                  structure = structure,
                                                  data = dat[indices,, drop = FALSE],
                                                  sample_weights = sample_weights[indices],
                                                  path_estimation = path_estimation,
                                                  max_iterations = max_iterations,
                                                  convergence = convergence))
    return(coef(fit_results, include_loadings = include_loadings))
  }
  # bootstrapping
  results <- boot::boot(PLS_result$input$data,
                        statistic = bootstrap_pls,
                        R = R,
                        measurement = PLS_result$input$measurement,
                        structure = PLS_result$input$structure,
                        sample_weights = PLS_result$input$sample_weights,
                        path_estimation = PLS_result$input$path_estimation,
                        max_iterations = PLS_result$input$max_iterations,
                        convergence = PLS_result$input$convergence,
                        include_loadings = include_loadings)
  colnames(results$t) <- names(results$t0)
  lower_ci <- apply(results$t, 2, quantile, alpha)
  upper_ci <- apply(results$t, 2, quantile, 1-alpha)

  return(list(confidence_intervals = data.frame(Parameter = names(results$t0),
                                                Estimate = unname(results$t0),
                                                type = ifelse(grepl("=~", names(results$t0)),
                                                              "loading",
                                                              ifelse(grepl("<~", names(results$t0)),
                                                                     "weight",
                                                                     "regression")),
                                                lower_ci = unname(lower_ci),
                                                upper_ci = unname(upper_ci)),
              full_results = results))
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


