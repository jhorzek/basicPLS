#' PLS
#'
#' Basic function estimating a PLS-SEM with formative measurement model.
#' The algorithm is adapted from https://github.com/gastonstat/plspm.
#'
#' @param measurement alist with formulas defining the measurements
#' @param structure alist with formulas defining structural model
#' @param data data set (data.frame)
#' @param max_iterations maximal number of iterations for the estimation
#' algorithm
#' @param convergence convergence criterion. If the maximal change in the
#' weights falls below this value, the estimation is stopped.
#' @return list with weights, effects, and components
#' @importFrom stats runif
#' @importFrom stats cor
#' @importFrom stats coef
#' @importFrom stats lm
#' @export
#' @examples
#' library(cSEM)
#' data_set <- 10*cSEM::threecommonfactors + 2
#' PLS_result <- PLS(measurement = alist(eta1 ~ y11 + y12 + y13,
#'                                       eta2 ~ y21 + y22 + y23,
#'                                       eta3 ~ y31 + y32 + y33),
#'                   structure = alist(eta2 ~ eta1,
#'                                     eta3 ~ eta1 + eta2),
#'                   data = data_set)
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
#'                  .PLS_modes = "modeB")
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
                max_iterations = 1000,
                convergence = 1e-5){
  # save input; this will make computing the standard errors easier
  input <- list(
    measurement = measurement,
    structure = structure,
    data = data,
    max_iterations = max_iterations,
    convergence = convergence
  )

  # We will need the sample size to correct the standard deviations from R.
  N <- nrow(data)
  data_std <- scale(data)*sqrt((N-1)/N)

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
  if(any(!unlist(observables) %in% colnames(data_std)))
    stop("The following observables were not found in the data set: ",
         paste0(unlist(observables)[!unlist(observables) %in% colnames(data_std)], collapse = ", "), ".")

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
                      if(length(x) == 1){
                        w <- c(1)}
                      else{
                        w <- runif(length(x), min = .1, max = .5)
                      }
                      names(w) <- x
                      return(w)
                    },
                    simplify = FALSE)

  # Now that we have everything prepared, we can get started with the
  # estimation process.
  for(i in seq_len(max_iterations)){

    #### Outer Model ####
    # First, we predict the components using the weights and the data.
    component_values <- sapply(weights,
                               function(x) data_std[, names(x), drop = FALSE] %*% matrix(x, ncol = 1))
    # Next, we scale the components
    component_values <- scale(component_values)*sqrt((N-1)/N)

    #### Inner Model ####
    # We compute the correlations of the components.
    component_correlation <- cor(component_values)
    # Now, we take the component structure into account by creating a structure
    # matrix that has 1 if adjacent components are positively correlated, -1
    # if they are negatively correlated, and 0 otherwise.
    effects_matrix <- sign(component_correlation) *
      adjacency_matrix[colnames(component_values), colnames(component_values)]

    # With this effects matrix, we update the components:
    component_values <- component_values %*% effects_matrix

    # Finally, we update the weights by predicting the component values with the
    # observed data using linear regressions.
    weights_upd <- sapply(measurement,
                          function(x) coef(lm(x, data = as.data.frame(cbind(data_std, component_values))))[-1],
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
  component_values <- scale(component_values)*sqrt((N-1)/N)

  # Compute the final weights
  weights_upd <- sapply(measurement,
                        function(x) coef(lm(x, data = as.data.frame(cbind(data_std, component_values))))[-1],
                        simplify = FALSE)
  names(weights_upd) <- names(weights)
  weights <- weights_upd

  # unstandardized weights
  weights_unstandardized <-  sapply(measurement,
                                    function(x) coef(lm(x, data = as.data.frame(cbind(data, component_values))))[-1],
                                    simplify = FALSE)
  names(weights_unstandardized) <- names(weights)

  # Compute effects as linear regressions between the components:
  effects <- sapply(structure,
                    function(x) coef(lm(x, data = as.data.frame(component_values)))[-1],
                    simplify = FALSE)
  names(effects) <- sapply(structure, function(x) all.vars(x)[1])

  # Return results:
  result <- list(components = component_values,
                 weights = weights,
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
#' @param ... not used
#' @returns vector with parameters
#' @importFrom stats coef
#' @export
coef.PLS_SEM <- function(object, ...){
  par_values <- c()
  par_names  <- c()
  for(comp in names(object$weights)){
    par_values <- c(par_values, object$weights[[comp]])
    par_names <- c(par_names, paste0(comp, " <~ ", names(object$weights[[comp]])))
  }
  for(comp in names(object$effects)){
    par_values <- c(par_values, object$effects[[comp]])
    par_names <- c(par_names, paste0(comp, " ~ ", names(object$effects[[comp]])))
  }
  names(par_values) <- par_names
  return(par_values)
}

#' standard_errors
#'
#' Compute standard errors using bootstrap samples
#' @param PLS_result results from PLS
#' @param R number of repetitions
#' @returns boostrap results
#' @importFrom boot boot
#' @importFrom stats sd
#' @export
standard_errors <- function(PLS_result,
                            R = 1000){

  bootstrap_pls <- function(dat,
                            indices,
                            measurement,
                            structure,
                            max_iterations,
                            convergence){
    fit_results <- suppressMessages(basicPLS::PLS(measurement = measurement,
                                                  structure = structure,
                                                  data = dat[indices,, drop = FALSE],
                                                  max_iterations = max_iterations,
                                                  convergence = convergence))
    return(coef(fit_results))
  }
  # bootstrapping with 1000 replications
  results <- boot::boot(PLS_result$input$data,
                        statistic = bootstrap_pls,
                        R = R,
                        measurement = PLS_result$input$measurement,
                        structure = PLS_result$input$structure,
                        max_iterations = PLS_result$input$max_iterations,
                        convergence = PLS_result$input$convergence)
  colnames(results$t) <- names(results$t0)
  se <- apply(results$t, 2, sd)
  return(se)
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
  if(!is(PLS_result, "PLS_SEM"))
    stop("PLS_result must be of class PLS_SEM")
  r2 <- sapply(PLS_result$input$structure,
               function(x) summary(lm(x, data = as.data.frame(PLS_result$components)))$r.squared,
               simplify = FALSE)
  names(r2) <- names(PLS_result$effects)
  return(r2)
}


