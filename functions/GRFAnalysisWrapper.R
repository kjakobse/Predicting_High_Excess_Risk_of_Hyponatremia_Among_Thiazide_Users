#' Wrapper around generalized linear forests
#'
#' Covariate data can contains both binary, discrete (with more
#' than two levels), and continuous columns. Discrete covariates with more than
#' two levels is prepared using one-hot encoding. The input data frame is 
#' prepared so it can be plugged into a generalized random forest.
#'
#' @param data A data frame or data frame like object (e.g. tibble) containing
#'   the covariates, exposure and outcome.
#' @param type string with the type of generalized random forest. 
#'   Can be regression, survival, causal, or causal_survival.
#' @param continuous <[dplyr::dplyr_tidy_select]> One unquoted expression
#'   specifying the continuous (and binary) covariates in data.
#' @param discrete <[dplyr::dplyr_tidy_select]> One unquoted expression
#'   specifying the discrete (more than two levels) covariates in data.
#' @param Y string with the name of the outcome.
#' @param W string with the name of the exposure.
#' @param D string with the name of the event indicator.
#' @param ... Additional arguments for the generalized random forest.

GRFAnalysisWrapper <- function (
    data,
    type = c("regression", "survival", "causal", "causal_survival"),
    continuous,
    discrete,
    Y = "Y",
    W = "W",
    D = "D",
    ...
) {
  type <- match.arg(type)
  if (!exists("DiscreteCovariatesToOneHot")) {
    stop ("The function 'DiscreteCovariatesToOneHot' must be available.")
  }
  X <- data |>
    select({{ continuous }}) |>
    bind_cols(
      data |>
        select({{ discrete }}) |>
        DiscreteCovariatesToOneHot()
    )
  if(type == "regression") {
    out <- regression_forest(
      X = X,
      Y = as.vector(data[[Y]]),
      ...
    )
  } else if (type == "survival") {
    out <- survival_forest(
      X = X,
      Y = as.vector(data[[Y]]),
      D = as.vector(data[[D]]),
      ...
    )
  } else if (type == "causal") {
    out <- causal_forest(
      X = X,
      Y = as.vector(data[[Y]]),
      W = as.vector(data[[W]]),
      ...
    )
  } else if (type == "causal_survival") {
    out <- causal_survival_forest(
      X = X,
      Y = as.vector(data[[Y]]),
      W = as.vector(data[[W]]),
      D = as.vector(data[[D]]),
      ...
    )
  }
  return(out)
}