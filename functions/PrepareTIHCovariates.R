#' Prepare data frame with covariates for causal survival forest
#'
#' Covariate data in the TIH cohort contains both binary, discrete (with more
#' than two levels), and continuous columns. Discrete covariates with more than
#' two levels is prepared using one-hot encoding. The output is a data frame
#' with covariates encoded so they can be plugged into a causal survival forest.
#'
#' @param data A data frame or data frame like object (e.g. tibble) containing
#'   the covariates.
#' @param continuous <[dplyr::dplyr_tidy_select]> One unquoted expression
#'   specifying the continuous (and binary) covariates in data.
#' @param discrete <[dplyr::dplyr_tidy_select]> One unquoted expression
#'   specifying the discrete (more than two levels) covariates in data.
#' @param simplify Boolean. If FALSE (default), returns a list with the
#'   covariate data frames. If TRUE, length one lists (e.g. if data was a data
#'   frame) will be simplified to a data frame.

PrepareTIHCovariates <- function(data, continuous, discrete, simplify = FALSE) {
  data <- list_flatten(list(data))
  out <- map(
    data,
    \(data, continuous, discrete) {
      X_c <- select(data, {{ continuous }})
      X_d <- data |> 
        select({{ discrete }}) |>
        DiscreteCovariatesToOneHot()
      bind_cols(X_c, X_d)
    },
    continuous = {{ continuous }},
    discrete = {{ discrete }}
  )
  if (simplify && length(out) == 1) {
    out <- out[[1]]
  }
  return(out)
}