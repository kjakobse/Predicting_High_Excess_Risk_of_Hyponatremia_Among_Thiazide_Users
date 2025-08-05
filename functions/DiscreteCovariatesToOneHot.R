DiscreteCovariatesToOneHot <- function(covariates) {
  if (!inherits(covariates, "data.frame")) {
    stop("'covariates' must be a data.frame or data.frame like object.")
  }
  if (ncol(covariates) == 0) {
    return(covariates)
  }
  for (i in seq_along(covariates)) {
    if (!is.factor(covariates[[i]])) {
      stop("Each covariate must be a factor.")
    }
  }
  X <- suppressMessages(
    map(
      names(covariates),
      ~ model.matrix(~ 0 + covariates[[.x]])
    ) |>
      map(~ as_tibble(.x)) |> 
      list_cbind()
  )
  X_names <- c()
  for (i in seq_along(names(covariates))) {
    for (j in seq_along(levels(covariates[[names(covariates)[i]]]))) {
      X_names <- c(
        X_names,
        paste0(
          names(covariates)[i],
          "_",
          levels(covariates[[names(covariates)[i]]])[j]
        )
      )
    }
  }
  names(X) <- X_names
  return(X)
}