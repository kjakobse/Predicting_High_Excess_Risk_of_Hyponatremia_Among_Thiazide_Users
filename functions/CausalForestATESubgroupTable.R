CausalForestATESubgroupTable <- function(x, y, cf, level, subset = NULL) {
  ci_names <- c(
    paste0(100 * level, "% CI - lower"),
    paste0(100 * level, "% CI - upper")
  )
  if (is.null(x) && is.null(y) && is.null(subset)) {
    ate <- grf::average_treatment_effect(cf)
    return(
      tibble(
        subgroup = "Full population",
        n = length(cf$predictions),
        estimate = ate[1],
        "{ci_names[1]}" := ate[1] + qnorm(0.5 - level / 2) * ate[2],
        "{ci_names[2]}" := ate[1] + qnorm(0.5 + level / 2) * ate[2]
      )
    )
  } else if (is.null(subset)) {
    subset <- cf[["X.orig"]][[y]] %in% x
    ate <- grf::average_treatment_effect(cf, subset = subset)
  } else {
    if (length(y) != 2) {
      stop(glue::glue("y must be a character vector of length 2, ",
                      "not {length(y)}, containing a column name ",
                      "and a name to use in the column"))
    }
    y <- c("custom", y)
    if (sum(subset) > 1) {
      ate <- grf::average_treatment_effect(cf, subset = subset)
    }
  }
  if (sum(subset) > 1) {
    tbl <- dplyr::tibble(
      estimate = ate[1],
      std_err = ate[2],
      n = sum(subset),
      "{ci_names[1]}" := ate[1] + qnorm(0.5 - level / 2) * ate[2],
      "{ci_names[2]}" := ate[1] + qnorm(0.5 + level / 2) * ate[2]
    )
  } else {
    tbl <- dplyr::tibble(
      estimate = NA_real_,
      std_err = NA_real_,
      n = sum(subset),
      "{ci_names[1]}" := NA_real_,
      "{ci_names[2]}" := NA_real_
    )
  }
  switch(
    EXPR = y[1],
    custom = {
      tbl <- dplyr::bind_cols(
        dplyr::tibble("{y[2]}" := y[3]),
        tbl
      )
    }
  )
  return(tbl)
}