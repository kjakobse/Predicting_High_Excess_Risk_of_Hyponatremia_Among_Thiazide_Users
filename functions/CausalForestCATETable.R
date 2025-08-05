#' @param cf An object of class causal_forest
#' @param cov_list A tibble (or tibble-like object) with covariates used to define axis-aligned subsets of the covariate space.
#' Each row specifies the part of each covariate to be included in the subset. 
#' For continuous covariates, a list-column with either length two character vectors specifying an interval, or a 
#' numeric vector with values to include (appropriate if covariate is discretized, e.g. only takes integer values).
#' For discrete covariates, a list column with a character vector specifying the levels to include in subset. These 
#' must be specified using the format 'covariate_level'.
#' The list-columns can be named, in which case these names are returned in the corresponding row of the output.
#' @param subset A logical vector, specifying the subset of the cohort to use in the CATE table
#' @param level The expected coverage of the provided confidence intervals

CausalForestCATETable <- function (cf,
                                   cov_list,
                                   subset = rep(TRUE, length(cf[["Y.orig"]])),
                                   level = 0.95) {
  if (!inherits(cf, c("causal_survival_forest", "causal_forest"))) {
    stop(
      glue::glue(
        "'cf' has class {class(cf)} but must be an object of class ",
        "'causal_survival_forest' or 'causal_forest'."
      )
    )
  }
  if (
    !(
      (is.list(subset) && 
       length(subset) == 1 && 
       is.logical(subset[[1]]) &&
       length(subset[[1]]) == length(cf[["Y.orig"]])) | 
      (is.logical(subset) && length(subset) == length(cf[["Y.orig"]]))
    )
  ) {
    stop(
      glue::glue(
        "subset must be a named one element list with a logical vector of ",
        "length {length(cf[['Y.orig']])}, or a logical vector of ",
        "length {length(cf[['Y.orig']])}."
      )
    )
  }
  str_helper <- function (...) {
    string = ""
    for (i in seq_along(cov_list)) {
      if (is.null(names(cov_list[[i]]))) {
        name <- str_remove(str_remove(...elt(i+1), names(cov_list)[i]), "^_")
      } else {
        name <- names(cov_list[[i]])[...elt(1)]
      }
      if (i == 1) {
        string <- paste0(string, name)
      } else {
        string <- paste0(string, "_", name)
      }
    }
    if (!is.null(subset_name)) {
      string <- paste0(string, "_", subset_name)
    }
    return(string)
  }
  subset_helper <- function (...) {
    for (i in seq_along(cov_list)) {
      name <- names(cov_list)[i]
      if (i == 1) {
        if (is.character(...elt(2)) && length(...elt(2)) == 1) {
          sub <- cf[["X.orig"]][[...elt(2)]] == 1
        } else if (is.character(...elt(2)) && length(...elt(2)) == 2) {
          if (suppressWarnings(name %in% names(cf[["X.orig"]]) & !any(is.na(as.numeric(...elt(2)))))) {
            sub <- cf[["X.orig"]][[name]] >= as.numeric(...elt(2)[1]) &
              cf[["X.orig"]][[name]] < as.numeric(...elt(2)[2])
          } else {
          sub <- cf[["X.orig"]][[...elt(2)[1]]] == 1 |
            cf[["X.orig"]][[...elt(2)[2]]] == 1
          }
        } else {
          if (suppressWarnings(name %in% names(cf[["X.orig"]]) & !any(is.na(as.numeric(...elt(2)))))) {
          sub <- cf[["X.orig"]][[name]] %in% ...elt(2)
          } else {
            sub <- (select(cf[["X.orig"]], all_of(...elt(2))) |>
              mutate(a = rowSums(across(everything()))) |>
              pull(a)) == 1
          }
        }
      } else {
        if (is.character(...elt(i+1)) && length(...elt(i+1)) == 1) {
          sub <- sub & cf[["X.orig"]][[...elt(i+1)]] == 1
        } else if (is.character(...elt(i+1)) && length(...elt(i+1)) == 2) {
          if (suppressWarnings(name %in% names(cf[["X.orig"]]) & !any(is.na(as.numeric(...elt(i+1)))))) {
            sub <- sub & cf[["X.orig"]][[name]] >= as.numeric(...elt(i+1)[1]) &
              cf[["X.orig"]][[name]] < as.numeric(...elt(i+1)[2])
          } else {
            sub <- sub & (cf[["X.orig"]][[...elt(i+1)[1]]] == 1 |
              cf[["X.orig"]][[...elt(i+1)[2]]] == 1)
          }
        } else {
          if (suppressWarnings(name %in% names(cf[["X.orig"]]) & !any(is.na(as.numeric(...elt(i+1)))))) {
            sub <- sub & cf[["X.orig"]][[name]] %in% ...elt(i+1)
          } else {
            sub <- sub & (select(cf[["X.orig"]], all_of(...elt(i+1))) |>
                            mutate(a = rowSums(across(everything()))) |>
                            pull(a)) == 1
          }
        }
      }
    }
    return(sub)
  }
  subset_name <- names(subset)
  if (is.list(subset)) subset <- subset[[1]]
  ci_names <- c(
    paste0(100 * level, "% CI - lower"),
    paste0(100 * level, "% CI - upper")
  )
  ci_names_pct <- c(
    paste0(100 * level, "% CI - lower (%)"),
    paste0(100 * level, "% CI - upper (%)")
  )
  cate_table <- purrr::pmap(
    c(
      list(id = seq_along(cov_list[[1]])),
      cov_list
    ),
    \(...) {
      string <- str_helper(...)
      sub <- subset_helper(...)
      CausalForestATESubgroupTable(
        NULL,
        c(paste0(names(cov_list), collapse = "_"), string),
        cf,
        level,
        sub & subset
      )
    }
  ) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      `estimate (%)` = sprintf("%.1f", 100 * estimate),
      !!sym(ci_names_pct[1]) := sprintf("%.1f", 100 * !!sym(ci_names[1])),
      !!sym(ci_names_pct[2]) := sprintf("%.1f", 100 * !!sym(ci_names[2]))
    )
  return(cate_table)
}