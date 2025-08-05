CausalForestCATEAllTable <- function (X.orig,
                                      Y.orig, 
                                      Y.hat,
                                      W.orig, 
                                      W.hat,
                                      tau.hat.pointwise,
                                      sample.weights,
                                      cov_list,
                                      subset = rep(TRUE, length(Y.orig)),
                                      level = 0.95,
                                      clusters = NULL,
                                      target.sample = c("all", "treated", "control"),
                                      num.trees.for.weights = 500) {
  target.sample <- match.arg(target.sample)
  if (
    !(
      (is.list(subset) && 
       length(subset) == 1 && 
       is.logical(subset[[1]]) &&
       length(subset[[1]]) == length(Y.orig)) | 
      (is.logical(subset) && length(subset) == length(Y.orig))
    )
  ) {
    stop(
      glue::glue(
        "subset must be a named one element list with a logical vector of ",
        "length {length(Y.orig)}, or a logical vector of ",
        "length {length(Y.orig)}."
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
          sub <- X.orig[[...elt(2)]] == 1
        } else if (is.character(...elt(2)) && length(...elt(2)) == 2) {
          if (suppressWarnings(name %in% names(X.orig) & !any(is.na(as.numeric(...elt(2)))))) {
            sub <- X.orig[[name]] >= as.numeric(...elt(2)[1]) &
              X.orig[[name]] < as.numeric(...elt(2)[2])
          } else {
            sub <- X.orig[[...elt(2)[1]]] == 1 |
              X.orig[[...elt(2)[2]]] == 1
          }
        } else {
          if (suppressWarnings(name %in% names(X.orig) & !any(is.na(as.numeric(...elt(2)))))) {
            sub <- X.orig[[name]] %in% ...elt(2)
          } else {
            sub <- (select(X.orig, all_of(...elt(2))) |>
                      mutate(a = rowSums(across(everything()))) |>
                      pull(a)) == 1
          }
        }
      } else {
        if (is.character(...elt(i+1)) && length(...elt(i+1)) == 1) {
          sub <- sub & X.orig[[...elt(i+1)]] == 1
        } else if (is.character(...elt(i+1)) && length(...elt(i+1)) == 2) {
          if (suppressWarnings(name %in% names(X.orig) & !any(is.na(as.numeric(...elt(i+1)))))) {
            sub <- sub & X.orig[[name]] >= as.numeric(...elt(i+1)[1]) &
              X.orig[[name]] < as.numeric(...elt(i+1)[2])
          } else {
            sub <- sub & (X.orig[[...elt(i+1)[1]]] == 1 |
                            X.orig[[...elt(i+1)[2]]] == 1)
          }
        } else {
          if (suppressWarnings(name %in% names(X.orig) & !any(is.na(as.numeric(...elt(i+1)))))) {
            sub <- sub & X.orig[[name]] %in% ...elt(i+1)
          } else {
            sub <- sub & (select(X.orig, all_of(...elt(i+1))) |>
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
      CausalForestATEAllSubgroupTable(
        NULL,
        c(paste0(names(cov_list), collapse = "_"), string),
        Y.orig, 
        Y.hat,
        W.orig, 
        W.hat,
        tau.hat.pointwise,
        sample.weights,
        level,
        sub & subset,
        clusters,
        target.sample,
        num.trees.for.weights
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