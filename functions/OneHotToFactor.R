OneHotToFactor <- function(X) {
  nm <- names(X)
  tofactor <- str_subset(nm, "^X_\\d{2}_.{1,}")|> 
    str_sub(1,4) |> 
    unique()
  for (cov in tofactor) {
    levels <- str_subset(nm, glue::glue("^{cov}"))
    str <- vector("character", nrow(X))
    for (lev in levels) {
      str[X[[lev]] == 1] <- str_remove(lev, paste0(cov, "_"))
    }
    factor <- factor(str)
    X <- X |>
      mutate("{cov}" := factor) |>
      select(!all_of(levels))
  }
  return(select(X, order(colnames(X))))
}