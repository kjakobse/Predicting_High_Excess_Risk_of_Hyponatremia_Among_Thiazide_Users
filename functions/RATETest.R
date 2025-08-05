RATETest <- function(forest,
                     priorities,
                     level = 0.95,
                     cov_type = c("continuous", "discrete"),
                     target = c("AUTOC", "QINI"),
                     q = seq(0.1, 1, by = 0.1),
                     R = 500,
                     subset = NULL,
                     debiasing.weights = NULL,
                     compliance.score = NULL,
                     num.trees.for.weights = 500) {
  if (is.null(priorities)) return(NULL)
  cov_type <- match.arg(cov_type)
  target <- match.arg(target)
  if(
    !(
      is.character(priorities) && length(priorities) == 1 && !is.null(names(forest$X.orig)) ||
      is.numeric(priorities) && length(priorities) == 1 ||
      is.numeric(priorities) && length(priorities) == length(forest$Y.orig)
    )
  ) {
    stop(
      "'priorities' must be a length one character vector, a length one ",
      "integer vector,\nor a numeric vector with the same length as ",
      "forest$Y.orig. If 'priorities' is a\ncharacter vector, forest$X.orig",
      "must have named columns."
    )
  }
  if (length(priorities) == 1) {
    if (is.matrix(forest$X.orig)) {
      priorities <- forest[["X.orig"]][, priorities]
    } else {
      priorities <- forest[["X.orig"]][[priorities]]
    }
  }
  if (!hasArg(q) && cov_type == "discrete") {
    q <- cumsum(rev(table(priorities))) /
      length(priorities)
    if (min(q) > 0.001) {
      q <- c(0.001, q)
    }
  }
  rate <- grf::rank_average_treatment_effect(
    forest = forest,
    priorities = priorities,
    target = target,
    q = q,
    R = R,
    subset = subset,
    debiasing.weights = debiasing.weights,
    compliance.score = compliance.score,
    num.trees.for.weights = num.trees.for.weights
  )
  confint <- rate$estimate + 
    dplyr::tibble(
      estimate = 0,
      lower = qnorm(0.5 - level / 2) * rate$std.err,
      upper = qnorm(0.5 + level / 2) * rate$std.err
    )
  pval <- 2 * pnorm(-abs(rate$estimate) / rate$std.err)
  out <- c(
    rate,
    list(
      confint = confint,
      pval = pval
    )
  )
  class(out) <- "rank_average_treatment_effect"
  return(out)
}