mbcal <- function(mu_0_hat,
                  mu_1_hat,
                  Y,
                  W,
                  p0 = NULL,
                  p1 = NULL,
                  type = c("both", "treated", "control")) {
  type <- match.arg(type)
  Z <- ifelse(W == 1, 1, -1)
  deltahatlp <- stats::qlogis(mu_1_hat) - stats::qlogis(mu_0_hat)
  deltahatlpstar <- ifelse(W == 1, deltahatlp, -deltahatlp)
  if (is.null(p0) && is.null(p1)) {
    if (type == "treated") {
      ind.B <- which(W == 1)
      mod <- stats::glm(
        Y[ind.B] ~ deltahatlp[ind.B],
        family = "quasibinomial",
        offset = stats::qlogis(mu_0_hat[ind.B])
      )
    } else if (type == "control") {
      ind.A <- which(W == 0)
      minus_one <- rep(-1, length(ind.A))
      mod <- stats::glm(
        Y[ind.A] ~ 0 + minus_one + I(-deltahatlp[ind.A]),
        family = "quasibinomial",
        offset = stats::qlogis(mu_1_hat[ind.A])
      )
    } else if (type == "both") {
      os <- ifelse(W == 1, stats::qlogis(mu_0_hat), stats::qlogis(mu_1_hat))
      mod <- stats::glm(
        Y ~ 0 + Z + deltahatlpstar,
        family = "quasibinomial",
        offset = os
      )
    }
  } else if (!is.null(p0) && !is.null(p1)) {
    if (type == "treated") {
      ind.B <- which(W == 1)
      mod <- stats::glm(
        p1[ind.B] ~ deltahatlp[ind.B],
        family = "quasibinomial",
        offset = stats::qlogis(p0[ind.B])
      )
    } else if (type == "control") {
      ind.A <- which(W == 0)
      minus_one <- rep(-1, length(ind.A))
      mod <- stats::glm(
        p0[ind.A] ~ 0 + minus_one + I(-deltahatlp[ind.A]),
        family = "quasibinomial",
        offset = stats::qlogis(p1[ind.A])
      )
    } else if (type == "both") {
      os <- ifelse(W == 1, stats::qlogis(p0), stats::qlogis(p1))
      p01 <- ifelse(W == 1, p1, p0)
      mod <- stats::glm(
        p01 ~ 0 + Z + deltahatlpstar,
        family = "quasibinomial",
        offset = os
      )
    }
  } else {
    stop("supply both p0 and p1 (for estimand under simulation) or neither (empirical)")
  }
  return(mod)
}