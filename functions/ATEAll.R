ATEAll <- function (Y.orig, Y.hat,
                    W.orig, W.hat,
                    tau.hat.pointwise,
                    sample.weights,
                    clusters = NULL,
                    target.sample = c("all", "treated", "control"),
                    num.trees.for.weights = 500) {
  target.sample <- match.arg(target.sample)
  cluster.se <- TRUE
  if(is.null(clusters)) {
    clusters <- 1:NROW(Y.orig)
    cluster.se <- FALSE
  }
  if (all(W.orig %in% c(0, 1))) {
    if (min(W.hat) <= 0.05 && max(W.hat) >= 0.95) {
      rng <- range(W.hat)
      warning(paste0("Estimated treatment propensities take values between ", 
                     round(rng[1], 3), " and ", round(rng[2], 3), 
                     " and in particular get very close to 0 and 1. "))
    }
    else if (min(W.hat) <= 0.05) {
      warning(paste0("Estimated treatment propensities go as low as ", 
                     round(min(W.hat), 3), " which means that treatment ", 
                     "effects for some controls may not be well identified. "))
      if (min(W.hat) == 0) {
        W.hat[W.hat == 0] <- 0.00001
        warning(paste0("An estimated treatment propensity is 0, breaking the overlap assumption. ",
                       "The value has been shifted away from 0. However, the result is of dubious quality."))
      }
    }
    else if (max(W.hat) >= 0.95) {
      warning(paste0("Estimated treatment propensities go as high as ", 
                     round(max(W.hat), 3), " which means that treatment ", 
                     "effects for some treated units may not be well identified. "))
    }
  }
  .sigma2.hat <- function(DR.scores, tau.hat) {
    correction <- 
      Matrix::sparse.model.matrix(~factor(clusters) + 0, transpose = TRUE) %*% 
      (sweep(as.matrix(DR.scores), 2, tau.hat, "-") * sample.weights)
    n.adj <- sum(rowsum(sample.weights, clusters) > 0)
    Matrix::colSums(correction^2) / sum(sample.weights)^2 * 
      n.adj / (n.adj - 1)
  }
  if (all(W.orig %in% c(0, 1))) {
    debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
  }
  else {
    clusters <- 1:length(Y.orig)
    variance_forest <- grf::regression_forest(
      X.orig, 
      (W.orig - W.hat)^2, 
      clusters = clusters, 
      sample.weights = sample.weights, 
      num.trees = num.trees.for.weights, 
      ci.group.size = 1
    )
    V.hat <- grf::predict.regression_forest(variance_forest)$predictions
    debiasing.weights <- (W.orig - W.hat) / V.hat
  }
  control.idx <- which(W.orig == 0)
  treated.idx <- which(W.orig == 1)
  Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
  Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise
  if (target.sample == "all") {
    Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))
    DR.scores <- tau.hat.pointwise + debiasing.weights * Y.residual
    tau.hat <- weighted.mean(DR.scores, sample.weights)
    sigma2.hat <- .sigma2.hat(DR.scores, tau.hat)
    if (is.nan(tau.hat) | is.nan(sigma2.hat)) stop("Invalid sample provided. Are there any units?")
    return(c(estimate = tau.hat, std.err = sqrt(sigma2.hat)))
  } else if (target.sample == "treated") {
    tau.avg.raw <- weighted.mean(tau.hat.pointwise[treated.idx], sample.weights[treated.idx])
    tau.avg.var <- sum(sample.weights[treated.idx]^2 * (tau.hat.pointwise[treated.idx] - tau.avg.raw)^2) / 
      sum(sample.weights[treated.idx])^2
    gamma.control.raw <- W.hat[control.idx] / (1 - W.hat[control.idx])
    gamma.treated.raw <- rep(1, length(treated.idx))
  }
  else if (target.sample == "control") {
    tau.avg.raw <- weighted.mean(tau.hat.pointwise[control.idx], sample.weights[control.idx])
    tau.avg.var <- sum(sample.weights[control.idx]^2 * (tau.hat.pointwise[control.idx] - tau.avg.raw)^2) / 
      sum(sample.weights[control.idx])^2
    gamma.control.raw <- rep(1, length(control.idx))
    gamma.treated.raw <- (1 - W.hat[treated.idx]) / W.hat[treated.idx]
  }
  gamma <- rep(0, length(W.orig))
  gamma[control.idx] <- gamma.control.raw / 
    sum(sample.weights[control.idx] * gamma.control.raw) * sum(sample.weights)
  gamma[treated.idx] <- gamma.treated.raw / 
    sum(sample.weights[treated.idx] * gamma.treated.raw) * sum(sample.weights)
  dr.correction.all <- W.orig * gamma * (Y.orig - Y.hat.1) - (1 - W.orig) * gamma * (Y.orig - Y.hat.0)
  dr.correction <- weighted.mean(dr.correction.all, sample.weights)
  if (cluster.se) {
    correction.clust <- Matrix::sparse.model.matrix(~factor(clusters) + 0, transpose = TRUE) %*% 
      (dr.correction.all * sample.weights)
    n.adj <- sum(rowsum(sample.weights, clusters) > 0)
    sigma2.hat <- sum(correction.clust^2) / sum(sample.weights)^2 * n.adj / (n.adj - 1)
  }
  else {
    sigma2.hat <- sum(sample.weights^2 * dr.correction.all^2 / sum(sample.weights)^2) * 
      length(sample.weights[sample.weights > 0]) / (length(sample.weights[sample.weights > 0]) - 1)
  }
  tau.avg <- tau.avg.raw + dr.correction
  tau.se <- sqrt(tau.avg.var + sigma2.hat)
  if (is.nan(unname(tau.avg)) | is.nan(tau.se)) stop("Invalid sample provided. Are there any units in target.sample?")
  return(c(estimate = unname(tau.avg), std.err = tau.se))
}