TestCalibrationData <- function(forest) {
  observation.weight <- observation_weights(forest)
  clusters <- if (length(forest$clusters) > 0) {
    forest$clusters
  }
  else {
    1:length(observation.weight)
  }
  if ("regression_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- dplyr::tibble(
      target = unname(forest$Y.orig), 
      mean.forest.prediction = mean.pred, 
      differential.forest.prediction = preds - mean.pred
    )
  }
  else if ("causal_forest" %in% class(forest)) {
    preds <- predict(forest)$predictions
    mean.pred <- weighted.mean(preds, observation.weight)
    DF <- dplyr::tibble(
      target = unname(forest$Y.orig - forest$Y.hat), 
      mean.forest.prediction = 
        unname(forest$W.orig - forest$W.hat) * mean.pred, 
      differential.forest.prediction = 
        unname(forest$W.orig - forest$W.hat) * (preds - mean.pred)
    )
  }
  else {
    stop("Calibration check not supported for this type of forest.")
  }
  best.linear.predictor <- lm(
    target ~ mean.forest.prediction + differential.forest.prediction + 0, 
    weights = observation.weight, 
    data = DF
  )
  DF <- DF |>
    dplyr::mutate(
      target.fitted = unname(predict(best.linear.predictor)),
      residuals = target - target.fitted
    )
  
  return(DF)
}

observation_weights <- function (forest) 
{
  if (is.null(forest$sample.weights)) {
    if (length(forest$clusters) == 0 || !forest$equalize.cluster.weights) {
      raw.weights <- rep(1, NROW(forest$Y.orig))
    }
    else {
      clust.factor <- factor(forest$clusters)
      inverse.counts <- 1 / 
        as.numeric(
          Matrix::colSums(Matrix::sparse.model.matrix(~clust.factor + 0))
        )
      raw.weights <- inverse.counts[as.numeric(clust.factor)]
    }
  }
  if (!is.null(forest$sample.weights)) {
    if (length(forest$clusters) == 0 || !forest$equalize.cluster.weights) {
      raw.weights <- forest$sample.weights
    }
    else {
      stop("Specifying non-null sample.weights is not allowed when equalize.cluster.weights = TRUE")
    }
  }
  return(raw.weights / sum(raw.weights))
}