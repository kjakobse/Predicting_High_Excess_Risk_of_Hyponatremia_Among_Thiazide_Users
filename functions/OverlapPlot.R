OverlapPlot <- function(csf, draw = TRUE) {
  if(!inherits(csf, c("causal_forest", "causal_survival_forest"))) {
    stop("csf must be a causal_survival_forest or causal_forest object.")
  }
  plot <- 
    ggplot(tibble(W_hat = csf$W.hat)) +
    geom_density(aes(x = W_hat, y = after_stat(density)),
                 color = "white", fill = "darkgray") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    labs(title = "Estimated propensity scores for exposure") +
    xlab("Propensity") + ylab("")
  if(draw) print(plot)
  return(invisible(plot))
}