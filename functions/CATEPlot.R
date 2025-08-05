CATEPlot <- function(csf, draw = TRUE, xlim = NULL, ylim = NULL) {
  if(!inherits(csf, c("causal_forest", "causal_survival_forest"))) {
    stop("csf must be a causal_survival_forest or causal_forest object.")
  }
  xintercept <- 100 * average_treatment_effect(csf)[1]
  x <- 100 * csf$predictions
  if (inherits(csf, c("causal_survival_forest"))) {
    xintercept <- -xintercept
    x <- -x
  }
  plot <- 
    tibble(x = x) |>
    ggplot() +
    geom_density(
      aes(x = x, y = after_stat(density)),
      color = "white", 
      fill = "darkgray"
    ) +
    geom_vline(
      xintercept = xintercept, 
      linetype = 2
    ) +
    xlab("cumulative incidence difference (percentage points)") + ylab("density") +
    scale_x_continuous(limits = xlim) + 
    scale_y_continuous(limits = ylim)
  if(draw) print(plot)
  return(invisible(plot))
}