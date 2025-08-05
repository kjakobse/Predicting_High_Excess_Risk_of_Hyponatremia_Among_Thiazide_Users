CATEPlots <- function(object, ...) {
  UseMethod("CATEPlots")
}

CATEPlots.default <- function(object,
                              X,
                              W,
                              predictions,
                              ate,
                              return = FALSE,
                              smooth_method = "lm",
                              hex_bins = 100L,
                              linewidth = 0.8,
                              path = NULL,
                              filename = NULL,
                              width = 12,
                              height = 8,
                              unit = "cm",
                              scale = 1,
                              res = 144,
                              quality = 75,
                              draw = FALSE,
                              name_pairs = NULL, # tibble with columns short and full (replace short with full)
                              ...) {
  # Convert one-hot encodings to factors
  X_fct <- OneHotToFactor(X)
  
  # Rename covariates if name_pairs is provided
  if (!is.null(name_pairs)) {
    X_fct <- X_fct |>
      dplyr::rename_with(
        \(nm, name_pairs) {
          dplyr::left_join(
            dplyr::tibble(short = nm),
            name_pairs,
            by = "short"
          ) |>
            dplyr::pull(full)
        },
        name_pairs = name_pairs
      )
  }
  
  # combine covariates and predictions as plot data
  plot_data <- dplyr::mutate(X_fct, W = W, predictions = predictions)
  
  # make a plot for each column in X_fct
  plots <- purrr::map(
    names(X_fct),
    \(nm) {
      # check for NAs in covariate and remove
      if (any(is.na(plot_data[[nm]]))) {
        data <- dplyr::filter(plot_data, !is.na(.data[[nm]]))
      } else {
        data <- plot_data
      }
      
      # create boxplot if factor and hexplot if continuous
      if(nrow(data) > 0) {
        plot <- ggplot2::ggplot(
          data, 
          ggplot2::aes(x = .data[[nm]], y = predictions)
        )
        if (is.factor(X_fct[[nm]]) || length(unique(X_fct[[nm]])) == 2) {
          plot <- plot + 
            ggplot2::geom_violin(
              ggplot2::aes(fill = ggplot2::cut_width(.data[[nm]], 1)), 
              width = 0.8
            ) +
            ggsci::scale_fill_jama() +
            ggplot2::ylab("Increased risk of hyponatremia") +
            ggplot2::theme(
              legend.position = "none"
            ) 
          if (!is.factor(X_fct[[nm]])) {
            plot <- plot + 
              ggplot2::geom_boxplot(
                ggplot2::aes(group = ggplot2::cut_width(.data[[nm]], 1)), 
                width = 0.15, 
                color = "grey", 
                alpha = 0.2,
                outlier.size = 0.8
              ) +
              ggplot2::scale_x_continuous(breaks = c(0, 1), labels = c("0", "1"))
          } else {
            plot <- plot + 
              ggplot2::geom_boxplot(
                width = 0.15, 
                color = "grey", 
                alpha = 0.2,
                outlier.size = 0.8
              )
          }
        } else {
          plot <- plot +
            ggplot2::geom_hex(bins = hex_bins) +
            ggplot2::geom_hline(yintercept = ate, linetype = 2, linewidth = linewidth) + 
            ggplot2::geom_smooth(method = smooth_method, linetype = 1, linewidth = linewidth, color = "red") + 
            ggplot2::ylab("Increased risk of hyponatremia")
        }
      } else {
        plot <- NULL
      }
      # print plots if user selected
      if(draw && !is.null(plot)) print(plot)
      if (!is.null(path) && !is.null(plot)) {
        if (is.null(filename)) filename <- "plot"
        ggplot2::ggsave(
          filename = paste0(path, filename, "_", nm, ".jpg"),
          plot = plot,
          device = "jpeg",
          width = width,
          height = height,
          units = unit,
          dpi = res,
          quality = quality,
          scale = scale
        )
      }
      return(plot)
    }
  )
  names(plots) <- names(X_fct)
  
  # create CATE density plot
  plot_density <- 
    plot_data |>
    ggplot2::ggplot() +
    ggplot2::geom_density(
      ggplot2::aes(x = predictions, y = ggplot2::after_stat(density)),
      color = "white", 
      fill = "darkgray"
    ) +
    ggplot2::geom_vline(
      xintercept = ate, 
      linetype = 2,
      linewidth = linewidth
    ) + 
    ggplot2::coord_cartesian(xlim = c(-0.01, 0.15), ylim = c(0, 105)) +
    ggplot2::xlab("Increased risk of hyponatremia") + 
    ggplot2::ylab("density")
  if(draw) print(plot_density)
  if (!is.null(path)) {
    if (is.null(filename)) filename <- "plot"
    ggplot2::ggsave(
      filename = paste0(path, filename, "_density", ".jpg"),
      plot = plot_density,
      device = "jpeg",
      width = width,
      height = height,
      units = unit,
      dpi = res,
      quality = quality,
      scale = scale
    )
  }
  
  # create CATE cdf plot (overall + among exposed)
  plot_cdf <- 
    plot_data |>
    ggplot2::ggplot() +
    ggplot2::stat_ecdf(
      ggplot2::aes(x = predictions, y = ggplot2::after_stat(ecdf)),
      geom = "step",
      color = "black", 
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = ate, 
      linetype = 2,
      linewidth = linewidth
    ) +
    ggplot2::xlab("Increased risk of hyponatremia") + 
    ggplot2::ylab("cumulative distribution")
  if(draw) print(plot_cdf)
  if (!is.null(path)) {
    if (is.null(filename)) filename <- "plot"
    ggplot2::ggsave(
      filename = paste0(path, filename, "_cdf_all", ".jpg"),
      plot = plot_cdf,
      device = "jpeg",
      width = width,
      height = height,
      units = unit,
      dpi = res,
      quality = quality,
      scale = scale
    )
  }
  plot_cdf_exposed <- 
    plot_data |>
    filter(W == 1) |>
    ggplot2::ggplot() +
    ggplot2::stat_ecdf(
      ggplot2::aes(x = predictions, y = ggplot2::after_stat(ecdf)),
      geom = "step",
      color = "black", 
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = ate, 
      linetype = 2,
      linewidth = linewidth
    ) +
    ggplot2::xlab("Increased risk of hyponatremia") + 
    ggplot2::ylab("cumulative distribution")
  if(draw) print(plot_cdf_exposed)
  if (!is.null(path)) {
    if (is.null(filename)) filename <- "plot"
    ggplot2::ggsave(
      filename = paste0(path, filename, "_cdf_exposed", ".jpg"),
      plot = plot_cdf_exposed,
      device = "jpeg",
      width = width,
      height = height,
      units = unit,
      dpi = res,
      quality = quality,
      scale = scale
    )
  }
  
  # return plots or invisible NULL depending on user choice
  if (return) {
    return(
      list(
        density = plot_density, 
        cdf_all = plot_cdf,
        cdf_exposed = plot_cdf_exposed,
        covariates = plots
      )
    )
  } else {
    return(invisible(NULL))
  }
}

CATEPlots.causal_forest <- function(
    object,
    X = object$X, 
    predictions = object$predictions, 
    ate = grf::average_treatment_effect(object)[1], 
    draw = TRUE,
    xlim = NULL,
    ylim = NULL
) {
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

GetScale <- function(plot,
                     width,
                     height,
                     unit = "in") {
  h <- grid::convertHeight(sum(plot$heights), unit, TRUE)
  w <- grid::convertWidth(sum(plot$width), unit, TRUE)
  max(c(w/width_wanted, h/height_wanted))
}