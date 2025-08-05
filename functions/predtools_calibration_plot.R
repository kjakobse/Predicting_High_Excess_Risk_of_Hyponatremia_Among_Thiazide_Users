predtools_calibration_plot <- function (data, obs, follow_up = NULL, pred, group = NULL, nTiles = 10, 
                                        legendPosition = "right", title = NULL, x_lim = NULL, y_lim = NULL, 
                                        xlab = "Prediction", ylab = "Observation", points_col_list = NULL, 
                                        data_summary = FALSE) {
  if (!exists("obs") | !exists("pred")) 
    stop("obs and pred can not be null.")
  n_groups <- ifelse(is.null(group), 0, nrow(unique(data[, group])))
  if (is.null(follow_up)) 
    data$follow_up <- 1
  if (!is.null(group)) {
    dataDec_mods <- data %>% 
      group_by(!!sym(group)) %>% 
      mutate(decile = ntile(!!sym(pred), nTiles)) %>% 
      group_by(.data$decile, !!sym(group)) %>% 
      summarise(
        obsRate = mean(!!sym(obs)/follow_up, na.rm = T), 
        obsRate_SE = sd(!!sym(obs)/follow_up, na.rm = T) / sqrt(n()), 
        obsNo = n(), 
        predRate = mean(!!sym(pred), na.rm = T)
      )
    colnames(dataDec_mods)[colnames(dataDec_mods) == "group"] <- group
  }
  else {
    dataDec_mods <- data %>% 
      mutate(decile = ntile(!!sym(pred), nTiles)) %>% 
      group_by(.data$decile) %>% 
      summarise(
        obsRate = mean(!!sym(obs)/follow_up, na.rm = T), 
        obsRate_SE = sd(!!sym(obs)/follow_up, na.rm = T)/sqrt(n()), 
        obsNo = n(), 
        predRate = mean(!!sym(pred), na.rm = T)
      )
  }
  dataDec_mods$obsRate_UCL <- dataDec_mods$obsRate + 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods$obsRate_LCL <- dataDec_mods$obsRate - 1.96 * dataDec_mods$obsRate_SE
  dataDec_mods <- as.data.frame(dataDec_mods)
  if (!is.null(group)) {
    dataDec_mods[, group] <- factor(dataDec_mods[, group])
    calibPlot_obj <- ggplot(
      data = dataDec_mods, 
      aes(y = .data$obsRate, x = .data$predRate, group = !!sym(group), color = !!sym(group))) + 
      geom_point() + 
      coord_cartesian(
        xlim = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim), 
        ylim = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)
      ) + 
      geom_errorbar(aes(ymax = .data$obsRate_UCL, ymin = .data$obsRate_LCL)) + 
      geom_abline(intercept = 0, slope = 1) + 
      scale_color_manual(
        values = ifelse(
          rep(is.null(points_col_list), n_groups), 
          (ggplot2::scale_colour_brewer(palette = "Set3")$palette(8)[c(4:8, 1:3)])[c(1:n_groups)], 
          points_col_list
        )
      ) + 
      labs(
        x = ifelse(is.null(xlab), pred, xlab), 
        y = ifelse(is.null(ylab), obs, ylab), 
        title = title
      ) + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        legend.key = element_rect(fill = "white"), 
        axis.text = element_text(colour = "black", size = 12), 
        legend.position = legendPosition
      )
  }
  else {
    calibPlot_obj <- ggplot(
      data = dataDec_mods, aes(y = .data$obsRate, x = .data$predRate)) + 
      geom_point(color = ifelse(is.null(points_col_list), "black", points_col_list)) + 
      coord_cartesian(
        xlim = ifelse(rep(is.null(x_lim), 2), c(min(dataDec_mods$predRate), max(dataDec_mods$predRate)), x_lim), 
        ylim = ifelse(rep(is.null(y_lim), 2), c(min(dataDec_mods$obsRate_LCL), max(dataDec_mods$obsRate_UCL)), y_lim)
      ) +
      geom_errorbar(
        aes(ymax = .data$obsRate_UCL, ymin = .data$obsRate_LCL), 
        color = ifelse(is.null(points_col_list), "black", points_col_list)
      ) + 
      geom_abline(intercept = 0, slope = 1) + 
      scale_color_manual(
        values = ifelse(is.null(points_col_list), "black", points_col_list)
      ) + 
      labs(
        x = ifelse(is.null(xlab), 
                   pred, xlab), 
        y = ifelse(is.null(ylab), obs, ylab), 
        title = title) + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(colour = "black", size = 12), 
        legend.position = legendPosition
      )
  }
  res_list <- list(calibration_plot = calibPlot_obj)
  if (data_summary) 
    res_list$data_summary <- dataDec_mods
  return(res_list)
}