CovariateBalance <- function(cf,
                             plots = c("all", "Love", "density", "ecdf"),
                             balance_table = TRUE,
                             covariates = NULL,
                             names = NULL,
                             factor = NULL,
                             treatment_name = "W",
                             love_breaks = NULL,
                             love_xlim = NULL,
                             love_scale_color = NULL,
                             cd_ncol = NULL,
                             cd_nrow = NULL,
                             cd_x_scale_width = NULL,
                             cd_bar_width = NULL,
                             cd_scale_fill = NULL,
                             ec_nrow = NULL,
                             ec_ncol = NULL,
                             ec_x_scale_width = NULL,
                             ec_scale_color = NULL,
                             plot_data_only = FALSE, 
                             plots_only = FALSE) {
  if (plots_only) {
    love_plot_data <- cf$love_data
    plot_type <- cf$cd_type
    plot_data <- cf$cd_data
    X_type <- cf$X_type
    X_orig <- plot_type
    which_fct <- purrr::map_lgl(
      dplyr::select(plot_data, -c(dplyr::all_of(treatment_name), "IPW")), 
      is.factor
    )
    if (any(names(cf) == "X_orig_old")) {
      X_orig_old <- cf$X_orig_old
    }
  }
  # check input
  plots <- match.arg(plots)
  love <- plots == "all" | plots == "Love"
  density <- plots == "all" | plots == "density"
  ecdf <- plots == "all" | plots == "ecdf"
  bal_tab <- balance_table
  if (!plots_only) {
    stopifnot(
      `cf must be a ´causal_forest´ object` = 
        inherits(cf, c("causal_forest", "causal_survival_forest"))
    )
    if (!inherits(cf, "grf")) {
      warning(
        "cf doesn't inherit from class ´grf´. Make sure your ´causal_forest'\nobject complies with the class structure used by grf::causal_forest()."
      )
    }
  }
  if (!plots_only) {
    if(is.null(names)) {
      if(!inherits(cf$X.orig, "tbl_df")) {
        X_orig <- dplyr::as_tibble(cf$X.orig, .name_repair = "unique")
      } else {
        X_orig <- cf$X.orig
      }
      if (!is.null(covariates)) {
        if (is.null(colnames(cf$X.orig))) {
          X_orig <- dplyr::select(X_orig, {{ covariates }})
        } else {
          X_orig <- dplyr::select(X_orig, dplyr::all_of(covariates))
        }
      }
    } else if (is.character(names)) {
      if(!inherits(cf$X.orig, "tbl_df")) {
        X_orig <- tibble::as_tibble(cf$X.orig, .name_repair = "minimal")
      } else {
        X_orig <- cf$X.orig
      }
      if (is.null(names(names))) {
        X_orig <- purrr::set_names(X_orig, names)
        if (!is.null(covariates)) {
          X_orig <- dplyr::select(X_orig, dplyr::all_of(covariates))
        }
      } else {
        which_fct <- purrr::map_lgl(
          names,
          function(names) {
            sum(grepl(paste0("^", names), names(X_orig))) > 1
          }
        )
        if (sum(which_fct) > 0) {
          X_orig_fct <- X_orig
          for (i in which(which_fct)) {
            which_var <- which(grepl(paste0("^", names[i]), names(X_orig_fct)))
            tmp <- vector("character", nrow(X_orig_fct))
            for (j in which_var) {
              tmp[dplyr::pull(X_orig_fct, j) == 1] <-
                sub(paste0("^", names[i], "_"), "", names(X_orig_fct)[j])
            }
            X_orig_fct <- X_orig_fct |>
              dplyr::select(!dplyr::all_of(which_var)) |>
              dplyr::mutate(
                "{names(names)[i]}" := factor(tmp)
              )
          }
          X_orig_fct <- dplyr::rename(X_orig_fct, dplyr::all_of(names[!which_fct]))
          X_orig <- dplyr::rename(X_orig, dplyr::all_of(names[!which_fct]))
          for (i in which(which_fct)) {
            names(X_orig)[grepl(paste0("^", names[i]), names(X_orig))] <-
              sub(
                paste0(names[i], "_"),
                paste0(names(names)[i], " - "),
                names(X_orig)[grepl(paste0("^", names[i]), names(X_orig))]
              )
            if (!is.null(factor[[names(names)[i]]])) {
              for (j in which(grepl(paste0("^", names(names)[i]), names(X_orig)))) {
                names(X_orig)[j] <-
                  sub(
                    factor[[names(names)[i]]][
                      factor[[names(names)[i]]] == 
                        sub(paste(names(names)[i], "- "), "", names(X_orig)[j])
                    ],
                    names(factor[[names(names)[i]]])[
                      factor[[names(names)[i]]] == 
                        sub(paste(names(names)[i], "- "), "", names(X_orig)[j])
                    ],
                    names(X_orig)[j]
                  )
              }
            }
          }
          if (!is.null(covariates)) {
            covariates <- unique(c(covariates, names(names)))
            X_orig_fct <- dplyr::select(X_orig_fct, dplyr::all_of(covariates))
            X_orig <- dplyr::select(X_orig, dplyr::starts_with(covariates))
          }
        } else {
          X_orig <- dplyr::rename(X_orig, dplyr::all_of(names))
          if (!is.null(covariates)) {
            covariates <- unique(c(covariates, names(names)))
            X_orig <- dplyr::select(X_orig, dplyr::all_of(covariates))
          }
        }
      }
    } else {
      stop(
        "names must either be NULL to use the original names form cf$X.orig,\na character vector with covariate names that replace names(cf$X.orig), or\na named character vector with the original names and the replacement names\nin the names attribute."
      )
    }
    
    if (!is.null(factor)) {
      if (!exists("X_orig_fct")) {
        X_orig_fct <- X_orig
      }
      factor_name <- intersect(names(X_orig_fct), names(factor))
      for (i in factor_name) {
        X_orig_fct[[i]] <- factor(
          X_orig_fct[[i]],
          levels = unlist(factor[[i]], use.names = FALSE),
          labels = names(factor[[i]])
        )
      }
      which_fct <- purrr::map_lgl(X_orig_fct, is.factor)
    }
  }
  if (!(is.null(cd_x_scale_width) || is.numeric(cd_x_scale_width))) {
    stop(
      glue::glue(
        "'cd_x_scale_width' must be null or a numeric vector of length 1 or ",
        "{nol(X_orig)},"
      )
    )
  } 
  if (!(length(cd_x_scale_width) %in% c(0, 1, ncol(X_orig)))) {
    stop(
      glue::glue(
        "'cd_x_scale_width' has the wrong length. Must be NULL or have ",
        "length 1 or {ncol(X_orig)}."
      )
    )
  }
  if (!(is.null(ec_x_scale_width) || is.numeric(ec_x_scale_width))) {
    stop(
      glue::glue(
        "'ec_x_scale_width' must be null or a numeric vector of length 1 or ",
        "{nol(X_orig)},"
      )
    )
  } 
  if (!(length(ec_x_scale_width) %in% c(0, 1, ncol(X_orig)))) {
    stop(
      glue::glue(
        "'ec_x_scale_width' has the wrong length. Must be NULL or have ",
        "length 1 or {ncol(X_orig)}."
      )
    )
  }
  if (!(is.character(treatment_name) && length(treatment_name) == 1)) {
    stop(
      "'treatment_name' must be a single string with the treatment name to
display in the covariate distribution plot."
    )
  }
  if (!(is.null(cd_bar_width) || is.numeric(cd_bar_width))) {
    stop(
      glue::glue(
        "'cd_bar_width' must be null or a numeric vectr of length 1 or ",
        "{nol(X_orig)},"
      )
    )
  } 
  if (!(length(cd_bar_width) %in% c(0, 1, ncol(X_orig)))) {
    stop(
      glue::glue(
        "'cd_bar_width' has the wrong length. Must be NULL or have ",
        "length 1 or {ncol(X_orig)}."
      )
    )
  }
  if (!(is.null(cd_nrow) || (is.numeric(cd_nrow) && cd_nrow > 0.5))) {
    stop("'cd_nrow' must be NULL or a positive integer.")
  }
  if (!(is.null(cd_ncol) || (is.numeric(cd_ncol) && cd_ncol > 0.5))) {
    stop("'ec_ncol' must be NULL or a positive integer.")
  }
  if (!(is.null(ec_nrow) || (is.numeric(ec_nrow) && ec_nrow > 0.5))) {
    stop("'ec_nrow' must be NULL or a positive integer.")
  }
  if (!(is.null(ec_ncol) || (is.numeric(ec_ncol) && ec_ncol > 0.5))) {
    stop("'ec_ncol' must be NULL or a positive integer.")
  }
  
  # initialize balance table
  if (bal_tab) {
    if (!plots_only) {
      balance_table <- list(unadjusted = NA, adjusted = NA)
    } else {
      balance_table <- cf$balance_table
    }
  }
  
  # create love plot
  if (love) {
    if (!plots_only) {
      funs <- c(mean = mean, var = var)
      df_0 <- X_orig
      df_0[cf$W.orig == 1,] <- NA
      df_1 <- X_orig
      df_1[cf$W.orig == 0,] <- NA
      df_ori <- dplyr::as_tibble(cbind(df_0, df_1), .name_repair = "minimal") |>
        purrr::map(\(x) tibble::as_tibble(purrr::map(funs, rlang::exec, x, na.rm = TRUE))) |>
        purrr::list_rbind() |>
        dplyr::mutate(name = rep(names({{ X_orig }}), 2)) |>
        (\(df) {
          out <- dplyr::summarise(
            df,
            mean_control = dplyr::first(.data$mean),
            mean_treated = dplyr::last(.data$mean),
            std_abs_mean_diff = abs(diff(.data$mean)) / sqrt(sum(.data$var)),
            var_ratio = dplyr::last(.data$var) / dplyr::first(.data$var),
            .by = "name"
          )
          if (any(stringr::str_detect(df$name, " - "))) {
            out <- rbind(
              out, 
              df |>
                dplyr::mutate(name_comb = stringr::str_remove(name, " - .{1,}$")) |>
                dplyr::mutate(
                  W = dplyr::case_when(dplyr::first(.data$mean) == .data$mean ~ 0, dplyr::last(.data$mean) == .data$mean ~ 1),
                  .by = "name"
                ) |>
                dplyr::summarise(
                  mean_control = mean(.data$mean[.data$W == 0]),
                  mean_treated = mean(.data$mean[.data$W == 1]),
                  std_abs_mean_diff = mean(abs(.data$mean[.data$W == 1] - .data$mean[.data$W == 0])) / 
                    sqrt(sum(.data$mean[.data$W == 1] + .data$mean[.data$W == 0]) / length(.data$mean[.data$W == 1])^2),
                  var_ratio = sum(.data$var[.data$W == 1]) / sum(.data$var[.data$W == 0]),
                  .by = "name_comb"
                ) |>
                rename("name" = name_comb)
            )
          }
          return(out)
        })()
        
      weighted.var <- function(x, w, na.rm = TRUE) {
        na_index <- is.na(w) | is.na(x)
        if (any(na_index)) {
          if(!na.rm)
            return(NA_real_)
          ok <- !na_index
          x <- x[ok]
          w <- w[ok]
        }
        cov.wt(matrix(x, ncol = 1), w)$cov[1, 1]
      }
      funs <- c(mean = weighted.mean, var = weighted.var)
      df_adj <- dplyr::as_tibble(cbind(df_0, df_1), .name_repair = "minimal") |>
        purrr::map(
          \(x) {
            dplyr::as_tibble(
              purrr::map(
                funs, 
                purrr::exec, 
                x,
                cf$W.orig / cf$W.hat + (1 - cf$W.orig) / (1 - cf$W.hat),
                na.rm = TRUE)
            )
          }
        ) |>
        purrr::list_rbind() |>
        dplyr::mutate(name = rep(names({{ X_orig }}), 2)) |>
        (\(df) {
          out <- dplyr::summarise(
            df,
            mean_control = dplyr::first(.data$mean),
            mean_treated = dplyr::last(.data$mean),
            std_abs_mean_diff_adj = abs(diff(.data$mean)) / sqrt(sum(.data$var)),
            var_ratio = dplyr::last(.data$var) / dplyr::first(.data$var),
            .by = "name"
          )
          if (any(stringr::str_detect(df$name, " - "))) {
            out <- rbind(
              out, 
              df |>
                dplyr::mutate(name_comb = stringr::str_remove(name, " - .{1,}$")) |>
                dplyr::mutate(
                  W = dplyr::case_when(dplyr::first(.data$mean) == .data$mean ~ 0, dplyr::last(.data$mean) == .data$mean ~ 1),
                  .by = "name"
                ) |>
                dplyr::summarise(
                  mean_control = mean(.data$mean[.data$W == 0]),
                  mean_treated = mean(.data$mean[.data$W == 1]),
                  std_abs_mean_diff_adj = mean(abs(.data$mean[.data$W == 1] - .data$mean[.data$W == 0])) / 
                    sqrt(sum(.data$mean[.data$W == 1] + .data$mean[.data$W == 0]) / length(.data$mean[.data$W == 1])^2),
                  var_ratio = sum(.data$var[.data$W == 1]) / sum(.data$var[.data$W == 0]),
                  .by = "name_comb"
                ) |>
                rename("name" = name_comb)
            )
          }
          return(out)
        })()
      love_plot_data <- dplyr::inner_join(
        dplyr::select(df_ori, "name", "std_abs_mean_diff"),
        dplyr::select(df_adj, "name", "std_abs_mean_diff_adj"),
        by = "name",
        relationship = "one-to-one"
      ) |>
        dplyr::arrange(.data$std_abs_mean_diff) |>
        dplyr::rename(
          covariate_name = "name",
          `Before adjustment` = "std_abs_mean_diff",
          `After adjustment` = "std_abs_mean_diff_adj"
        ) |>
        tidyr::pivot_longer(
          cols = c("Before adjustment", "After adjustment"),
          names_to = "Cohort",
          values_to = "value"
        )
      covariate_name_ordered <- love_plot_data |>
        dplyr::filter(.data$Cohort == "Before adjustment") |>
        dplyr::arrange(.data$value)
      love_plot_data <- love_plot_data |>
        dplyr::mutate(
          covariate_name = factor(
            .data$covariate_name,
            levels = covariate_name_ordered$covariate_name
          )
        )
      
      if (bal_tab) {
        balance_table$unadjusted <- df_ori |>
          mutate(
            var_ratio = ifelse(grepl(" - ", name), NA_real_, var_ratio)
          )
        balance_table$adjusted <- df_adj |>
          mutate(
            var_ratio = ifelse(grepl(" - ", name), NA_real_, var_ratio)
          )
      }
    }
    
    if (!plot_data_only) {
      if (is.null(love_scale_color)) {
        if (requireNamespace("ggsci", quitely = TRUE)) {
          love_scale_color <- ggsci::scale_color_jama()
        } else {
          love_scale_color <- ggplot2::scale_color_discrete()
        }
      }
      love_plot <- love_plot_data |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data$value, 
            y = .data$covariate_name, 
            color = .data$Cohort
          )
        ) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_line(
          ggplot2::aes(group = .data$Cohort), 
          orientation = "y"
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = 1) +
        ggplot2::geom_vline(xintercept = 0.1, linetype = 2) +
        ggplot2::xlab("Absolute standardized mean difference") + 
        ggplot2::ylab("") + 
        love_scale_color
      if (!(is.null(love_breaks) && is.null(love_xlim))) {
        love_plot <- love_plot +
          ggplot2::scale_x_continuous(
            breaks = love_breaks,
            limits = love_xlim
          )
      } else if (!(is.null(love_breaks))) {
        love_plot <- love_plot +
          ggplot2::scale_x_continuous(
            breaks = love_breaks
          )
      } else if (!(is.null(love_xlim))) {
        love_plot <- love_plot +
          ggplot2::scale_x_continuous(
            limits = love_xlim
          )
      }
    }
  }
  
  # create covariate distribution plot
  if(density | ecdf) {
    if (exists("X_orig_fct")) {
      X_orig_old <- X_orig
      X_orig <- X_orig_fct
    }
    
    if (!plots_only) {
      plot_data <- X_orig |>
        dplyr::mutate(
          dplyr::across(dplyr::where(is.integer), as.numeric),
          "{treatment_name}" := as.factor(cf$W.orig),
          IPW = ifelse(
            .data[[treatment_name]] == 1, 
            1 / cf$W.hat, 
            1 / (1 - cf$W.hat)
          )
        )
    }
    
    if (density) {
      if (!plots_only) {
        plot_type <- X_orig |>
          dplyr::summarise(
            dplyr::across(
              dplyr::everything(), \(x) ifelse(is.factor(x), "bar", "hist")
            )
          )
      }
      
      if (!plot_data_only) {
        if (is.null(cd_x_scale_width)) {
          nm_num <- names(plot_type)
          if (exists("which_fct")) {
            nm_num <- nm_num[!(names(plot_type) %in% names(names)[which_fct])]
          }
          if (length(nm_num) == 0) {
            cd_x_scale_width <- vector("numeric", if(is.null(dim(plot_type))) length(plot_type) else ncol(plot_type))
          } else {
            scale_width_dt <- data.table::as.data.table(plot_data)
            scale_width_dt <- scale_width_dt[
              , 
              (nm_num) := lapply(.SD, \(x) signif(diff(range(x)), 1) / 5), 
              .SDcols = nm_num
            ]
            cd_x_scale_width <- suppressWarnings(
              as.numeric(scale_width_dt[1])[seq_len(if(is.null(dim(plot_type))) length(plot_type) else ncol(plot_type))]
            )
          }
        } else if (length(cd_x_scale_width) == 1) {
          cd_x_scale_width <- rep(
            cd_x_scale_width,
            if(is.null(dim(plot_type))) length(plot_type) else ncol(plot_type)
          )
        }
        
        if (is.null(cd_bar_width)) {
          bar_width_dt <- plot_data |>
            dplyr::mutate(dplyr::across(dplyr::where(is.factor), as.numeric)) |>
            data.table::as.data.table()
          bar_width_dt <- bar_width_dt[
            , 
            (names(plot_type)) := purrr::map2(
              .SD, 
              plot_type,
              \(x, y) {
                if (any(y == "bar")) {
                  0.9 * ggplot2::resolution(x)
                } else {
                  diff(range(x)) / 50
                }
              }
            ), 
            .SDcols = names(X_orig)
          ]
          cd_bar_width <- as.numeric(bar_width_dt[1])[seq_len(if(is.null(dim(X_orig))) length(X_orig) else ncol(X_orig))]
        } else if (length(cd_bar_width) == 1) {
          cd_bar_width <- rep(
            cd_bar_width,
            if(is.null(dim(X_orig))) length(X_orig) else ncol(X_orig)
          )
        }
        
        # helper functions
        floor_dec <- function(x, digits = 1) {
          round(x - 5 * 10^(-digits - 1), digits)
        }
        ceiling_dec <- function(x, digits = 1) {
          round(x + 5 * 10^(-digits - 1), digits)
        }
        decimalplaces <- function(x) {
          if (abs(x - round(x)) > .Machine$double.eps^0.5) {
            nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
          } else {
            0
          }
        }
        
        if (is.null(cd_scale_fill)) {
          if (requireNamespace("ggsci", quitely = TRUE)) {
            cd_scale_fill <- ggsci::scale_fill_jama()
          } else {
            cd_scale_fill <- ggplot2::scale_fill_discrete()
          }
        }
        
        cov_plots_unadjusted <- purrr::pmap(
          list(
            names = names(plot_type),
            type = plot_type,
            cd_x_scale_width = cd_x_scale_width,
            cd_bar_width = cd_bar_width
          ),
          \(names, type, cd_x_scale_width, cd_bar_width) {
            plot_data <- plot_data |>
              dplyr::select(
                dplyr::all_of(treatment_name), 
                "IPW", 
                "covariate_values" = dplyr::all_of(names)
              ) |>
              dplyr::mutate(covariate_name = {{ names }})
            covariate_values <- plot_data[["covariate_values"]]
            p <- plot_data |>
              ggplot2::ggplot() + 
              ggplot2::facet_grid(~ .data$covariate_name)
            if (type == "bar") {
              p <- p +
                ggplot2::geom_bar(
                  ggplot2::aes(
                    x = .data$covariate_values,
                    fill = .data[[treatment_name]]
                  ),
                  alpha = 0.5,
                  position = "dodge",
                  width = cd_bar_width
                )
            } else {
              p <- p +
                ggplot2::geom_histogram(
                  ggplot2::aes(
                    x = .data$covariate_values,
                    y = ggplot2::after_stat(density),
                    fill = .data[[treatment_name]]
                  ),
                  alpha = 0.5,
                  position = "identity",
                  binwidth = cd_bar_width
                ) + 
                ggplot2::scale_x_continuous(
                  breaks = seq(
                    floor_dec(
                      min(covariate_values),
                      digits = decimalplaces(cd_x_scale_width)
                    ),
                    ceiling_dec(
                      max(covariate_values),
                      digits = decimalplaces(cd_x_scale_width)
                    ),
                    cd_x_scale_width
                  )
                )
            }
            p <- p +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              cd_scale_fill
          }
        )
        cov_plots_adjusted <- purrr::pmap(
          list(
            names = names(plot_type),
            type = plot_type,
            cd_x_scale_width = cd_x_scale_width,
            cd_bar_width = cd_bar_width
          ),
          \(names, type, cd_x_scale_width, cd_bar_width) {
            plot_data <- plot_data |>
              dplyr::select(
                dplyr::all_of(treatment_name), 
                "IPW", 
                covariate_values = dplyr::all_of(names)
              ) |>
              dplyr::mutate(covariate_name = {{ names }})
            covariate_values <- plot_data[["covariate_values"]]
            p <- plot_data |>
              ggplot2::ggplot() + 
              ggplot2::facet_grid(~ .data$covariate_name)
            if (type == "bar") {
              p <- p +
                ggplot2::geom_bar(
                  ggplot2::aes(
                    x = .data$covariate_values,
                    weight = .data$IPW,
                    fill = .data[[treatment_name]]
                  ),
                  alpha = 0.5,
                  position = "dodge",
                  width = cd_bar_width
                )
            } else {
              p <- p +
                ggplot2::geom_histogram(
                  ggplot2::aes(
                    x = .data$covariate_values,
                    y = after_stat(density),
                    weight = .data$IPW,
                    fill = .data[[treatment_name]]
                  ),
                  alpha = 0.5,
                  position = "identity",
                  binwidth = cd_bar_width
                ) +
                ggplot2::scale_x_continuous(
                  breaks = seq(
                    floor_dec(
                      min(covariate_values),
                      digits = decimalplaces(cd_x_scale_width)
                    ),
                    ceiling_dec(
                      max(covariate_values),
                      digits = decimalplaces(cd_x_scale_width)
                    ),
                    cd_x_scale_width
                  )
                )
            }
            p <- p +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              cd_scale_fill
          }
        )
        cd_plot_unadjusted <- cowplot::ggdraw(
          gridExtra::arrangeGrob(
            patchwork::patchworkGrob(
              (
                Reduce("+", cov_plots_unadjusted) +
                  patchwork::plot_layout(
                    nrow = cd_nrow,
                    ncol = cd_ncol,
                    guides = "collect"
                  ) &
                  ggplot2::theme(legend.position = "top")
              ) +
                patchwork::plot_annotation(title = "Covariate plots (before adjustment)")
            ),
            left = "density/count", 
            bottom = "covariate values"
          )
        ) +
          ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
        cd_plot_adjusted <- cowplot::ggdraw(
          gridExtra::arrangeGrob(
            patchwork::patchworkGrob(
              (
                Reduce("+", cov_plots_adjusted) +
                  patchwork::plot_layout(
                    nrow = cd_nrow,
                    ncol = cd_ncol,
                    guides = "collect"
                  ) &
                  ggplot2::theme(legend.position = "top")
              ) +
                patchwork::plot_annotation(title = "Covariate plots (after adjustment)")
            ),
            left = "density/count", 
            bottom = "covariate values"
          )
        ) +
          ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", color = NA)
          )
      }
    }
    
    if (ecdf) {
      if (!plots_only) {
        X_type <- X_orig |>
          dplyr::summarise(
            dplyr::across(
              dplyr::everything(), \(x) ifelse(is.factor(x), "dis", "con")
            )
          )
      }
      
      if (!plot_data_only) {
        if (is.null(ec_scale_color)) {
          if (requireNamespace("ggsci", quietly = TRUE)) {
            ec_scale_color <- ggsci::scale_color_jama()
          } else {
            ec_scale_color <- ggplot2::scale_color_discrete()
          }
        }
        if (is.null(ec_x_scale_width)) {
          nm_num <- names(X_type)
          if (exists("which_fct")) {
            nm_num <- nm_num[!(names(X_type) %in% names(names)[which_fct])]
          }
          if (length(nm_num) == 0) {
            ec_x_scale_width <- vector("numeric", if(is.null(dim(X_type))) length(X_type) else ncol(X_type))
          } else {
            scale_width_dt <- data.table::as.data.table(plot_data)
            scale_width_dt <- scale_width_dt[
              ,
              (nm_num) := lapply(.SD, \(x) signif(diff(range(x)), 1) / 5),
              .SDcols = nm_num
            ]
            ec_x_scale_width <- suppressWarnings(
              as.numeric(scale_width_dt[1])[seq_len(if(is.null(dim(X_type))) length(X_type) else ncol(X_type))]
            )
          }
        } else if (length(ec_x_scale_width) == 1) {
          ec_x_scale_width <- rep(
            ec_x_scale_width,
            if(is.null(dim(X_type))) length(X_type) else ncol(X_type)
          )
        }
        # initialise mean and max eCDF table
        if (bal_tab) {
          balance_table$unadjusted[["eCDF_mean"]] <- NA
          balance_table$unadjusted[["eCDF_max"]] <- NA
          balance_table$adjusted[["eCDF_mean"]] <- NA
          balance_table$adjusted[["eCDF_max"]] <- NA
        }
        # generate eCDF plots
        ec_plots <- purrr::pmap(
          list(
            names = names(X_type),
            ec_x_scale_width = ec_x_scale_width,
            X_type = X_type
          ),
          \(names, ec_x_scale_width, X_type) {
            plot_data <- plot_data |>
              dplyr::select(
                dplyr::all_of(treatment_name),
                "IPW",
                covariate_values = dplyr::all_of(names)
              ) |>
              dplyr::mutate(covariate_name = {{ names }}) |>
              dplyr::group_by(dplyr::across(dplyr::all_of(treatment_name))) |>
              dplyr::arrange(.data$covariate_values) |>
              dplyr::mutate(
                cum_pct_ori = seq_len(n()) / n(),
                cum_pct_wei = cumsum(.data$IPW) / sum(.data$IPW)
              ) |>
              dplyr::ungroup()
            ## calculate mean diff and max diff of eCDF
            if (bal_tab) {
              if (!is.factor(plot_data$covariate_values)) {
                cum_trt_ori <- cum_trt_wei <- cum_ctr_ori <- cum_ctr_wei <- 0
                cov_val <- -Inf
                covariate_values <- plot_data$covariate_values
                eCDF_mean_ori <- eCDF_mean_wei <- eCDF_max_ori <- eCDF_max_wei <- 0
                for (i in seq_len(nrow(plot_data))) {
                  if (!is.infinite(cov_val)) {
                    cum_diff_ori <- abs(cum_trt_ori - cum_ctr_ori)
                    cum_diff_wei <- abs(cum_trt_wei - cum_ctr_wei)
                    eCDF_max_ori <- max(eCDF_max_ori, cum_diff_ori)
                    eCDF_max_wei <- max(eCDF_max_wei, cum_diff_wei)
                    eCDF_mean_ori <- eCDF_mean_ori +
                      (covariate_values[i] - cov_val) * cum_diff_ori
                    eCDF_mean_wei <- eCDF_mean_wei +
                      (covariate_values[i] - cov_val) * cum_diff_wei
                  }
                  cov_val <- covariate_values[i]
                  if (plot_data[[treatment_name]][i] == 0) {
                    cum_ctr_ori <- plot_data$cum_pct_ori[i]
                    cum_ctr_wei <- plot_data$cum_pct_wei[i]
                  } else {
                    cum_trt_ori <- plot_data$cum_pct_ori[i]
                    cum_trt_wei <- plot_data$cum_pct_wei[i]
                  }
                }
                balance_table$unadjusted$eCDF_mean[
                  balance_table$unadjusted$name == names
                ] <<- eCDF_mean_ori / diff(range(covariate_values))
                balance_table$adjusted$eCDF_mean[
                  balance_table$adjusted$name == names
                ] <<- eCDF_mean_wei / diff(range(covariate_values))
                balance_table$unadjusted$eCDF_max[
                  balance_table$unadjusted$name == names
                ] <<- eCDF_max_ori
                balance_table$adjusted$eCDF_max[
                  balance_table$adjusted$name == names
                ] <<- eCDF_max_wei
              } else if (sum(grepl(names, names(X_orig_old))) == 1) {
                eCDF_data <- plot_data |>
                  dplyr::group_by(
                    .data$covariate_values,
                    dplyr::across(dplyr::all_of(treatment_name))
                  ) |>
                  dplyr::summarise(
                    eCDF_mean_ori = n(),
                    eCDF_mean_wei = sum(.data$IPW),
                    .groups = "drop"
                  ) |>
                  dplyr::group_by(dplyr::across(dplyr::all_of(treatment_name))) |>
                  dplyr::mutate(
                    eCDF_mean_ori = .data$eCDF_mean_ori / sum(.data$eCDF_mean_ori),
                    eCDF_mean_wei = .data$eCDF_mean_wei / sum(.data$eCDF_mean_wei)
                  ) |>
                  dplyr::group_by(.data$covariate_values) |>
                  dplyr::summarise(
                    dplyr::across(eCDF_mean_ori:eCDF_mean_wei, \(x) abs(diff(x)))
                  ) |>
                  dplyr::mutate(
                    covariate = names,
                    eCDF_max_ori = .data$eCDF_mean_ori,
                    eCDF_max_wei = .data$eCDF_mean_wei
                  )
                balance_table$unadjusted$eCDF_mean[
                  balance_table$unadjusted$name == names
                ] <<- eCDF_data$eCDF_mean_ori[1]
                balance_table$adjusted$eCDF_mean[
                  balance_table$adjusted$name == names
                ] <<- eCDF_data$eCDF_mean_wei[1]
                balance_table$unadjusted$eCDF_max[
                  balance_table$unadjusted$name == names
                ] <<- eCDF_data$eCDF_max_ori[1]
                balance_table$adjusted$eCDF_max[
                  balance_table$adjusted$name == names
                ] <<- eCDF_data$eCDF_max_wei[1]
              } else {
                eCDF_data <- plot_data |>
                  dplyr::group_by(
                    covariate_values,
                    dplyr::across(dplyr::all_of(treatment_name))
                  ) |>
                  dplyr::summarise(
                    eCDF_mean_ori = n(),
                    eCDF_mean_wei = sum(.data$IPW),
                    .groups = "drop"
                  ) |>
                  dplyr::group_by(dplyr::across(dplyr::all_of(treatment_name))) |>
                  dplyr::mutate(
                    eCDF_mean_ori = .data$eCDF_mean_ori / sum(.data$eCDF_mean_ori),
                    eCDF_mean_wei = .data$eCDF_mean_wei / sum(.data$eCDF_mean_wei)
                  ) |>
                  dplyr::group_by(.data$covariate_values) |>
                  dplyr::summarise(
                    dplyr::across(eCDF_mean_ori:eCDF_mean_wei, \(x) abs(diff(x)))
                  ) |>
                  dplyr::mutate(
                    covariate = paste(names, "-", covariate_values),
                    eCDF_max_ori = .data$eCDF_mean_ori,
                    eCDF_max_wei = .data$eCDF_mean_wei
                  )
                balance_table$unadjusted[
                  which(balance_table$unadjusted$name %in% eCDF_data$covariate),
                  "eCDF_mean"
                ] <<-
                  eCDF_data$eCDF_mean_ori[
                    na.omit(match(balance_table$unadjusted$name, eCDF_data$covariate))
                  ]
                balance_table$unadjusted[
                  which(balance_table$unadjusted$name %in% eCDF_data$covariate),
                  "eCDF_max"
                ] <<-
                  eCDF_data$eCDF_max_ori[
                    na.omit(match(balance_table$unadjusted$name, eCDF_data$covariate))
                  ]
                balance_table$adjusted[
                  which(balance_table$adjusted$name %in% eCDF_data$covariate),
                  "eCDF_mean"
                ] <<-
                  eCDF_data$eCDF_mean_wei[
                    na.omit(match(balance_table$unadjusted$name, eCDF_data$covariate))
                  ]
                balance_table$adjusted[
                  which(balance_table$adjusted$name %in% eCDF_data$covariate),
                  "eCDF_max"
                ] <<-
                  eCDF_data$eCDF_max_wei[
                    na.omit(match(balance_table$adjusted$name, eCDF_data$covariate))
                  ]
              }
            }
            ## plot eCDF
            if (X_type == "dis") {
              plot_data <- plot_data |>
                dplyr::mutate(
                  covariate_values_fct = .data$covariate_values,
                  covariate_values = as.numeric(.data$covariate_values)
                )
              covariate_levels <- levels(plot_data[["covariate_values_fct"]])
            }
            covariate_values <- plot_data[["covariate_values"]]
            p_ori <- plot_data |>
              ggplot2::ggplot() +
              ggplot2::facet_grid(~ .data$covariate_name) +
              ggplot2::geom_step(
                ggplot2::aes(
                  x = .data$covariate_values,
                  y = .data$cum_pct_ori,
                  color = .data[[treatment_name]]
                )
              ) +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              ec_scale_color
            p_wei <- plot_data |>
              ggplot2::ggplot() +
              ggplot2::facet_grid(~ .data$covariate_name) +
              ggplot2::geom_step(
                ggplot2::aes(
                  x = .data$covariate_values,
                  y = .data$cum_pct_wei,
                  color = .data[[treatment_name]]
                )
              ) +
              ggplot2::xlab("") +
              ggplot2::ylab("") +
              ec_scale_color
            if (X_type == "con") {
              p_ori <- p_ori +
                ggplot2::scale_x_continuous(
                  breaks = seq(
                    floor_dec(
                      min(covariate_values),
                      digits = decimalplaces(ec_x_scale_width)
                    ),
                    ceiling_dec(
                      max(covariate_values),
                      digits = decimalplaces(ec_x_scale_width)
                    ),
                    ec_x_scale_width
                  )
                )
              p_wei <- p_wei +
                ggplot2::scale_x_continuous(
                  breaks = seq(
                    floor_dec(
                      min(covariate_values),
                      digits = decimalplaces(ec_x_scale_width)
                    ),
                    ceiling_dec(
                      max(covariate_values),
                      digits = decimalplaces(ec_x_scale_width)
                    ),
                    ec_x_scale_width
                  )
                )
            } else {
              p_ori <- p_ori +
                ggplot2::scale_x_continuous(
                  breaks = seq_along(covariate_levels),
                  labels = covariate_levels
                )
              p_wei <- p_wei +
                ggplot2::scale_x_continuous(
                  breaks = seq_along(covariate_levels),
                  labels = covariate_levels
                )
            }
            return(list(ori = p_ori, wei = p_wei))
          }
        )
        ec_plots_unadjusted <- purrr::list_transpose(ec_plots)[[1]]
        ec_plots_adjusted <- purrr::list_transpose(ec_plots)[[2]]
        ecdf_plot_unadjusted <- cowplot::ggdraw(
          gridExtra::arrangeGrob(
            patchwork::patchworkGrob(
              (
                Reduce("+", ec_plots_unadjusted) +
                  patchwork::plot_layout(
                    nrow = ec_nrow,
                    ncol = ec_ncol,
                    guides = "collect"
                  ) &
                  ggplot2::theme(legend.position = "top")
              ) +
                patchwork::plot_annotation(title = "eCDF plots (before adjustment)")
            ), 
            left = "cumulative probability", 
            bottom = "covariate values"
          )
        ) +
          ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
        ecdf_plot_adjusted <- cowplot::ggdraw(
          gridExtra::arrangeGrob(
            patchwork::patchworkGrob(
              (
                Reduce("+", ec_plots_adjusted) +
                  patchwork::plot_layout(
                    nrow = ec_nrow,
                    ncol = ec_ncol,
                    guides = "collect"
                  ) &
                  ggplot2::theme(legend.position = "top")
              ) +
                patchwork::plot_annotation(title = "eCDF plots (after adjustment)")
            ), 
            left = "cumulative probability", 
            bottom = "covariate values"
          )
        ) +
          ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
      }
    }
  }
  
  # return list with plots
  out <- c()
  if (love) {
    out <- c(
      out,
      list(
        love_data = love_plot_data
      )
    )
    if (!plot_data_only) {
      out <- c(
        out,
        list(
          love = love_plot
        )
      )
    }
  }
  if (density) {
    out <- c(
      out,
      list(
        cd_data = plot_data,
        cd_type = plot_type
      )
    )
    if (!plot_data_only) {
      out <- c(
        out,
        list(
          cd_unadjusted = cd_plot_unadjusted,
          cd_adjusted = cd_plot_adjusted
        )
      )
    }
  }
  if (ecdf) {
    out <- c(
      out,
      list(
        X_type = X_type
      )
    )
    if (!plot_data_only) {
      out <- c(
        out,
        list(
          ecdf_unadjusted = ecdf_plot_unadjusted,
          ecdf_adjusted = ecdf_plot_adjusted
        )
      )
    }
  }
  if (bal_tab) {
    out <- c(
      out,
      list(
        balance_table = balance_table
      )
    )
  }
  if (plot_data_only) {
    if (exists("X_orig_old")) {
      out <- c(
        out,
        list(
          X_orig_old = X_orig_old
        )
      )
    }
  }
  return(out)
}
