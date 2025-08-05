MapCovariateBalance <- function(data, 
                                nm_ind, 
                                nm = dplyr::everything(), 
                                plot_data_only = FALSE, 
                                plots_only = FALSE,
                                ...) {
  if (plots_only) {
    data_list <- data
    grepl_string <- data_list$grepl_string
    names_string <- data_list$names_string
  } else {
    grepl_string <- data$X.orig |>
      dplyr::select(dplyr::all_of({{ nm }})) |>
      names() |>
      stringr::str_subset("^X_\\d{2}_", negate = TRUE) |>
      (\(x) {
        suppressWarnings(
          split(
            x, 
            factor(
              rep(
                seq_len(length(x) %/% 16 + 1), 
                each = ceiling(length(x) / (length(x) %/% 16 + 1))
              )
            )
          )
        )
      })() |>
      purrr::map_chr(\(string) paste0(paste0("^", string, "$"), collapse = "|"))
    
    names_string <- data$X.orig |>
      dplyr::select(dplyr::all_of({{ nm }})) |>
      names() |>
      stringr::str_subset("^X_\\d{2}_") |>
      stringr::str_extract("^X_\\d{2}") |>
      unique() |>
      (\(names) {
        nm_ind |>
          dplyr::filter(.data$short %in% .env$names) |>
          (\(df) {
            out <- df$short
            names(out) <- df$full
            return(out)
          })()
      })() |>
      (\(x) if (length(x) > 0) list(x) else x)()
  }
  
  out <- list()
  
  for (i in seq_along(grepl_string)) {
    if (plots_only) {
      data <- data_list[[paste0("numeric", i)]]
      if (!is.null(data$X_orig_old)) {
        data$X.orig <- data$X_orig_old
        names(data$X.orig) <- nm_ind |>
          dplyr::filter(full %in% names(data$X_orig_old)) |>
          dplyr::select("short", "full") |>
          dplyr::left_join(
            dplyr::tibble(
              full = names(data$X_orig_old),
              id = seq_along(names(data$X_orig_old))
            ),
            by = "full"
          ) |>
          dplyr::arrange(id) |>
          dplyr::pull("short")
      }
    }
    names <- nm_ind |>
      dplyr::filter(stringr::str_detect(short, "^X")) |>
      tidyr::unnest(levels, keep_empty = TRUE) |>
      dplyr::mutate(
        out = ifelse(is.na(levels), short, paste0(short, "_", levels)),
        nm = ifelse(is.na(levels), full, paste0(full, " - ", levels)),
        .keep = "none"
      ) |>
      (\(df) structure(df$out, names = df$nm))() |>
      (\(x) subset(x, grepl(grepl_string[i], x)))()
    
    factor <- dplyr::select(data$X.orig, dplyr::all_of(names)) |>
      dplyr::select(dplyr::where(is.logical)) |>
      purrr::map(~list("FALSE" = FALSE, "TRUE" = TRUE))
    factor <- c(
      factor,
      dplyr::select(data$X.orig, dplyr::all_of(names)) |>
        dplyr::select(dplyr::where(\(x) is.numeric(x) & all(unique(x) %in% c(0, 1)))) |>
        purrr::map(~list("FALSE" = 0, "TRUE" = 1))
    )
    out[[paste0("numeric", i)]] <- CovariateBalance(
      data,
      covariates = character(0),
      names = names,
      factor = factor,
      cd_ncol = 4,
      plot_data_only = plot_data_only,
      plots_only = plots_only,
      ...
    )
  }
  
  for (i in seq_along(names_string)) {
    if (plots_only) data <- data_list[[paste0("discrete", i)]]
    out[[paste0("discrete", i)]] <- CovariateBalance(
      data, 
      covariates = character(0),
      names = names_string[[i]],
      cd_ncol = 3, 
      plot_data_only = plot_data_only,
      plots_only = plots_only,
      ...
    )
  }
  
  if (plot_data_only) {
    out$grepl_string <- grepl_string
    out$names_string <- names_string
  }
  return(out)
}