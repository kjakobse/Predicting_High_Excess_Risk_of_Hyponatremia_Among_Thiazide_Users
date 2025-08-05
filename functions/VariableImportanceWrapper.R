VariableImportanceWrapper <- function(csf, name_index) {
  tibble(
    variable_name = names(csf$X.orig),
    variable_importance = as.numeric(variable_importance(csf))
  ) |>
    arrange(desc(variable_importance)) |>
    mutate(
      var_name = str_sub(variable_name, 1, 4),
      variable_importance_num = variable_importance,
      variable_importance = sprintf("%.3f", variable_importance)
    ) |>
    left_join(
      name_index,
      by = c("var_name" = "short")
    ) |>
    mutate(
      variable_name = trimws(paste(full, str_sub(variable_name, 6, -1)))
    ) |>
    select(variable_name, variable_importance_num, variable_importance)
}