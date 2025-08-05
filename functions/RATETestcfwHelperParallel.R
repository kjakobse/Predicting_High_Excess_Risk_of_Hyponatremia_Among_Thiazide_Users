RATETestcfwHelperParallel <- function(cf, parallel = FALSE, n_cores = 1, seed = NULL) {
  RATETest <- function(forest,
                       priorities,
                       level = 0.95,
                       cov_type = c("continuous", "discrete"),
                       target = c("AUTOC", "QINI"),
                       q = seq(0.1, 1, by = 0.1),
                       R = 500,
                       subset = NULL,
                       debiasing.weights = NULL,
                       compliance.score = NULL,
                       num.trees.for.weights = 500) {
    if (is.null(priorities)) return(NULL)
    cov_type <- match.arg(cov_type)
    target <- match.arg(target)
    if(
      !(
        is.character(priorities) && length(priorities) == 1 && !is.null(names(forest$X.orig)) ||
        is.numeric(priorities) && length(priorities) == 1 ||
        is.numeric(priorities) && length(priorities) == length(forest$Y.orig)
      )
    ) {
      stop(
        "'priorities' must be a length one character vector, a length one ",
        "integer vector,\nor a numeric vector with the same length as ",
        "forest$Y.orig. If 'priorities' is a\ncharacter vector, forest$X.orig ",
        "must have named columns."
      )
    }
    if (length(priorities) == 1) {
      if (is.matrix(forest$X.orig)) {
        priorities <- forest[["X.orig"]][, priorities]
      } else {
        priorities <- forest[["X.orig"]][[priorities]]
      }
    }
    if (!hasArg(q) && cov_type == "discrete") {
      q <- cumsum(rev(table(priorities))) /
        length(priorities)
      if (min(q) > 0.001) {
        q <- c(0.001, q)
      }
    }
    rate <- rank_average_treatment_effect(
      forest = forest,
      priorities = priorities,
      target = target,
      q = q,
      R = R,
      subset = subset,
      debiasing.weights = debiasing.weights,
      compliance.score = compliance.score,
      num.trees.for.weights = num.trees.for.weights
    )
    confint <- rate$estimate + 
      tibble(
        estimate = 0,
        lower = qnorm(0.5 - level / 2) * rate$std.err,
        upper = qnorm(0.5 + level / 2) * rate$std.err
      )
    pval <- 2 * pnorm(-abs(rate$estimate) / rate$std.err)
    out <- c(
      rate,
      list(
        confint = confint,
        pval = pval
      )
    )
    class(out) <- "rank_average_treatment_effect"
    return(out)
  }
  
  cov_nm <- c(
    names(cf$X.orig),
    if (any(grepl("X_01", names(cf$X.orig)))) "quantile_X_01",
    if (any(grepl("X_04", names(cf$X.orig)))) "comb_X_04",
    if (any(grepl("X_05", names(cf$X.orig)))) "comb_X_05",
    if (any(grepl("X_06", names(cf$X.orig)))) "comb_X_06",
    if (any(grepl("X_32", names(cf$X.orig)))) "comb_X_32",
    if (any(grepl("X_51", names(cf$X.orig)))) "comb_X_51",
    if (any(grepl("X_52", names(cf$X.orig)))) "comb_X_52",
    if (any(grepl("X_53", names(cf$X.orig)))) "comb_X_53",
    if (any(grepl("X_54", names(cf$X.orig)))) "comb_X_54",
    if (any(grepl("X_55", names(cf$X.orig)))) "comb_X_55",
    if (any(grepl("X_56", names(cf$X.orig)))) "comb_X_56",
    if (any(grepl("X_57", names(cf$X.orig)))) "comb_X_57",
    if (any(grepl("X_58", names(cf$X.orig)))) "comb_X_58",
    if (any(grepl("X_59", names(cf$X.orig)))) "comb_X_59",
    if (any(grepl("X_60", names(cf$X.orig)))) "comb_X_60",
    if (any(grepl("X_61", names(cf$X.orig)))) "comb_X_61",
    if (any(grepl("X_62", names(cf$X.orig)))) "comb_X_62",
    if (any(grepl("X_63", names(cf$X.orig)))) "comb_X_63",
    if (any(grepl("X_64", names(cf$X.orig)))) "comb_X_64",
    if (any(grepl("X_65", names(cf$X.orig)))) "comb_X_65",
    if (any(grepl("X_66", names(cf$X.orig)))) "comb_X_66",
    if (any(grepl("X_67", names(cf$X.orig)))) "comb_X_67",
    if (any(grepl("X_68", names(cf$X.orig)))) "comb_X_68",
    if (any(grepl("X_69", names(cf$X.orig)))) "quantile_X_69",
    if (any(grepl("X_70", names(cf$X.orig)))) "quantile_X_70",
    if (any(grepl("X_71", names(cf$X.orig)))) "quantile_X_71",
    if (any(grepl("X_72", names(cf$X.orig)))) "quantile_X_72",
    if (any(grepl("X_73", names(cf$X.orig)))) "quantile_X_73",
    if (any(grepl("X_74", names(cf$X.orig)))) "quantile_X_74",
    if (any(grepl("X_75", names(cf$X.orig)))) "quantile_X_75",
    if (any(grepl("X_76", names(cf$X.orig)))) "quantile_X_76",
    if (any(grepl("X_77", names(cf$X.orig)))) "quantile_X_77",
    if (any(grepl("X_78", names(cf$X.orig)))) "quantile_X_78",
    if (any(grepl("X_79", names(cf$X.orig)))) "quantile_X_79",
    if (any(grepl("X_80", names(cf$X.orig)))) "quantile_X_80",
    if (any(grepl("X_81", names(cf$X.orig)))) "quantile_X_81",
    if (any(grepl("X_82", names(cf$X.orig)))) "quantile_X_82"
  )
  if (parallel) {
    plan(multisession, workers = n_cores)
    on.exit(plan(sequential))
  } else {
    plan(sequential)
  }
  rate_table <- future_map(
    .x = cov_nm,
    .options = furrr_options(
      packages = c("grf", "dplyr", "purrr", "stringr", "glue", "stats"),
      globals = c("DiscreteCovariatesToOneHot", "cf"),
      seed = seed
    ),
    .f = \(cov_nm) {
      if (str_detect(cov_nm, "^X_(01|69|7[0-9]|8[0-2])$")) {
        rate <- RATETest(cf, cov_nm, cov_type = "continuous")
        as_tibble(c(
          covariate = cov_nm,
          rate$confint,
          std.err = rate$std.err, 
          pval = rate$pval
        ))
      } else if (str_detect(cov_nm, "^X_(0[2-9]|[1-5]\\d|6[0-8])")) { # no $ because discrete variables have levels appended
        rate <- RATETest(cf, cov_nm, cov_type = "discrete")
        as_tibble(c(
          covariate = cov_nm,
          rate$confint,
          std.err = rate$std.err, 
          pval = rate$pval
        ))
      } else if (cov_nm == "quantile_X_01") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_01, quantile(X_01), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Age - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "comb_X_04") {
        mutate(
          cf$X.orig,
          priority = case_when(
            X_04_Q1 == 1 ~ 1, `X_04_Q2` == 1 ~ 2, `X_04_Q3` == 1 ~ 3,
            `X_04_Q4` == 1 ~ 4, X_04_Q5 == 1 ~ 5
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "House hold income", y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_05") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_05_<10` == 1 ~ 1, `X_05_10-12` == 1 ~ 2, 
            `X_05_13-15` == 1 ~ 3, `X_05_>15` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Years of education", y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_06") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_06_2018` == 1 ~ 1, `X_06_2017` == 1 ~ 2, `X_06_2016` == 1 ~ 3,
            `X_06_2015` == 1 ~ 4, `X_06_2014` == 1 ~ 5
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Calendar year of inclusion", y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_32") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_32_never` == 1 ~ 1, `X_32_<4 months` == 1 ~ 2, 
            `X_32_≥4 months` == 1 ~ 3
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "History of hyponatremia", y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_51") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_51_0` == 1 ~ 1, `X_51_1-3` == 1 ~ 2, `X_51_4-6` == 1 ~ 3,
            `X_51_7-9` == 1 ~ 4, `X_51_≥10` == 1 ~ 5
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) {
            as_tibble(c(
              covariate = "No. of different prescription drugs used", 
              y$confint, 
              std.err = y$std.err, 
              pval = y$pval
            ))
          })()
      } else if (cov_nm == "comb_X_52") {
        mutate(
          cf$X.orig,
          priority = case_when(
            X_52_0 == 1 ~ 1, `X_52_1-7` == 1 ~ 2, `X_52_8-14` == 1 ~ 3,
            `X_52_≥15` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) {
            as_tibble(c(
              covariate = "Days of hospitalization in the past year prior to index data", 
              y$confint, 
              std.err = y$std.err, 
              pval = y$pval
            ))
          })()
      } else if (cov_nm == "comb_X_53") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_53_0` == 1 ~ 1, `X_53_1-3` == 1 ~ 2,
            `X_53_4-6` == 1 ~ 3, `X_53_≥7` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) {
            as_tibble(c(
              covariate = "No. of outpatient contacts in the last year prior to index data", 
              y$confint, 
              std.err = y$std.err, 
              pval = y$pval
            ))
          })()
      } else if (cov_nm == "comb_X_54") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_54_0` == 1 ~ 1, `X_54_1-3` == 1 ~ 2,
            `X_54_4-6` == 1 ~ 3, `X_54_≥7` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) {
            as_tibble(c(
              covariate = "No. of primary care contacts in the last 2 years", 
              y$confint,
              std.err = y$std.err, 
              pval = y$pval
            ))
          })()
      } else if (cov_nm == "comb_X_55") {
        # set NA level for blood test measurements to have average priority
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_55_NA` == 1 ~ 2.5, `X_55_Q1` == 1 ~ 1,
            `X_55_Q2` == 1 ~ 2, `X_55_Q3` == 1 ~ 3,
            `X_55_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Sodium",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_56") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_56_NA` == 1 ~ 2.5, `X_56_Q1` == 1 ~ 1,
            `X_56_Q2` == 1 ~ 2, `X_56_Q3` == 1 ~ 3,
            `X_56_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "eGRF",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_57") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_57_NA` == 1 ~ 2.5, `X_57_Q1` == 1 ~ 1,
            `X_57_Q2` == 1 ~ 2, `X_57_Q3` == 1 ~ 3,
            `X_57_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Potassium",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_58") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_58_NA` == 1 ~ 2.5, `X_58_Q1` == 1 ~ 1,
            `X_58_Q2` == 1 ~ 2, `X_58_Q3` == 1 ~ 3,
            `X_58_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Hemoglobin",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_59") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_59_NA` == 1 ~ 2.5, `X_59_Q1` == 1 ~ 1,
            `X_59_Q2` == 1 ~ 2, `X_59_Q3` == 1 ~ 3,
            `X_59_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "ALAT",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_60") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_60_NA` == 1 ~ 2.5, `X_60_Q1` == 1 ~ 1,
            `X_60_Q2` == 1 ~ 2, `X_60_Q3` == 1 ~ 3,
            `X_60_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Alkaline phosphatase",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_61") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_61_NA` == 1 ~ 2.5, `X_61_Q1` == 1 ~ 1,
            `X_61_Q2` == 1 ~ 2, `X_61_Q3` == 1 ~ 3,
            `X_61_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Albumin",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_62") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_62_NA` == 1 ~ 2.5, `X_62_Q1` == 1 ~ 1,
            `X_62_Q2` == 1 ~ 2, `X_62_Q3` == 1 ~ 3,
            `X_62_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "LDH",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_63") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_63_NA` == 1 ~ 2.5, `X_63_Q1` == 1 ~ 1,
            `X_63_Q2` == 1 ~ 2, `X_63_Q3` == 1 ~ 3,
            `X_63_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Carbamide",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_64") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_64_NA` == 1 ~ 2.5, `X_64_Q1` == 1 ~ 1,
            `X_64_Q2` == 1 ~ 2, `X_64_Q3` == 1 ~ 3,
            `X_64_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Thrombocytes",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_65") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_65_NA` == 1 ~ 2.5, `X_65_Q1` == 1 ~ 1,
            `X_65_Q2` == 1 ~ 2, `X_65_Q3` == 1 ~ 3,
            `X_65_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Leukocytes",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_66") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_66_NA` == 1 ~ 2.5, `X_66_Q1` == 1 ~ 1,
            `X_66_Q2` == 1 ~ 2, `X_66_Q3` == 1 ~ 3,
            `X_66_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "CRP",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_67") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_67_NA` == 1 ~ 2.5, `X_67_Q1` == 1 ~ 1,
            `X_67_Q2` == 1 ~ 2, `X_67_Q3` == 1 ~ 3,
            `X_67_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "Cholesterol",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "comb_X_68") {
        mutate(
          cf$X.orig,
          priority = case_when(
            `X_68_NA` == 1 ~ 2.5, `X_68_Q1` == 1 ~ 1,
            `X_68_Q2` == 1 ~ 2, `X_68_Q3` == 1 ~ 3,
            `X_68_Q4` == 1 ~ 4
          ),
          .keep = "none"
        ) |>
          (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
          (\(y) as_tibble(c(covariate = "LDL-C",  y$confint, std.err = y$std.err, pval = y$pval)))()
      } else if (cov_nm == "quantile_X_69") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_69, quantile(X_69), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Sodium - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_70") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_70, quantile(X_70), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("eGRF - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_71") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_71, quantile(X_71), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Potassium - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_72") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_72, quantile(X_72), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Hemoglobin - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_73") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_73, quantile(X_73), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("ALAT - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_74") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_74, quantile(X_74), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Alkaline phosphatase - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_75") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_75, quantile(X_75), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Albumin - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_76") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_76, quantile(X_76), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Lactate dehydrogenase - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_77") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_77, quantile(X_77), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Carbamide - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_78") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_78, quantile(X_78), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Thrombocytes - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_79") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_79, quantile(X_79), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Leukocytes - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_80") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_80, quantile(X_80), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("C-reactive protein - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_81") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_81, quantile(X_81), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("Cholesterol - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else if (cov_nm == "quantile_X_82") {
        cf$X.orig |>
          transmute(
            quartile = cut(X_82, quantile(X_82), labels = paste0("Q", 1:4), include.lowest = TRUE)
          ) |>
          DiscreteCovariatesToOneHot() |>
          (\(y) 
           {
             map2(
               y,
               str_extract(names(y), "Q\\d"),
               \(priority, nm) {
                 RATETest(cf, priority, cov_type = "discrete") |>
                   (\(rate) as_tibble(c(covariate = paste0("LDL-C - ", nm),  rate$confint, std.err = rate$std.err, pval = rate$pval)))()
               }
             ) |>
               list_rbind()
          })()
      } else {
        stop("cov_nm incorrect value: ", cov_nm)
      }
    }
  )
  
  return(
    rate_table |>
      list_rbind() |>
      arrange(pval) |>
      mutate(
        pval_bonf = p.adjust(pval, "bonferroni"),
        pval_bh = p.adjust(pval, "BH")
      )
  )
}