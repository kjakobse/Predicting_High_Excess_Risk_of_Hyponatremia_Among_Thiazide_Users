RATETestcfwHelper <- function(cf) {
  map(
    names(cf$X.orig),
    \(cov_nm) {
      if (str_detect(cov_nm, "^X_(01|69|7[0-9]|8[0-2])")) {
        rate <- RATETest(cf, cov_nm, cov_type = "continuous")
      } else {
        rate <- RATETest(cf, cov_nm, cov_type = "discrete")
      }
      print(glue::glue("Time: {Sys.time()}, covariate: {cov_nm}"))
      as_tibble(c(
        covariate = cov_nm,
        rate$confint,
        pval = rate$pval
      ))
    }
  ) |> 
    bind_rows() |>
    (\(x) {
      if (any(grepl("X_04", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
          mutate(
            cf$X.orig,
            priority = case_when(
              X_04_Q1 == 1 ~ 1, `X_04_Q2` == 1 ~ 2, `X_04_Q3` == 1 ~ 3,
              `X_04_Q4` == 1 ~ 4, X_04_Q5 == 1 ~ 5
            ),
            .keep = "none"
          ) |>
            (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
            (\(y) as_tibble(c(covariate = "House hold income", y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: House hold income"))
        return(out)
      } else {
        return(x)
      }
    })() |>
    (\(x) {
      if (any(grepl("X_05", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
          mutate(
            cf$X.orig,
            priority = case_when(
              `X_05_<10` == 1 ~ 1, `X_05_10-12` == 1 ~ 2, 
              `X_05_13-15` == 1 ~ 3, `X_05_>15` == 1 ~ 4
            ),
            .keep = "none"
          ) |>
            (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
            (\(y) as_tibble(c(covariate = "Years of education", y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Years of education"))
        return(out)
      } else {
        return(x)
      }
    })() |>
    (\(x) {
      if (any(grepl("X_06", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
          mutate(
            cf$X.orig,
            priority = case_when(
              `X_06_2018` == 1 ~ 1, `X_06_2017` == 1 ~ 2, `X_06_2016` == 1 ~ 3,
              `X_06_2015` == 1 ~ 4, `X_06_2014` == 1 ~ 5
            ),
            .keep = "none"
          ) |>
            (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
            (\(y) as_tibble(c(covariate = "Calendar year of inclusion", y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Calendar year of inclusion"))
        return(out)
      } else {
        return(x)
      } 
    })()|>
    (\(x) {
      if (any(grepl("X_32", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
          mutate(
            cf$X.orig,
            priority = case_when(
              `X_32_never` == 1 ~ 1, `X_32_<4 months` == 1 ~ 2, 
              `X_32_≥4 months` == 1 ~ 3
            ),
            .keep = "none"
          ) |>
            (\(y) RATETest(cf, y$priority, cov_type = "discrete"))() |>
            (\(y) as_tibble(c(covariate = "History of hyponatremia", y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: History of hyponatremia"))
        return(out)
      } else {
        return(x)
      } 
    })()|>
    (\(x) {
      if (any(grepl("X_51", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
                pval = y$pval
              ))
            })()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: No. of different prescription drugs used"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_52", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
                pval = y$pval
              ))
            })()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Days of hospitalization in the past year prior to index data"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_53", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
                pval = y$pval
              ))
            })()
        ) 
        print(glue::glue("Time: {Sys.time()}, covariate: No. of outpatient contacts in the last year prior to index data"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_54", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
                pval = y$pval
              ))
            })()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: No. of primary care contacts in the last 2 years"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    arrange(pval) |>
    mutate(
      pval_bonf = p.adjust(pval, "bonferroni"),
      pval_bh = p.adjust(pval, "BH")
    ) |>
    # set NA level for blood test measurements to have average priority
    (\(x) {
      if (any(grepl("X_55", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Sodium",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Sodium"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_56", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "eGRF",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: eGRF"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_57", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Potassium",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Potassium"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_58", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Hemoglobin",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Hemoglobin"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_59", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "ALAT",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: ALAT"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_60", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Alkaline phosphatase",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Alkaline phosphatase"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_61", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Albumin",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Albumin"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_62", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "LDH",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: LDH"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_63", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Carbamide",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Carbamide"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_64", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Thrombocytes",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Thrombocytes"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_65", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Leukocytes",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Leukocytes"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_66", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "CRP",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: CRP"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_67", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "Cholesterol",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: Cholesterol"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    (\(x) {
      if (any(grepl("X_68", names(cf$X.orig)))) {
        out <- bind_rows(
          x,
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
            (\(y) as_tibble(c(covariate = "LDL-C",  y$confint, pval = y$pval)))()
        )
        print(glue::glue("Time: {Sys.time()}, covariate: LDL-C"))
        return(out)
      } else {
        return(x)
      } 
    })() |>
    arrange(pval) |>
    mutate(
      pval_bonf = p.adjust(pval, "bonferroni"),
      pval_bh = p.adjust(pval, "BH")
    )
}