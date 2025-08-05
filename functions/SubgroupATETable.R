SubgroupATETable <- function(x, y, cf, level, names_index = NULL, subset = NULL) {
  ci_names <- c(
    paste0(100 * level, "% CI - lower"),
    paste0(100 * level, "% CI - upper")
  )
  if (is.null(x) && is.null(y) && is.null(subset)) {
    ate <- grf::average_treatment_effect(cf)
    return(
      dplyr::tibble(
        subgroup = "Full population",
        n = length(cf$predictions),
        estimate = ate[1],
        std_err = ate[2],
        "{ci_names[1]}" := ate[1] - qnorm(1 - (1 - level) / 2) * ate[2],
        "{ci_names[2]}" := ate[1] + qnorm(1 - (1 - level) / 2) * ate[2]
      )
    )
  } else if (is.null(subset)) {
    if (is.list(x) && length(x) == 2) {
      subset <- cf[["X.orig"]][[y]] >= x[[1]] & cf[["X.orig"]][[y]] < x[[2]]
    } else if (is.atomic(x) && length(x) > 0) {
      subset <- cf[["X.orig"]][[y]] %in% x
    } else {
      stop("x must be a two element list with lower and upper bounds on a continuous interval ",
           "or an atomic vector of elements of y to include in the subset.")
    }
    ate <- grf::average_treatment_effect(cf, subset = subset)
  } else {
    if (length(y) != 2) {
      stop(
        glue::glue(
          "y must be a character vector of length 2, not {length(y)}, containing",
          "a column name and a name to use in the column."
        )
      )
    }
    y <- c("custom", y)
    ate <- grf::average_treatment_effect(cf, subset = subset)
  }
  tbl <- dplyr::tibble(
    n = sum(subset),
    estimate = ate[1],
    std_err = ate[2],
    "{ci_names[1]}" := ate[1] - qnorm(1 - (1 - level) / 2) * ate[2],
    "{ci_names[2]}" := ate[1] + qnorm(1 - (1 - level) / 2) * ate[2]
  )
  if (!is.null(names_index) && length(y) != 2) {
    subgroup <- names_index[[2]][names_index[[1]] == y]
    if (str_detect(y, "^X_\\d{2}$")) {
      if (str_detect(y, "X_(01|69|7[0-9]|8[012])")) {
        if (is.list(x) && length(x) == 2) {
          subgroup <- paste0(subgroup, " - ", x[[1]], "-", x[[2]])
        } else if (is.atomic(x) && length(x) > 0) {
          subgroup <- paste0(subgroup, " - ", x[1], "-", x[length(x)])
        }
      } else if (y == "X_02") {
        subgroup <- paste0(subgroup, " - ", ifelse(x == 0, "female", "male"))
      } else {
        subgroup <- paste0(subgroup, " - ", ifelse(x == 0, "no", "yes"))
      }
    }
    tbl <- dplyr::bind_cols(
      dplyr::tibble(subgroup = subgroup),
      tbl
    )
  } else {
    switch(
      EXPR = y[1],
      X_01 = {# age
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = paste0(y, " - ", x[1], "-", x[length(x)])),
          tbl
        )
      },
      X_02 = {# sex
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = ifelse(x == 0, "Sex - female", "Sex - male")), 
          tbl
        ) 
      },
      X_03_Hovedstaden = {# region_of_residence
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Region of residence - Hovedstaden"),
          tbl
        ) 
      },
      X_03_Sjælland = {# region_of_residence
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Region of residence - Sjælland"),
          tbl
        ) 
      },
      X_03_Syddanmark = {# region_of_residence
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Region of residence - Syddanmark"),
          tbl
        ) 
      },
      X_03_Nordjylland = {# region_of_residence
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "region of residence - Nordjylland"),
          tbl
        ) 
      },
      X_03_Midtjylland = {# region_of_residence
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "region of residence - Midtjylland"),
          tbl
        ) 
      },
      X_04_Q1 = {# house_hold_income
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "House hold income - Q1"),
          tbl
        ) 
      },
      X_04_Q2 = {# house_hold_income
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "House hold income - Q2"),
          tbl
        ) 
      },
      X_04_Q3 = {# house_hold_income
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "House hold income - Q3"),
          tbl
        ) 
      },
      X_04_Q4 = {# house_hold_income
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "House hold income - Q4"),
          tbl
        ) 
      },
      X_04_Q5 = {# house_hold_income
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "House hold income - Q5"),
          tbl
        ) 
      },
      `X_05_<10` = {# years_of_education
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Years of education - <10"),
          tbl
        ) 
      },
      `X_05_10-12` = {# years_of_education
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Years of education - 10-12"),
          tbl
        ) 
      },
      `X_05_13-15` = {# years_of_education
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Years of education - 13-15"),
          tbl
        ) 
      },
      `X_05_>15` = {# years_of_education
        tbl <- dplyr::bind_cols(
          dplyr::tibble(subgroup = "Years of education - >15"),
          tbl
        ) 
      },
      X_07 = {# heart_failure
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Heart failure - no", 
              "Heart failure - yes"
            )
          ),
          tbl
        ) 
      },
      X_08 = {# essential_hypertension
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Essential hypertension - no", 
              "Essential hypertension - yes"
            )
          ),
          tbl
        ) 
      },
      X_09 = {# secondary_hypertension
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Secondary hypertension - no", 
              "Secondary hypertension - yes"
            )
          ),
          tbl
        ) 
      },
      X_10 = {# ischemic_heart_disease 
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Ischemic heart disease - no", 
              "Ischemic heart disease - yes"
            )
          ),
          tbl
        ) 
      },
      X_11 = {# cerebrovascular_disease
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Cerebrovascular disease - no", 
              "Cerebrovascular disease - yes"
            )
          ),
          tbl
        ) 
      },
      X_12 = {# other_CNS_disorders
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Other CNS disorders - no", 
              "Other CNS disorders - yes"
            )
          ),
          tbl
        ) 
      },
      X_13 = {# arrythmias 
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Arrythmias - no", 
              "Arrythmias - yes"
            )
          ),
          tbl
        ) 
      },
      X_14 = {# any_malignancy
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Any malignancy - no", 
              "Any malignancy - yes"
            )
          ),
          tbl
        ) 
      },
      X_15 = {# malignancy_associated_with_hyponatremia
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Malignancy associated with hyponatremia - no", 
              "Malignancy associated with hyponatremia - yes"
            )
          ),
          tbl
        ) 
      },
      X_16 = {# other_malignancy
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Other malignancy - no", 
              "Other malignancy - yes"
            )
          ),
          tbl
        ) 
      },
      X_17 = {# renal_disorders
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Renal disorders - no", 
              "Renal disorders - yes"
            )
          ),
          tbl
        ) 
      },
      X_18 = {# liver_disease_and_peritonitis
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Liver disease and peritonitis - no", 
              "Liver disease and peritonitis - yes"
            )
          ),
          tbl
        ) 
      },
      X_19 = {# Pancreatitis
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Pancreatitis - no", 
              "Pancreatitis - yes"
            )
          ),
          tbl
        ) 
      },
      X_20 = {# COPD
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "COPD - no", 
              "COPD - yes"
            )
          ),
          tbl
        ) 
      },
      X_21 = {# diabetes
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Diabetes - no", 
              "Diabetes - yes"
            )
          ),
          tbl
        ) 
      },
      X_22 = {# dehydration
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Dehydration - no", 
              "Dehydration - yes"
            )
          ),
          tbl
        ) 
      },
      X_23 = {# frail_general_health
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Frail general health - no", 
              "Frail general health - yes"
            )
          ),
          tbl
        ) 
      },
      X_24 = {# physical_impairment
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Physical impairment - no", 
              "Physical impairment - yes"
            )
          ),
          tbl
        ) 
      },
      X_25 = {# mental_impairment
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Mental impairment - no", 
              "Mental impairment - yes"
            )
          ),
          tbl
        ) 
      },
      X_26 = {# rehabilitation_contacts
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Rehabilitation contacts - no", 
              "Rehabilitation contacts - yes"
            )
          ),
          tbl
        ) 
      },
      X_27 = {# podiatric_contacts
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Podiatric contacts - no", 
              "Podiatric contacts - yes"
            )
          ),
          tbl
        ) 
      },
      X_28 = {# alcohol_abuse
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Alcohol abuse - no", 
              "Alcohol abuse - yes"
            )
          ),
          tbl
        ) 
      },
      X_29 = {# drug_abuse
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Drug abuse - no", 
              "Drug abuse - yes"
            )
          ),
          tbl
        ) 
      },
      X_30 = {# HIV
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "HIV - no", 
              "HIV - yes"
            )
          ),
          tbl
        ) 
      },
      X_31 = {# anorexia_and_primary_polydipsia
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = ifelse(
              x == 0, 
              "Anorexia and primary polydipsia - no", 
              "Anorexia and primary polydipsia - yes"
            )
          ),
          tbl
        ) 
      },
      X_32_never = {# history_of_hyponatremia_never
        tbl <- dplyr::bind_cols(
          dplyr::tibble(
            subgroup = "History of hyponatremia - never"
          ),
          tbl
        ) 
      },
      custom = {
        tbl <- dplyr::bind_cols(
          tibble("{y[2]}" := y[3]),
          tbl
        )
      }
    )
  }
  return(tbl)
}