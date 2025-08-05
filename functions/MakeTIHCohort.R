MakeTIHCohort <- function(data, 
                          h = 120,
                          test_prop = 0, 
                          seed = NULL) {
  stopifnot(
    "test_prop must be between 0 and 1." = 
      test_prop >= 0 & test_prop <= 1
  )
  df <- data |>
    mutate(
      pnr = pnr,
      # time of event (either censoring or hyponatremia). Outcome used with
      # causal survival forest
      Y_130 = t_130,
      # event indicator (0 = censoring, 1 = hyponatremia). Used with causal
      # survival forest
      D_130 = event_130,
      Y_125 = t_125,
      D_125 = event_125,
      .keep = "none"
    )
  for (hh in h) {
    df <- df |>
      mutate(
        # indicator of hyponatremia before horizon. Outcome used with censoring
        # weighted causal forest. Only considers observed events, defined as
        # hyponatremia status at horizon. Observations censored before the horizon
        # are removed (represented by NA),
        "Y_130_{hh}" := (\() {
          out <- data.table::fcase(
            Y_130 > hh, 0L,
            Y_130 <= hh & D_130 == 1, 1L,
            default = NA_integer_
          )
          attr(out, "label") <- paste("Hyponatremia (<130 mM) event inside", hh, "days")
          return(out)
        })(),
        # restricted event indicator (1 = hyponatremia or survived until horizon) 
        "D_130_{hh}" := (\() {
          out <- pmax(D_130, Y_130 >= hh)
          attr(out, "label") <- paste("Hyponatremia (<130 mM) event or >", hh, "days")
          return(out)
        })(),
        "Y_125_{hh}" := (\() {
          out <- data.table::fcase(
            Y_125 > hh, 0L,
            Y_125 <= hh & D_125 == 1, 1L,
            default = NA_integer_
          )
          attr(out, "label") <- paste("Hyponatremia (<125 mM) event inside", hh, "days")
          return(out)
        })(),
        # restricted event indicator (1 = hyponatremia or survived until horizon) 
        "D_125_{hh}" := (\() {
          out <- pmax(D_125, Y_125 >= hh)
          attr(out, "label") <- paste("Hyponatremia (<125 mM) event or >", hh, "days")
          return(out)
        })()
      )
  }
  df <- bind_cols(
    df,
    data |>
      mutate(
        # treatment indicator BFZ vs. CCB
        W_bvc = (\() {
          out <- data.table::fcase(
            indexdrug == "ccb", 0L,
            indexdrug == "bfz", 1L,
            default = NA_integer_
          )
          attr(out, "label") <- "BFZ vs. CCB"
          return(out)
        })(),
        # treatment indicator HCTZ vs. RAS
        W_hvr = (\() {
          out <- data.table::fcase(
            indexdrug == "ras", 0L,
            indexdrug == "hctz", 1L,
            default = NA_integer_
          )
          attr(out, "label") <- "HCTZ vs. RAS"
          return(out)
        })(),
        # treatment indicator Tiazide (BFZ+HCTZ) vs. non-tiazide (CCB+RAS)
        W = (\() {
          out <- data.table::fcase(
            indexdrug %in% c("ccb", "ras"), 0L,
            indexdrug %in% c("bfz", "hctz"), 1L,
            default = NA_integer_
          )
          attr(out, "label") <- "tiazide vs. non-tiazide"
          return(out)
        })(),
        # Age at index date
        X_01 = age,
        # Sex (female = 0, male = 1)
        X_02 = (\() {
          out <- ifelse(sex == "Female", 0, 1)
          attr(out, "label") <- "Sex"
          return(out)
        })(),
        # region of residence (middle Denmark excluded)
        X_03 = (\() {
          out <- factor(
            str_sub(str_replace_all(reg, "\UFFFD", "æ"), 3, -1),
            levels = c("Hovedstaden", "Nordjylland", "Sjælland", "Syddanmark")
          )
          attr(out, "label") <- "Region of residence"
          return(out)
        })(),
        # household income quintiles
        X_04 = (\() {
          out <- factor(
            ifelse(inc_missing == "missing", NA, inc), 
            levels = 1:5, 
            labels = paste0("Q", 1:5),
            ordered = TRUE
          )
          attr(out, "label") <- "Income quintile"
          return(out)
        })(),
        # years of education
        X_05 = (\() {
          out <- factor(
            data.table::fcase(
              edu_missing == "missing", NA_character_,
              str_detect(edu, "01"), "<10",
              str_detect(edu, "02"), "10-12",
              str_detect(edu, "03"), "13-15",
              str_detect(edu, "04"), ">15",
              default = NA_character_
            ),
            levels = c("<10", "10-12", "13-15", ">15"),
            ordered = TRUE
          )
          attr(out, "label") <- "Education, years"
          return(out)
        })(),
        # calendar year of inclusion
        X_06 = (\() {
          out <- factor(
            data.table::fcase(
              y01 > 0.1, "2014",
              y02 > 0.1, "2015",
              y03 > 0.1, "2016",
              y04 > 0.1, "2017",
              y05 > 0.1, "2018",
              y06 > 0.1, "2019",
              y07 > 0.1, "2020",
              default = NA_character_
            ),
            levels = c("2014", "2015", "2016", "2017", "2018", "2019", "2020")
          )
          attr(out, "label") <- "Calendar year of inclusion"
          return(out)
        })(),
        # Heart failure
        X_07 = d01,
        # Essential hypertension
        X_08 = d02,
        # Secondary hypertension
        X_09 = d03,
        # Ischemic heart disease
        X_10 = d04,
        # Cerebrovascular disease
        X_11 = d05,
        # Other CNS disorders
        X_12 = d06,
        # Arrhythmias
        X_13 = d07,
        # Any malignancy,
        X_14 = d08,
        # Malignancy associated with hyponatremia
        X_15 = d09,
        # Other malignancy
        X_16 = d10,
        # Renal disorders
        X_17 = d11,
        # Liver disease and peritonitis
        X_18 = d12,
        # Pancreatitis
        X_19 = d13,
        # Chronic obstructive airway disorder
        X_20 = d14,
        # Diabetes
        X_21 = d15,
        # Dehydration
        X_22 = d16,
        # Frail general health
        X_23 = d17,
        # Physical impairment
        X_24 = d18,
        # Mental impairment
        X_25 = d19,
        # Rehabilitation contacts
        X_26 = d20,
        # Podiatric contacts
        X_27 = d21,
        # Alcohol abuse
        X_28 = d22,
        # Drug abuse
        X_29 = d23,
        # HIV
        X_30 = d24,
        # Anorexia and primary polydipsia
        X_31 = d25,
        # History of hyponatremia
        X_32 = (\() {
          out <- factor(
            data.table::fcase(
              d26 == "never", "never",
              d26 == "lt_4mths", "<4 months",
              d26 == "mt_4mths", "\u22654 months",
              default = NA_character_
            ),
            # never is longer since last hyponatremia event than any recorded event
            levels = c("<4 months", "\u22654 months", "never") 
          )
          attr(out, "label") <- "History of hyponatremia"
          return(out)
        })(),
        # PPI and anticides
        X_33 = m01,
        # Obstipation/diarrhea
        X_34 = m02,
        # Antidiabetics (not insulin)
        X_35 = m03,
        # Insulin
        X_36 = m04,
        # Oral anticoagulants
        X_37 = m05,
        # Aspirin
        X_38 = m06,
        # ADPi
        X_39 = m07,
        # Nitrates
        X_40 = m08,
        # Beta-blockers
        X_41 = m09,
        # Lipid lowering drugs
        X_42 = m10,
        # Desmopressin
        X_43 = m11,
        # Eltroxin
        X_44 = m12,
        # NSAIDs
        X_45 = m13,
        # Antiepileptics
        X_46 = m14,
        # Opioids
        X_47 = m15,
        # Antidepressives
        X_48 = m16,
        # Antipsychotics
        X_49 = m17,
        # Agents for COPD
        X_50 = m18,
        # No. of different prescription drugs used (in past year to index date)
        X_51 = (\() {
          out <- factor(
            data.table::fcase(
              n_drugs == 0, "0",
              n_drugs %in% 1:3, "1-3",
              n_drugs %in% 4:6, "4-6",
              n_drugs %in% 7:9, "7-9",
              n_drugs > 9, paste0("\U2265", 10),
              default =  NA_character_
            ),
            levels = c("0", "1-3", "4-6", "7-9", paste0("\U2265", "10"))
          )
          attr(out, "label") <- "No. of different prescription drugs used"
          return(out)
        })(),
        # Days of hospitalization (in past year from index date)
        X_52 = h01,
        # No. of outpatient contacts (in past year from index date)
        X_53 = h02,
        # No of primary care contacts (in past 2 years from index date)
        X_54 = h03,
        # Sodium
        X_69 = l01,
        X_69_tbi = indexdate - date_l01, # time between measurement and index date
        X_55 = (\() {
          out <- cut(
            x = l01, 
            breaks = quantile(l01, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Sodium"
          return(out)
        })(),
        # eGRF
        X_70 = structure(egfr, label = "eGFR"),
        X_70_tbi = indexdate - date_egfr,
        X_56 = (\() {
          out <- cut(
            x = egfr, 
            breaks = quantile(egfr, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "eGFR"
          return(out)
        })(),
        # Potassium
        X_71 = l03,
        X_71_tbi = indexdate - date_l03,
        X_57 = (\() {
          out <- cut(
            x = l03, 
            breaks = quantile(l03, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Potassium"
          return(out)
        })(),
        # Hemoglobin
        X_72 = l04,
        X_72_tbi = indexdate - date_l04,
        X_58 = (\() {
          out <- cut(
            x = l04, 
            breaks = quantile(l04, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Hemoglobin"
          return(out)
        })(),
        # ALAT
        X_73 = l05,
        X_73_tbi = indexdate - date_l05,
        X_59 = (\() {
          out <- cut(
            x = l05, 
            breaks = quantile(l05, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "ALAT"
          return(out)
        })(),
        # Alkaline phosphatase
        X_74 = l06,
        X_74_tbi = indexdate - date_l06,
        X_60 = (\() {
          out <- cut(
            x = l06, 
            breaks = quantile(l06, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Alkaline phosphatase"
          return(out)
        })(),
        # Albumin
        X_75 = l07,
        X_75_tbi = indexdate - date_l07,
        X_61 = (\() {
          out <- cut(
            x = l07, 
            breaks = quantile(l07, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Albumin"
          return(out)
        })(),
        # LDH
        X_76 = l08,
        X_76_tbi = indexdate - date_l08,
        X_62 = (\() {
          out <- cut(
            x = l08, 
            breaks = quantile(l08, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "LDH"
          return(out)
        })(),
        # Carbamide
        X_77 = l09,
        X_77_tbi = indexdate - date_l09,
        X_63 = (\() {
          out <- cut(
            x = l09, 
            breaks = quantile(l09, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Carbamide"
          return(out)
        })(),
        # Thrombocytes
        X_78 = l10,
        X_78_tbi = indexdate - date_l10,
        X_64 = (\() {
          out <- cut(
            x = l10, 
            breaks = quantile(l10, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Thrombocytes"
          return(out)
        })(),
        # Leukocytes
        X_79 = l11,
        X_79_tbi = indexdate - date_l11,
        X_65 = (\() {
          out <- cut(
            x = l11, 
            breaks = quantile(l11, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Leukocytes"
          return(out)
        })(),
        # CRP
        X_80 = l12,
        X_80_tbi = indexdate - date_l12,
        X_66 = (\() {
          out <- cut(
            x = l12, 
            breaks = quantile(l12, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "CRP"
          return(out)
        })(),
        # Cholesterol
        X_81 = l13,
        X_81_tbi = indexdate - date_l13,
        X_67 = (\() {
          out <- cut(
            x = l13, 
            breaks = quantile(l13, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "Cholesterol"
          return(out)
        })(),
        # LDL-C
        X_82 = l14,
        X_82_tbi = indexdate - date_l14,
        X_68 = (\() {
          out <- cut(
            x = l14, 
            breaks = quantile(l14, seq(0, 1, 0.25), na.rm = TRUE),
            labels = paste0("Q", 1:4)
          ) |>
            addNA()
          attr(out, "label") <- "LDL-C"
          return(out)
        })(),
        .keep = "none"
      )
  )
  if (test_prop > 0 & test_prop < 1) {
    if(!is.null(seed)) set.seed(seed)
    test_index <- sample(seq_len(nrow(df)), floor(test_prop * nrow(df)))
    test <- df[test_index,]
    train <- df[-test_index,]
    return(list(original_train = train, original_test = test))
  }
  return(df)
}