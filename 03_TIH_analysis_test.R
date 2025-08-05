# ===   1. Setup   ================================================================================================ ====
## load and attach packages ....................................................................................... ----
library(tidyverse)
library(rlang)
library(lattice)
library(glue)
library(grf)
library(ggpubr)
library(ggsci)
library(MatchIt)
library(mice)
library(furrr)
library(future.callr)
library(missForest)
library(xlsx2dfs)
library(tableone)
library(DescTools)
library(ggrepel)
library(cowplot)

## source helper functions ........................................................................................ ----
source("functions/ATEAll.R", encoding = 'UTF-8')
source("functions/CATEPlot.R", encoding = 'UTF-8')
source("functions/CATEPlots.R", encoding = 'UTF-8')
source("functions/CausalForestATEAllSubgroupTable.R", encoding = 'UTF-8')
source("functions/CausalForestCATEAllTable.R", encoding = 'UTF-8')
source("functions/CausalForestDynamicSubgroups.R", encoding = 'UTF-8')
source("functions/CForBenefit.R", encoding = 'UTF-8')
source("functions/CovariateBalance.R", encoding = 'UTF-8')
source("functions/DiscreteCovariatesToOneHot.R", encoding = 'UTF-8')
source("functions/futuremice.R", encoding = 'UTF-8')
source("functions/GRFAnalysisWrapper.R", encoding = 'UTF-8')
source("functions/MakeTIHCohort.R", encoding = 'UTF-8')
source("functions/MapCovariateBalance.R", encoding = 'UTF-8')
source("functions/mbcal.R", encoding = 'UTF-8')
source("functions/OneHotToFactor.R", encoding = 'UTF-8')
source("functions/OverlapPlot.R", encoding = 'UTF-8')
source("functions/predtools_calibration_plot.R", encoding = 'UTF-8')
source("functions/PrepareTIHCovariates.R", encoding = 'UTF-8')
source("functions/RATETest.R", encoding = 'UTF-8')
source("functions/RATETestcfwHelper.R", encoding = 'UTF-8')
source("functions/RATETestcfwHelperParallel.R", encoding = 'UTF-8')
source("functions/SubgroupATETable.R", encoding = 'UTF-8')
source("functions/TestCalibrationData.R", encoding = 'UTF-8')
source("functions/VariableImportanceWrapper.R", encoding = 'UTF-8')
source("functions/vcovHC.R", encoding = 'UTF-8')

Rcpp::sourceCpp("functions/mbcb.cpp")
Rcpp::sourceCpp("functions/rcorr.cpp")

## set ggplot theme ............................................................................................... ----
theme_set(theme_pubr())

## set training flags ............................................................................................. ----
# should data be imputed:
impute_data <- TRUE
# should causal forest training happen:
# Note that training only progresses if no existing trained model object exists in the provided folder
train_cf <- TRUE
# should results based on trained models be computed:
compute_results <- TRUE
# should forest objects used to predict sample weights, expected outcome, and propensity to treat be saved on disc:
save_sampleweight_forest <- FALSE
save_exposure_forest <- FALSE
save_outcome_forest <- FALSE
# should models on folds to dynamically determine subgroups with differential CATE estimates be trained:
dynamic_subgroups <- TRUE
# Which setups should be run:
setup_index <- 1:10
# Which number of most most important variables by variable importance should be run:
vi_models <- c(1:5, 20)
# cutoffs for excess risk strategy
excess_risk_cutoff <- c(0, 0.01, 0.05)
excess_risk_quantiles <- c(0, 0.5, 0.9, 0.95)
# how many folds and rank groups to use for dynamic subgroups:
n_folds <- 5L
n_rankings <- c(5L, 10L, 20L)
# how many threads and clusters to use during training:
num_clusters <- 30L
num_threads <- 120L
num_threads_impute <- 90L
# Name of clustering variable
cluster_name <- "X_06"
# how far from 0 and 1 should propensity score estimates at least be:
tolerence <- 0.01
# the horizon where excess risk is predicted:
horizons <- c(90, 120)
# how many imputations to run with MICE
n_imputations <- 30L
# how many iterations to run MICE
mice_maxit <- 30L

## path to save analysis objects .................................................................................. ----
training_path <- "PATH"
results_path <- "PATH"
tables_path <- "PATH"
figures_path <- "PATH"
## list of seeds used for reproducability ......................................................................... ----
seeds <- c(
  -735087720, -10699101, 663442021, 449770544, 617909963, -413161926, 976055038, -72918859, -485090593, -596908609, 
  570210709, 661278320, -411759050, 563817149, 480611678, -580082426, -858973366, -710513816, 577107071, -338327084, 
  978676888, -926563065, -734251329, -426467832, -20145358, -615079854, -889653913, -679644580, -775810847, 30566385, 
  813796033, 323467803, -195989677, 655610715, 870182586, -208929002, 907365007, -298395367, -334384224, -835655092, 
  15342222, -450686552, 457400328, 951463918, -876456446, 156816919, -699829902, -330489018, -365759238, 281839434, 
  -598677075, 229848288, 263689993, -355309037, 902892950, -453387783, 514213700, 438382152, 324725423, 209974196, 
  -972541913, 565574127, 747476937, -180374364, -458221123, 828721326, -772133844, -982699607, -282218739, 843858245, 
  715296968, -797070600, 842088228, 816781881, -7661028, 859660052, -931862659, -708166124, -50800671, -191032638, 
  -851990868, 212070311, -767906417, 815197927, -952137283, 764419573, -30791771, -509241984, -337012672, 991684812, 
  -301811156, -193930035, 662946266, -94388585, -223694519, -183020550, -680821091, 934320028, -951060409, 948854802
)

## load thiazide cohort ............................................................................................ ----
study_cohort <- readRDS("data/cohort_raw_data/study_cohort.rds")

# ===   2. Data wrangling   ======================================================================================= ====
## list of covariate names connected to a short-hand (X_xx) ....................................................... ----
names_index <- tibble(
  short = c(
    "Y_130", "D_130", "Y_125", "D_125", 
    "Y_130_90", "D_130_90", "Y_125_90", "D_125_90", 
    "Y_130_120", "D_130_120", "Y_125_120", "D_125_120", 
    "W_bvc", "W_hvr", "W", 
    paste0(
      "X_", 
      c(rep(0, 9), rep(1:7, each = 10), rep(8, 3)), 
      c(1:9, rep(0:9, times = 7), 0:2)
    )
  ),
  full = c( # event indicator: (0: censored, 1: hyponatraemia)
    "Follow-up for <130mM event", "hyponatraemia (<130mM) event", 
    "Follow-up for <125mM event", "hyponatraemia (<125mM) event", 
    "hyponatraemia (<130 mM) event inside 90 days", "hyponatraemia (<130 mM) event or >90 days",
    "hyponatraemia (<125 mM) event inside 90 days", "hyponatraemia (<125 mM) event or >90 days",
    "hyponatraemia (<130 mM) event inside 120 days", "hyponatraemia (<130 mM) event or >120 days",
    "hyponatraemia (<125 mM) event inside 120 days", "hyponatraemia (<125 mM) event or >120 days",
    "BFZ initiation", "HCTZ initiation", "thiazide drug initiation", 
    "Age", "Sex", "Region of residence", "House hold income", "Years of education", "Calendar year of inclusion", 
    "Heart failure", "Essential hypertension", "Secondary hypertension", "Ischemic heart disease", 
    "Cerebrovascular disease", "Other CNS disorders", "Arrhythmias", "Any malignancy", 
    "Malignancy associated with hyponatraemia",  "Other malignancy", "Renal disorders", "Liver disease and peritonitis", 
    "Pancreatitis", "Chronic obstructive pulmonary disease", "Diabetes", "Dehydration", "Frail general health", 
    "Physical impairment", "Mental impairment", "Rehabilitation contacts", "Podiatric contacts", "Alcohol abuse", 
    "drug abuse", "HIV", "Anorexia and primary polydipsia", "History of hyponatraemia", 
    "Drugs for peptic ulcer and gastroesophageal reflux", "Drugs used for obstipation or diarrhea", 
    "Antidiabetic (not insulin)", "Insulin", "Oral anticoagulants", "Aspirin", "ADPi", "Nitrates", "Beta-blockers", 
    "Lipid lowering drugs", "Desmopressin", "Eltroxin", "NSAIDs", "Antiepileptics", "Opioids", "Antidepressives", 
    "Antipsychotics", "Agents for obstructive pulmonary disease", "No. of different prescription drugs used",
    "Days of hospitalization in the past year prior to index data", 
    "No. of outpatient contacts in the past year prior to index data", 
    "No. of primary care contacts in the last 2 years", "cSodium", "ceGRF", "cPotassium", "cHemoglobin", "cALAT", 
    "cAlkaline phosphatase", "cAlbumin", "cLactate dehydrogenase", "cCarbamide", "cThrombocytes", "cLeukocytes", 
    "cC-reactive protein", "cCholesterol", "cLDL-C", "Sodium", "eGFR", "Potassium", "Hemoglobin", "ALAT", 
    "Alkaline phosphatase", "Albumin", "Lactate dehydrogenase", "Carbamide", "Thrombocytes", "Leukocytes", 
    "C-reactive protein", "Cholesterol", "LDL-C"
  ),
  levels = c(
    rep(list(character(0)), 17), 
    list(c("Hovedstaden", "Nordjylland", "Sjælland", "Syddanmark")),
    list(c("Q1", "Q2", "Q3", "Q4", "Q5")),
    list(c("<10", "10-12", "13-15", ">15")),
    list(c("2014", "2015", "2016", "2017", "2018")),
    rep(list(character(0)), 25),
    list(c("never", "<4 months", "\u22654 months")),
    rep(list(character(0)), 18),
    list(c("0", "1-3", "4-6", "7-9", paste0("\U2265", "10"))),
    list(c("0", "1-7", "8-14", paste0("\U2265", "15"))),
    list(c("0", "1-3", "4-6", paste0("\U2265", "7"))),
    list(c("0", "1-3", "4-6", paste0("\U2265", "7"))),
    rep(list(c(paste0("Q", 1:4), "NA")), 14),
    rep(list(character(0)), 14)
  )
)
names_index_comb <- names_index |>
  unnest(levels, keep_empty = TRUE) |>
  mutate(
    short = ifelse(is.na(levels), short, paste0(short, "_", levels)),
    full = ifelse(is.na(levels), full, paste(full, "-", levels)),
    .keep = "none"
  )

## Create cohort dataset .......................................................................................... ----
# Remove individuals without an indication of hypertension.
# cohort with all covariates included. Selection of different subsets of the
# covariates happens later.
tih_cohort_all <- study_cohort |>
  filter(grepl("5(4|7)$", indo)) |>
  MakeTIHCohort(h = horizons)

## Split into training and test data .............................................................................. ----
# split temporally. test data 2019-2020
tih_cohort_test <- tih_cohort_all |>
  filter(!(X_06 %in% c("2014", "2015", "2016", "2017", "2018")))
tih_cohort_test$X_06 <- droplevels(tih_cohort_test$X_06)
attr(tih_cohort_test$X_06, "label") <- "Calendar year of inclusion"

# filter on complete cases and four least missing blood test results
tih_cohort_test_complete_case <- 
  tih_cohort_test |>
  select(
    all_of(
      c("Y_130", "D_130", "Y_125", "D_125", 
        "Y_130_90", "D_130_90", "Y_125_90", "D_125_90", 
        "Y_130_120", "D_130_120", "Y_125_120", "D_125_120", 
        "W_bvc", "W_hvr", "W", 
        paste0("X_0", c(1:9)),
        paste0("X_", c(10:57, 59, 69, 70, 71, 73)),
        paste0("X_", c(69, 70, 71, 73), "_tbi"))
    )
  ) |>
  filter(!is.na(rowSums(select(tih_cohort_test, "X_69", "X_70", "X_71", "X_73")))) |>
  replace_na(
    list(
      X_04 = "Q4",
      X_05 = "10-12"
    )
  )

# systematic temporal shift from thiazide use to non-thiazide use
# any thiazide
tih_cohort_all |>
  count(W, X_06) |>
  mutate(percent = 100 * n / sum(n), .by = X_06) |>
  arrange(X_06, W)

# bfz
tih_cohort_all |>
  filter(!is.na(W_bvc)) |>
  count(W_bvc, X_06) |>
  mutate(percent = 100 * n / sum(n), .by = X_06) |>
  arrange(X_06, W_bvc)

# hctz
tih_cohort_all |>
  filter(!is.na(W_hvr)) |>
  count(W_hvr, X_06) |>
  mutate(percent = 100 * n / sum(n), .by = X_06) |>
  arrange(X_06, W_hvr)

## Construct a table 1 from development and validation cohorts .................................................... ----
table_one_vars <- tih_cohort_test |>
  select(starts_with("X")) |>
  select(!ends_with("tbi")) |>
  select(matches("[0-4]\\d|5[0-4]|69|7\\d|8[0-2]")) |>
  map_chr(\(x) attr(x, "label"))
table_one_vars_cc <- tih_cohort_test_complete_case |>
  select(starts_with("X")) |>
  select(!ends_with("tbi")) |>
  select(matches("[0-4]\\d|5[0-4]|69|7\\d|8[0-2]")) |>
  map_chr(\(x) attr(x, "label"))

table_one_val <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort_test |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_val_cc <- CreateTableOne(
  vars = as.character(table_one_vars_cc), 
  strata = "W", 
  data = tih_cohort_test_complete_case |>
    rename(
      !!!setNames(as.list(names(table_one_vars_cc)), table_one_vars_cc)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_val_bfz <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort_test |>
    filter(!is.na(W_bvc)) |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_val_hctz <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort_test |>
    filter(!is.na(W_hvr)) |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)

table_one_val_exp <- print(
  table_one_val,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_val_cc_exp <- print(
  table_one_val_cc,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_val_bfz_exp <- print(
  table_one_val_bfz,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_val_hctz_exp <- print(
  table_one_val_hctz,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)

for (col in seq_len(ncol(table_one_val_exp))) {
  for (row in seq_len(nrow(table_one_val_exp))) {
    cell <- table_one_val_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_val_exp[row, col] <- ifelse(
        as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, 
        "    <5 ( 0.0) ", 
        cell
      )
    }
  }
}
for (col in seq_len(ncol(table_one_val_cc_exp))) {
  for (row in seq_len(nrow(table_one_val_cc_exp))) {
    cell <- table_one_val_cc_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_val_cc_exp[row, col] <- ifelse(
        as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, 
        "    <5 ( 0.0) ", 
        cell
      )
    }
  }
}
for (col in seq_len(ncol(table_one_val_bfz_exp))) {
  for (row in seq_len(nrow(table_one_val_bfz_exp))) {
    cell <- table_one_val_bfz_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_val_bfz_exp[row, col] <- ifelse(
        as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, 
        "    <5 ( 0.0) ", 
        cell
      )
    }
  }
}
for (col in seq_len(ncol(table_one_val_hctz_exp))) {
  for (row in seq_len(nrow(table_one_val_hctz_exp))) {
    cell <- table_one_val_hctz_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_val_hctz_exp[row, col] <- ifelse(
        as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, 
        "    <5 ( 0.0) ", 
        cell
      )
    }
  }
}

write.csv2(table_one_val_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_val.csv"))
write.csv2(table_one_val_cc_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_val_cc.csv"))
write.csv2(table_one_val_bfz_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_val_bfz.csv"))
write.csv2(table_one_val_hctz_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_val_hctz.csv"))

## count of censoring by thiazide group ----
censtab <- left_join(
  tih_cohort_test,
  study_cohort |> 
    select("pnr", "dead_exit", "emigrated_exit", "hyp_stop_add_on_exit", "reg_midt_censor_date_exit"), 
  by = "pnr"
)
censtab |> filter(Y < 120) |> count(W, Y_h)
censtab |> filter(is.na(Y_h) & Y < 120) |> count(W, dead_exit)
censtab |> filter(is.na(Y_h) & Y < 120) |> count(W, emigrated_exit)
censtab |> filter(is.na(Y_h) & Y < 120) |> count(W, hyp_stop_add_on_exit)
censtab |> filter(is.na(Y_h) & Y < 120) |> count(W, reg_midt_censor_date_exit)

# Imputation of missing blood test measurements ------------------------------------------------------------------- ----
## Imputation using multiple imputation ........................................................................... ----
# use multiple imputation to impute m complete data matrices on validation set 
if (impute_data) {
  plan(multisession, workers = 3)
  future_pmap(
    .l = list( 
      cohort =
        list(
          tih_cohort_test, 
          tih_cohort_test |> filter(!is.na(W_bvc)),
          tih_cohort_test |> filter(!is.na(W_hvr))
        ), 
      exposure = 
        c("thiazide", "bfz", "hctz"),
      tbi = rep(c(years(1)), each = 3),
      parallelseed = seeds[1:3]
    ),
    .options = furrr_options(
      packages = c("future", "mice", "dplyr", "glue"),
      globals = c("num_threads_impute", "n_imputations", "mice_maxit", "futuremice"),
      seed = TRUE
    ),
    .f = \(cohort, exposure, tbi, parallelseed) {
      # set option used to retrive available cores for futuremice
      options("parallelly.availableCores.methods" = "system")
      
      # replace any measurement not within the lookback of 1 year with NA
      for (i in 69:82) {
        cohort <- dplyr::mutate(
          cohort,
          "X_{i}" := ifelse(
            cohort[[glue("X_{i}_tbi")]] < tbi,
            cohort[[glue("X_{i}")]],
            NA_real_
          )
        )
      }
      
      # define matrix used in mice to determine variables used to impute
      pred_matrix <- mice(
        cohort, 
        # impute income quintiles, education level, and lab test results
        # Income quintiles and education levels are ordered multilevel factors, 
        # using a proportional odds logistic regression model.
        # lab test results are numeric variables, using predictive mean matching
        method = c(rep("", 19), "polr", "polr", rep("", 49), rep(c("pmm", "", ""), times = 14)), 
        m = 1L, 
        maxit = 0L
      )$predictorMatrix
      # Remove outcome, exposure and unrelated variables from being used for imputation
      pred_matrix[, c("Y_130", "D_130", "Y_125", "D_125")] <- 0
      pred_matrix[, c("Y_130_90", "D_130_90", "Y_125_90", "D_125_90")] <- 0
      pred_matrix[, c("Y_130_120", "D_130_120", "Y_125_120", "D_125_120")] <- 0
      pred_matrix[, c("W_bvc", "W_hvr", "W")] <- 0
      pred_matrix[, which(colnames(pred_matrix) %in% paste0("X_", 55:68))] <- 0
      pred_matrix[, which(colnames(pred_matrix) %in% paste0("X_", 69:82, "_tbi"))] <- 0
      # Only impute on relevant variables
      pred_matrix[-(which(rownames(pred_matrix) %in% c("X_04", "X_05", paste0("X_", 69:82)))), ] <- 0
      mice_imp <- futuremice(
        cohort,
        method = c(rep("", 19), "polr", "polr", rep("", 49), rep(c("pmm", "", ""), times = 14)),
        predictorMatrix = pred_matrix,
        m = n_imputations, 
        maxit = mice_maxit,
        parallelseed = parallelseed,
        n.core = min(n_imputations, num_threads_impute %/% 3L),
        available = 80L,
        # arguments passed down to pmm():
        blots = list(
          X_04 = list(
            nnet.maxit = 100, # max number of iterations
            nnet.trace = FALSE, # switch for tracing optimization
            nnet.MaxNWts = 1500 # maximum allowed number of weights.
          ),
          X_05 = list(
            nnet.maxit = 100,
            nnet.trace = FALSE,
            nnet.MaxNWts = 1500
          ),
          X_69 = list(
            donors = 5L, # size of donor pool. consider 5L (5L default)
            exclude = -99999999 # values to exclude from donor pool
          ),
          X_70 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_71 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_72 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_73 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_74 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_75 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_76 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_77 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_78 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_79 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_80 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_81 = list(
            donors = 5L,
            exclude = -99999999
          ),
          X_82 = list(
            donors = 5L,
            exclude = -99999999
          )
        )
      )
      saveRDS(
        mice_imp, 
        paste0("data/cohorts/", exposure, "_tbi_", tbi@year, "_test_cohort_imputed_mice.rds")
      )
    } 
  )
  plan(sequential)
}

## read in imputed data from multiple imputation .................................................................. ----
tih_mice_imp_test <- readRDS("data/cohorts/thiazide_tbi_1_test_cohort_imputed_mice.rds")
tih_tbi_1_test_cohort_imputed_mice <- list()
for (i in seq_len(tih_mice_imp_test$m)) {
  tih_tbi_1_test_cohort_imputed_mice[[i]] <- complete(tih_mice_imp_test, i)
}
bfz_mice_imp_test <- readRDS("data/cohorts/bfz_tbi_1_test_cohort_imputed_mice.rds")
bfz_tbi_1_test_cohort_imputed_mice <- list()
for (i in seq_len(bfz_mice_imp_test$m)) {
  bfz_tbi_1_test_cohort_imputed_mice[[i]] <- complete(bfz_mice_imp_test, i)
}
hctz_mice_imp_test <- readRDS("data/cohorts/hctz_tbi_1_test_cohort_imputed_mice.rds")
hctz_tbi_1_test_cohort_imputed_mice <- list()
for (i in seq_len(hctz_mice_imp_test$m)) {
  hctz_tbi_1_test_cohort_imputed_mice[[i]] <- complete(hctz_mice_imp_test, i)
}

tih_tbi_1_test_cohort_imputed <- tih_tbi_1_test_cohort_imputed_mice
bfz_tbi_1_test_cohort_imputed <- bfz_tbi_1_test_cohort_imputed_mice
hctz_tbi_1_test_cohort_imputed <- hctz_tbi_1_test_cohort_imputed_mice

# filter on complete cases and four least missing blood test results
tih_cohort_test_cc_comp <-  map(
  tih_tbi_1_test_cohort_imputed,
  \(cohort) {
    cohort |>
      select(
        all_of(
          c("Y_130", "D_130", "Y_125", "D_125", 
            "Y_130_90", "D_130_90", "Y_125_90", "D_125_90", 
            "Y_130_120", "D_130_120", "Y_125_120", "D_125_120", 
            "W_bvc", "W_hvr", "W",
            paste0("X_0", c(1:9)),
            paste0("X_", c(10:57, 59, 69, 70, 71, 73)),
            paste0("X_", c(69, 70, 71, 73), "_tbi"))
        )
      )
  }
)

# ===   3. Analysis   ============================================================================================= ====
# The test dataset is multiply imputed separately from the training dataset.

### Train survival forest for sample weights, regression forests for expected outcome and propensity to treat ..... ---- 
if(train_cf) {
  # set seed to insure consistent indices for dynamic subgroups
  set.seed(seeds[4])
  plan(sequential)
  # map over different analyses (e.g. thiazide combined, BFZ, and HCTZ)
  future_pmap(
    .l = list(
      name = list(
        # tih = all thiazide drugs, 4mo = events within 4 months of index date, abt(n) = all blood tests n years back
        # nbt = no blood tests, kbt(n) kidney related blood tests n years back, 
        # rds = removed Denmark specific covariates
        # bfz = only bfz's, hctz = only hctz's
        "tih_4mo_abt1",
        "tih_3mo_abt1_severe",
        "tih_4mo_nbt",
        "tih_4mo_kbt1",
        "tih_4mo_kbt1_rds",
        "bfz_4mo_abt1",
        "bfz_4mo_kbt1_rds",
        "hctz_4mo_abt1",
        "hctz_4mo_kbt1_rds",
        "tih_4mo_cc",
        "tih_4mo_cc_comp"
      ),
      horizon = list(
        120,
        90,
        120,
        120,
        120,
        120,
        120,
        120,
        120,
        120,
        120
      ),
      thiazide_limit = list(
        130,
        125,
        130,
        130,
        130,
        130,
        130,
        130,
        130,
        130,
        130
      ),
      cohort = list(
        tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
        list(tih_cohort_test), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
        tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
        tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        list(tih_cohort_test_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
        tih_cohort_test_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
      ),
      seed = as.list(seeds[5:15]),
      path = list(
        paste0(training_path, "4mo_thiazide_all/1year/"),
        paste0(training_path, "3mo_tiazide_all_severe/1year/"),
        paste0(training_path, "4mo_thiazide_no_lab/"),
        paste0(training_path, "4mo_thiazide_lab_fluid/1year/"),
        paste0(training_path, "4mo_thiazide_international/1year/"),
        paste0(training_path, "4mo_bfz_all/1year/"),
        paste0(training_path, "4mo_bfz_international/1year/"),
        paste0(training_path, "4mo_hctz_all/1year/"),
        paste0(training_path, "4mo_hctz_international/1year/"),
        paste0(training_path, "4mo_thiazide_complete_case/"),
        paste0(training_path, "4mo_thiazide_cc_comp/1year/")
      ),
      continuous = list(
        list(# thiazide all
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# thiazide all severe
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# thiazide no lab
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89])$|"
        ),
        list(# thiazide only fluid lab results
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# thiazide only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# bfz all
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# bfz only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# hctz all
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# hctz only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# thiazide complete case
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        ),
        list(# thiazide complete case comparison
          sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
          apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
        )
      ),
      discrete = list(
        list(# thiazide all
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# thiazide all severe
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# thiazide no lab
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# thiazide only fluid lab results
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# thiazide only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(32|5[1-4])$")),
          apriori = "(^X_32$|"
        ),
        list(# bfz all
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# bfz only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(32|5[1-4])$")),
          apriori = "(^X_32$|"
        ),
        list(# hctz all
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# hctz only fluid lab results and no Denmark specific covariates
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(32|5[1-4])$")),
          apriori = "(^X_32$|"
        ),
        list(# thiazide complete case
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        ),
        list(# thiazide complete case comparison
          sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
          exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
          cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
          apriori = "(^X_(0[45]|32)$|"
        )
      ),
      tune_parameters = list(
        list(# thiazide all
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide all severe
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide no lab
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide only fluid lab results
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide only fluid lab results and no Denmark specific covariates
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# bfz all
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# bfz only fluid lab results and no Denmark specific covariates
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# hctz all
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# hctz only fluid lab results and no Denmark specific covariates
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide complete case
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        ),
        list(# thiazide complete case comparison
          exp_out = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves"),
          cf = c("sample.fraction", "mtry", "min.node.size", "honesty.fraction", "honesty.prune.leaves")
        )
      )
    ) |> map(\(x) x[setup_index]),
    .options = furrr_options(
      packages = c("grf", "rlang", "dplyr", "stringr", "tidyr", "purrr", "glue", "stats", "writexl", "future", 
                   "furrr", "future.callr"),
      globals = c("num_threads", "num_clusters", "names_index", "names_index_comb", "horizon", "GRFAnalysisWrapper", 
                  "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel", "VariableImportanceWrapper", "tolerence",
                  "save_sampleweight_forest", "save_exposure_forest", "save_outcome_forest", "n_folds", "n_rankings",
                  "dynamic_subgroups", "vcovHC", "vi_models", "cluster_name"),
      seed = TRUE
    ),
    .f = \(name, horizon, thiazide_limit, cohort, seed, path, continuous, discrete, tune_parameters) {
      # limit number of clusters to number of imputations in cohort
      if (length(cohort) < num_clusters) num_clusters <- length(cohort)
      # set up future plan for analysis of each imputed dataset
      if (num_clusters == 1) {
        plan(sequential)
      } else {
        plan(multisession, workers = num_clusters)
      }
      # map over each imputation
      future_pmap(
        .l = list(
          cohort = cohort,
          data_index = seq_along(cohort)
        ),
        .options = furrr_options(
          packages = c("grf", "rlang", "dplyr", "stringr", "tidyr", "purrr", "glue", "stats", "writexl", "future", 
                       "furrr", "future.callr", "DescTools"),
          globals = c("path", "num_threads", "num_clusters", "names_index", "names_index_comb", "horizon", 
                      "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel", "tolerence",
                      "VariableImportanceWrapper", "save_sampleweight_forest", "save_exposure_forest", 
                      "save_outcome_forest", "continuous", "discrete", "tune_parameters", "indices_dyn", 
                      "n_folds", "n_rankings", "name", "vcovHC", "vi_models", "dynamic_subgroups",
                      "thiazide_limit", "cluster_name"),
          seed = seed
        ),
        .f = \(cohort, data_index) {
          # create seeds for grf forest generation
          rng_seed <- sample(-999999999:999999999, 4)
          
          if (file.exists(paste0(path, "cfw_full_test_mod_", data_index, ".rds"))) {
            # loading existing trained test model
            cfw_full_test <- readRDS(paste0(path, "cfw_full_test_mod_", data_index, ".rds"))
          } else {
            ### Survival forest model to predict sample weights ..................................................... ----
            # train survival forest to obtain sample weights
            # Note: automatic tuning not implemented for survival forests
            cfw_sample_weight_test <-
              (function(data) {
                # set up variables based on the covariates used
                continuous <- continuous$sample_weight
                discrete <- discrete$sample_weight
                # calculate sample weights based on censoring process
                # The censoring probabilities are P(C > T|X, W), therefore the 
                # variables used in the survival forest include both X and W.
                forest_C <- GRFAnalysisWrapper(
                  data |> mutate(D_C = 1 - !!sym(paste0("D_", thiazide_limit))),
                  type = "survival",
                  continuous = c(!!continuous, "W"),
                  discrete = !!discrete,
                  Y = paste0("Y_", thiazide_limit), D = "D_C",
                  clusters = data[[cluster_name]], # cluster on calendar year of inclusion
                  num.trees = 2000, alpha = 0.0002,
                  seed = rng_seed[1],
                  num.threads = floor(num_threads / num_clusters)
                )
                # we use the survival probability of the censoring process at the "outcome time",
                # which is the observed time when considering a hyponatraemia event and at the 
                # horizon when observing no hyponatraemia up to the horizon.
                # If one makes it past the horizon, we know they survived up to this time.
                predict_C <- predict(
                  forest_C, 
                  failure.times = pmin(data[[paste0("Y_", thiazide_limit)]], horizon), 
                  prediction.times = "time"
                )
                censoring_prob <- predict_C$predictions
                observed_events <- !is.na(data[[paste0("Y_", thiazide_limit, "_", horizon)]])
                sample_weights <- 1 / censoring_prob[observed_events]
                return(
                  list(
                    survival_forest = forest_C,
                    sample_weights = sample_weights
                  )
                )
              })(data = cohort)
            
            if (save_sampleweight_forest) {
              saveRDS(cfw_sample_weight_test, paste0(path, "cfw_sample_weight_test_", data_index, ".rds"))
            }
            
            sample_weights_test <- cfw_sample_weight_test$sample_weights
            
            # remove survival forest from memory
            rm("cfw_sample_weight_test")
            gc()
            
            # create data without censored observations (complete-case)
            data_obs_test <- filter(cohort, !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
            
            # specify tuning parameters
            tune.num.trees <- if (thiazide_limit == 125) 400 else 200
            tune.num.reps <- 100
            tune.num.draws <- if (thiazide_limit == 125) 2000 else 1000
            
            ### Regression forest models to predict expected outcome and propensity ............................... ----
            cfw_exp_out_test <- (function(data, sample_weights, tolerence) {
                # set up variables based on the covariates used
                continuous <- continuous$exp_out
                discrete <- discrete$exp_out
                tune_parameters <- tune_parameters$exp_out
                
                # calculate expected outcome and treatment propensity
                forest_Y <- GRFAnalysisWrapper(
                  data,
                  type = "regression",
                  continuous = !!continuous,
                  discrete = !!discrete,
                  Y = paste0("Y_", thiazide_limit, "_", horizon),
                  sample.weights = sample_weights,
                  clusters = data[[cluster_name]], # cluster on year of inclusion
                  num.trees = 2000, alpha = 0.0002,
                  tune.parameters = tune_parameters,
                  tune.num.trees = tune.num.trees, tune.num.reps = tune.num.reps, tune.num.draws = tune.num.draws,
                  seed = rng_seed[2],
                  num.threads = floor(num_threads / num_clusters)
                )
                Sys.sleep(1)
                forest_W <- GRFAnalysisWrapper(
                  data,
                  type = "regression",
                  continuous = !!continuous,
                  discrete = !!discrete,
                  Y = "W",
                  sample.weights = sample_weights,
                  clusters = data[[cluster_name]], # cluster on year of inclusion
                  num.trees = 2000, alpha = 0.0002,
                  tune.parameters = tune_parameters,
                  tune.num.trees = tune.num.trees, tune.num.reps = tune.num.reps, tune.num.draws = tune.num.draws,
                  seed = rng_seed[3],
                  num.threads = floor(num_threads / num_clusters)
                )
                Y_hat <- predict(forest_Y, num.threads = floor(num_threads / num_clusters))$predictions
                W_hat <- predict(forest_W, num.threads = floor(num_threads / num_clusters))$predictions
                # expected outcome and treatment propensity in censored
                Y_hat_cens <- predict(
                  forest_Y,
                  newdata = cohort |>
                    filter(is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
                    select(!!continuous) |>
                    bind_cols(
                      cohort |>
                        filter(is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
                        select(!!discrete) |>
                        DiscreteCovariatesToOneHot()
                    ) |>
                    select(all_of(names(forest_Y$X.orig))),
                  num.threads = floor(num_threads / num_clusters)
                )$predictions
                W_hat_cens <- predict(
                  forest_W,
                  newdata = cohort |>
                    filter(is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
                    select(!!continuous) |>
                    bind_cols(
                      cohort |>
                        filter(is.na(.data[[paste0("Y_", thiazide_limit, "_", horizon)]])) |>
                        select(!!discrete) |>
                        DiscreteCovariatesToOneHot()
                    ) |>
                    select(all_of(names(forest_W$X.orig))),
                  num.threads = floor(num_threads / num_clusters)
                )$predictions
                # return list with forest objects as well as expected outcome and treatment propensity in full cohort
                out <- list(
                  outcome_forest = forest_Y,
                  exposure_forest = forest_W,
                  exp_out = tibble(
                    id = c(
                      which(!is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]])),
                      which(is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]]))
                    ),
                    observed = id %in% which(!is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]])),
                    Y_hat = c(Y_hat, Y_hat_cens),
                    W_hat = c(W_hat, W_hat_cens)
                  ) |>
                    arrange(id)
                )
                if (any(out$exp_out$W_hat < tolerence | 1 - tolerence < out$exp_out$W_hat)) {
                  out$exp_out$W_hat[which(out$exp_out$W_hat < tolerence)] <- tolerence
                  out$exp_out$W_hat[which(out$exp_out$W_hat > 1 - tolerence)] <- 1 - tolerence
                  warning("Some units have propensity scores close to 0 or 1. ",
                          "Scores are set to be at least ", tolerence, " away from 0 and 1.")
                }
                return(out)
              })(data = data_obs_test, sample_weights = sample_weights_test, tolerence = tolerence)
            
            if (save_exposure_forest) {
              saveRDS(cfw_exp_out_test$exposure_forest, paste0(path, "cfw_exposure_test_", data_index, ".rds"))
            }
            if (save_outcome_forest) {
              saveRDS(cfw_exp_out_test$outcome_forest, paste0(path, "cfw_outcome_test_", data_index, ".rds"))
            }
            
            # save table with expected outcome and propensity scores in the full cohort
            saveRDS(cfw_exp_out_test$exp_out, paste0(path, "cfw_exp_out_test_", data_index, ".rds"))
            
            # extract expected outcome and propensity
            Y_hat_test <- filter(cfw_exp_out_test$exp_out, observed)$Y_hat
            W_hat_test <- filter(cfw_exp_out_test$exp_out, observed)$W_hat
            
            # remove outcome and exposure forests from memory
            rm("cfw_exp_out_test")
            gc()
            
            ### Causal forest models with all covariates ............................................................ ----
            cfw_full_test <-
              (function(data, sample_weights, Y_hat, W_hat) {
                # set up variables
                continuous <- continuous$cf
                discrete <- discrete$cf
                tune_parameters <- tune_parameters$cf
                
                out <- GRFAnalysisWrapper(
                  data,
                  type = "causal",
                  continuous = !!continuous,
                  discrete = !!discrete,
                  Y = paste0("Y_", thiazide_limit, "_", horizon), W = "W",
                  Y.hat = Y_hat, W.hat = W_hat,
                  sample.weights = sample_weights,
                  clusters = data[[cluster_name]], # cluster on calendar year of inclusion
                  num.trees = 2000, alpha = 0.0002,
                  tune.parameters = tune_parameters,
                  tune.num.trees = tune.num.trees, tune.num.reps = tune.num.reps, tune.num.draws = tune.num.draws,
                  seed = rng_seed[4],
                  num.threads = floor(num_threads / num_clusters)
                )
                if("X_03_SjÃ¦lland" %in% names(out$X.orig)) {
                  names(out$X.orig)[
                    which(names(out$X.orig) == "X_03_SjÃ¦lland")
                  ] <- "X_03_Sjælland"
                }
                return(out)
              })(
                data = data_obs_test,
                sample_weights = sample_weights_test,
                Y_hat = Y_hat_test,
                W_hat = W_hat_test
              )
            
            saveRDS(cfw_full_test, paste0(path, "cfw_full_test_mod_", data_index, ".rds"))
          }
          
          # calculate variable importance from model with all covariates
          cfw_full_test_vi <- VariableImportanceWrapper(cfw_full_test, names_index)
          saveRDS(cfw_full_test_vi, paste0(path, "cfw_full_test_vi_", data_index, ".rds"))
          writexl::write_xlsx(cfw_full_test_vi, paste0(path, "cfw_full_test_vi_", data_index, ".xlsx"))
          
          # remove causal forest from memory
          rm("cfw_full_test")
          
          # combine variable importance from different levels of categorical variables
          cfw_full_test_vi_comb <- (\(vi, nm_index) {
            for (nm in names_index$full) {
              if (
                nrow(
                  filter(
                    cfw_full_test_vi,
                    str_detect(cfw_full_test_vi$variable_name, str_c("^", nm, " "))
                  )
                ) > 1
              ) {
                vi <- bind_rows(
                  vi,
                  cfw_full_test_vi |>
                    filter(str_detect(variable_name, str_c("^", nm, " "))) |>
                    summarise(
                      variable_name = nm,
                      variable_importance_num = sum(variable_importance_num)
                    )
                ) |>
                  filter(str_detect(variable_name, str_c("^", nm, " "), negate = TRUE))
              }
            } 
            
            return(
              vi |>
                mutate(
                  variable_name = trimws(variable_name),
                  variable_importance = sprintf("%.3f", variable_importance_num)
                ) |>
                arrange(desc(variable_importance_num)) |>
                left_join(
                  select(nm_index, "short", "full"),
                  by = c("variable_name" = "full")
                )
            )
          })(
            vi = cfw_full_test_vi[, 1:2],
            nm_index = filter(names_index, !grepl("^(X_(5[5-9]|6[0-8])$|[^X])", short))
          )
          writexl::write_xlsx(cfw_full_test_vi_comb, paste0(path, "cfw_full_test_vi_comb_", data_index, ".xlsx"))
          
          return(NULL)
        }
      )
      
      # aggregate variable importance from each imputed full model:
      if (file.exists(paste0(path, "cfw_full_test_vi_agg.rds"))) {
        # load existing aggregated variable importance values
        cfw_full_test_vi_agg <- readRDS(paste0(path, "cfw_full_test_vi_agg.rds"))
        cfw_full_test_vi_comb <- readRDS(paste0(path, "cfw_full_test_vi_comb_agg.rds"))
      } else {
        # load variable importance values from each imputed model
        cfw_full_test_vi <- list()
        for (i in seq_along(list.files(path, pattern = "^cfw_full_test_vi_\\d{1,}.rds"))) {
          cfw_full_test_vi[[i]] <- readRDS(paste0(path, "cfw_full_test_vi_", i, ".rds"))
        }
        # pool variable importance estimates
        cfw_full_test_vi_agg <- cfw_full_test_vi |>
          map(\(df) df |> dplyr::arrange(variable_name, .locale = "en")) |>
          purrr::list_transpose(simplify = FALSE) |>
          (\(list) {
            list <- list[c("variable_name", "variable_importance_num")]
            i <- 1
            repeat {
              var_name <- c()
              for (j in seq_along(list$variable_name)) {
                var_name[j] <- list$variable_name[[j]][i]
              }
              if (length(unique(var_name)) > 1) {
                name_to_append <- sort(var_name)[1]
              } else if (!is.na(unique(var_name))){
                i <- i + 1
                next
              } else {
                break
              }
              for (j in seq_along(list$variable_name)) {
                if (is.na(list$variable_name[[j]][i]) || list$variable_name[[j]][i] != name_to_append) {
                  list$variable_name[[j]] <- append(list$variable_name[[j]], name_to_append, i - 1)
                  list$variable_importance_num[[j]] <- append(list$variable_importance_num[[j]], 0, i - 1)
                }
              }
              i <- i + 1
            }
            return(list)
          })() |>
          map(\(x) list_transpose(x, simplify = TRUE)) |>
          list_transpose(simplify = FALSE) |>
          map(
            \(x) {
              if (length(x$variable_importance_num) > 1) {
                pool <- pool.scalar(Q = x$variable_importance_num, U = rep(0, 5), n = 1, k = 1, rule = "rubin1987")
                out <- tibble(
                  variable_name = x$variable_name[1],
                  variable_importance_num = pool$qbar,
                  variable_importance = sprintf("%.5f", variable_importance_num)
                )
              } else {
                pool <- NULL
                out <- x |>
                  as_tibble() |>
                  mutate(
                    variable_importance = sprintf("%.5f", variable_importance_num)
                  )
              }
              return(list(out, pool))
            }
          ) |> 
          (\(model) {
            tbl <- map(model,\(x) x[[1]]) |> 
              list_rbind() |>
              arrange(desc(variable_importance_num))
            pool <- map(model, \(x) x[[2]])
            names(pool) <- map_chr(model, \(x) x[[1]]$variable_name)
            return(list(table = tbl, pool_results = pool))
          })()
        saveRDS(cfw_full_test_vi_agg, paste0(path, "cfw_full_test_vi_agg.rds"))
        # combine variable importance from different levels of categorical variables
        cfw_full_test_vi_comb <- (\(vi, nm_index) {
          for (nm in names_index$full) {
            if (
              nrow(
                filter(
                  cfw_full_test_vi_agg$table,
                  str_detect(cfw_full_test_vi_agg$table$variable_name, str_c("^", nm, " "))
                )
              ) > 1
            ) {
              vi <- bind_rows(
                vi,
                cfw_full_test_vi_agg$table |>
                  filter(str_detect(variable_name, str_c("^", nm, " "))) |>
                  summarise(
                    variable_name = nm,
                    variable_importance_num = sum(variable_importance_num)
                  )
              ) |>
                filter(str_detect(variable_name, str_c("^", nm, " "), negate = TRUE))
            }
          } 
          return(
            vi |>
              mutate(
                variable_name = trimws(variable_name),
                variable_importance = sprintf("%.3f", variable_importance_num)
              ) |>
              arrange(desc(variable_importance_num)) |>
              left_join(
                select(nm_index, "short", "full"),
                by = c("variable_name" = "full")
              )
          )
        })(
          vi = cfw_full_test_vi_agg[["table"]][, 1:2],
          nm_index = filter(names_index, !grepl("^(X_(5[5-9]|6[0-8])$|[^X])", short))
        )
        saveRDS(cfw_full_test_vi_comb, paste0(path, "cfw_full_test_vi_comb_agg.rds"))
        writexl::write_xlsx(cfw_full_test_vi_comb, paste0(path, "cfw_full_test_vi_comb_agg.xlsx"))
      }
      cfw_full_vi_comb <- readRDS(paste0(path, "cfw_full_vi_comb_agg.rds"))
      
      future_pmap(
        .l = list(
          cohort = cohort,
          data_index = seq_along(cohort)
        ),
        .options = furrr_options(
          packages = c("grf", "rlang", "dplyr", "stringr", "tidyr", "purrr", "glue", "stats", "writexl", "future", 
                       "furrr", "future.callr", "DescTools"),
          globals = c("path", "num_threads", "num_clusters", "names_index", "names_index_comb", "horizon", 
                      "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel",
                      "VariableImportanceWrapper", "save_sampleweight_forest", "save_exposure_forest", 
                      "save_outcome_forest", "continuous", "discrete", "tune_parameters", "indices_dyn", 
                      "n_folds", "n_rankings", "name", "vcovHC", "vi_models", "thiazide_limit", "cluster_name",
                      "dynamic_subgroups", "cfw_full_test_vi_comb", "cfw_full_vi_comb"),
          seed = seed
        ),
        .f = \(cohort, data_index) {
          ### read in data from model with all covariates ......................................................... ----
          # loading existing trained model
          cfw_full_test <- readRDS(paste0(path, "cfw_full_test_mod_", data_index, ".rds"))
          sample_weights_test <- cfw_full_test$sample.weights
          W_hat_test <- cfw_full_test$W.hat
          Y_hat_test <- cfw_full_test$Y.hat
          # create data without censored observations (complete-case)
          data_obs_test <- filter(cohort, !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
          # set flag for further tuning if successful for all covariates
          tune <- cfw_full_test$tuning.output$status != "default"
          tune_param <- if (tune)  tune_parameters$cf else "none"
          # remove causal forest
          rm("cfw_full_test")
          gc()
          
          # specify tuning parameters
          tune.num.trees <- if (thiazide_limit == 125) 400 else 200
          tune.num.reps <- 100
          tune.num.draws <- if (thiazide_limit == 125) 2000 else 1000
          
          ### Models with highest variable importance from full_cont model ........................................ ----
          # create seeds for grf forest generation
          rng_seed <- sample(-999999999:999999999, length(vi_models))
          # run causal forest analysis for n most important covariates
          for (n in vi_models) {
            (function(data, Y_hat, W_hat, sample_weights, continuous, discrete) {
              if (!file.exists(paste0(path, "cfw_vi_test_mod_", n, "_", data_index, ".rds"))) {
                # n most important covariates
                continuous$cf <- cfw_full_vi_comb$short[seq_len(n)] |>
                  intersect(
                    names_index |>
                      unnest(levels, keep_empty = TRUE) |> 
                      filter(is.na(levels) & grepl("^X", short)) |> 
                      pull(short)
                  ) |>
                  (\(x) paste0("^", x, "$"))() |>
                  paste0(collapse = "|") |>
                  (\(x) paste0("(", x, ")"))() |>
                  (\(x) expr(matches(!!x)))() |>
                  (\(x) {
                    if (identical(x, expr(matches("(^)")))) {
                      # match nothing if no continuous covariates in top n VI
                      return(expr(matches("^ABC$"))) 
                    }
                    return(x)
                  })()
                discrete$cf <- cfw_full_vi_comb$short[seq_len(n)] |>
                  intersect(
                    names_index |>
                      unnest(levels, keep_empty = TRUE) |> 
                      filter(!is.na(levels) & grepl("^X", short)) |> 
                      pull(short)
                  ) |>
                  (\(x) paste0("^", x, "$"))() |>
                  paste0(collapse = "|") |>
                  (\(x) paste0("(", x, ")"))() |>
                  (\(x) expr(matches(!!x)))() |>
                  (\(x) {
                    if (identical(x, expr(matches("(^)")))) {
                      # match nothing if no discrete covariates in top n VI
                      return(expr(matches("^ABC$"))) 
                    }
                    return(x)
                  })()
                
                out <- GRFAnalysisWrapper(
                  data,
                  type = "causal",
                  continuous = !!(continuous$cf),
                  discrete = !!(discrete$cf),
                  Y = paste0("Y_", thiazide_limit, "_", horizon), W = "W",
                  Y.hat = Y_hat, W.hat = W_hat,
                  sample.weights = sample_weights,
                  clusters = data[[cluster_name]], # cluster on year of inclusion
                  num.trees = 2000, alpha = 0.0002,
                  tune.parameters = tune_param,
                  tune.num.trees = tune.num.trees, tune.num.reps = tune.num.reps, tune.num.draws = tune.num.draws,
                  seed = rng_seed[which(vi_models == n)],
                  num.threads = floor(num_threads / num_clusters)
                )
                if("X_03_SjÃ¦lland" %in% names(out$X.orig)) {
                  names(out$X.orig)[
                    which(names(out$X.orig) == "X_03_SjÃ¦lland")
                  ] <- "X_03_Sjælland"
                }
                saveRDS(out, paste0(path, "cfw_vi_test_mod_", n, "_", data_index, ".rds"))
                
                # remove causal forest from memory
                rm(list = "out")
              }
              
              return(NULL)
            })(data = data_obs_test,
               Y_hat = Y_hat_test,
               W_hat = W_hat_test,
               sample_weights = sample_weights_test,
               continuous = continuous,
               discrete = discrete)
          }
          
          # empty return
          return(NULL)
        }
      )
    }
  )
}


### Obtain predictions from training set models on imputed validation samples ..................................... ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    seed = as.list(seeds[16:25])
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
      "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
      "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
      "vi_models", "excess_risk_cutoff"
    ),
    seed = TRUE
  ),
  .f = \(path, result_path, seed) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    future_map(
      .x = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c(
          "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
          "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
          "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
          "excess_risk_cutoff", "path", "result_path", "vi_models"
        ),
        seed = seed
      ),
      .f = \(i)  {
        options("parallelly.availableCores.methods" = "system")
        ### Training model CATE estimates on test dataset --------------------------------------------------
        cf_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
        cfw_cate_test <- map(
          .x = seq_along(list.files(path, pattern = "^cfw_full_mod_\\d{1,}")),
          .f = \(j) {
            workers <- min(floor(num_threads / num_clusters), length(vi_models) + 1)
            threads_per_worker <- floor(num_threads / num_clusters / workers)
            if(workers == 1) {
              plan(sequential)
            } else {
              plan(multisession, workers = workers)
            }
            out <- furrr::future_map(
              .x = c(0, vi_models),
              .options = furrr_options(
                packages = c(
                  "dplyr", "stringr", "tibble", "grf", "glue", "stats", "future", "furrr"
                ),
                globals = c("i", "j", "path", "cf_test", "num_threads", "num_clusters", "threads_per_worker"),
                seed = TRUE
              ),
              .f = \(vi_model) {
                if (vi_model == 0) {
                  cf <- readRDS(paste0(path, "cfw_full_mod_", j, ".rds"))
                } else {
                  cf <- readRDS(paste0(path, "cfw_vi_mod_", vi_model, "_", j, ".rds"))
                }
                X_test <- select(cf_test[["X.orig"]], all_of(names(cf[["X.orig"]])))
                predictions <- predict(
                  cf, 
                  newdata = X_test,
                  estimate.variance = TRUE, 
                  num.threads = threads_per_worker
                )
                return(predictions)
              }
            )
            plan(sequential)
            names(out) <- c(
              "full_mod",
              paste0("vi_mod_", vi_models)
            )
            return(out)
          }
        )
        
        # pool CATE estimates from each development model
        cfw_cate_test_pooled <- cfw_cate_test |>
          list_transpose(simplify = FALSE) |>
          future_map(
            .options = furrr_options(
              packages = c("mice", "dplyr", "purrr")
            ),
            .f = \(model) {
              model |>
                list_transpose(simplify = FALSE) |>
                map(\(x) list_transpose(x, simplify = TRUE)) |>
                list_transpose(simplify = FALSE) |>
                (\(x) {
                  n <- length(x)
                  map(
                    x,
                    \(x, n) {
                      if (is.list(x) && length(x[[1]]) > 1) {
                        pool <- pool.scalar(Q = x$predictions, U = x$variance.estimates, n = n, k = 1, rule = "rubin1987")
                        out <- tibble(
                          predictions = pool$qbar,
                          variance_estimates = pool$t
                        )
                      } else {
                        pool <- NULL
                        out <- tibble(
                          predictions = x$predictions,
                          variance_estimates = x$variance.estimates
                        )
                      }
                      return(list(table = out, pool_results = pool))
                    }, 
                    n = n
                  )
                })() |>
                list_transpose()
            }
          )
        
        saveRDS(cfw_cate_test_pooled, paste0(result_path, "cfw_cate_test_", i, ".rds"))
        
        return(NULL)
      }
    )
  }
)

### Obtain sample weights, expected outcome and propensity to treat on imputed validation samples ................. ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    seed = as.list(seeds[26:36])
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
      "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
      "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
      "vi_models", "excess_risk_cutoff"
    ),
    seed = TRUE
  ),
  .f = \(path, result_path, seed) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    future_map(
      .x = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c(
          "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
          "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
          "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
          "excess_risk_cutoff", "path", "result_path", "vi_models"
        ),
        seed = seed
      ),
      .f = \(i)  {
        cfw_full_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
        cfw_sample_weights_test <- cfw_full_test[["sample.weights"]]
        cfw_exp_out_test <- cfw_full_test[["Y.hat"]]
        cfw_exp_trt_test <- cfw_full_test[["W.hat"]]
        # If zero estimate of W.hat, add small amount to avoid degenerate case. 
        # NOTE: Estimated propensities so close to zero indicates problems with the overlap assumption,
        #       so care must be taken when drawing any conclusions!
        cfw_exp_trt_test[cfw_exp_trt_test == 0] <- 1e-6
        
        saveRDS(cfw_sample_weights_test, paste0(result_path, "cfw_sample_weights_test_", i, ".rds"))
        saveRDS(cfw_exp_out_test, paste0(result_path, "cfw_exp_out_test_", i, ".rds"))
        saveRDS(cfw_exp_trt_test, paste0(result_path, "cfw_exp_trt_test_", i, ".rds"))
        
        return(NULL)
      }
    )
  }
)

### Obtain expected potential outcomes in each treatment group ----------------------------------
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    seed = as.list(seeds[37:47])
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c("num_threads", "num_clusters", "names_index", "compute_results",
                "names_index_comb", "mice", "vi_models", "excess_risk_cutoff"),
    seed = TRUE
  ),
  .f = \(path, result_path, seed) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    future_map(
      .x = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra",
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c("num_threads", "num_clusters", "names_index", "compute_results",
                    "names_index_comb", "mice", "excess_risk_cutoff", "path", "result_path", "vi_models"),
        seed = seed
      ),
      .f = \(i)  {
        options("parallelly.availableCores.methods" = "system")
        cfw_outcome_trt_grp_test <- map(
          .x = c(0, vi_models),
          .f = \(vi_model, i, path) {
            if (vi_model == 0) {
              cf_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
            } else {
              cf_test <- readRDS(paste0(path, "cfw_vi_test_mod_", vi_model, "_", i, ".rds"))
            }
            cf_test$W.hat[cf_test$W.hat == 0] <- 1e-6
            p_0_test <- as.numeric(cf_test$Y.hat - cf_test$W.hat * cf_test$predictions)
            p_1_test <- as.numeric(cf_test$Y.hat + (1 - cf_test$W.hat) * cf_test$predictions)
            
            return(
              list(
                p_0_test = p_0_test,
                p_1_test = p_1_test
              )
            )
          },
          i = i,
          path = path
        )
        names(cfw_outcome_trt_grp_test) <- c(
          "full_mod",
          paste0("vi_mod_", vi_models)
        )
        
        saveRDS(cfw_outcome_trt_grp_test, paste0(result_path, "cfw_outcome_trt_grp_test_", i, ".rds"))
        
        return(NULL)
      }
    )
  }
)

### Obtain aipw scores -------------------------------------------------------------------
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    seed = as.list(seeds[48:58])
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c("num_threads", "num_clusters", "names_index", "compute_results",
                "names_index_comb", "mice", "vi_models", "excess_risk_cutoff"),
    seed = TRUE
  ),
  .f = \(path, result_path, seed) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    future_map(
      .x = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c("num_threads", "num_clusters", "names_index", "compute_results",
                    "names_index_comb", "mice", "excess_risk_cutoff", "path", "result_path", "vi_models"),
        seed = seed
      ),
      .f = \(i)  {
        options("parallelly.availableCores.methods" = "system")
        cfw_cate_test <- readRDS(paste0(result_path,"cfw_cate_test_", i, ".rds"))
        names(cfw_cate_test) <- c(
          "full_mod",
          paste0("vi_mod_", vi_models)
        )
        cfw_outcome_trt_grp_test <- readRDS(paste0(result_path, "cfw_outcome_trt_grp_test_", i, ".rds"))
        cfw_aipw_scores_test <- pmap(
          .l = list(
            vi_model = c(0, vi_models),
            tau_hat = map(cfw_cate_test, \(x) x[["table"]][["predictions"]]),
            mu_hat_0 = map(cfw_outcome_trt_grp_test, \(x) x[["p_0_test"]]),
            mu_hat_1 = map(cfw_outcome_trt_grp_test, \(x) x[["p_1_test"]]) 
          ),
          .f = \(vi_model, tau_hat, mu_hat_0, mu_hat_1, i, path, result_path) {
            if (vi_model == 0) {
              cf_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
            } else {
              cf_test <- readRDS(paste0(path, "cfw_vi_test_mod_", vi_model, "_", i, ".rds"))
            }
            cf_test$W.hat[cf_test$W.hat == 0] <- 1e-6
            aipw_scores <- 
              tau_hat +
              cf_test[["W.orig"]] / cf_test[["W.hat"]] * (cf_test[["Y.orig"]] - mu_hat_1) -
              (1 - cf_test[["W.orig"]]) / (1 - cf_test[["W.hat"]]) * (cf_test[["Y.orig"]] - mu_hat_0)
            
            return(aipw_scores)
          },
          i = i,
          path = path,
          result_path = result_path
        )
        names(cfw_aipw_scores_test) <- c(
          "full_mod",
          paste0("vi_mod_", vi_models)
        )
        
        saveRDS(cfw_aipw_scores_test, paste0(result_path, "cfw_aipw_scores_test_", i, ".rds"))
        
        return(NULL)
      }
    )
  }
)

### pool sample weights, expected outcome, propensity to treat, and treatment effect predictions .................. ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    W_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W)
    ),
    Y_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(Y_125_90),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120)
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli","future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c("tih_cohort_test", "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results",
                "names_index_comb", "mice", "vi_models", "excess_risk_cutoff"),
    seed = TRUE
  ),
  .f = \(path, result_path, W_orig, Y_orig) {
    cfw_cate_test <- list()
    cfw_sample_weights_test <- list()
    cfw_exp_out_test <- list()
    cfw_exp_trt_test <- list()
    cfw_outcome_trt_grp_test <- list()
    cfw_aipw_scores_test <- list()
    for (i in seq_along(list.files(result_path, pattern = "^cfw_cate_test_\\d{1,}"))) {
      cfw_cate_test[[i]] <- readRDS(paste0(result_path, "cfw_cate_test_", i, ".rds"))
      names(cfw_cate_test[[i]]) <- c("full_mod", paste0("vi_mod_", vi_models))
      cfw_sample_weights_test[[i]] <- readRDS(paste0(result_path, "cfw_sample_weights_test_", i, ".rds"))
      cfw_exp_out_test[[i]] <- readRDS(paste0(result_path, "cfw_exp_out_test_", i, ".rds"))
      cfw_exp_trt_test[[i]] <- readRDS(paste0(result_path, "cfw_exp_trt_test_", i, ".rds"))
      cfw_outcome_trt_grp_test[[i]] <- readRDS(paste0(result_path, "cfw_outcome_trt_grp_test_", i, ".rds"))
      cfw_aipw_scores_test[[i]] <- readRDS(paste0(result_path, "cfw_aipw_scores_test_", i, ".rds"))
    }
    plan(multisession, workers = min(num_clusters, length(vi_models) + 1))
    cfw_validation_data_agg <- list()
    cfw_validation_data_agg$cate <- cfw_cate_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          map(model, \(x) x$table) |>
            purrr::list_transpose(simplify = FALSE) |>
            map(\(x) purrr::list_transpose(x, simplify = TRUE)) |>
            purrr::list_transpose(simplify = FALSE) |>
            (\(x) {
              n <- length(x)
              map(
                x,
                \(x, n) {
                  if (is.list(x) && length(x[[1]]) > 1) {
                    pool <- pool.scalar(Q = x$predictions, U = x$variance_estimates, n = n, k = 1, rule = "rubin1987")
                    out <- tibble(
                      estimate = pool$qbar,
                      std_err = sqrt(pool$t),
                      lower = estimate - qnorm(0.975) * std_err,
                      upper = estimate + qnorm(0.975) * std_err
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      estimate = x[["predictions"]],
                      std_err = sqrt(x[["variance_estimates"]]),
                      lower = estimate - qnorm(0.975) * std_err,
                      upper = estimate + qnorm(0.975) * std_err
                    )
                  }
                  return(list(table = out, pool_results = pool))
                }, 
                n = 1
              )
            })() |>
            list_transpose()
        }
      )
    cfw_validation_data_agg$sample_weights <- cfw_sample_weights_test |>
      list_transpose(simplify = TRUE) |>
      (\(x) {
        n <- length(x)
        map(
          x,
          \(x, n) {
            if (length(x[1]) > 1) {
              pool <- pool.scalar(Q = x, U = 0, n = n, k = 1, rule = "rubin1987")
              out <- pool$qbar
            } else {
              pool <- NULL
              out <- x[1]
            }
            return(list(table = out, pool_results = pool))
          }, 
          n = n
        )
      })() |>
      list_transpose()
    cfw_validation_data_agg$W_hat <- cfw_exp_trt_test |>
      list_transpose(simplify = TRUE) |>
      (\(x) {
        n <- length(x)
        map(
          x,
          \(x, n) {
            if (length(x[1]) > 1) {
              pool <- pool.scalar(Q = x, U = 0, n = n, k = 1, rule = "rubin1987")
              out <- pool$qbar
            } else {
              pool <- NULL
              out <- x[1]
            }
            return(list(table = out, pool_results = pool))
          }, 
          n = n
        )
      })() |>
      list_transpose()
    cfw_validation_data_agg$Y_hat <- cfw_exp_out_test |>
      list_transpose(simplify = TRUE) |>
      (\(x) {
        n <- length(x)
        map(
          x,
          \(x, n) {
            if (length(x[1]) > 1) {
              pool <- pool.scalar(Q = x, U = 0, n = n, k = 1, rule = "rubin1987")
              out <- pool$qbar
            } else {
              pool <- NULL
              out <- x[1]
            }
            return(list(table = out, pool_results = pool))
          }, 
          n = n
        )
      })() |>
      list_transpose()
    cfw_validation_data_agg$out_trt_grp <- cfw_outcome_trt_grp_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          model |>
            list_transpose(simplify = FALSE) |>
            map(\(x) list_transpose(x, simplify = TRUE)) |>
            list_transpose(simplify = FALSE) |>
            (\(x) {
              n <- length(x)
              map(
                x,
                \(x, n) {
                  if (is.list(x) && length(x[[1]]) > 1) {
                    pool <- list(
                      p_0_test = pool.scalar(Q = x$p_0_test, U = 0, n = n, k = 1, rule = "rubin1987"),
                      p_1_test = pool.scalar(Q = x$p_1_test, U = 0, n = n, k = 1, rule = "rubin1987")
                    )
                    out <- tibble(
                      p_0_test = pool$p_0_test$qbar,
                      p_1_test = pool$p_1_test$qbar
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      p_0_test = x[[1]],
                      p_1_test = x[[2]]
                    )
                  }
                  return(list(table = out, pool_results = pool))
                }, 
                n = n
              )
            })() |>
            list_transpose()
        }
      )
    cfw_validation_data_agg$aipw_scores <- cfw_aipw_scores_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          model |>
            list_transpose(simplify = TRUE) |>
            (\(x) {
              n <- length(x)
              map(
                x,
                \(x, n) {
                  if (length(x[1]) > 1) {
                    pool <- pool.scalar(Q = x, U = 0, n = n, k = 1, rule = "rubin1987")
                    out <- pool$qbar
                  } else {
                    pool <- NULL
                    out <- x[1]
                  }
                  return(list(table = out, pool_results = pool))
                }, 
                n = n
              )
            })() |>
            list_transpose()
        }
      )
    plan(sequential)
    cfw_validation_data <- map(
      seq_along(cfw_validation_data_agg$cate),
      \(i) {
        tibble(
          sample_weights = cfw_validation_data_agg$sample_weights$table, # obtained from test model
          W_obs = W_orig, # test data
          W_hat = cfw_validation_data_agg$W_hat$table, # obtained from test model
          Y_obs = Y_orig, # test data
          Y_hat = cfw_validation_data_agg$Y_hat$table, # obtained from test model
          aipw_scores = (cfw_validation_data_agg$aipw_scores)[[i]]$table, # tau_hat from training mod, rest from test mod
          tau_hat = (cfw_validation_data_agg$cate)[[i]]$table$estimate, # obtained from training model
          tau_hat_var = ((cfw_validation_data_agg$cate)[[i]]$table$std_err)^2, # obtained from training model
          p_0 = (cfw_validation_data_agg$out_trt_grp)[[i]]$table$p_0_test, # obtained from test model
          p_1 = (cfw_validation_data_agg$out_trt_grp)[[i]]$table$p_1_test # obtained from test model
        )
      }
    ) |>
      structure(names = names(cfw_validation_data_agg$cate))
    saveRDS(cfw_validation_data, paste0(result_path, "cfw_validation_data.rds"))
  }
)

### expected loss using R-loss criterion .......................................................................... ----
# This error measure is suggested in Shuler et.al. (2018) as a way of model selection by choosing the 
# model with the lowest mean R-loss. They call this measure \tau-risk_R. Note that in the validation sample,
# we use estimates of mean outcome and propensity to treat estimated on the validation sample. Thus the only
# part of the R-loss derived from the training sample are the treatment effect estimates.
plan(sequential)
pmap(
  list(
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    )
  )|> map(\(x) x[setup_index]),
  \(result_path, table_path) {
    cfw_validation_data <- readRDS(paste0(result_path, "cfw_validation_data.rds"))
    cfw_tau_risk_R <- map_dbl(
      cfw_validation_data,
      \(tbl) {
        pred <- tbl$tau_hat * (tbl$W_obs - tbl$W_hat)
        obs <- tbl$Y_obs - tbl$Y_hat
        sq_err <- (pred - obs)^2
        return(weighted.mean(sq_err, tbl$sample_weights))
      }
    )
    saveRDS(cfw_tau_risk_R, paste0(result_path, "cfw_tau_risk_R.rds"))
    
    # write table with tau_risk_R
    cfw_tau_risk_R |>
      (\(x) {
        tibble(
          model = names(x),
          tau_risk_R = as.numeric(x),
          relative_error = tau_risk_R / min(tau_risk_R)
        ) 
      }
      )() |>
      (
        \(x) {
          dfs2xlsx(
            list(`tau_risk_R` = x), 
            paste0(table_path, "cfw_tau_risk_R.xlsx"), 
            rowNames = FALSE
          )
        }
      )()
    
    return(NULL)
  }
)

### discrimination performance on validation set .................................................................. ----
# calculate C-for-benefit using aggregated numbers
plan(sequential)
pmap(
  .l = list(
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
  )|> map(\(x) x[setup_index]),
  .f = \(result_path, table_path) {
    cfw_validation_data <- readRDS(paste0(result_path, "cfw_validation_data.rds"))
    plan(multisession, workers = min(num_clusters, length(cfw_validation_data)))
    
    # Model-based c-for-benefit
    cfw_mbcfb_test <- future_pmap(
      .l = list(
        val_dat = cfw_validation_data
      ),
      .options = furrr_options(
        packages = c("grf", "dplyr", "rlang", "glue", "stats", "Rcpp"),
        globals = c("num_threads", "num_clusters"),
        seed = TRUE
      ),
      .f = \(val_dat) {
        Rcpp::sourceCpp("functions/mbcb.cpp")
        cindex <- mbcb(val_dat$p_0, val_dat$p_1, val_dat$tau_hat, val_dat$sample_weights)
        c_for_benefit <- cindex$`C Index`
        c_for_benefit_se <- cindex$S.D. / 2
        lower_CI <- c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se
        upper_CI <- c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
        return(
          list(
            c_for_benefit = c_for_benefit,
            c_for_benefit_se = c_for_benefit_se,
            lower_CI = lower_CI,
            upper_CI = upper_CI,
            cfb_tbl = dplyr::tibble(
              "ci_{0.5 - 0.95 / 2}" := lower_CI,
              "estimate" = c_for_benefit,
              "ci_{0.5 + 0.95 / 2}" := upper_CI
            )
          )
        )
      }
    )
    saveRDS(cfw_mbcfb_test, paste0(result_path, "cfw_mbcfb_test.rds"))
    
    # CATE matched c-for-benefit
    cfw_cfb_cate_test <- future_pmap(
      .l = list(
        val_dat = cfw_validation_data
      ),
      .options = furrr_options(
        packages = c("grf", "dplyr", "rlang", "glue", "stats", "Rcpp"),
        globals = c("num_threads", "num_clusters"),
        seed = TRUE
      ),
      .f = \(val_dat) {
        Rcpp::sourceCpp("functions/rcorr.cpp")
        matched <- MatchIt::matchit(
          W_obs ~ tau_hat,
          data = val_dat,
          method = "nearest",
          distance = "mahalanobis",
          estimand = "ATT"
        )
        matched_patients <- MatchIt::match.data(matched)
        matched_patients$subclass <- as.numeric(matched_patients$subclass)
        matched_patients <- dplyr::as_tibble(matched_patients)
        matched_patients <- matched_patients |>
          dplyr::arrange(subclass, W_obs)
        observed_te <- matched_patients |>
          dplyr::select(subclass, Y_obs, W_obs, weights) |>
          dplyr::summarise(
            Y = sum(weights * Y_obs) / sum(weights),
            .by = c(subclass, W_obs)
          ) |>
          dplyr::summarise(
            matched_tau_obs = sum((W_obs == 1) * Y - (W_obs == 0) * Y),
            .by = subclass
          )
        matched_patients <- matched_patients |>
          inner_join(observed_te, by = "subclass")
        matched_patients <- matched_patients |>
          dplyr::mutate(
            matched_p_0 = sum((W_obs == 0) * weights * p_0) / sum((W_obs == 0) * weights),
            matched_p_1 = sum((W_obs == 1) * weights * p_1) / sum((W_obs == 1) * weights),
            matched_tau_hat = matched_p_1 - matched_p_0,
            .by = "subclass"
          )
        cindex <- rcorr(
          matched_patients$matched_tau_hat[seq(1, nrow(matched_patients), 2)],
          matched_patients$matched_tau_obs[seq(1, nrow(matched_patients), 2)],
          matched_patients$sample_weights[seq(1, nrow(matched_patients), 2)]
        )
        c_for_benefit <- cindex$`C Index`
        c_for_benefit_se <- cindex$S.D. / 2
        lower_CI <- c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se
        upper_CI <- c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
        return(
          list(
            c_for_benefit = c_for_benefit,
            c_for_benefit_se = c_for_benefit_se,
            lower_CI = lower_CI,
            upper_CI = upper_CI,
            cfb_tbl = dplyr::tibble(
              "ci_{0.5 - 0.95 / 2}" := lower_CI,
              "estimate" = c_for_benefit,
              "ci_{0.5 + 0.95 / 2}" := upper_CI
            )
          )
        )
      }
    )
    saveRDS(cfw_cfb_cate_test, paste0(result_path, "cfw_cfb_cate_test.rds"))
    
    # Control-risk matched c-for-benefit
    cfw_cfb_control_risk_test <- future_pmap(
      .l = list(
        val_dat = cfw_validation_data
      ),
      .options = furrr_options(
        packages = c("grf", "dplyr", "rlang", "glue", "stats", "Rcpp"),
        globals = c("num_threads", "num_clusters"),
        seed = TRUE
      ),
      .f = \(val_dat) {
        Rcpp::sourceCpp("functions/rcorr.cpp")
        matched <- MatchIt::matchit(
          W_obs ~ p_0,
          data = val_dat,
          method = "nearest",
          distance = "mahalanobis",
          estimand = "ATT"
        )
        matched_patients <- MatchIt::match.data(matched)
        matched_patients$subclass <- as.numeric(matched_patients$subclass)
        matched_patients <- dplyr::as_tibble(matched_patients)
        matched_patients <- matched_patients |>
          dplyr::arrange(subclass, W_obs)
        observed_te <- matched_patients |>
          dplyr::select(subclass, Y_obs, W_obs, weights) |>
          dplyr::summarise(
            Y = sum(weights * Y_obs) / sum(weights),
            .by = c(subclass, W_obs)
          ) |>
          dplyr::summarise(
            matched_tau_obs = sum((W_obs == 1) * Y - (W_obs == 0) * Y),
            .by = subclass
          )
        matched_patients <- matched_patients |>
          inner_join(observed_te, by = "subclass")
        matched_patients <- matched_patients |>
          dplyr::mutate(
            matched_tau_hat = sum((W_obs == 1) * weights * tau_hat) / sum((W_obs == 1) * weights),
            .by = "subclass"
          )
        cindex <- rcorr(
          matched_patients$matched_tau_hat[seq(1, nrow(matched_patients), 2)],
          matched_patients$matched_tau_obs[seq(1, nrow(matched_patients), 2)],
          matched_patients$sample_weights[seq(1, nrow(matched_patients), 2)]
        )
        c_for_benefit <- cindex$`C Index`
        c_for_benefit_se <- cindex$S.D. / 2
        lower_CI <- c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se
        upper_CI <- c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
        return(
          list(
            c_for_benefit = c_for_benefit,
            c_for_benefit_se = c_for_benefit_se,
            lower_CI = lower_CI,
            upper_CI = upper_CI,
            cfb_tbl = dplyr::tibble(
              "ci_{0.5 - 0.95 / 2}" := lower_CI,
              "estimate" = c_for_benefit,
              "ci_{0.5 + 0.95 / 2}" := upper_CI
            )
          )
        )
      }
    )
    saveRDS(cfw_cfb_control_risk_test, paste0(result_path, "cfw_cfb_control_risk_test.rds"))
    
    # combined control-risk and exposed-risk matched c-for-benefit
    cfw_cfb_combined_test <- future_pmap(
      .l = list(
        val_dat = cfw_validation_data
      ),
      .options = furrr_options(
        packages = c("grf", "dplyr", "rlang", "glue", "stats", "Rcpp"),
        globals = c("num_threads", "num_clusters"),
        seed = TRUE
      ),
      .f = \(val_dat) {
        Rcpp::sourceCpp("functions/rcorr.cpp")
        matched_trt <- MatchIt::matchit(
          W_obs ~ p_0,
          data = val_dat,
          method = "nearest",
          distance = "mahalanobis",
          estimand = "ATT",
          replace = TRUE
        )
        matched_ctr <- MatchIt::matchit(
          W_obs ~ p_1,
          data = val_dat,
          method = "nearest",
          distance = "mahalanobis",
          estimand = "ATC",
          replace = TRUE
        )
        matched_patients_trt <- val_dat |>
          filter(matched_trt$treat == 1) |>
          mutate(subclass = seq_len(n())) |>
          bind_rows(
            val_dat |>
              slice(as.integer(matched_trt$match.matrix[, 1])) |>
              mutate(subclass = seq_len(n()))
          ) |>
          arrange(subclass, W_obs) |>
          mutate(
            sample_weights = sum((W_obs == 1) * sample_weights),
            .by = subclass
          ) |>
          mutate(
            Y_obs = Y_obs - p_0,
            weights = 1
          )
        matched_patients_ctr <- val_dat |>
          filter(matched_ctr$treat == 0) |>
          mutate(subclass = seq_len(n()) + sum(val_dat$W_obs == 1)) |>
          bind_rows(
            val_dat |>
              slice(as.integer(matched_ctr$match.matrix[, 1])) |>
              mutate(subclass = seq_len(n()) + sum(val_dat$W_obs == 1))
          ) |>
          arrange(subclass, W_obs) |>
          mutate(
            sample_weights = sum((W_obs == 0) * sample_weights),
            .by = subclass
          ) |>
          mutate(
            Y_obs = Y_obs - p_1,
            weights = 1
          )
        matched_patients <- bind_rows(
          matched_patients_trt,
          matched_patients_ctr
        ) |>
          dplyr::arrange(subclass, W_obs)
        observed_te <- matched_patients |>
          dplyr::select(subclass, Y_obs, W_obs, weights) |>
          dplyr::summarise(
            Y_obs = sum(weights * Y_obs) / sum(weights),
            .by = c(subclass, W_obs)
          ) |>
          dplyr::summarise(
            matched_tau_obs = sum((W_obs == 1) * Y_obs - (W_obs == 0) * Y_obs),
            .by = subclass
          )
        matched_patients <- matched_patients |>
          inner_join(observed_te, by = "subclass")
        matched_patients <- matched_patients |>
          dplyr::mutate(
            matched_p_0 = sum((W_obs == 0) * weights * p_0) / sum((W_obs == 0) * weights),
            matched_p_1 = sum((W_obs == 1) * weights * p_1) / sum((W_obs == 1) * weights),
            matched_tau_hat = matched_p_1 - matched_p_0,
            .by = "subclass"
          )
        matched_obs_pred <- matched_patients |> 
          slice_head(n = 1, by = subclass) |> 
          select("matched_tau_hat", "matched_tau_obs", "sample_weights")
        cindex <- rcorr(
          matched_obs_pred |> pull("matched_tau_hat"),
          matched_obs_pred |> pull("matched_tau_obs"),
          matched_obs_pred |> pull("sample_weights")
        )
        c_for_benefit <- cindex$`C Index`
        c_for_benefit_se <- cindex$S.D. / 2
        lower_CI <- c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se
        upper_CI <- c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
        return(
          list(
            c_for_benefit = c_for_benefit,
            c_for_benefit_se = c_for_benefit_se,
            lower_CI = lower_CI,
            upper_CI = upper_CI,
            cfb_tbl = dplyr::tibble(
              "ci_{0.5 - 0.95 / 2}" := lower_CI,
              "estimate" = c_for_benefit,
              "ci_{0.5 + 0.95 / 2}" := upper_CI
            )
          )
        )
      }
    )
    saveRDS(cfw_cfb_combined_test, paste0(result_path, "cfw_cfb_combined_test.rds"))
    
    plan(sequential)
    
    # write table with c-for-benefit
    list(
      `CATE matched` = cfw_cfb_cate_test |>
        map(\(x) as_tibble(x[1:4])) |>
        list_rbind() |>
        mutate(model = names(cfw_cfb_cate_test), .before = 1),
      `Control risk matched` = cfw_cfb_control_risk_test |>
        map(\(x) as_tibble(x[1:4])) |>
        list_rbind() |>
        mutate(model = names(cfw_cfb_control_risk_test), .before = 1),
      `Combined` = cfw_cfb_combined_test |>
        map(\(x) as_tibble(x[1:4])) |>
        list_rbind() |>
        mutate(model = names(cfw_cfb_combined_test), .before = 1),
      `Model based` = cfw_mbcfb_test |>
        map(\(x) as_tibble(x[1:4])) |>
        list_rbind() |>
        mutate(model = names(cfw_mbcfb_test), .before = 1)
    ) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_c_for_benefit_test.xlsx"
        ), 
        rowNames = FALSE
      )
    
    return(NULL)
  }
)

# Calibration performance on validation set ....................................................................... ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_tiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_tiazide_no_lab/"),
      paste0(figures_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_tiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_tiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_tiazide_cc_comp/1year/")
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid",
      "future", "furrr", "future.callr", "callr", "MatchIt"
    ),
    globals = c("n_rankings", "vcovHC"),
    seed = TRUE
  ),
  .f = \(path, result_path, table_path, figure_path) {
    cfw_validation_data <- readRDS(paste0(result_path, "cfw_validation_data.rds"))
    plan(multisession, workers = min(length(cfw_validation_data), num_clusters))
    cfw_dsg_test <- future_imap(
      .x = cfw_validation_data,
      .options = furrr_options(
        packages = c("dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid", "MatchIt"),
        globals = c("n_rankings", "vcovHC"),
        seed = TRUE
      ),
      .f = \(data_tbl, nm) {
        out <- list()
        for (i in n_rankings) {
          # Ranking by treatment effect size -----
          tau_hat_quantiles <- DescTools::Quantile(
            x = data_tbl$tau_hat,
            weights = data_tbl$sample_weights,
            probs = seq(0, 1, by = 1 / i)
          )
          
          if (!any(abs(tau_hat_quantiles - lag(tau_hat_quantiles)) < sqrt(.Machine$double.eps), na.rm = TRUE)) {
            data_tbl <- data_tbl |>
              mutate(
                "ranking_{i}" := cut(
                  tau_hat, 
                  breaks = tau_hat_quantiles,
                  include.lowest = TRUE,
                  labels = seq_len(i)
                )
              )
          } else {
            len <- length(data_tbl$tau_hat)
            ranking <- data_tbl$tau_hat |>
              (\(x) {
                dplyr::tibble(
                  id = seq_along(x),
                  tau_hat = x
                )
              }
              )() |>
              dplyr::arrange(.data$tau_hat) |>
              dplyr::mutate(
                id_2 = seq_along(.data$tau_hat),
                rank = cut(
                  x = .data$id_2,
                  breaks = c(
                    seq(0, len %% i, by = 1) *
                      (len %/% i + 1),
                    seq(
                      len %% i + 1,
                      i,
                      length.out = i - len %% i
                    ) *
                      len %/% i + len %% i
                  ),
                  include.lowest = TRUE,
                  labels = seq_len(i)
                )
              ) |>
              dplyr::arrange(.data$id) |>
              dplyr::pull(.data$rank)
            data_tbl <- mutate(data_tbl, "ranking_{i}" := ranking)
          }
          
          # ATE in subgroups ----
          # fit linear model of aipw scores to find average in each rank
          ranking <- paste0("ranking_", i)
          aipw <- lm(data_tbl$aipw_scores ~ 0 + factor(data_tbl[[ranking]]), weights = data_tbl$sample_weights)
          ols <- lm(data_tbl$tau_hat ~ 0 + factor(data_tbl[[ranking]]), weights = data_tbl$sample_weights)
          forest_rank_ate <- dplyr::tibble(
            method = rep(c("aipw", "ols"), each = i),
            ranking = rep(paste0("Q", seq_len(i)), 2),
            estimate = c(coef(aipw), coef(ols)),
            std_err = c(sqrt(diag(vcov(aipw))), sqrt(diag(vcov(aipw))))
          )
          
          # plot with estimates and 95% confidence intervals within each ranking:
          forest_rank_ate_plot <- list(
            data = forest_rank_ate,
            plot = expr(
              ggplot2::ggplot(
                !!forest_rank_ate,
                ggplot2::aes(
                  x = factor(.data$ranking, levels = paste0("Q", seq_len(!!i))), 
                  y = .data$estimate,
                  color = .data$method
                )
              ) +
                ggplot2::geom_point(position = position_dodge(width = 0.3)) +
                ggplot2::geom_errorbar(
                  ggplot2::aes(
                    ymin = .data$estimate + qnorm(0.025) * .data$std_err,
                    ymax = .data$estimate + qnorm(0.975) * .data$std_err
                  ),
                  width = 0.2,
                  position = position_dodge(width = 0.3)
                ) + 
                ggplot2::xlab("") + 
                ggplot2::ylab("") + 
                ggplot2::ggtitle(
                  "OLS/AIPW score within each ranking (defined by predicted CATE)"
                )+
                ggpubr::theme_pubr() +
                ggplot2::theme(
                  plot.title = element_text(size = 11)
                ) 
            )
          )
          
          # table with tests for differences between ranking groups
          forest_rank_diff_test <- dplyr::tibble()
          for (j in seq_len(i - 1)) {
            lev <- seq_len(i)
            aipw <- lm(
              data_tbl$aipw_scores ~
                1 + factor(data_tbl[[ranking]], levels = c(lev[j], lev[-j]))
            )
            forest_rank_diff_test <- coef(summary(aipw))[
              seq(j + 1, i),
              c(1, 2, 4),
              drop = FALSE
            ] |>
              dplyr::as_tibble() |>
              dplyr::mutate(
                id = paste("Rank", seq(j + 1, i), "- Rank ", j)
              ) |>
              (\(x) rbind(forest_rank_diff_test, x))()
          }
          
          # Adjust for multiple testing using the Benjamini-Hockberg procedure
          forest_rank_diff_test <- forest_rank_diff_test |>
            dplyr::rename("Orig. p-value" = "Pr(>|t|)") |>
            dplyr::mutate(
              `95% CI` = paste0(
                "(",
                sprintf(
                  fmt = "%.3f",
                  .data$Estimate + qnorm(0.025) * .data$`Std. Error`
                ),
                ", ",
                sprintf(
                  fmt = "%.3f",
                  .data$Estimate + qnorm(0.975) * .data$`Std. Error`
                ),
                ")"
              ),
              `Orig. p-value` = sprintf(
                fmt = "%.3f",
                .data$`Orig. p-value`
              ),
              `Adj. p-value` = sprintf(
                fmt = "%.3f",
                p.adjust(.data$`Orig. p-value`, method = "BH")
              )
            ) |>
            dplyr::select(
              "id",
              "Estimate",
              "Std. Error",
              "95% CI",
              dplyr::everything()
            )
          
          # calibration for benefit ----
          # match patients in each ranking group
          data_tbl_matched <- list(
            CATE = map(
              levels(data_tbl[[ranking]]),
              \(lvl) {
                matched <- MatchIt::matchit(
                  W_obs ~ tau_hat,
                  data = filter(data_tbl, !!sym(ranking) == lvl),
                  method = "quick",
                  distance = "mahalanobis",
                  estimand = "ATE",
                  s.weights = ~sample_weights
                )
                matched_patients <- MatchIt::match.data(matched)
                matched_patients$subclass <- as.numeric(matched_patients$subclass)
                return(dplyr::as_tibble(matched_patients))
              }
            ) |>
              list_rbind() |>
              arrange(subclass, ranking),
            control_risk = map(
              levels(data_tbl[[ranking]]),
              \(lvl) {
                matched <- MatchIt::matchit(
                  W_obs ~ p_0,
                  data = filter(data_tbl, !!sym(ranking) == lvl),
                  method = "quick",
                  distance = "mahalanobis",
                  estimand = "ATE",
                  s.weights = ~sample_weights
                )
                matched_patients <- MatchIt::match.data(matched)
                matched_patients$subclass <- as.numeric(matched_patients$subclass)
                return(dplyr::as_tibble(matched_patients))
              }
            ) |>
              list_rbind() |>
              mutate(Y_adj = Y_obs - p_0) |>
              arrange(subclass, ranking),
            combined = map(
              levels(data_tbl[[ranking]]),
              \(lvl) {
                matched_trt <- MatchIt::matchit(
                  W_obs ~ p_0,
                  data = filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 0),
                  method = "nearest",
                  distance = "mahalanobis",
                  estimand = "ATT",
                  replace = TRUE
                )
                matched_ctr <- MatchIt::matchit(
                  W_obs ~ p_1,
                  data = filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 1),
                  method = "nearest",
                  distance = "mahalanobis",
                  estimand = "ATC",
                  replace = TRUE
                )
                matched_patients_trt <- filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 0) |>
                  filter(matched_trt$treat == 1) |>
                  mutate(subclass = seq_len(n())) |>
                  bind_rows(
                    filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 0) |>
                      slice(as.integer(matched_trt$match.matrix[, 1])) |>
                      mutate(subclass = seq_len(n()))
                  ) |>
                  arrange(subclass, W_obs) |>
                  mutate(
                    tau_hat = sum((W_obs == 1) * tau_hat),
                    aipw_scores = sum((W_obs == 1) * aipw_scores),
                    weights = sum((W_obs == 1) * sample_weights),
                    .by = subclass
                  ) |>
                  mutate(Y_adj = Y_obs - p_0)
                matched_patients_ctr <- filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 1) |>
                  filter(matched_ctr$treat == 0) |>
                  mutate(subclass = seq_len(n()) + sum(filter(data_tbl, !!sym(ranking) == lvl)$W_obs == 1)) |>
                  bind_rows(
                    filter(data_tbl, (!!sym(ranking) == lvl) | W_obs == 1) |>
                      slice(as.integer(matched_ctr$match.matrix[, 1])) |>
                      mutate(subclass = seq_len(n()) + sum(filter(data_tbl, !!sym(ranking) == lvl)$W_obs == 1))
                  ) |>
                  arrange(subclass, W_obs) |>
                  mutate(
                    tau_hat = sum((W_obs == 0) * tau_hat),
                    aipw_scores = sum((W_obs == 0) * aipw_scores),
                    weights = sum((W_obs == 0) * sample_weights),
                    .by = subclass
                  ) |>
                  mutate(Y_adj = Y_obs - p_1)
                matched_patients <- bind_rows(
                  matched_patients_trt,
                  matched_patients_ctr
                )
                return(matched_patients)
              }
            ) |>
              list_rbind() |>
              arrange(subclass, ranking)
          )
          
          # Mean bias
          # expected difference between observed effect and predicted effect in matched sample
          mean_bias <- map(
            data_tbl_matched,
            \(tbl) {
              tbl |> 
                summarise(
                  observed_treatment_effect = 
                    sum((W_obs == 1) * weights * Y_obs) / sum((W_obs == 1) * weights) - 
                    sum((W_obs == 0) * weights * Y_obs) / sum((W_obs == 0) * weights),
                  expected_treatment_effect = sum(
                    (W_obs == 1) * weights * p_1 / sum((W_obs == 1) * weights) -
                      (W_obs == 0) * weights * p_0 / sum((W_obs == 0) * weights)
                  ),
                  mean_bias = observed_treatment_effect - expected_treatment_effect
                )
            }
          )
          
          # calculate adjusted expected effect (expected outcome under treatment of treated minus expected outcome under control for non-treated)
          # Partition AIPW score into components for treated and untreated individual in each pair:
          # Treated: mu_1(X) + e(X)^-1 * (Y - mu_1(X)) 
          # Untreated: mu_0(X) + (1 - e(X))^-1 * (Y - mu_0(X)) 
          # 
          # The expected effect is adjusted to better match the observed benefit, which for each matched pair
          # consists of the treatment effect plus the difference in untreated risk.
          # The standard error of expected effect is based on the standard deviation of predicted effect from
          # the causal forest model within the matched samples (by ranking). 
          # Control risk matched:
          # When matched on control risk, mu_hat_0 will be similar within a matched pair, so the expected effect
          # should be a good representation of the true effect. This is in contrast to CATE matching, where the
          # difference in control risk will make up a component of the expected effect. 
          expected_treatment_effect <- list()
          expected_treatment_effect$ols <- list(
            CATE = data_tbl_matched$CATE |> 
              summarise(
                expected_treatment_effect = 
                  sum((W_obs == 1) * weights * p_1 / sum((W_obs == 1) * weights)) -
                  sum((W_obs == 0) * weights * p_0 / sum((W_obs == 0) * weights)),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              ),
            control_risk = data_tbl_matched$control_risk |> 
              summarise(
                expected_treatment_effect = weighted.mean(tau_hat, weights),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              ),
            combined = data_tbl_matched$combined |> 
              summarise(
                tau_hat = mean(tau_hat),
                aipw_scores = mean(aipw_scores),
                weights = mean(weights),
                .by = c(all_of(ranking), subclass)
              ) |>
              summarise(
                expected_treatment_effect = weighted.mean(tau_hat, weights),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              )
          )
          expected_treatment_effect$aipw <- list(
            CATE = data_tbl_matched$CATE |> 
              summarise(
                expected_treatment_effect = 
                  sum((W_obs == 1) * weights * p_1 / sum((W_obs == 1) * weights)) -
                  sum((W_obs == 0) * weights * p_0 / sum((W_obs == 0) * weights)),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              ),
            control_risk = data_tbl_matched$control_risk |> 
              summarise(
                expected_treatment_effect = weighted.mean(aipw_scores, weights),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              ),
            combined = data_tbl_matched$combined |> 
              summarise(
                tau_hat = mean(tau_hat),
                aipw_scores = mean(aipw_scores),
                weights = mean(weights),
                .by = c(all_of(ranking), subclass)
              ) |>
              summarise(
                expected_treatment_effect = weighted.mean(aipw_scores, weights),
                expected_treatment_effect_se = sqrt(sum(weights^2 * (aipw_scores - tau_hat)^2) / sum(weights)^2),
                .by = all_of(ranking)
              )
          )
          
          # calculate observed effect (average of observed effect for each matched pair within a ranking group)
          # 
          # CATE matched:
          # The observed effect is calculated as the mean difference between the two treatment 
          # arms in the matched sample. The standard error is based on Y having a bernoulli distribution,
          # such that sum((W==1) * Y) is binomially distributed with estimated parameter 
          # sum((W_obs == 1) * Y_obs) / sum(W_obs == 1).
          # Then the variance is 
          # sum(W_obs == 1) * 
          # sum((W_obs == 1) * Y_obs) / sum(W_obs == 1) * 
          # sum((W_obs == 1) * (1 - Y_obs)) / sum(W_obs == 1)
          # Combining the two treatment arms the standard error of the estimate of the observed effect is
          # sqrt(sum((W == 1) * (Y == 1)) * sum((W == 1) * (Y == 0)) / sum(W == 1)^3 +
          #        sum((W == 0) * (Y == 1)) * sum((W == 0) * (Y == 0)) / sum(W == 0)^3)
          # 
          # Control risk matched:
          # The observed effect is adjusted by the control risk for each individual in a matched pair. Since we have
          # matched the pair based on control risk, these adjustments will (approximately) cancel out, giving us the
          # same observed effects as with CATE matching.
          observed_treatment_effect <- list(
            CATE = data_tbl_matched$CATE |>
              summarise(
                observed_treatment_effect = 
                  sum((W_obs == 1) * weights * Y_obs) / sum((W_obs == 1) * weights) - 
                  sum((W_obs == 0) * weights * Y_obs) / sum((W_obs == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W_obs == 1) * weights * (Y_obs == 1)) * sum((W_obs == 1) * weights * (Y_obs == 0)) / sum((W_obs == 1) * weights)^3 +
                    sum((W_obs == 0) * weights * (Y_obs == 1)) * sum((W_obs == 0) * weights * (Y_obs == 0)) / sum((W_obs == 0) * weights)^3
                ),
                .by = all_of(ranking)
              ),
            control_risk = data_tbl_matched$control_risk |>
              summarise(
                observed_treatment_effect = 
                  sum((W_obs == 1) * weights * Y_adj) / sum((W_obs == 1) * weights) - 
                  sum((W_obs == 0) * weights * Y_adj) / sum((W_obs == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W_obs == 1) * weights * (Y_obs == 1)) * sum((W_obs == 1) * weights * (Y_obs == 0)) / sum((W_obs == 1) * weights)^3 +
                    sum((W_obs == 0) * weights * (Y_obs == 1)) * sum((W_obs == 0) * weights * (Y_obs == 0)) / sum((W_obs == 0) * weights)^3
                ),
                .by = all_of(ranking)
              ),
            combined = data_tbl_matched$combined |>
              summarise(
                observed_treatment_effect = 
                  sum((W_obs == 1) * weights * Y_adj) / sum((W_obs == 1) * weights) - 
                  sum((W_obs == 0) * weights * Y_adj) / sum((W_obs == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W_obs == 1) * weights * (Y_obs == 1)) * sum((W_obs == 1) * weights * (Y_obs == 0)) / sum((W_obs == 1) * weights)^3 +
                    sum((W_obs == 0) * weights * (Y_obs == 1)) * sum((W_obs == 0) * weights * (Y_obs == 0)) / sum((W_obs == 0) * weights)^3
                ),
                .by = all_of(ranking)
              )
          )
          
          # calibration plot
          calibration_data <- map(
            expected_treatment_effect,
            \(data) {
              map2(
                data,
                observed_treatment_effect,
                \(exp, obs) inner_join(exp, obs, by = ranking)
              )
            }
          )
          
          calibration_coef <- map(
            expected_treatment_effect,
            \(data) {
              map2(
                data,
                observed_treatment_effect,
                \(exp, obs) {
                  coef(lm(obs$observed_treatment_effect ~ exp$expected_treatment_effect)) |>
                    structure(names = c("\U03B2_0", "\U03B2_1"))
                }
              )
            }
          )
          
          # confidence intervals on parameters from calibration linear regression
          calibration_confint <- map(
            expected_treatment_effect,
            \(data) {
              map2(
                data,
                observed_treatment_effect,
                \(exp, obs) {
                  confint(lm(obs$observed_treatment_effect ~ exp$expected_treatment_effect), level = 0.95) |>
                    structure(names = c("\U03B2_0", "\U03B2_1"))
                }
              )
            }
          )
          
          calibration_plot <- map2(
            calibration_data,
            calibration_coef,
            \(data, coef) {
              map2(
                data,
                coef,
                \(data, coef) {
                  expr(
                    !!data |>
                      ggplot(aes(x = expected_treatment_effect, y = observed_treatment_effect)) +
                      geom_point() +
                      geom_errorbar(
                        aes(
                          ymin = observed_treatment_effect - observed_treatment_effect_se,
                          ymax = observed_treatment_effect + observed_treatment_effect_se,
                          width = 0.001
                        )
                      ) +
                      geom_errorbarh(
                        aes(
                          xmin = expected_treatment_effect - expected_treatment_effect_se,
                          xmax = expected_treatment_effect + expected_treatment_effect_se,
                          height = 0.001
                        )
                      ) +
                      geom_abline(slope = 1, intercept = 0, linetype = 2) +
                      scale_x_continuous(
                        breaks = seq(-0.03, 0.15, 0.03),
                        labels = paste0(seq(-3, 15, 3), "%")
                      ) +
                      scale_y_continuous(
                        breaks = seq(-0.03, 0.15, 0.03),
                        labels = paste0(seq(-3, 15, 3), "%")
                      ) +
                      xlab("Predicted excess risk") + ylab("Observed excess risk") +
                      coord_cartesian(xlim = c(-0.05, 0.15), ylim = c(-0.05, 0.15))
                  )
                }
              )
            }
          )
          
          out[[ranking]] <- list(
            forest_rank_ate = forest_rank_ate,
            forest_rank_ate_plot = forest_rank_ate_plot,
            forest_rank_diff_test = forest_rank_diff_test,
            mean_bias = mean_bias,
            calibration_data = calibration_data,
            calibration_plot = calibration_plot,
            calibration_coef = calibration_coef,
            calibration_confint = calibration_confint
          )
        }
        
        # calibration with individual observed effects
        data_matched <- list(
          combined =
            (\() {
              matched_trt <- MatchIt::matchit(
                W_obs ~ p_0,
                data = data_tbl,
                method = "nearest",
                distance = "mahalanobis",
                estimand = "ATT",
                replace = TRUE
              )
              matched_ctr <- MatchIt::matchit(
                W_obs ~ p_1,
                data = data_tbl,
                method = "nearest",
                distance = "mahalanobis",
                estimand = "ATC",
                replace = TRUE
              )
              matched_patients_trt <- data_tbl |>
                filter(matched_trt$treat == 1) |>
                mutate(subclass = seq_len(n())) |>
                bind_rows(
                  data_tbl |>
                    slice(as.integer(matched_trt$match.matrix[, 1])) |>
                    mutate(subclass = seq_len(n()))
                ) |>
                arrange(subclass, W_obs) |>
                mutate(
                  tau_hat = sum((W_obs == 1) * tau_hat),
                  aipw_scores = sum((W_obs == 1) * aipw_scores),
                  weights = sum((W_obs == 1) * sample_weights),
                  .by = subclass
                ) |>
                mutate(Y_adj = Y_obs - p_0)
              matched_patients_ctr <- data_tbl |>
                filter(matched_ctr$treat == 0) |>
                mutate(subclass = seq_len(n()) + sum(data_tbl$W_obs == 1)) |>
                bind_rows(
                  data_tbl |>
                    slice(as.integer(matched_ctr$match.matrix[, 1])) |>
                    mutate(subclass = seq_len(n()) + sum(data_tbl$W_obs == 1))
                ) |>
                arrange(subclass, W_obs) |>
                mutate(
                  tau_hat = sum((W_obs == 0) * tau_hat),
                  aipw_scores = sum((W_obs == 0) * aipw_scores),
                  weights = sum((W_obs == 0) * sample_weights),
                  .by = subclass
                ) |>
                mutate(Y_adj = Y_obs - p_1)
              matched_patients <- bind_rows(
                matched_patients_trt,
                matched_patients_ctr
              )
              return(matched_patients)
            })()
        )
        expected_treatment_effect <- list(
          combined = data_matched$combined |> 
            summarise(
              weights = mean(weights),
              expected_treatment_effect = mean(tau_hat),
              .by = "subclass"
            )
        )
        observed_treatment_effect <- list(
          combined = data_matched$combined |>
            summarise(
              weights = mean(weights),
              observed_treatment_effect = sum((W_obs == 1) * Y_adj - (W_obs == 0) * Y_adj),
              .by = "subclass"
            )
        )
        calibration_data <- map2(
          expected_treatment_effect,
          observed_treatment_effect,
          \(exp, obs) {
            inner_join(exp, obs, by = c("subclass", "weights"))
          }
        )
        calibration_coef <-
          map(
            calibration_data,
            \(data) {
              lm_mod <- with(data, lm(observed_treatment_effect ~ expected_treatment_effect), weights = weights)
              coef(lm_mod) |>
                as_tibble() |>
                structure(names = "est") |>
                bind_cols(
                  confint(lm_mod) |>
                    as_tibble() |>
                    mutate(coef = c("\U03B2_0", "\U03B2_1"))
                ) |>
                select(4, 2, 1, 3)
            }
          )
        
        calibration_plot <- map2(
          calibration_data,
          calibration_coef,
          \(data, coef) {
            expr(
              !!data |>
                ggplot(aes(x = expected_treatment_effect, y = observed_treatment_effect)) +
                geom_hex(binwidth = 0.01) +
                geom_smooth(method = "lm") +
                geom_label(
                  aes(x = x, y = y, label = label),
                  data = tibble(
                    x = -0.2,
                    y = 0.7,
                    label = c(
                      paste0("\u03B2\u2080 = ", sprintf('%.3f', !!coef$est[1]), "\n\u03B2\u2081 = ", sprintf('%.3f', !!coef$est[2]))
                    )
                  )
                ) +
                geom_abline(slope = 1, intercept = 0, linetype = 2) +
                scale_x_continuous(
                  breaks = seq(-0.03, 0.12, 0.03),
                  labels = paste0(seq(-3, 12, 3), "%")
                ) +
                scale_y_continuous(
                  breaks = seq(-1, 1, 0.2),
                  labels = paste0(seq(-100, 100, 20), "%")
                ) +
                xlab("Predicted excess risk") + ylab("Observed excess risk") +
                coord_cartesian(xlim = c(-0.03, 0.12), ylim = c(-1, 1))
            )
          }
        )
        out$individual <- list(
          calibration_data = calibration_data,
          calibration_plot = calibration_plot,
          calibration_coef = calibration_coef
        )
        
        return(out)
      }
    )
    plan(sequential)
    
    # save dynamic subgroups results
    saveRDS(
      map(
        cfw_dsg_test, # map over model, e.g. full, vi_1, etc.
        \(mod) {
          out <- imap(
            mod, # map over n_rankings
            \(mod, nm) {
              if (nm == "individual") {
                mod$calibration_plot$combined <- deparse(mod$calibration_plot$combined)
              } else {
                mod$forest_rank_ate_plot$plot <- deparse(mod$forest_rank_ate_plot$plot)
                mod$calibration_plot$ols$CATE <- deparse(mod$calibration_plot$ols$CATE)
                mod$calibration_plot$ols$control_risk <- deparse(mod$calibration_plot$ols$control_risk)
                mod$calibration_plot$ols$combined <- deparse(mod$calibration_plot$ols$combined)
                mod$calibration_plot$aipw$CATE <- deparse(mod$calibration_plot$aipw$CATE)
                mod$calibration_plot$aipw$control_risk <- deparse(mod$calibration_plot$aipw$control_risk)
                mod$calibration_plot$aipw$combined <- deparse(mod$calibration_plot$aipw$combined)
              }
              return(mod)
            }
          )
          return(out)
        }
      ),
      paste0(result_path, "cfw_dsg_test.rds")
    )
    
    # save calibration plots using dynamic subgroups
    cal_plot_combined_aipw_comb_test <- imap(
      cfw_dsg_test,
      \(mod, nm) {
        if (!dir.exists(paste0(figure_path, nm, "/"))) dir.create(paste0(figure_path, nm, "/"))
        data_comb_ind <- mod$individual$calibration_data$combined
        coef_comb_ind <- mod$individual$calibration_coef$combined
        cal_plot_combined_aipw <- imap(
          mod[str_detect(names(mod), "ranking")],
          \(mod, rank, nm) {
            ate_plot <- eval(mod$forest_rank_ate_plot$plot)
            cal_plot_cate_ols <- eval(mod$calibration_plot$ols$CATE)
            cal_plot_control_risk_ols <- eval(mod$calibration_plot$ols$control_risk)
            cal_plot_combined_ols <- eval(mod$calibration_plot$ols$combined) +
              geom_label(
                aes(x = x, y = y, label = label),
                data = tibble(
                  x = -0.02,
                  y = 0.13,
                  label = c(
                    paste0(
                      "\u03B2\u2080 = ", 
                      sprintf('%.3f', coef_comb_ind$est[1]), 
                      "\n\u03B2\u2081 = ",
                      sprintf('%.3f', coef_comb_ind$est[2])
                    )
                  )
                )
              )
            cal_plot_cate_aipw <- eval(mod$calibration_plot$aipw$CATE)
            cal_plot_control_risk_aipw <- eval(mod$calibration_plot$aipw$control_risk)
            cal_plot_combined_aipw <- eval(mod$calibration_plot$aipw$combined) +
              geom_label(
                aes(x = x, y = y, label = label),
                data = tibble(
                  x = -0.02,
                  y = 0.13,
                  label = c(
                    paste0(
                      "\u03B2\u2080 = ", 
                      sprintf('%.3f', coef_comb_ind$est[1]), 
                      "\n\u03B2\u2081 = ", 
                      sprintf('%.3f', coef_comb_ind$est[2])
                    )
                  )
                )
              )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/cal_ate_plot_test_", rank, ".jpg"),
              plot = ate_plot,
              device = "jpeg",
              width = 16,
              height = 7,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/ols/cal_plot_cate_test_", rank, ".jpg"),
              plot = cal_plot_cate_ols,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/aipw/cal_plot_cate_test_", rank, ".jpg"),
              plot = cal_plot_cate_aipw,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/ols/cal_plot_control_risk_test_", rank, ".jpg"),
              plot = cal_plot_control_risk_ols,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/aipw/cal_plot_control_risk_test_", rank, ".jpg"),
              plot = cal_plot_control_risk_aipw,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/ols/cal_plot_combined_test_", rank, ".jpg"),
              plot = cal_plot_combined_ols,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/aipw/cal_plot_combined_test_", rank, ".jpg"),
              plot = cal_plot_combined_aipw,
              device = "jpeg",
              width = 12,
              height = 12,
              unit = "cm",
              scale = 1,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            return(cal_plot_combined_aipw)
          },
          nm = nm
        )
        return(cal_plot_combined_aipw)
      }
    )
    
    # save combined figure with calibration plots from all models
    dsg_plot_comb_5_test <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb_test, 
        \(x, nm) {
          if (nm == "full_mod") {
            label <- "66-cov (full) model"
          } else {
            label <- paste0(str_extract(nm, "[:digit:]{1,}"), "-cov model")
          }
          x <- x$ranking_5 +
            ggtitle(label = label) +
            theme(
              plot.title = element_text(size = 14, hjust = -0.2)
            )
        }
      )[c(2:7, 1)],
      nrow = 4,
      ncol = 2
    )
    dsg_plot_comb_10_test <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb_test, 
        \(x, nm) {
          if (nm == "full_mod") {
            label <- "66-cov (full) model"
          } else {
            label <- paste0(str_extract(nm, "[:digit:]{1,}"), "-cov model")
          }
          x <- x$ranking_10 +
            ggtitle(label = label) +
            theme(
              plot.title = element_text(size = 14, hjust = -0.2)
            )
        }
      )[c(2:7, 1)],
      nrow = 4,
      ncol = 2
    )
    dsg_plot_comb_20_test <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb_test, 
        \(x, nm) {
          if (nm == "full_mod") {
            label <- "66-cov (full) model"
          } else {
            label <- paste0(str_extract(nm, "[:digit:]{1,}"), "-cov model")
          }
          x <- x$ranking_20 +
            ggtitle(label = label) +
            theme(
              plot.title = element_text(size = 14, hjust = -0.2)
            )
        }
      )[c(2:7, 1)],
      nrow = 4,
      ncol = 2
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_5_test.jpg"),
      plot = dsg_plot_comb_5_test,
      device = "jpeg",
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      dpi = 240,
      quality = 90,
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_5_test.pdf"),
      plot = dsg_plot_comb_5_test,
      device = Cairo::CairoPDF,
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      family = "Arial Unicode MS",
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_10_test.jpg"),
      plot = dsg_plot_comb_10_test,
      device = "jpeg",
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      dpi = 240,
      quality = 90,
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_10_test.pdf"),
      plot = dsg_plot_comb_10_test,
      device = Cairo::CairoPDF,
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      family = "Arial Unicode MS",
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_20_test.jpg"),
      plot = dsg_plot_comb_20_test,
      device = "jpeg",
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      dpi = 240,
      quality = 90,
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_20_test.pdf"),
      plot = dsg_plot_comb_20_test,
      device = Cairo::CairoPDF,
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      family = "Arial Unicode MS",
      create.dir = TRUE
    )
    
    # save tables with coefficients from calibration plots
    cfw_dsg_test |> 
      (\(mod) {
        imap(
          mod, 
          \(mod, nm_mod) {
            imap(
              mod,
              \(mod, nm_grp, nm_mod) {
                if (nm_grp == "individual") {
                  list(
                    ols = list(
                      combined = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_coef$combined$`2.5 %`[1], 
                          mod$calibration_coef$combined$est[1],
                          mod$calibration_coef$combined$`97.5 %`[1]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_coef$combined$`2.5 %`[2], 
                          mod$calibration_coef$combined$est[2],
                          mod$calibration_coef$combined$`97.5 %`[2]
                        )
                      )
                    )
                  )
                } else {
                  list(
                    ols = list(
                      CATE = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$ols$CATE[1, 1], 
                          mod$calibration_coef$ols$CATE[1],
                          mod$calibration_confint$ols$CATE[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$ols$CATE[2, 1], 
                          mod$calibration_coef$ols$CATE[2],
                          mod$calibration_confint$ols$CATE[2, 2] 
                        )
                      ),
                      control_risk = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$ols$control_risk[1, 1], 
                          mod$calibration_coef$ols$control_risk[1],
                          mod$calibration_confint$ols$control_risk[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$ols$control_risk[2, 1], 
                          mod$calibration_coef$ols$control_risk[2],
                          mod$calibration_confint$ols$control_risk[2, 2] 
                        )
                      ),
                      combined = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$ols$combined[1, 1], 
                          mod$calibration_coef$ols$combined[1],
                          mod$calibration_confint$ols$combined[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$ols$combined[2, 1], 
                          mod$calibration_coef$ols$combined[2],
                          mod$calibration_confint$ols$combined[2, 2] 
                        )
                      )
                    ),
                    aipw = list(
                      CATE = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$aipw$CATE[1, 1], 
                          mod$calibration_coef$aipw$CATE[1],
                          mod$calibration_confint$aipw$CATE[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$aipw$CATE[2, 1], 
                          mod$calibration_coef$aipw$CATE[2],
                          mod$calibration_confint$aipw$CATE[2, 2] 
                        )
                      ),
                      control_risk = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$aipw$control_risk[1, 1], 
                          mod$calibration_coef$aipw$control_risk[1],
                          mod$calibration_confint$aipw$control_risk[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$aipw$control_risk[2, 1], 
                          mod$calibration_coef$aipw$control_risk[2],
                          mod$calibration_confint$aipw$control_risk[2, 2] 
                        )
                      ),
                      combined = tibble(
                        "model" = rep(nm_mod, 3),
                        "type" = c("lower", "estimate", "upper"),
                        "\U03B2_0" := c(
                          mod$calibration_confint$aipw$combined[1, 1], 
                          mod$calibration_coef$aipw$combined[1],
                          mod$calibration_confint$aipw$combined[1, 2]
                        ),
                        "\U03B2_1" := c(
                          mod$calibration_confint$aipw$combined[2, 1], 
                          mod$calibration_coef$aipw$combined[2],
                          mod$calibration_confint$aipw$combined[2, 2] 
                        )
                      )
                    )
                  )
                }
              },
              nm_mod = nm_mod
            )
          }
        )
      })() |> 
      list_transpose(simplify = FALSE) |>
      map(\(x) map(x, \(y) {list_transpose(y, simplify = FALSE) |> map(\(z) list_rbind(z, names_to = "exp_type"))})) |>
      map(\(x) list_transpose(x, simplify = FALSE) |> map(\(y) list_rbind(y, names_to = "model"))) |>
      list_transpose(simplify = FALSE) |>
      map(\(x) {
        list_rbind(x, names_to = "num_ranks") |> 
          mutate(
            num_ranks = replace_na(
              as.integer(str_extract(num_ranks, "\\d{1,}$")), 
              nrow(cfw_dsg_test$full_mod$individual$calibration_data$combined)
            )
          ) |>
          arrange(num_ranks, exp_type, model)
      }) |>
      dfs2xlsx(paste0(table_path, "cfw_dyn_cal_coef_test.xlsx"), rowNames = FALSE)
    
    return(NULL)
  }
)

### calibration plot bfz and hctz combined ------------------------------------------------------------------------ ----
cfw_dsg_bfz <- readRDS(paste0("results/data/mice/4mo_bfz_all/1year/cfw_dsg.rds"))
cfw_dsg_hctz <- readRDS(paste0("results/data/mice/4mo_hctz_all/1year/cfw_dsg.rds"))
cfw_dsg_test_bfz <- readRDS(paste0("results/data/mice/4mo_bfz_all/1year/cfw_dsg_test.rds"))
cfw_dsg_test_hctz <- readRDS(paste0("results/data/mice/4mo_hctz_all/1year/cfw_dsg_test.rds"))

cfw_bfz_hctz_dsg_plots_10 <- pmap(
  list(
    list(
      bfz_dev = cfw_dsg_bfz$vi_mod_4$ranking_10$calibration_data$aipw$combined,
      hctz_dev = cfw_dsg_hctz$vi_mod_4$ranking_10$calibration_data$aipw$combined,
      bfz_val = cfw_dsg_test_bfz$vi_mod_4$ranking_10$calibration_data$aipw$combined,
      hctz_val = cfw_dsg_test_hctz$vi_mod_4$ranking_10$calibration_data$aipw$combined
    ),
    list(
      bfz_dev = cfw_dsg_bfz$vi_mod_4$ranking_10$calibration_coef$aipw$combined,
      hctz_dev = cfw_dsg_hctz$vi_mod_4$ranking_10$calibration_coef$aipw$combined,
      bfz_val = cfw_dsg_test_bfz$vi_mod_4$ranking_10$calibration_coef$aipw$combined,
      hctz_val = cfw_dsg_test_hctz$vi_mod_4$ranking_10$calibration_coef$aipw$combined
    ),
    nm = c("bfz_dev", "hctz_dev", "bfz_val", "hctz_val")
  ),
  \(data, coef, nm) {
    data |>
      ggplot(aes(x = expected_treatment_effect, y = observed_treatment_effect)) +
      geom_point() +
      geom_errorbar(
        aes(
          ymin = observed_treatment_effect - observed_treatment_effect_se,
          ymax = observed_treatment_effect + observed_treatment_effect_se,
          width = 0.001
        )
      ) +
      geom_errorbarh(
        aes(
          xmin = expected_treatment_effect - expected_treatment_effect_se,
          xmax = expected_treatment_effect + expected_treatment_effect_se,
          height = 0.001
        )
      ) +
      geom_label(
        aes(x = x, y = y, label = label),
        data = tibble(
          x = -0.02,
          y = 0.13,
          label = c(
            paste0("\u03B2\u2080 = ", sprintf('%.3f', !!coef[1]), "\n\u03B2\u2081 = ", sprintf('%.3f', !!coef[2]))
          )
        )
      ) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      scale_x_continuous(
        breaks = seq(-0.03, 0.15, 0.03),
        labels = paste0(seq(-3, 15, 3), "%")
      ) +
      scale_y_continuous(
        breaks = seq(-0.03, 0.15, 0.03),
        labels = paste0(seq(-3, 15, 3), "%")
      ) +
      xlab(
        if (str_detect(nm, "dev")) {
          ""
        } else {
          "Predicted excess risk"
        }
      ) + 
      ylab(
        if (str_detect(nm, "bfz")) {
          "Observed excess risk"
        } else {
          ""
        }
      ) +
      coord_cartesian(xlim = c(-0.05, 0.15), ylim = c(-0.05, 0.15)) +
      ggplot2::ggtitle(
        if (str_detect(nm, "dev")) {
          "Development cohort"
        } else {
          "Validation cohort"
        }
      ) +
      ggplot2::theme(
        plot.title = element_text(size = 10, hjust = -0.3)
      ) 
  }
)
dsg_plot_10_bfz_hctz_titles <- list(
  ggdraw() +
    draw_label(
      "Bendroflumethiazide model",
      fontface = "bold",
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    ),
  ggdraw() +
    draw_label(
      "Hydrochlorothiazide model",
      fontface = "bold",
      x = 0.5,
      hjust = 0.5,
      size = 12
    ) +
    theme(
      plot.margin = margin(0, 0, 0, 7)
    )
)
dsg_plot_10_bfz_hctz <- cowplot::plot_grid(
  plotlist = c(dsg_plot_10_bfz_hctz_titles, cfw_bfz_hctz_dsg_plots_10),
  nrow = 3,
  ncol = 2,
  rel_heights = c(0.1, 1, 1)
)
ggplot2::ggsave(
  filename = paste0("results/figures/mice/4mo_bfz_hctz_comp/1year/cal_plot_aipw_bfz_hctz_10.jpg"),
  plot = dsg_plot_10_bfz_hctz,
  device = "jpeg",
  width = 15,
  height = 15,
  unit = "cm",
  scale = 1.25,
  dpi = 240,
  quality = 90,
  create.dir = TRUE
)
ggplot2::ggsave(
  filename = paste0("results/figures/mice/4mo_bfz_hctz_comp/1year/cal_plot_aipw_bfz_hctz_10.pdf"),
  plot = dsg_plot_10_bfz_hctz,
  device = Cairo::CairoPDF,
  width = 15,
  height = 15,
  unit = "cm",
  scale = 1.25,
  family = "Arial Unicode MS",
  create.dir = TRUE
)

### Distribution of CATEs ----------------------------------------------------------------------------------------- ----
plan(sequential)
future_pmap(
  .l = list(
    cohort = list(
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
      list(tih_cohort_test), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      list(tih_cohort_test_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
      tih_cohort_test_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
    ),
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_tiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_tiazide_no_lab/"),
      paste0(figures_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_tiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_tiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_tiazide_cc_comp/1year/")
    ),
    continuous = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89])$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      )
    ),
    discrete = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      )
    ),
    seed = as.list(seeds[59:69])
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
      "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
      "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
      "vi_models", "excess_risk_cutoff"
    ),
    seed = TRUE
  ),
  .f = \(cohort, path, result_path, figure_path, continuous, discrete, seed) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    cfw_cate_test <- future_pmap(
      .l = list(
        i = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
        cohort = cohort
      ),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c(
          "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
          "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
          "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
          "excess_risk_cutoff", "path", "result_path", "vi_models", "continuous", "discrete"
        ),
        seed = seed
      ),
      .f = \(i, cohort)  {
        options("parallelly.availableCores.methods" = "system")
        cfw_cate_test <- map(
          .x = seq_along(list.files(path, pattern = "^cfw_full_mod_\\d{1,}")),
          .f = \(j) {
            workers <- min(floor(num_threads / num_clusters), length(vi_models) + 1)
            threads_per_worker <- floor(num_threads / num_clusters / workers)
            if(workers == 1) {
              plan(sequential)
            } else {
              plan(multisession, workers = workers)
            }
            out <- furrr::future_map(
              .x = c(0, vi_models),
              .options = furrr_options(
                packages = c(
                  "dplyr", "stringr", "tibble", "grf", "glue", "stats", "future", "furrr"
                ),
                globals = c("i", "j", "path", "num_threads", "num_clusters", "threads_per_worker",
                            "continuous", "discrete", "cohort", "DiscreteCovariatesToOneHot"),
                seed = TRUE
              ),
              .f = \(vi_model) {
                if (vi_model == 0) {
                  cf <- readRDS(paste0(path, "cfw_full_mod_", j, ".rds"))
                } else {
                  cf <- readRDS(paste0(path, "cfw_vi_mod_", vi_model, "_", j, ".rds"))
                }
                continuous_cf <- continuous$cf
                discrete_cf <- discrete$cf
                predictions <- predict(
                  cf, 
                  newdata = cohort |>
                    select({{ continuous_cf }}) |>
                    bind_cols(
                      cohort |>
                        select({{ discrete_cf }}) |>
                        DiscreteCovariatesToOneHot()
                    ) |>
                    select(all_of(names(cf$X.orig))),
                  estimate.variance = TRUE, 
                  num.threads = threads_per_worker
                )
                return(predictions)
              }
            )
            plan(sequential)
            names(out) <- c(
              "full_mod",
              paste0("vi_mod_", vi_models)
            )
            return(out)
          }
        )
        
        # pool CATE estimates from each development model
        plan(multisession, workers = min(num_clusters, length(vi_models) + 1))
        cfw_cate_test_pooled <- cfw_cate_test |>
          list_transpose(simplify = FALSE) |>
          future_map(
            .options = furrr_options(
              packages = c("mice", "dplyr", "purrr")
            ),
            .f = \(model) {
              model |>
                list_transpose(simplify = FALSE) |>
                map(\(x) list_transpose(x, simplify = TRUE)) |>
                list_transpose(simplify = FALSE) |>
                (\(x) {
                  n <- length(x)
                  map(
                    x,
                    \(x, n) {
                      if (is.list(x) && length(x[[1]]) > 1) {
                        pool <- pool.scalar(Q = x$predictions, U = x$variance.estimates, n = n, k = 1, rule = "rubin1987")
                        out <- tibble(
                          predictions = pool$qbar,
                          variance_estimates = pool$t
                        )
                      } else {
                        pool <- NULL
                        out <- tibble(
                          predictions = x$predictions,
                          variance_estimates = x$variance.estimates
                        )
                      }
                      return(list(table = out, pool_results = pool))
                    }, 
                    n = n
                  )
                })() |>
                list_transpose()
            }
          )
        plan(sequential)
        
        saveRDS(cfw_cate_test_pooled, paste0(result_path, "cfw_cate_test_with_censored_units_", i, ".rds"))
        
        return(map(cfw_cate_test_pooled, \(x) x$table))
      }
    )
    
    # aggregate cate estimates
    cfw_cate_test_agg <- cfw_cate_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          model |>
            purrr::list_transpose(simplify = FALSE) |>
            map(\(x) purrr::list_transpose(x, simplify = TRUE)) |>
            purrr::list_transpose(simplify = FALSE) |>
            (\(x) {
              n <- length(x)
              map(
                x,
                \(x, n) {
                  if (is.list(x) && length(x[[1]]) > 1) {
                    pool <- pool.scalar(Q = x$predictions, U = x$variance_estimates, n = n, k = 1, rule = "rubin1987")
                    out <- tibble(
                      estimate = pool$qbar,
                      std_err = sqrt(pool$t),
                      lower = estimate - qnorm(0.975) * std_err,
                      upper = estimate + qnorm(0.975) * std_err
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      estimate = x[["predictions"]],
                      std_err = sqrt(x[["variance_estimates"]]),
                      lower = estimate - qnorm(0.975) * std_err,
                      upper = estimate + qnorm(0.975) * std_err
                    )
                  }
                  return(list(table = out, pool_results = pool))
                }, 
                n = 1
              )
            })() |>
            list_transpose()
        }
      )
    plan(sequential)
    
    saveRDS(cfw_cate_test_agg, paste0(result_path, "cfw_cate_test_with_censored_units_agg.rds"))
    
    # plot density and rank 
    imap(
      cfw_cate_test_agg,
      \(x, nm) {
        cate_table <- x$table
        cate_table <- cate_table |> 
          arrange(estimate) |>
          mutate(
            order = seq_len(n()),
            quantiles = factor(case_when(
              order / n() <= 0.5 ~ "0%-50%",
              order / n() <= 0.9 ~ "50%-90%",
              order / n() <= 0.95 ~ "90%-95%",
              TRUE ~ "95%-100%",
            ))
          )
        plot_a <- ggplot(cate_table) +
          geom_histogram(aes(x = estimate, fill = quantiles), binwidth = 0.002) + 
          scale_x_continuous(breaks = seq(
            floor(min(cate_table$estimate) * 40) / 40, 
            ceiling(max(cate_table$estimate) * 40) / 40, 
            0.025
          ),
          labels = paste0(
            seq(floor(min(cate_table$estimate)*40)/40, ceiling(max(cate_table$estimate)*40)/40, 0.025)*100, "%"
          )
          ) +
          scale_y_continuous(labels = \(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE)) +
          scale_fill_jama() +
          coord_cartesian(xlim = c(floor(min(cate_table$estimate)*100)/100, ceiling(max(cate_table$estimate)*100)/100)) +
          xlab("Excess risk") + 
          ylab("Number of patients") +
          theme(legend.title = element_blank())
        plot_b <- ggplot(cate_table) +
          geom_line(aes(x = order, y = estimate, color = quantiles), linewidth = 1.1) + 
          scale_y_continuous(breaks = seq(
            floor(min(cate_table$estimate) * 40) / 40, 
            ceiling(max(cate_table$estimate) * 40) / 40, 
            0.025
          ),
          labels = paste0(
            seq(floor(min(cate_table$estimate)*40)/40, ceiling(max(cate_table$estimate)*40)/40, 0.025)*100, "%"
          )) +
          scale_x_continuous(
            breaks = round(nrow(cate_table) * c(0, 0.25, 0.5, 0.75, 1)),
            labels = \(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE)
          ) +
          scale_color_jama() +
          coord_cartesian(ylim = c(floor(min(cate_table$estimate)*100)/100, ceiling(max(cate_table$estimate)*100)/100)) +
          xlab("Patients ranked by excess risk") + 
          ylab("Excess risk") +
          theme(legend.title = element_blank())
        plot_comb <- cowplot::plot_grid(
          plot_a,
          plot_b,
          labels = c("A", "B"),
          label_size = 10,
          nrow = 2,
          ncol = 1
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_density_plot_test", ".jpg"),
          plot = plot_a,
          device = "jpeg",
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_density_plot_test", ".pdf"),
          plot = plot_a,
          device = Cairo::CairoPDF,
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          family = "Arial Unicode MS",
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_rank_plot_test", ".jpg"),
          plot = plot_b,
          device = "jpeg",
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_rank_plot_test", ".pdf"),
          plot = plot_b,
          device = Cairo::CairoPDF,
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          family = "Arial Unicode MS",
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_plot_comb_test", ".jpg"),
          plot = plot_comb,
          device = "jpeg",
          width = 10,
          height = 14,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_plot_comb_test", ".pdf"),
          plot = plot_comb,
          device = Cairo::CairoPDF,
          width = 10,
          height = 14,
          unit = "cm",
          scale = 1.3,
          family = "Arial Unicode MS",
          create.dir = TRUE
        )
        
        return(invisible(NULL))
      }
    )
  }
)

### excess risk strategies ---------------------------------------------------------------------------------------- ----
plan(sequential)
future_pmap(
  .l = list(
    cohort = list(
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
      list(tih_cohort_test), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      list(tih_cohort_test_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
      tih_cohort_test_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
    ),
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    W_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W)
    ),
    Y_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(Y_125_90),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120)
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "tih_cohort_test", "n_rankings", "ATEAll", "num_threads", "num_clusters", "names_index", 
      "compute_results", "names_index_comb", "mice", "excess_risk_quantiles", "vi_models", "excess_risk_cutoff"
    ),
    seed = TRUE
  ),
  .f = \(cohort, path, result_path, table_path, W_orig, Y_orig) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    
    # Look at the expected excess risk among those with a CATE above certain
    # quantile cutoffs (0%, 50%, 90%, 95%). Also look at population excess risk when applying 
    # treatment strategy where treatment is avoided based on these cut-offs.
    # we consider a strategy of avoiding thiazides above a particular excess risk,
    # i.e. treatment pattern among low excess risk is the same, but everyone 
    # above cutoff get non-thiazide.
    future_pmap(
      .l = list(
        i = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
        cohort = cohort
      ),
      .options = furrr_options(
        packages = c("stringr", "purrr", "dplyr", "grf", "future", "furrr"),
        globals = c("path", "result_path", "excess_risk_cutoff", "num_threads", "num_clusters",
                    "ATEAll", "Y_orig", "W_orig", "mice", "vi_models", "excess_risk_quantiles"),
        seed = TRUE
      ),
      .f = \(i, cohort) {
        options("parallelly.availableCores.methods" = "system")
        cfw_full_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
        sample_weights <- cfw_full_test$sample.weights
        Y_hat <- cfw_full_test$Y.hat
        W_hat <- cfw_full_test$W.hat
        cfw_cate_test <- readRDS(paste0(result_path, "cfw_cate_test_", i, ".rds"))
        if (is.null(names(cfw_cate_test))) {
          names(cfw_cate_test) <- c("full_mod", paste0("vi_mod_", vi_models))
        }
        tau_hat_pointwise <- map(cfw_cate_test, \(x) x$table$predictions)
        cfw_cate_all_test <- readRDS(paste0(result_path, "cfw_cate_test_with_censored_units_", i, ".rds"))
        if (is.null(names(cfw_cate_all_test))) {
          names(cfw_cate_all_test) <- c("full_mod", paste0("vi_mod_", vi_models))
        }
        tau_hat_all <- map(cfw_cate_all_test, \(x) x$table$predictions)
        
        plan(sequential)
        
        cfw_risk_strategy_test <- future_map2(
          .x = tau_hat_pointwise,
          .y = tau_hat_all,
          .options = furrr_options(
            packages = c("stringr", "purrr", "dplyr", "grf"),
            globals = c("excess_risk_cutoff", "ATEAll", "Y_orig", "W_orig", "mice", "excess_risk_quantiles",
                        "sample_weights", "Y_hat", "W_hat", "i", "vi_models"),
            seed = TRUE
          ),
          .f = \(tau_hat_pointwise, tau_hat_all, W) {
            CFRiskStrategy(
              list(
                Y.orig = Y_orig, 
                Y.hat = Y_hat,
                W.orig = W_orig,
                W.hat = W_hat,
                predictions = tau_hat_pointwise,
                sample.weights = sample_weights
              ), 
              tau_hat_all,
              W,
              type = "test",
              level = 0.95,
              excess_risk_cutoff = excess_risk_cutoff,
              excess_risk_quantiles = excess_risk_quantiles
            )
          },
          W = cohort$W
        )
        saveRDS(cfw_risk_strategy_test, paste0(result_path, "cfw_risk_strategy_test_", i, ".rds"))
        rm(cfw_risk_strategy_test)
      }
    )
    
    # aggregate results from excess risk strategies 
    cfw_risk_strategy_test <- list()
    for (i in seq_along(list.files(result_path, pattern = "^cfw_risk_strategy_test_\\d{1,}"))) {
      cfw_risk_strategy_test[[i]] <- readRDS(paste0(result_path, "cfw_risk_strategy_test_", i, ".rds"))
    }
    cfw_risk_strategy_test_agg <- cfw_risk_strategy_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose(simplify = FALSE) |>
              (\(y) {
                map2(
                  y,
                  names(y),
                  \(list, nm) {
                    list |>
                      (\(list) {
                        out <- list()
                        for (i in seq_along(list)) {
                          for (j in seq_along(list[[i]])) {
                            if (i == 1) {
                              if (nm == "name") out[[j]] <- vector("character", length(list))
                              else out[[j]] <- vector("double", length(list))
                            }
                            out[[j]][i] <- list[[i]][j]
                          }
                        }
                        return(out)
                      })() 
                  }
                )
              })() |>
              list_transpose(simplify = FALSE) |>
              map(
                \(x) {
                  if (any(str_detect(x$name, "^units_nt"))) {
                    pool <- list(
                      "lower_0.025" = pool.scalar(Q = 1 / x[["lower_0.025"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                      "estimate" = pool.scalar(Q = 1 / x$estimate, U = 0, n = 1, k = 1, rule = "rubin1987"),
                      "upper_0.975" = pool.scalar(Q = 1 / x[["upper_0.975"]], U = 0, n = 1, k = 1, rule = "rubin1987")
                    )
                    out <- tibble(
                      "name" = x$name[1],
                      "lower_0.025" = 1 / pool[["lower_0.025"]]$qbar,
                      "estimate" = 1 / pool$estimate$qbar,
                      "upper_0.975" = 1 / pool[["upper_0.975"]]$qbar,
                      "std_err" = NA
                    )
                    return(list(table = out, pool = pool))
                  } else if (any(is.na(x$std_err))) {
                    pool <- list(
                      "lower_0.025" = pool.scalar(Q = x[["lower_0.025"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                      "estimate" = pool.scalar(Q = x$estimate, U = 0, n = 1, k = 1, rule = "rubin1987"),
                      "upper_0.975" = pool.scalar(Q = x[["upper_0.975"]], U = 0, n = 1, k = 1, rule = "rubin1987")
                    )
                    out <- tibble(
                      "name" = x$name[1],
                      "lower_0.025" = pool[["lower_0.025"]]$qbar,
                      "estimate" = pool$estimate$qbar,
                      "upper_0.975" = pool[["upper_0.975"]]$qbar,
                      "std_err" = NA
                    )
                    return(list(table = out, pool = pool))
                  } else {
                    pool <- pool.scalar(Q = x$estimate, U = x$std_err, n = 1, k = 1, rule = "rubin1987")
                    out <- tibble(
                      "name" = x$name[1],
                      "lower_0.025" = pool$qbar + qnorm(0.5 - 0.95 / 2) * pool$t,
                      "estimate" = pool$qbar,
                      "upper_0.975" = pool$qbar + qnorm(0.5 + 0.95 / 2) * pool$t,
                      "std_err" = pool$t
                    )
                    # for cases avoided, multiply with total number of treated above risk threshold
                    if (identical(x[["lower_0.025"]], x[["upper_0.975"]])) {
                      out <- mutate(out, across(2:5, \(y) y * mean(x[["lower_0.025"]])))
                    }
                    return(list(table = out, pool = pool))
                  }
                }
              ) |>
              (\(x) {
                out <- list_transpose(x)
                names(out$pool) <- out$table$name
                return(out)
              })()
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                (\(x) {
                  tibble(
                    "name" = x$name,
                    "lower_0.025" = x$`lower_0.025`, 
                    "estimate" = x$estimate, 
                    "upper_0.975" = x$`upper_0.975`,
                    "std_err" = x$std_err
                  )
                })()
            )
            for (i in seq_len(nrow(out$table))) {
              tbl <- out$table[i,]
              if (!is.na(tbl[["lower_0.025"]]) && !is.na(tbl[["upper_0.975"]]) && identical(tbl[["lower_0.025"]], tbl[["upper_0.975"]])) {
                count <- tbl[["lower_0.025"]]
                out$table[i,] <- mutate(
                  tbl, 
                  "lower_0.025" = estimate + qnorm(0.5 - 0.95 / 2) * std_err, 
                  "upper_0.975" = estimate + qnorm(0.5 + 0.95 / 2) * std_err, 
                  across(2:5, \(y) y * count)
                )
              }
            }
          }
          return(out)
        }
      )
    saveRDS(cfw_risk_strategy_test_agg, paste0(result_path, "cfw_risk_strategy_test_agg.rds"))
    # Risk strategy table
    cfw_risk_strategy_test_agg |> 
      map(\(x) {
        x$table
      }) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_strategy_test.xlsx"
        ), 
        rowNames = FALSE
      )
    plan(sequential)
  }
)

### excess risk reduction plots ----------------------------------------------------------------------------------- ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    W_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(W_bvc),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(W_hvr),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(W),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(W)
    ),
    Y_orig = list(
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_125_90)) |> pull(Y_125_90),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_bvc)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120) & !is.na(W_hvr)) |> pull(Y_130_120),
      filter(tih_cohort_test_complete_case, !is.na(Y_130_120)) |> pull(Y_130_120),
      filter(tih_cohort_test, !is.na(Y_130_120)) |> pull(Y_130_120)
    ),
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "tih_cohort_test", "n_rankings", "ATEAll", "num_threads", "num_clusters", "names_index", "compute_results",
      "names_index_comb", "mice", "excess_risk_quantiles", "vi_models", "excess_risk_cutoff"
    ),
    seed = TRUE
  ),
  .f = \(path, result_path, table_path, W_orig, Y_orig) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    
    # excess risk reduction
    future_map(
      .x = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
      .options = furrr_options(
        packages = c("stringr", "purrr", "dplyr", "grf", "future", "furrr"),
        globals = c("path", "result_path", "excess_risk_cutoff", "num_threads", "num_clusters",
                    "ATEAll", "Y_orig", "W_orig", "mice", "vi_models", "excess_risk_quantiles"),
        seed = TRUE
      ),
      .f = \(i) {
        options("parallelly.availableCores.methods" = "system")
        cfw_full_test <- readRDS(paste0(path, "cfw_full_test_mod_", i, ".rds"))
        sample_weights <- cfw_full_test$sample.weights
        Y_hat <- cfw_full_test$Y.hat
        W_hat <- cfw_full_test$W.hat
        cfw_cate_test <- readRDS(paste0(result_path, "cfw_cate_test_", i, ".rds"))
        if (is.null(names(cfw_cate_test))) {
          names(cfw_cate_test) <- c("full_mod", paste0("vi_mod_", vi_models))
        }
        tau_hat_pointwise <- map(cfw_cate_test, \(x) x$table$predictions)
        
        plan(sequential)
        cfw_risk_reduction_plot_test <- future_map(
          .x = tau_hat_pointwise,
          .options = furrr_options(
            packages = c("stringr", "purrr", "dplyr", "grf"),
            globals = c("workers", "excess_risk_cutoff", "excess_risk_quantiles",
                        "W_orig", "Y_orig", "W_hat", "Y_hat", "sample_weights"),
            seed = TRUE
          ),
          .f = \(predictions) {
            aipw <- map(
              seq(0, 1, 0.02),
              \(percentile) {
                # determine the cutoff in excess risk among treated units to treat with new strategy
                cutoff <- tibble(
                  predictions = predictions,
                  W = W_orig,
                  sample_weight = sample_weights
                ) |>
                  filter(W == 1) |>
                  arrange(predictions) |>
                  mutate(cum_weight = cumsum(sample_weight)) |>
                  (\(x) slice_head(x, n = sum(x$cum_weight <= percentile * sum(x$sample_weight))))() |>
                  pull(predictions) |>
                  max()
                
                # extract rows of data with risk above cutoff of treated with new strategy
                data_nt <- tibble(
                  predictions = predictions,
                  W = W_orig,
                  Y = Y_orig,
                  W_hat = W_hat,
                  Y_hat = Y_hat,
                  sample_weight = sample_weights,
                  id = seq_along(predictions)
                ) |>
                  filter(predictions > cutoff)
                # ATE among treated in low percentile of risk
                ate_treated_nt <- with(
                  data_nt,
                  tryCatch(
                    ATEAll(
                      Y.orig = Y, 
                      Y.hat = Y_hat,
                      W.orig = W,
                      W.hat = W_hat,
                      tau.hat.pointwise = predictions,
                      sample.weights = sample_weight,
                      target.sample = "treated"
                    ) |>
                      (\(x) {
                        tibble(
                          estimate = x["estimate"],
                          std_err = x["std.err"]
                        )
                      })(),
                    error = \(e) {
                      tibble(
                        estimate = NA_real_,
                        std_err = NA_real_
                      )
                    }
                  ) 
                )
                # Difference in risk between full population and part with low percentile of risk
                risk_avoid <- tryCatch(
                  tibble(
                    percentile = percentile,
                    lower_0.025 = (1 - percentile) * ate_treated_nt[["estimate"]] + 
                      qnorm(0.025) * (1 - percentile) * ate_treated_nt[["std_err"]],
                    estimate = (1 - percentile) * ate_treated_nt[["estimate"]],
                    upper_0.975 = (1 - percentile) * ate_treated_nt[["estimate"]] + 
                      qnorm(0.975) * (1 - percentile) * ate_treated_nt[["std_err"]],
                    std_err = (1 - percentile) * ate_treated_nt[["std_err"]],
                    estimate_nt = ate_treated_nt[["estimate"]]
                  ),
                  error = \(e) {
                    tibble(
                      percentile = percentile,
                      lower_0.025 = NA_real_,
                      estimate = (1 - percentile) * ate_treated_nt[["estimate"]],
                      upper_0.975 = NA_real_,
                      std_err = (1 - percentile) * ate_treated_nt[["std_err"]],
                      estimate_nt = ate_treated_nt[["estimate"]]
                    )
                  }
                )
              }
            ) |>
              list_rbind()
            aipw <- add_row(
              aipw, 
              percentile = 1, 
              lower_0.025 = 0,
              estimate = 0,
              upper_0.975 = 0,
              std_err = 0,
              .before = 51L
            ) |>
              slice(1:51)
            
            ols <- map(
              seq(0, 1, 0.02),
              \(percentile) {
                pred <- tibble(
                  predictions = predictions,
                  W = W_orig,
                  sample_weight = sample_weights
                ) |>
                  filter(W == 1) |>
                  arrange(predictions) |>
                  mutate(cum_weight = cumsum(sample_weight))
                n_keep <- sum(pred$cum_weight <= percentile * sum(pred$sample_weight))
                ate_treated_cutoff <- tryCatch(
                  tibble(
                    estimate = weighted.mean(head(pred$predictions, n_keep), head(pred$sample_weight, n_keep))
                  ),
                  error = \(e) {
                    tibble(
                      estimate = NA_real_
                    )
                  }
                ) |>
                  mutate(estimate = estimate * sum(head(pred$sample_weight, n_keep)) / sum(pred$sample_weight))
                n_keep <- sum(pred$cum_weight > percentile * sum(pred$sample_weight))
                ate_treated_nt <- tryCatch(
                  tibble(
                    estimate = weighted.mean(tail(pred$predictions, n_keep), tail(pred$sample_weight, n_keep))
                  ),
                  error = \(e) {
                    tibble(
                      estimate = NA_real_
                    )
                  }
                ) 
                ate_treated <- tibble(estimate = weighted.mean(pred$predictions, pred$sample_weight))
                risk_avoid <- tryCatch(
                  tibble(
                    percentile = percentile,
                    estimate = ate_treated[["estimate"]] - ate_treated_cutoff[["estimate"]],
                    estimate_nt = ate_treated_nt[["estimate"]]
                  ),
                  error = \(e) {
                    tibble(
                      percentile = percentile,
                      estimate = ate_treated[["estimate"]],
                      estimate_nt = ate_treated_nt[["estimate"]]
                    )
                  }
                )
                return(risk_avoid)
              }
            ) |>
              list_rbind()
            ols <- add_row(
              ols, 
              percentile = 0,
              estimate = with(
                tibble(
                  predictions = predictions,
                  W = W_orig,
                  sample_weight = sample_weights
                ) |>
                  filter(W == 1),
                weighted.mean(predictions, sample_weight)
              ),
              estimate_nt = with(
                tibble(
                  predictions = predictions,
                  W = W_orig,
                  sample_weight = sample_weights
                ) |>
                  filter(W == 1),
                weighted.mean(predictions, sample_weight)
              ),
              .before = 2L
            ) |>
              slice(2:n())
            
            # ATE among treated in full population
            ate_all <- ATEAll(
              Y.orig = Y_orig, 
              Y.hat = Y_hat,
              W.orig = W_orig,
              W.hat = W_hat,
              tau.hat.pointwise = predictions,
              sample.weights = sample_weights
            ) |>
              (\(x) {
                tibble(
                  estimate = x["estimate"],
                  std_err = x["std.err"]
                )
              })()
            
            out <- list(aipw = aipw, ols = ols, ate_all = ate_all)
            return(out)
          }
        )
        saveRDS(
          map(cfw_risk_reduction_plot_test, \(x) x$ate_all), 
          paste0(result_path, "cfw_ate_test_", i, ".rds")
        )
        saveRDS(
          map(cfw_risk_reduction_plot_test, \(x) x[c("aipw", "ols")]), 
          paste0(result_path, "cfw_risk_reduction_plot_test_", i, ".rds")
        )
        cfw_risk_reduction_plot_test_zp <- future_map(
          .x = tau_hat_pointwise,
          .options = furrr_options(
            packages = c("stringr", "purrr", "dplyr", "grf"),
            globals = c("workers", "excess_risk_cutoff", "excess_risk_quantiles"),
            seed = TRUE
          ),
          .f = \(predictions) {
            mean(predictions <= 0)
          }
        )
        saveRDS(
          cfw_risk_reduction_plot_test_zp, 
          paste0(result_path, "cfw_risk_reduction_plot_test_zp_", i, ".rds")
        )
      }
    )
    
    # pool results
    cfw_risk_reduction_plot_test <- list()
    cfw_risk_reduction_plot_test_zp <- list()
    cfw_ate_test <- list()
    for (i in seq_along(list.files(result_path, pattern = "^cfw_risk_reduction_plot_test_\\d{1,}"))) {
      cfw_risk_reduction_plot_test[[i]] <- readRDS(paste0(result_path, "cfw_risk_reduction_plot_test_", i, ".rds"))
      cfw_risk_reduction_plot_test_zp[[i]] <- readRDS(paste0(result_path, "cfw_risk_reduction_plot_test_zp_", i, ".rds"))
      cfw_ate_test[[i]] <- readRDS(paste0(result_path, "cfw_ate_test_", i, ".rds"))
    }
    cfw_risk_reduction_plot_test_agg <- cfw_risk_reduction_plot_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose() |>
              map(
                \(model) {
                  model |>
                    list_transpose(simplify = FALSE) |>
                    (\(y) {
                      map2(
                        y,
                        names(y),
                        \(list, nm) {
                          list |>
                            (\(list) {
                              out <- list()
                              for (i in seq_along(list)) {
                                for (j in seq_along(list[[i]])) {
                                  if (i == 1) {
                                    out[[j]] <- vector("double", length(list))
                                  }
                                  out[[j]][i] <- list[[i]][j]
                                }
                              }
                              return(out)
                            })() 
                        }
                      )
                    })() |>
                    list_transpose(simplify = FALSE) |>
                    map(
                      \(x) {
                        if (any(is.na(x$std_err))) {
                          pool <- list(
                            "lower_0.025" = pool.scalar(Q = x[["lower_0.025"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "estimate" = pool.scalar(Q = x$estimate, U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "upper_0.975" = pool.scalar(Q = x[["upper_0.975"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "estimate_nt" = pool.scalar(Q = x$estimate_nt, U = 0, n = 1, k = 1, rule = "rubin1987")
                          )
                          out <- tibble(
                            "percentile" = x$percentile[1],
                            "lower_0.025" = pool[["lower_0.025"]]$qbar,
                            "estimate" = pool$estimate$qbar,
                            "upper_0.975" = pool[["upper_0.975"]]$qbar,
                            "std_err" = NA,
                            "estimate_nt" = pool$estimate_nt$qbar,
                            "unt_estimate" = 1 / estimate_nt
                          )
                          return(list(table = out, pool = pool))
                        } else {
                          pool <-  list(
                            "estimate" = pool.scalar(Q = x$estimate, U = (x$std_err)^2, n = 1, k = 1, rule = "rubin1987"),
                            "estimate_nt" = pool.scalar(Q = x$estimate_nt, U = 0, n = 1, k = 1, rule = "rubin1987")
                          )
                          out <- tibble(
                            "percentile" = x$percentile[1],
                            "lower_0.025" = pool$estimate$qbar + qnorm(0.5 - 0.95 / 2) * sqrt(pool$estimate$t),
                            "estimate" = pool$estimate$qbar,
                            "upper_0.975" = pool$estimate$qbar + qnorm(0.5 + 0.95 / 2) * sqrt(pool$estimate$t),
                            "std_err" = sqrt(pool$estimate$t),
                            "estimate_nt" = pool$estimate_nt$qbar,
                            "unt_estimate" = 1 / estimate_nt
                          )
                          return(list(table = out, pool = pool))
                        }
                      }
                    ) |>
                    (\(x) {
                      out <- list_transpose(x)
                      names(out$pool) <- out$table$percentile
                      return(out)
                    })()
                }
              )
          } else {
            out <- map(
              model[[1]],
              \(model) {
                list(
                  pool = NULL,
                  table = model |>
                    (\(x) {
                      tibble(
                        "percentile" = x$percentile,
                        "lower_0.025" = x$`lower_0.025`, 
                        "estimate" = x$estimate, 
                        "upper_0.975" = x$`upper_0.975`,
                        "std_err" = x$std_err,
                        "estimate_nt" = x$estimate_nt,
                        "unt_estimate" = 1 / estimate_nt
                      )
                    })()
                )
              }
            )
          }
          return(out)
        }
      )
    saveRDS(
      cfw_risk_reduction_plot_test_agg, 
      paste0(result_path, "cfw_risk_reduction_plot_test_agg.rds")
    )
    cfw_risk_reduction_plot_test_zp_agg <- cfw_risk_reduction_plot_test_zp |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            pool <- pool.scalar(Q = unlist(model), U = 0, n = 1, k = 1, rule = "rubin1987")
            return(list(quantile = pool$qbar, pool = pool))
          } else {
            return(list(quantile = unlist(model), pool = NULL))
          }
        }
      )
    saveRDS(
      cfw_risk_reduction_plot_test_zp_agg, 
      paste0(result_path, "cfw_risk_reduction_plot_test_zp_agg.rds")
    )
    
    cfw_ate_test_agg <- cfw_ate_test |>
      list_transpose(simplify = FALSE) |>
      map(\(x) list_transpose(x, simplify = TRUE)) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            pool <- pool.scalar(Q = model$estimate, U = (model$std_err)^2, n = 1, k = 1, rule = "rubin1987")
            out <- tibble(
              "lower_0.025" = pool$qbar + qnorm(0.5 - 0.95 / 2) * sqrt(pool$t),
              "estimate" = pool$qbar,
              "upper_0.975" = pool$qbar + qnorm(0.5 + 0.95 / 2) * sqrt(pool$t),
              "std_err" = sqrt(pool$t)
            )
            return(list(table = out, pool = pool))
          } else {
            pool <- NULL
            out <- tibble(
              "lower_0.025" = model["estimate"] + qnorm(0.5 - 0.95 / 2) * model["std_err"],
              "estimate" = model["estimate"],
              "upper_0.975" = model["estimate"] + qnorm(0.5 + 0.95 / 2) * model["std_err"],
              "std_err" = model["std_err"]
            )
            return(list(table = out, pool = pool))
          }
        }
      )
    saveRDS(
      cfw_ate_test_agg, 
      paste0(result_path, "cfw_ate_test_agg.rds")
    )
    
    # excess risk reduction table
    list_transpose(cfw_risk_reduction_plot_test_agg)$aipw |>
      map(\(x) mutate(x$table, percentile = 1 - percentile)) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_reduction_aipw_test.xlsx"
        ), 
        rowNames = FALSE
      )
    
    list_transpose(cfw_risk_reduction_plot_test_agg)$ols |>
      map(\(x) mutate(x$table, percentile = 1 - percentile)) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_reduction_ols_test.xlsx"
        ), 
        rowNames = FALSE
      )
    
    # excess risk reduction figure
    pmap(
      list(
        cfw_risk_reduction_plot_test_agg,
        cfw_risk_reduction_plot_test_zp_agg,
        names(cfw_risk_reduction_plot_test_agg)
      ),
      \(x, zp, nm) {
        plot_aipw_exr <- x$aipw$table |>
          slice(2:n()) |>
          ggplot() + 
          geom_linerange(
            aes(ymin = ymin, ymax = ymax, x = x),
            data = tibble(
              ymin = -Inf,
              ymax = x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(estimate),
              x = 100 * 0.1
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_linerange(
            aes(xmin = xmin, xmax = xmax, y = y),
            data = tibble(
              xmin = -Inf,
              xmax = 100 * 0.1,
              y = x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(estimate)
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(xx = 10, yy = -max(x$aipw$table$estimate)/10, label = "10%"), 
            color = "red3",
            size = 4.25
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(
              xx = -11.6, 
              yy = x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(estimate), 
              label = paste0(sprintf("%.2f", 100 * pull(filter(x$aipw$table, percentile == (1 - 0.1)), estimate)), "%")
            ), 
            color = "red3",
            size = 4.25
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate), linewidth = 0.9) +
          geom_hline(yintercept = x$aipw$table$estimate[1], linetype = 2, linewidth = 0.6) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Proportion not treated with thiazide") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(x$aipw$table$estimate)), clip = "off") +
          scale_x_continuous(
            breaks = seq(0, 100, 20),
            labels = paste0(seq(0, 100, 20), "%")
          ) +
          scale_y_continuous(
            breaks = seq(
              0, 
              round(max(x$aipw$table$estimate), 3) + 0.001, 
              trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-3)
            ),
            labels = paste0(
              seq(
                0, 
                100 * round(max(x$aipw$table$estimate), 3) + 0.001, 
                trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-1)
              ), 
              "%"
            )
          )
        plot_aipw_unt <- x$aipw$table |>
          slice(2:n()) |>
          ggplot() + 
          geom_linerange(
            aes(ymin = ymin, ymax = ymax, x = x),
            data = tibble(
              ymin = -Inf,
              ymax = x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(unt_estimate),
              x = 100 * 0.1
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_linerange(
            aes(xmin = xmin, xmax = xmax, y = y),
            data = tibble(
              xmin = -Inf,
              xmax = 100 * 0.1,
              y = x$aipw$table |> filter(percentile == (1 - 0.1)) |> pull(unt_estimate)
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(xx = 10, yy = 5-(max(x$aipw$table$unt_estimate, na.rm = TRUE)-5)/10, label = "10%"), 
            color = "red3",
            size = 4.25
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(
              xx = -8.8, 
              yy = x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(unt_estimate), 
              label = paste0(sprintf("%.0f", pull(filter(x$aipw$table, percentile == (1 - 0.1)), unt_estimate)))
            ), 
            color = "red3",
            size = 4.25
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = unt_estimate), linewidth = 0.9) +
          geom_hline(yintercept = x$aipw$table$unt_estimate[1], linetype = 2, linewidth = 0.6) +
          ylab("Patients to not treat to avoid one TIH case") + xlab("Proportion not treated with thiazide") +
          coord_cartesian(xlim = c(0, 100), ylim = c(5, max(x$aipw$table$unt_estimate, na.rm = TRUE)), clip = "off") +
          scale_x_continuous(
            breaks = seq(0, 100, 20),
            labels = paste0(seq(0, 100, 20), "%")
          ) +
          scale_y_continuous(
            breaks = seq(
              0, 
              trunc(max(x$aipw$table$unt_estimate, na.rm = TRUE)/10+1)*10, 
              trunc(trunc(max(x$aipw$table$unt_estimate, na.rm = TRUE)/10+1)*10/60+0.99)*10
              ),
            labels = paste0(seq(
              0, 
              trunc(max(x$aipw$table$unt_estimate, na.rm = TRUE)/10+1)*10, 
              trunc(trunc(max(x$aipw$table$unt_estimate, na.rm = TRUE)/10+1)*10/60+0.99)*10
            ))
          )
        plot_aipw_comb <- cowplot::plot_grid(
          plot_aipw_exr + xlab(""),
          plot_aipw_unt,
          labels = c("A", "B"),
          label_size = 10,
          nrow = 2,
          ncol = 1,
          align = "v"
        )
        plot_ols <- x$ols$table |>
          slice(2:n()) |>
          ggplot() + 
          geom_linerange(
            aes(ymin = ymin, ymax = ymax, x = x),
            data = tibble(
              ymin = -Inf,
              ymax = x$ols$table|>filter(percentile == (1 - 0.1)) |> pull(estimate),
              x = 100 * 0.1
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_linerange(
            aes(xmin = xmin, xmax = xmax, y = y),
            data = tibble(
              xmin = -Inf,
              xmax = 100 * 0.1,
              y = x$ols$table|>filter(percentile == (1 - 0.1)) |> pull(estimate)
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.6
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(xx = 10, yy = -max(x$ols$table$estimate, na.rm = TRUE)/10, label = "10%"), 
            color = "red3",
            size = 4.25
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(
              xx = -11.6, 
              yy = x$ols$table|>filter(percentile == (1 - 0.1)) |> pull(estimate), 
              label = paste0(sprintf("%.2f", 100 * pull(filter(x$ols$table, percentile == (1 - 0.1)), estimate)), "%")
            ), 
            color = "red3",
            size = 4.25
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate), linewidth = 0.9) +
          geom_hline(yintercept = x$ols$table$estimate[1], linetype = 2, linewidth = 0.6) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Proportion not treated with thiazide") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(x$ols$table$estimate, na.rm = TRUE)), clip = "off") +
          scale_x_continuous(
            breaks = seq(0, 100, 20),
            labels = paste0(seq(0, 100, 20), "%")
          ) +
          scale_y_continuous(
            breaks = seq(
              0, 
              round(max(x$aipw$table$estimate), 3) + 0.001, 
              trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-3)
            ),
            labels = paste0(
              seq(
                0, 
                100 * round(max(x$aipw$table$estimate), 3) + 0.001, 
                trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-1)
              ), 
              "%"
            )
          )
        plot_combined <- bind_rows(
          slice(mutate(x$aipw$table, est = "AIPW"), 2:n()),
          slice(mutate(x$ols$table, est = "MER"), 2:n()) # mean excess risk
        ) |>
          ggplot() + 
          geom_linerange(
            aes(ymin = ymin, ymax = ymax, x = x),
            data = tibble(
              ymin = rep(-Inf, 2),
              ymax = c(
                x$aipw$table |> filter(percentile == (1 - 0.1)) |> pull(estimate),
                x$ols$table |> filter(percentile == (1 - 0.1)) |> pull(estimate)
              ),
              x = rep(100 * 0.1, 2)
            ),
            color = "red3",
            linewidth = 0.6,
            alpha = 0.8
          ) +
          geom_linerange(
            aes(xmin = xmin, xmax = xmax, y = y),
            data = tibble(
              xmin = rep(-Inf, 2),
              xmax = rep(100 * 0.1, 2),
              y = c(
                x$aipw$table |> filter(percentile == (1 - 0.1)) |> pull(estimate),
                x$ols$table |> filter(percentile == (1 - 0.1)) |> pull(estimate)
              )
            ),
            color = "red3",
            linewidth = 0.6,
            alpha = 0.8
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(xx = 10, yy = -max(x$aipw$table$estimate)/10, label = "10%"), 
            color = "red3",
            size = 4.25
          ) +
          geom_text(
            aes(x = xx, y = yy, label = label), 
            data = tibble(
              xx = rep(-11.6, 2), 
              yy = c(
                x$aipw$table|>filter(percentile == (1 - 0.1)) |> pull(estimate),
                x$ols$table|>filter(percentile == (1 - 0.1)) |> pull(estimate)
              ), 
              label = c(
                paste0(sprintf("%.2f", 100 * pull(filter(x$aipw$table, percentile == (1 - 0.1)), estimate)), "%"),
                paste0(sprintf("%.2f", 100 * pull(filter(x$ols$table, percentile == (1 - 0.1)), estimate)), "%")
              )
            ), 
            color = "red3",
            size = 4.25
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate, color = est), linewidth = 0.9) +
          geom_hline(yintercept = x$aipw$table$estimate[1], linetype = 2, linewidth = 0.6) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Proportion not treated with thiazide") +
          scale_color_jama() +
          labs(color = "Estimator") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(c(x$aipw$table$estimate, x$ols$table$estimate))), clip = "off") +
          scale_x_continuous(
            breaks = seq(0, 100, 20),
            labels = paste0(seq(0, 100, 20), "%")
          ) +
          scale_y_continuous(
            breaks = seq(
              0, 
              round(max(x$aipw$table$estimate), 3) + 0.001, 
              trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-3)
            ), #seq(0, 0.015, 0.005),
            labels = paste0(
              seq(
                0, 
                100 * round(max(x$aipw$table$estimate), 3) + 0.001, 
                trunc(10^3 * (round(max(x$aipw$table$estimate), 3) + 0.001) / 3) * 10^(-1)
              ), 
              "%"
            )
          )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_aipw_test", ".jpg"),
          plot = plot_aipw_exr,
          device = "jpeg",
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_aipw_test", ".pdf"),
          plot = plot_aipw_exr,
          device = Cairo::CairoPDF,
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.4,
          family = "Arial Unicode MS",
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_aipw_comb_test", ".jpg"),
          plot = plot_aipw_comb,
          device = "jpeg",
          width = 10,
          height = 15,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_aipw_comb_test", ".pdf"),
          plot = plot_aipw_comb,
          device = Cairo::CairoPDF,
          width = 10,
          height = 15,
          unit = "cm",
          scale = 1.4,
          family = "Arial Unicode MS",
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_ols_test", ".jpg"),
          plot = plot_ols,
          device = "jpeg",
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_comb_test", ".jpg"),
          plot = plot_combined,
          device = "jpeg",
          width = 10,
          height = 7,
          unit = "cm",
          scale = 1.3,
          dpi = 240,
          quality = 90,
          create.dir = TRUE
        )
        return(NULL)
      }
    )
  }
)


### gold standard comparison sensitivity analysis ----------------------------------------------------------------- ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_tiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_tiazide_no_lab/"),
      paste0(figures_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_tiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_tiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_tiazide_cc_comp/1year/")
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid",
      "future", "furrr", "future.callr", "callr", "MatchIt", "ggrepel"
    ),
    globals = c("n_rankings", "vcovHC", "excess_risk_quantiles"),
    seed = TRUE
  ),
  .f = \(path, result_path, table_path, figure_path) {
    cfw_cate_test_agg <- readRDS(
      paste0(
        result_path,
        "cfw_cate_test_with_censored_units_agg.rds"
      )
    )
    plan(multisession, workers = min(num_clusters, length(cfw_cate_test_agg)))
    cfw_gs_sens_test <- future_pmap(
      .l = list(
        gold_std = cfw_cate_test_agg, # gold standard model
        comp_list = map(seq_along(cfw_cate_test_agg), \(i) cfw_cate_test_agg[-i]), # list of comparison models
        gs_nm = names(cfw_cate_test_agg) # name of gold standard model
      ),
      .options = furrr_options(
        packages = c("dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid", "MatchIt", "ggrepel"),
        globals = c("W", "Y", "n_rankings", "vcovHC", "excess_risk_quantiles", "figure_path"),
        seed = TRUE
      ),
      .f = \(gold_std, comp_list, gs_nm) {
        # NOTE: 
        # Since size of positive groups are equal in gold standard and comparison model (pre-specified proportion),
        # the number of false positives is the same as the number of false negatives. Thus any difference in specificity
        # and sensitivity are solely a function of different denominators, i.e. difference between the proportion and
        # one minus the proportion. 
        theme_set(theme_pubr())
        map2(
          comp_list,
          names(comp_list),
          \(comp, comp_nm, gold_std) {
            gold_std_cutoffs <- as.numeric(quantile(gold_std$table$estimate, excess_risk_quantiles))
            comp_cutoffs <- as.numeric(quantile(comp$table$estimate, excess_risk_quantiles))
            
            plot <- tibble(
              gold_standard = gold_std$table$estimate,
              comparison = comp$table$estimate
            ) |>
              ggplot(aes(x = gold_standard, y = comparison)) +
              geom_hex(bins = 100) +
              geom_vline(xintercept = gold_std_cutoffs, linetype = 2) +
              geom_hline(yintercept = comp_cutoffs, linetype = 2) +
              geom_label_repel(
                aes(x = x, y = y, label = label),
                data = tibble(
                  x = gold_std_cutoffs,
                  y = range(comp$table$estimate)[2] - diff(range(comp$table$estimate)) * 0.05,
                  label = sprintf("%.2f", excess_risk_quantiles)
                )
              ) +
              geom_label_repel(
                aes(x = x, y = y, label = label),
                data = tibble(
                  x = range(gold_std$table$estimate)[2] - diff(range(gold_std$table$estimate)) * 0.05,
                  y = comp_cutoffs,
                  label = sprintf("%.2f", excess_risk_quantiles)
                )
              ) +
              geom_label(
                aes(x = x, y = y, label = label),
                data = tibble(
                  x = range(gold_std$table$estimate)[2] - diff(range(gold_std$table$estimate)) * 0.2,
                  y = range(comp$table$estimate)[2] - diff(range(comp$table$estimate)) * 0.05,
                  label = paste("Correlation:", sprintf("%.2f", cor(comp$table$estimate, gold_std$table$estimate)))
                )
              ) +
              xlab(gs_nm) + ylab(comp_nm) +
              theme(
                legend.position = "none"
              )
            ggplot2::ggsave(
              filename = paste0(figure_path, "/", gs_nm, "/gold_standard/", comp_nm, "/cate_scatter_plot_test", ".jpg"),
              plot = plot,
              device = "jpeg",
              width = 10,
              height = 7,
              unit = "cm",
              scale = 1.3,
              dpi = 240,
              quality = 90,
              create.dir = TRUE
            )
            
            out <- map(
              seq_along(excess_risk_quantiles),
              \(i) {
                quantile <- excess_risk_quantiles[i]
                gold_std_cutoff <- gold_std_cutoffs[i]
                comp_cutoff <- comp_cutoffs[i]
                # index of top x% of excess risk estimates in gold standard and comparison models
                gold_std_idx <- which(gold_std$table$estimate >= gold_std_cutoff)
                gold_std_eidx <- which(gold_std$table$estimate == gold_std_cutoff)
                gold_std_nidx <- which(gold_std$table$estimate < gold_std_cutoff)
                comp_idx <- which(comp$table$estimate > comp_cutoff)
                comp_eidx <- which(comp$table$estimate == comp_cutoff)
                comp_nidx <- which(comp$table$estimate < comp_cutoff)
                
                eidx_split <- if (quantile > 0) {
                  min(
                    ceiling(length(comp$table$estimate) * (1 - quantile)) - length(comp_idx),
                    length(comp_eidx)
                  )
                } else {
                  length(comp_eidx)
                }
                comp_idx <- c(comp_idx, comp_eidx[seq_len(eidx_split)])
                comp_nidx <- c(comp_nidx, comp_eidx[seq_len(length(comp_eidx) - eidx_split) + eidx_split])
                tibble(
                  gold_std = gs_nm,
                  comp = comp_nm,
                  quantile = paste0(sprintf("%.0f", 100 * quantile), "%"),
                  sensitivity = sum(gold_std_idx %in% comp_idx) / length(gold_std_idx),
                  specificity = sum(gold_std_nidx %in% comp_nidx) / length(gold_std_nidx)
                )
              }
            ) |>
              list_rbind()
            return(out)
          },
          gold_std = gold_std
        ) |>
          list_rbind()
      }
    )
    plan(sequential)
    
    # save gold standard model comparison results
    saveRDS(cfw_gs_sens_test, paste0(result_path, "cfw_gs_sens_test.rds"))
    
    # save tables with gold standard model comparison results
    cfw_gs_sens_test |>
      dfs2xlsx(paste0(table_path, "cfw_gs_sens_test.xlsx"), rowNames = FALSE)
    
    # list of plots for combined figure
    cfw_gs_plots_test <- pmap(
      .l = list(
        gold_std = cfw_cate_test_agg, # gold standard model
        comp_list = map(seq_along(cfw_cate_test_agg), \(i) cfw_cate_test_agg[-i]), # list of comparison models
        gs_nm = names(cfw_cate_test_agg) # name of gold standard model
      ),
      .f = \(gold_std, comp_list, gs_nm) {
        theme_set(theme_pubr())
        map2(
          comp_list,
          names(comp_list),
          \(comp, comp_nm, gold_std) {
            gold_std_cutoffs <- as.numeric(quantile(gold_std$table$estimate, excess_risk_quantiles))
            comp_cutoffs <- as.numeric(quantile(comp$table$estimate, excess_risk_quantiles))
            
            plot <- tibble(
              gold_standard = gold_std$table$estimate,
              comparison = comp$table$estimate
            ) |>
              ggplot(aes(x = gold_standard, y = comparison)) +
              geom_hex(bins = 100) +
              geom_vline(xintercept = gold_std_cutoffs, linetype = 2) +
              geom_hline(yintercept = comp_cutoffs, linetype = 2) +
              xlab(
                if(
                  gs_nm == "full_mod"  & comp_nm == "vi_mod_20" |
                  gs_nm == "vi_mod_20" & comp_nm == "vi_mod_5"  |
                  gs_nm == "vi_mod_5"  & comp_nm == "vi_mod_4"  |
                  gs_nm == "vi_mod_4"  & comp_nm == "vi_mod_3"  |
                  gs_nm == "vi_mod_3"  & comp_nm == "vi_mod_2"  |
                  gs_nm == "vi_mod_2"  & comp_nm == "vi_mod_1" 
                ) {
                  "Excess risk"
                } else {
                  ""
                }
              ) + 
              ylab(
                if(
                  gs_nm == "full_mod"  & comp_nm == "vi_mod_20" |
                  gs_nm == "vi_mod_20" & comp_nm == "vi_mod_5"  |
                  gs_nm == "vi_mod_5"  & comp_nm == "vi_mod_4"  |
                  gs_nm == "vi_mod_4"  & comp_nm == "vi_mod_3"  |
                  gs_nm == "vi_mod_3"  & comp_nm == "vi_mod_2"  |
                  gs_nm == "vi_mod_2"  & comp_nm == "vi_mod_1" 
                ) {
                  "Excess risk"
                } else {
                  ""
                }
              ) + 
              scale_x_continuous(
                breaks = seq(0, 0.12, 0.03),
                labels = paste0(seq(0, 12, 3), "%")
              ) +
              scale_y_continuous(
                breaks = seq(0, 0.12, 0.03),
                labels = paste0(seq(0, 12, 3), "%")
              ) +
              coord_cartesian(xlim = c(0, 0.12), ylim = c(0, 0.12)) +
              theme(
                legend.position = "none"
              )
            return(plot)
          },
          gold_std = gold_std
        )
      }
    )
    cfw_gs_top_title_test <- c(
      rep(list(ggplot()), 3), 
      list(ggdraw() + draw_label("Compared model", x = 0, hjust = 0, size = 14, fontface = "bold")),
      rep(list(ggplot()), 4)
    )
    cfw_gs_top_titles_test <- map(
      c("20-cov", "5-cov", "4-cov", "3-cov", "2-cov", "1-cov", "", ""),
      \(text) {
        ggdraw() +
          draw_label(
            text,
            x = 0.5,
            hjust = 0.5,
            size = 12
          )
      }
    )
    cfw_gs_right_title_test <- c(
      rep(list(ggplot()), 4), 
      list(
        ggdraw(clip = "off") + 
          draw_label("Reference model", y = 0.25, vjust = 1, size = 14, fontface = "bold", angle = 270)
      ),
      rep(list(ggplot()), 3)
    )
    cfw_gs_right_titles_test <- map(
      c("", "", "66-cov", "20-cov", "5-cov", "4-cov", "3-cov", "2-cov"),
      \(text) {
        ggdraw() +
          draw_label(
            text,
            x = 0.5,
            hjust = 0.5,
            size = 12
          )
      }
    )
    cfw_gs_plot_grid_test <- c(
      cfw_gs_top_title_test,
      cfw_gs_top_titles_test,
      cfw_gs_plots_test$full_mod[6:1], cfw_gs_right_titles_test[3], cfw_gs_right_title_test[3],
      list(ggplot()), cfw_gs_plots_test$vi_mod_20[6:2], cfw_gs_right_titles_test[4], cfw_gs_right_title_test[4],
      rep(list(ggplot()), 2), cfw_gs_plots_test$vi_mod_5[5:2], cfw_gs_right_titles_test[5], cfw_gs_right_title_test[5],
      rep(list(ggplot()), 3), cfw_gs_plots_test$vi_mod_4[4:2], cfw_gs_right_titles_test[6], cfw_gs_right_title_test[6],
      rep(list(ggplot()), 4), cfw_gs_plots_test$vi_mod_3[3:2], cfw_gs_right_titles_test[7], cfw_gs_right_title_test[7],
      rep(list(ggplot()), 5), cfw_gs_plots_test$vi_mod_2[2], cfw_gs_right_titles_test[8], cfw_gs_right_title_test[8]
    )
    cfw_gs_plot_test <- cowplot::plot_grid(
      plotlist = cfw_gs_plot_grid_test,
      nrow = 8,
      ncol = 8,
      rel_heights = c(0.1, 0.1, 1, 1, 1, 1, 1, 1),
      rel_widths = c(1, 1, 1, 1, 1, 1, 0.3, 0.1)
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/gs_comb_plot_test.jpg"),
      plot = cfw_gs_plot_test,
      device = "jpeg",
      width = 24,
      height = 24,
      unit = "cm",
      scale = 2,
      dpi = 180,
      quality = 80,
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/gs_comb_plot_test.pdf"),
      plot = cfw_gs_plot_test,
      device = Cairo::CairoPDF,
      width = 24,
      height = 24,
      unit = "cm",
      scale = 2,
      family = "Arial Unicode MS",
      create.dir = TRUE
    )
    
    return(NULL)
  }
)

### Excess risk distribution plots by covariate combinations ------------------------------------------------------ ----
plan(sequential)
future_pmap(
  .l = list(
    cohort = list(
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
      list(tih_cohort_test), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      list(tih_cohort_test_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
      tih_cohort_test_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
    ),
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_tiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_tiazide_no_lab/"),
      paste0(figures_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_tiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_tiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_tiazide_cc_comp/1year/")
    ),
    continuous = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89])$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      )
    ),
    discrete = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      )
    ),
    plot_model = list(
      "vi_test_mod_4", # tiazide all
      "vi_test_mod_4", # tiazide all severe
      "vi_test_mod_4", # tiazide no lab
      "vi_test_mod_4", # tiazide only fluid lab results
      "vi_test_mod_4", # tiazide only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # bfz all
      "vi_test_mod_4", # bfz only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # hctz all
      "vi_test_mod_4", # hctz only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # tiazide complete case
      "vi_test_mod_4" # tiazide complete case comparison
    ),
    plot_covs = list(
      c("X_01", "X_69", "X_72", "X_80"), # tiazide all
      c("X_01", "X_69", "X_72", "X_80"), # tiazide all severe
      c("X_01", "X_69", "X_72", "X_80"), # tiazide no lab
      c("X_01", "X_69", "X_72", "X_80"), # tiazide only fluid lab results
      c("X_01", "X_69", "X_72", "X_80"), # tiazide only fluid lab results and no Denmark specific covariates
      c("X_01", "X_69", "X_76", "X_80"), # bfz all
      c("X_01", "X_69", "X_72", "X_80"), # bfz only fluid lab results and no Denmark specific covariates
      c("X_01", "X_69", "X_72", "X_78"), # hctz all
      c("X_01", "X_69", "X_72", "X_80"), # hctz only fluid lab results and no Denmark specific covariates
      c("X_01", "X_05", "X_52", "X_69"), # tiazide complete case
      c("X_01", "X_05", "X_52", "X_69") # tiazide complete case comparison
    ),
    cutpoints = list(
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide all severe
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide no lab
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide only fluid lab results
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        `Lactate dehydrogenase` = c(150, 250),
        `C-reactive protein` = c(5, 20)
      ), # bfz all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        Thrombocytes = c(180, 350)
      ), # hctz all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = c(135, 138, 141)
      ), # tiazide complete case
      list(
        Age = c(55, 70),
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = c(135, 138, 141)
      ) # tiazide complete case comparison
    ),
    cut_labels = list(
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide all severe
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide no lab
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide only fluid lab results
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        `Lactate dehydrogenase` = c("< 150 U/L", "150-250 U/L", "\u2265 250 U/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # bfz all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        Thrombocytes = c("< 180 x 10\u2079/L", "180-350 x 10\u2079/L", "\u2265 350 x 10\u2079/L")
      ), # hctz all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        `Years of education` = c("< 10 years", "10-12 years", "\u2265 13 years"),
        `Days of hospitalization in the past year prior to index data` = c("0 days", "1-7 days", "\u2265 8 days"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L")
      ), # thiazide complete case
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        `Years of education` = c("< 10 years", "10-12 years", "\u2265 13 years"),
        `Days of hospitalization in the past year prior to index data` = c("0 days", "1-7 days", "\u2265 8 days"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L")
      ) # thiazide complete case comparison
    ),
    discrete_labels = list(
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide all severe
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide no lab
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide only fluid lab results
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # bfz all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # hctz all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = NA
      ), # tiazide complete case
      list(
        Age = NA,
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = NA
      ) # tiazide complete case comparison
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid",
      "future", "furrr", "future.callr", "callr", "MatchIt", "ggrepel"
    ),
    globals = c("n_rankings", "vcovHC", "excess_risk_quantiles"),
    seed = TRUE
  ),
  .f = \(path, result_path, table_path, figure_path, continuous, discrete, plot_model, plot_covs, cutpoints, cut_labels, discrete_labels) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    # prepare plot data tables for each imputation
    cfw_cate_histogram_data_test <- future_pmap(
      .l = list(
        i = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
        cohort = cohort
      ),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c(
          "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
          "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
          "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
          "excess_risk_cutoff", "path", "result_path", "vi_models", "continuous", "discrete", "plot_model", "plot_covs", 
          "cutpoints", "cut_labels"
        ),
        seed = TRUE
      ),
      .f = \(i, cohort)  {
        cfw_test_mod <- readRDS(paste0(path, "cfw_", plot_model, "_", i, ".rds"))
        cfw_cate_test <- readRDS(paste0(result_path, "cfw_cate_test_with_censored_units_", i, ".rds"))
        cfw_ate_test <- readRDS(paste0(result_path, "cfw_ate_test_", i, ".rds"))
        # covariate data for all (censored and uncensored) observations
        plot_model <- str_remove(plot_model, "_test")
        continuous_cf <- continuous$cf
        discrete_cf <- discrete$cf
        cov_data <- cohort |>
          select({{ continuous_cf }}) |>
          bind_cols(
            cohort |>
              select({{ discrete_cf }}) |>
              DiscreteCovariatesToOneHot()
          ) |>
          select(all_of(names(cfw_test_mod$X.orig))) |>
          as_tibble()
        
        # plot data 
        plot_data <- cov_data |>
          bind_cols(cfw_cate_test[[plot_model]]$table) |>
          mutate(
            ate_diff_signif = 
              (cfw_ate_test[[plot_model]] |> pull(estimate) < 
                 predictions - qnorm(0.975) * sqrt(variance_estimates)) |
              (cfw_ate_test[[plot_model]] |> pull(estimate) >
                 predictions + qnorm(0.975) * sqrt(variance_estimates))
          ) |>
          (\(table) {
            for (i in seq_along(plot_covs)) {
              cov <- plot_covs[i]
              cut <- cutpoints[i]
              lab <- cut_labels[i]
              if (ncol(select(table, starts_with(cov))) > 1) {
                table <- table |>
                  mutate(
                    "{names(cut)}" := (\(tbl) {
                      out <- vector("integer", nrow(tbl))
                      for (j in seq_along(cut[[1]])) {
                        sub <- select(tbl, paste0(cov, "_", cut[[1]][[j]])) |>
                          mutate(a = rowSums(across(everything()))) |>
                          pull(a)
                        out[sub == 1] <- j
                      }
                      out <- factor(out, labels = lab[[1]])
                      return(out)
                    })(tbl = table)
                  )
              } else {
                table <- table |>
                  mutate(
                    "{names(cut)}" := cut(!!sym(cov), c(-Inf, cut[[1]], Inf), labels = lab[[1]], right = FALSE)
                  )
              }
            }
            return(table)
          })()
        
        # center of each group
        group_median <- list(
          Age = if ("X_01" %in% plot_covs) plot_data |> summarise(Age_grp = median(X_01), .by = Age),
          `Years of education` = if ("X_05" %in% plot_covs) plot_data |> count(`Years of education`) |> select("Years of education") |> mutate(`Years of education_grp` = c(6, 11, 16)),
          `Days of hospitalization in the past year prior to index data` = if ("X_05" %in% plot_covs) plot_data |> count(`Days of hospitalization in the past year prior to index data`) |> select("Days of hospitalization in the past year prior to index data") |> mutate(`Days of hospitalization in the past year prior to index data_grp` = c(0, 4, 11)),
          Sodium = if ("X_69" %in% plot_covs) plot_data |> summarise(Sodium_grp = median(X_69), .by = Sodium),
          Hemoglobin = if ("X_72" %in% plot_covs) plot_data |> summarise(Hemoglobin_grp = median(X_72), .by = Hemoglobin),
          `Lactate dehydrogenase` = if ("X_76" %in% plot_covs) plot_data |> summarise(`Lactate dehydrogenase_grp` = median(X_76), .by = `Lactate dehydrogenase`),
          Thrombocytes = if ("X_78" %in% plot_covs) plot_data |> summarise(Thrombocytes_grp = median(X_78), .by = Thrombocytes),
          `C-reactive protein` = if ("X_80" %in% plot_covs) plot_data |> summarise(`C-reactive protein_grp` = median(X_80), .by = `C-reactive protein`)
        )
        group_median <- group_median[!sapply(group_median, is.null)]
        
        # save plot data
        saveRDS(plot_data, paste0(result_path, "cfw_cate_histogram_data_", plot_model, "_test_", i, ".rds"))
        saveRDS(
          group_median, 
          paste0(result_path, "cfw_cate_histogram_group_median_test_", plot_model, "_", i, ".rds"),
        )
        return(plot_data)
      }
    )
    plan(sequential)
    
    # pool tables from each imputation
    cfw_cate_histogram_data_test_agg <- list_rbind(cfw_cate_histogram_data_test) |>
      select(-starts_with("X"))
    saveRDS(
      cfw_cate_histogram_data_test_agg, 
      paste0(result_path, "cfw_cate_histogram_data_", plot_model, "_test_agg.rds")
    )
    
    # average excess risk in population
    cfw_ate_test_agg <- readRDS(paste0(result_path, "cfw_ate_test_agg.rds"))
    ate_all <- cfw_ate_test_agg[[str_remove(plot_model, "_test")]]$table$estimate
    
    label_data_oneway <- map(
      as.list(names(cutpoints)),
      \(var) {
        cfw_cate_histogram_data_test_agg |>
          summarise(
            count = n(),
            lower = sum(predictions < ate_all) / n(),
            higher = 1 - lower,
            .by = all_of(var)
          ) |>
          mutate(
            x_l = 100 * ate_all - 3.5,
            x_h = 100 * ate_all + 3.5,
            y = 2.82,
            lower = paste0(sprintf("%.1f", 100 * lower), "%"),
            higher = paste0(sprintf("%.1f", 100 * higher), "%")
          ) |>
          arrange(!!sym(var))
      }
    ) |>
      structure(names = names(cutpoints))
    label_data_twoway <- map(
      crossing(names(cutpoints), names(cutpoints)) |>
        t() |>
        (\(x) {
          x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
        })() |>
        as_tibble(.name_repair = "unique") |>
        as.list() |>
        structure(names = NULL),
      \(vars) {
        cfw_cate_histogram_data_test_agg |>
          summarise(
            count = n(),
            lower = sum(predictions < ate_all) / n(),
            higher = 1 - lower,
            .by = all_of(vars)
          ) |>
          mutate(
            x_l = 100 * ate_all - 3.5,
            x_h = 100 * ate_all + 3.5,
            y = 2.82,
            lower = paste0(sprintf("%.1f", 100 * lower), "%"),
            higher = paste0(sprintf("%.1f", 100 * higher), "%")
          ) |>
          arrange(!!!syms(vars))
      }
    ) |>
      structure(
        names = crossing(names(cutpoints), names(cutpoints)) |>
          t() |>
          (\(x) {
            x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
          })() |>
          as_tibble(.name_repair = "unique") |>
          as.list() |>
          structure(names = NULL) |>
          map_chr(\(x) paste0(x, collapse = "_"))
      )
    
    # create plots
    n_bins <- round(
      (
        round(100 * max(cfw_cate_histogram_data_test_agg$predictions) * 5) / 5 - 
          round(100 * min(cfw_cate_histogram_data_test_agg$predictions) * 5) / 5
      ) * 5
    ) + 
      1 +
      if (sign(min(label_data_oneway[[1]]$x_l) - 100 * min(cfw_cate_histogram_data_test_agg$predictions)) == -1) {
        round(
          (
            round(100 * min(cfw_cate_histogram_data_test_agg$predictions) * 5) / 5 - 
              round(min(label_data_oneway[[1]]$x_l) * 5) / 5
          ) * 5
        )
      } else {
        0
      }
    map(
      as.list(names(cutpoints)),
      \(var) {
        plot <- cfw_cate_histogram_data_test_agg |>
          arrange(!!sym(var)) |>
          ggplot() +
          geom_histogram(aes(x = 100 * predictions, 
                             y = after_stat(count) *
                               5 / 
                               imap(
                                 label_data_oneway[[var[1]]]$count,
                                 \(x, i) {
                                   if (summarise(cfw_cate_histogram_data_test_agg, sum(!!sym(var) == levels(!!sym(var))[i] & ate_diff_signif))[[1]] == 0) {
                                     out <- rep(x, n_bins)
                                   } else {
                                     out <- rep(x, 2 * n_bins)
                                   }
                                   return(out)
                                 } 
                               ) |>
                               unlist(), 
                             alpha = ate_diff_signif),
                         binwidth = 0.2,
                         fill = "red3") +
          geom_vline(xintercept = 100 * ate_all, linetype = 2, linewidth = 0.25) +
          geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5, nudge_y = 0,
                     data = label_data_oneway[[var[1]]],
                     label.size = 0.12, size = 3,
                     label.padding = unit(0.09, "lines"),
                     label.r = unit(0.06, "lines")) +
          geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5, nudge_y = 0,
                     data = label_data_oneway[[var[1]]],
                     label.size = 0.12, size = 3,
                     label.padding = unit(0.09, "lines"),
                     label.r = unit(0.06, "lines")) +
          facet_wrap(
            vars(!!sym(var[1])),
            nrow = 1,
            labeller = label_both
          ) +
          xlab("Excess risk") + ylab("density") +
          scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) + 
          scale_x_continuous(
            breaks = seq(-5, 15, 5),
            labels = paste0(seq(-5, 15, 5), "%")
          ) +
          scale_y_continuous(breaks = seq(0, 3, 1)) + 
          coord_cartesian(xlim = c(-5, 15), ylim = c(0,3)) +
          theme(
            legend.position = "none",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25),
            strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
            strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/density_", gsub(" ", "", var[1]), "_plot_test", ".jpg"),
          plot = plot,
          device = "jpeg",
          width = 18,
          height = 7,
          unit = "cm",
          scale = 1,
          dpi = 240,
          quality = 90
        )
        return(NULL)
      }
    )
    
    map(
      crossing(names(cutpoints), names(cutpoints)) |>
        t() |>
        (\(x) {
          x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
        })() |>
        as_tibble(.name_repair = "unique") |>
        as.list() |>
        structure(names = NULL),
      \(vars) {
        plot <- cfw_cate_histogram_data_test_agg |>
          arrange(!!!syms(vars)) |>
          ggplot() +
          geom_histogram(aes(x = 100 * predictions, 
                             y = after_stat(count) *
                               5 / 
                               pmap(
                                 list(
                                   x = label_data_twoway[[paste0(vars[1], "_", vars[2])]]$count,
                                   i = as.integer(label_data_twoway[[paste0(vars[1], "_", vars[2])]][[vars[1]]]),
                                   j = as.integer(label_data_twoway[[paste0(vars[1], "_", vars[2])]][[vars[2]]])
                                 ),
                                 \(x, i, j) {
                                   if (summarise(cfw_cate_histogram_data_test_agg, sum(!!sym(vars[1]) == levels(!!sym(vars[1]))[i] & !!sym(vars[2]) == levels(!!sym(vars[2]))[j] & ate_diff_signif))[[1]] == 0) {
                                     out <- rep(x, n_bins)
                                   } else {
                                     out <- rep(x, 2 * n_bins)
                                   }
                                   return(out)
                                 } 
                               ) |>
                               unlist(), 
                             alpha = ate_diff_signif),
                         binwidth = 0.2,
                         fill = "red3") +
          geom_vline(xintercept = 100 * ate_all, linetype = 2, linewidth = 0.25) +
          geom_label(aes(x = x_l, y = y, label = lower), vjust = 0.5, nudge_y = 0,
                     data = label_data_twoway[[paste0(vars[1], "_", vars[2])]],
                     label.size = 0.12, size = 3,
                     label.padding = unit(0.09, "lines"),
                     label.r = unit(0.06, "lines")) +
          geom_label(aes(x = x_h, y = y, label = higher), vjust = 0.5, nudge_y = 0,
                     data = label_data_twoway[[paste0(vars[1], "_", vars[2])]],
                     label.size = 0.12, size = 3,
                     label.padding = unit(0.09, "lines"),
                     label.r = unit(0.06, "lines")) +
          facet_grid(
            as.formula(
              glue(
                "
                {if(vars[1] == 'C-reactive protein') '`C-reactive protein`' else if (vars[1] == 'Lactate dehydrogenase') '`Lactate dehydrogenase`' else if (vars[1] == 'Days of hospitalization in the past year prior to index data') '`Days of hospitalization in the past year prior to index data`' else if (vars[1] == 'Years of education') '`Years of education`' else vars[1]} ~ 
                {if(vars[2] == 'C-reactive protein') '`C-reactive protein`' else if (vars[2] == 'Lactate dehydrogenase') '`Lactate dehydrogenase`' else if (vars[2] == 'Days of hospitalization in the past year prior to index data') '`Days of hospitalization in the past year prior to index data`' else if (vars[2] == 'Years of education') '`Years of education`' else vars[2]}
                "
              )
            ),
            labeller = label_both
          ) +
          xlab("Excess risk") + ylab("density") +
          scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5)) + 
          scale_x_continuous(
            breaks = seq(-5, 15, 5),
            labels = paste0(seq(-5, 15, 5), "%")
          ) +
          scale_y_continuous(breaks = seq(0, 3, 1)) + 
          coord_cartesian(xlim = c(-5, 15), ylim = c(0,3)) +
          theme(
            legend.position = "none",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25),
            strip.background = element_rect(linewidth = 0.25, fill = "#eff3f2"),
            strip.text = element_text(family = "sans", margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/density_", gsub(" ", "", vars[1]), "_", gsub(" ", "", vars[2]), "_plot_test", ".jpg"),
          plot = plot,
          device = "jpeg",
          width = 18,
          height = 14.5,
          unit = "cm",
          scale = 1,
          dpi = 240,
          quality = 90
        )
        return(NULL)
      }
    )
    return(NULL)
  }
)


### plots with average risk of thiazide-induced hyponatraemia in subgroups by covariate combinations -------------- ----
plan(sequential)
future_pmap(
  .l = list(
    cohort = list(
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
      list(tih_cohort_test), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
      tih_tbi_1_test_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      bfz_tbi_1_test_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
      hctz_tbi_1_test_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
      list(tih_cohort_test_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
      tih_cohort_test_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
    ),
    path = list(
      paste0(training_path, "4mo_tiazide_all/1year/"),
      paste0(training_path, "3mo_tiazide_all_severe/1year/"),
      paste0(training_path, "4mo_tiazide_no_lab/"),
      paste0(training_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_tiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_tiazide_complete_case/1year/"),
      paste0(training_path, "4mo_tiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_tiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_tiazide_no_lab/"),
      paste0(results_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_tiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_tiazide_complete_case/1year/"),
      paste0(results_path, "4mo_tiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_tiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_tiazide_no_lab/"),
      paste0(tables_path, "4mo_tiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_tiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_tiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_tiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_tiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_tiazide_no_lab/"),
      paste0(figures_path, "4mo_tiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_tiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_tiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_tiazide_cc_comp/1year/")
    ),
    continuous = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50)$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89])$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[01])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        exp_out = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        cf = expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[013])$")),
        apriori = "(^X_(0[127]|1[5789]|2[1-8]|3[456]|4[89]|69|70)$|"
      )
    ),
    discrete = list(
      list(# tiazide all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# thiazide all severe
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide no lab
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# bfz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# bfz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# hctz all
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# hctz only fluid lab results and no Denmark specific covariates
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(32|5[1-4])$")),
        apriori = "(^X_32$|"
      ),
      list(# tiazide complete case
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      ),
      list(# tiazide complete case comparison
        sample_weight = expr(matches("^X_(0[45]|32|5[1-4])$")),
        exp_out = expr(matches("^X_(0[45]|32|5[1-4])$")),
        cf = expr(matches("^X_(0[45]|32|5[1-4])$")),
        apriori = "(^X_(0[45]|32)$|"
      )
    ),
    plot_model = list(
      "vi_test_mod_4", # tiazide all
      "vi_test_mod_4", # tiazide all severe
      "vi_test_mod_4", # tiazide no lab
      "vi_test_mod_4", # tiazide only fluid lab results
      "vi_test_mod_4", # tiazide only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # bfz all
      "vi_test_mod_4", # bfz only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # hctz all
      "vi_test_mod_4", # hctz only fluid lab results and no Denmark specific covariates
      "vi_test_mod_4", # tiazide complete case
      "vi_test_mod_4" # tiazide complete case comparison
    ),
    plot_covs = list(
      c("X_01", "X_69", "X_72", "X_80"), # tiazide all
      c("X_01", "X_69", "X_72", "X_80"), # tiazide all severe
      c("X_01", "X_69", "X_72", "X_80"), # tiazide no lab
      c("X_01", "X_69", "X_72", "X_80"), # tiazide only fluid lab results
      c("X_01", "X_69", "X_72", "X_80"), # tiazide only fluid lab results and no Denmark specific covariates
      c("X_01", "X_69", "X_76", "X_80"), # bfz all
      c("X_01", "X_69", "X_72", "X_80"), # bfz only fluid lab results and no Denmark specific covariates
      c("X_01", "X_69", "X_72", "X_78"), # hctz all
      c("X_01", "X_69", "X_72", "X_80"), # hctz only fluid lab results and no Denmark specific covariates
      c("X_01", "X_05", "X_52", "X_69"), # tiazide complete case
      c("X_01", "X_05", "X_52", "X_69") # tiazide complete case comparison
    ),
    cutpoints = list(
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide all severe
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide no lab
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide only fluid lab results
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # tiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        `Lactate dehydrogenase` = c(150, 250),
        `C-reactive protein` = c(5, 20)
      ), # bfz all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        Thrombocytes = c(180, 350)
      ), # hctz all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = c(55, 70),
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = c(135, 138, 141)
      ), # tiazide complete case
      list(
        Age = c(55, 70),
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = c(135, 138, 141)
      ) # tiazide complete case comparison
    ),
    cut_labels = list(
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide all severe
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide no lab
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide only fluid lab results
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # thiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        `Lactate dehydrogenase` = c("< 150 U/L", "150-250 U/L", "\u2265 250 U/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # bfz all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        Thrombocytes = c("< 180 x 10\u2079/L", "180-350 x 10\u2079/L", "\u2265 350 x 10\u2079/L")
      ), # hctz all
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L"),
        Hemoglobin = c("< 8 mmol/L", "8-9 mmol/L", "9-10 mmol/L", "\u2265 10 mmol/L"),
        `C-reactive protein` = c("< 5 mg/L", "5-20 mg/L", "\u2265 20 mg/L")
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        `Years of education` = c("< 10 years", "10-12 years", "\u2265 13 years"),
        `Days of hospitalization in the past year prior to index data` = c("0 days", "1-7 days", "\u2265 8 days"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L")
      ), # thiazide complete case
      list(
        Age = c("40-55 years", "55-70 years", "\u2265 70 years"),
        `Years of education` = c("< 10 years", "10-12 years", "\u2265 13 years"),
        `Days of hospitalization in the past year prior to index data` = c("0 days", "1-7 days", "\u2265 8 days"),
        Sodium = c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L")
      ) # thiazide complete case comparison
    ),
    discrete_labels = list(
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide all severe
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide no lab
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide only fluid lab results
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # tiazide only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # bfz all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # bfz only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # hctz all
      list(
        Age = NA,
        Sodium = NA,
        Hemoglobin = NA,
        `C-reactive protein` = NA
      ), # hctz only fluid lab results and no Denmark specific covariates
      list(
        Age = NA,
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = NA
      ), # tiazide complete case
      list(
        Age = NA,
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = NA
      ) # tiazide complete case comparison
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid",
      "future", "furrr", "future.callr", "callr", "MatchIt", "ggrepel"
    ),
    globals = c(
      "n_rankings", "vcovHC", "excess_risk_quantiles", "vi_models", "num_threads", "num_clusters", "names_index", 
      "compute_results", "CausalForestCATEAllTable", "ATEAll", "CausalForestATEAllSubgroupTable"
    ),
    seed = TRUE
  ),
  .f = \(path, result_path, table_path, figure_path, continuous, discrete, plot_model, plot_covs, cutpoints, cut_labels, discrete_labels) {
    mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
    if (num_clusters == 1 || !mice) {
      plan(sequential)
    } else {
      plan(multisession, workers = num_clusters)
    }
    # prepare plot data tables for each imputation
    cfw_subgroup_cate_test <- future_pmap(
      .l = list(
        i = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}"))
      ),
      .options = furrr_options(
        packages = c(
          "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
          "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
          "callr", "fuzzyjoin"
        ),
        globals = c(
          "DiscreteCovariatesToOneHot", "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results",
          "names_index_comb", "mice", "excess_risk_cutoff", "path", "result_path", "vi_models", "continuous", 
          "discrete", "plot_model", "plot_covs", "cutpoints", "cut_labels", "discrete_labels",
          "CausalForestCATEAllTable", "ATEAll", "CausalForestATEAllSubgroupTable"
        ),
        seed = seed,
        stdout = FALSE
      ),
      .f = \(i)  {
        cfw_test_mod <- readRDS(paste0(path, "cfw_", plot_model, "_", i, ".rds"))
        cfw_cate_test <- readRDS(paste0(result_path, "cfw_cate_test_", i, ".rds"))[[str_remove(plot_model, "_test")]]$table
        group_median <- readRDS(
          paste0(
            result_path,
            "cfw_cate_histogram_group_median_test_",
            str_remove(plot_model, "_test"),
            "_",
            i,
            ".rds"
          )
        )
        
        cfw_subgroup_cate_test <- list()
        cfw_subgroup_cate_test$oneway <- map(
          as.list(names(cutpoints)),
          \(var) {
            cate_table <- CausalForestCATEAllTable(
              X.orig = (\(mod) {
                if (length(names(cutpoints)) == length(names(mod$X.orig))) {
                  names(mod$X.orig) <- names(cutpoints)
                } else {
                  nm_dis <- unique(str_extract(names(mod$X.orig[grep("^X_\\d{2}_", names(mod$X.orig))]), "^X_\\d{2}"))
                  nm <- names(mod$X.orig)
                  for (i in seq_along(plot_covs)) {
                    nm <- str_replace(nm, plot_covs[i], names(cutpoints)[i])
                  }
                  names(mod$X.orig) <- nm
                }
                return(mod$X.orig)
              })(mod = cfw_test_mod),
              Y.orig = cfw_test_mod$Y.orig,
              Y.hat = cfw_test_mod$Y.hat,
              W.orig = cfw_test_mod$W.orig,
              W.hat = cfw_test_mod$W.hat,
              tau.hat.pointwise = cfw_cate_test$predictions,
              sample.weights = cfw_test_mod$sample.weights,
              cov_list = if (sum(grepl(plot_covs[names(cutpoints) == var], names(cfw_test_mod$X.orig))) == 1) {
                list2(
                  "{var[1]}" := structure(map2(c(-Inf, cutpoints[[var[1]]]), c(cutpoints[[var[1]]], Inf), \(x, y) as.character(c(x, y))), names = cut_labels[[var[1]]])
                )
              } else {
                list2(
                  "{var[1]}" := structure(map(discrete_labels[[var]], \(x) paste0(var, "_", x)), names = cut_labels[[var]])
                )
              }
            ) |>
              select(1:4) |>
              left_join(group_median[[var]], by = var)
          }
        ) |>
          structure(names = names(cutpoints))
        
        cfw_subgroup_cate_test$twoway <-
          map(
            crossing(names(cutpoints), names(cutpoints)) |>
              t() |>
              (\(x) {
                x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
              })() |>
              as_tibble(.name_repair = "unique") |>
              as.list() |>
              structure(names = NULL),
            \(vars) {
              cate_table <- CausalForestCATEAllTable(
                X.orig = (\(mod) {
                  if (length(names(cutpoints)) == length(names(mod$X.orig))) {
                    names(mod$X.orig) <- names(cutpoints)
                  } else {
                    nm_dis <- unique(str_extract(names(mod$X.orig[grep("^X_\\d{2}_", names(mod$X.orig))]), "^X_\\d{2}"))
                    nm <- names(mod$X.orig)
                    for (i in seq_along(plot_covs)) {
                      nm <- str_replace(nm, plot_covs[i], names(cutpoints)[i])
                    }
                    names(mod$X.orig) <- nm
                  }
                  return(mod$X.orig)
                })(mod = cfw_test_mod),
                Y.orig = cfw_test_mod$Y.orig,
                Y.hat = cfw_test_mod$Y.hat,
                W.orig = cfw_test_mod$W.orig,
                W.hat = cfw_test_mod$W.hat,
                tau.hat.pointwise = cfw_cate_test$predictions,
                sample.weights = cfw_test_mod$sample.weights,
                cov_list = if (
                  sum(grepl(plot_covs[names(cutpoints) == vars[1]], names(cfw_test_mod$X.orig))) == 1 &
                  sum(grepl(plot_covs[names(cutpoints) == vars[2]], names(cfw_test_mod$X.orig))) == 1
                ) {
                  list2(
                    "{vars[1]}" := rep(structure(map2(c(-Inf, cutpoints[[vars[1]]]), c(cutpoints[[vars[1]]], Inf), \(x, y) as.character(c(x, y))), names = cut_labels[[vars[1]]]), times = length(cutpoints[[vars[2]]]) + 1),
                    "{vars[2]}" := rep(structure(map2(c(-Inf, cutpoints[[vars[2]]]), c(cutpoints[[vars[2]]], Inf), \(x, y) as.character(c(x, y))), names = cut_labels[[vars[2]]]), each = length(cutpoints[[vars[1]]]) + 1)
                  )
                } else if (
                  sum(grepl(plot_covs[names(cutpoints) == vars[1]], names(cfw_test_mod$X.orig))) == 1
                ) {
                  list2(
                    "{vars[1]}" := rep(structure(map2(c(-Inf, cutpoints[[vars[1]]]), c(cutpoints[[vars[1]]], Inf), \(x, y) as.character(c(x, y))), names = cut_labels[[vars[1]]]), times = length(cutpoints[[vars[2]]])),
                    "{vars[2]}" := rep(structure(map(discrete_labels[[vars[2]]], \(x) paste0(vars[2], "_", x)), names = cut_labels[[vars[2]]]), each = length(cutpoints[[vars[1]]]) + 1)
                  )
                } else if (
                  sum(grepl(plot_covs[names(cutpoints) == vars[2]], names(cfw_test_mod$X.orig))) == 1
                ) {
                  list2(
                    "{vars[1]}" := rep(structure(map(discrete_labels[[vars[1]]], \(x) paste0(vars[1], "_", x)), names = cut_labels[[vars[1]]]), times = length(cutpoints[[vars[2]]]) + 1),
                    "{vars[2]}" := rep(structure(map2(c(-Inf, cutpoints[[vars[2]]]), c(cutpoints[[vars[2]]], Inf), \(x, y) as.character(c(x, y))), names = cut_labels[[vars[2]]]), each = length(cutpoints[[vars[1]]]))
                  )
                } else {
                  list2(
                    "{vars[1]}" := rep(structure(map(discrete_labels[[vars[1]]], \(x) paste0(vars[1], "_", x)), names = cut_labels[[vars[1]]]), times = length(cutpoints[[vars[2]]])),
                    "{vars[2]}" :=rep(structure(map(discrete_labels[[vars[2]]], \(x) paste0(vars[2], "_", x)), names = cut_labels[[vars[2]]]), each = length(cutpoints[[vars[1]]]))
                  )
                }
              ) |>
                select(1:4) |>
                left_join(
                  bind_cols(
                    slice(arrange(group_median[[vars[1]]], !!sym(vars[1])), rep(seq_len(n()), times = nrow(group_median[[vars[2]]]))),
                    slice(arrange(group_median[[vars[2]]], !!sym(vars[2])), rep(seq_len(n()), each = nrow(group_median[[vars[1]]])))
                  ) |>
                    transmute(
                      "{paste0(vars[1], '_', vars[2])}" := paste0(!!sym(vars[1]), "_", !!sym(vars[2])),
                      "{paste0(vars[1], '_grp')}" := !!sym(paste0(vars[1], "_grp")),
                      "{paste0(vars[2], '_grp')}" := !!sym(paste0(vars[2], "_grp"))
                    ),
                  by = paste0(vars[1], '_', vars[2])
                )
            }
          ) |>
          structure(
            names = crossing(names(cutpoints), names(cutpoints)) |>
              t() |>
              (\(x) {
                x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
              })() |>
              as_tibble(.name_repair = "unique") |>
              as.list() |>
              structure(names = NULL) |>
              map_chr(\(x) paste0(x, collapse = "_"))
          )
        
        saveRDS(
          cfw_subgroup_cate_test, 
          paste0(
            result_path,
            "cfw_subgroup_cate_test_",
            plot_model,
            "_",
            i,
            ".rds"
          )
        )
        
        
        return(cfw_subgroup_cate_test)
      }
    )
    
    # pool tables from each imputation
    cfw_subgroup_cate_test_agg <- cfw_subgroup_cate_test |>
      list_transpose(simplify = FALSE) |>
      map(\(x) purrr::list_transpose(x, simplify = FALSE)) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(n_cov) {
          future_map(
            n_cov,
            .options = furrr_options(
              packages = c("mice", "dplyr", "purrr")
            ),
            .f = \(vars_data) {
              out <- vars_data |>
                purrr::list_transpose(simplify = FALSE) |>
                map(\(x) {
                  out <- list()
                  for(i in seq_along(x[[1]])) {
                    out[[i]] <- vector(typeof(x[[1]]), length = length(x))
                    names(out)[i] <- names(x[[1]])[i]
                    for (j in seq_along(x)) {
                      out[[i]][j] <- x[[j]][i]
                    }
                    names(out[[i]]) <- names(x)
                  }
                  return(out)
                }) |>
                purrr::list_transpose(simplify = FALSE) |>
                map(
                  \(x) {
                    if (str_detect(names(x)[1], "_")) {
                      var1 <- str_extract(names(x)[1], ".+?(?=_)")
                      var2 <- str_extract(names(x)[1], "(?<=_).+")
                    } else {
                      var <- names(x)[1]
                    }
                    if (is.list(x) && length(x[[1]]) > 1) {
                      pool <- pool.scalar(Q = x$estimate, U = (x$std_err)^2, n = x$n[1], k = 1, rule = "rubin1987")
                      out <- tibble(
                        "{names(x)[1]}" := x[[1]][1],
                        estimate = pool$qbar,
                        std_err = sqrt(pool$t),
                        "95% CI - lower" := estimate - qnorm(0.975) * std_err,
                        "95% CI - upper" := estimate + qnorm(0.975) * std_err,
                        n = x$n[1],
                        `estimate (%)` = sprintf("%.1f", 100 * estimate),
                        "95% CI - lower (%)" := sprintf("%.1f", 100 * (estimate - qnorm(0.975) * std_err)),
                        "95% CI - upper (%)" := sprintf("%.1f", 100 * (estimate + qnorm(0.975) * std_err))
                      )
                    } else {
                      pool <- NULL
                      out <- tibble(
                        "{names(x)[1]}" := x[[1]][1],
                        estimate = x$estimate,
                        std_err = x$std_err,
                        "95% CI - lower" := estimate - qnorm(0.975) * std_err,
                        "95% CI - upper" := estimate + qnorm(0.975) * std_err,
                        n = x$n[1],
                        `estimate (%)` = sprintf("%.1f", 100 * estimate),
                        "95% CI - lower (%)" := sprintf("%.1f", 100 * (estimate - qnorm(0.975) * std_err)),
                        "95% CI - upper (%)" := sprintf("%.1f", 100 * (estimate + qnorm(0.975) * std_err))
                      )
                    }
                    if (str_detect(names(x)[1], "_")) {
                      out <- out |>
                        mutate(
                          "{var1}" := str_extract(x[[1]][1], ".+?(?=_)"),
                          "{var2}" := str_extract(x[[1]][1], "(?<=_).+"),
                          "{var1}_grp" := x[[paste0(var1, "_grp")]][1],
                          "{var2}_grp" := x[[paste0(var2, "_grp")]][1]
                        )
                    } else {
                      out <- mutate(out, "{var}_grp" := x[[paste0(var, "_grp")]][1])
                    }
                    return(list(table = out, pool_results = pool))
                  }
                ) |>
                list_transpose()
              if ("Age" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  Age = factor(Age, levels = cut_labels$Age)
                )
              }
              if ("Sodium" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  Sodium = factor(Sodium, levels = cut_labels$Sodium)
                )
              }
              if ("Hemoglobin" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  Hemoglobin = factor(Hemoglobin, levels = cut_labels$Hemoglobin)
                )
              }
              if ("C-reactive protein" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  `C-reactive protein` = factor(`C-reactive protein`, levels = cut_labels$`C-reactive protein`)
                )
              }
              if ("Lactate dehydrogenase" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  `Lactate dehydrogenase` = factor(`Lactate dehydrogenase`, levels = cut_labels$`Lactate dehydrogenase`)
                )
              }
              if ("Thrombocytes" %in% names(out$table)) {
                out$table <- mutate(
                  out$table,
                  Thrombocytes = factor(Thrombocytes, levels = cut_labels$Thrombocytes)
                )
              }
              return(out)
            }
          )
        }
      )
    plan(sequential)
    saveRDS(
      cfw_subgroup_cate_test_agg, 
      paste0(
        result_path,
        "cfw_subgroup_cate_test_agg_",
        plot_model,
        ".rds"
      )
    )
    
    # average excess risk in population
    cfw_ate_test_agg <- readRDS(paste0(result_path, "cfw_ate_test_agg.rds"))
    ate_all <- cfw_ate_test_agg[[str_remove(plot_model, "_test")]]$table$estimate
    
    # plots with subgroup ATEs
    map(
      as.list(names(cutpoints)),
      \(var) {
        plot <- cfw_subgroup_cate_test_agg$oneway[[var]]$table |>
          ggplot() +
          geom_line(aes(x = !!sym(paste0(var, "_grp")), 
                        y = 100 * estimate),
                    linewidth = 0.7) +
          geom_errorbar(aes(x = !!sym(paste0(var, "_grp")),
                            ymin = 100 * `95% CI - lower`, 
                            ymax = 100 * `95% CI - upper`),
                        linewidth = 0.7,
                        width = 0) +
          geom_point(aes(x = !!sym(paste0(var, "_grp")), y = 100 * estimate), 
                     size = 1) +
          geom_hline(yintercept = 100 * ate_all, linetype = 2, linewidth = 0.3) +
          xlab(glue("{var} ({if(var %in% c('Age', 'Years of education')) 'years' else if (var == 'C-reactive protein') 'mg/L' else if (var == 'Lactate dehydrogenase') 'U/L' else if (var == 'Thrombocytes') 'x 10\u2079/L' else if (var == 'Sodium') 'mmol/L' else 'days'})")) +
          ylab("Excess risk") +
          scale_x_continuous(
            breaks = cfw_subgroup_cate_test_agg$oneway[[var]]$table[[paste0(var, "_grp")]],
            labels = str_extract(cfw_subgroup_cate_test_agg$oneway[[var]]$table[[var]], ".+?(?=\\s[a-zA-Z])"),
            guide = guide_axis(angle = 0)
          ) +
          scale_y_continuous(
            breaks = seq(-5, 20, 5),
            labels = paste0(seq(-5, 20, 5), "%")
          ) +
          coord_cartesian(xlim = if(var == 'Age') c(45, 80) else if (var == 'Sodium') c(131, 144) else if (var == 'Hemoglobin') c(7, 11) else if (var == 'Lactate dehydrogenase') c(130, 290) else if (var == "Thrombocytes") c(155, 390) else if (var == "C-reactive protein") c(0, 40) else if (var == "Years of education") c(4, 18) else c(-2, 13), ylim = c(-2, 10)) +
          theme(
            legend.title = element_blank(),
            legend.position = "top",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25)
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/subgroup_cate_", gsub(" ", "", var[1]), "_plot_test", ".jpg"),
          plot = plot,
          device = "jpeg",
          width = 9,
          height = 7,
          unit = "cm",
          scale = 1,
          dpi = 240,
          quality = 90
        )
        return(NULL)
      }
    )
    
    plot_list <- map(
      crossing(names(cutpoints), names(cutpoints)) |>
        t() |>
        (\(x) {
          x[, !duplicated(t(apply(x, 2, sort))) & apply(x, 2, \(x) length(unique(x)) > 1)]
        })() |>
        as_tibble(.name_repair = "unique") |>
        as.list() |>
        structure(names = NULL),
      \(vars) {
        bar_width <- switch(vars[1], Thrombocytes = 8, `Lactate dehydrogenase` = 6, `C-reactive protein` = , Age = 2, Sodium = 0.9, `Years of education` = , `Days of hospitalization in the past year prior to index data` = 0.7, Hemoglobin = 0.25)
        dodge_width <- 1.5 * switch(vars[1], Thrombocytes = 8, `Lactate dehydrogenase` = 6, `C-reactive protein` = , Age = 2, Sodium = 0.9, `Years of education` = , `Days of hospitalization in the past year prior to index data` = 0.7, Hemoglobin = 0.25)
        guide_nrow <- switch(vars[2], Thrombocytes = , Sodium = , Hemoglobin = , `Years of education` = , `Days of hospitalization in the past year prior to index data` = , `Lactate dehydrogenase` = , Age = , `C-reactive protein` = 2)
        guide_ncol <- switch(vars[2], Thrombocytes = , Sodium = , Hemoglobin = , `Years of education` = , `Days of hospitalization in the past year prior to index data` = , `Lactate dehydrogenase` = , Age = , `C-reactive protein` = 2)
        plot <- cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table |>
          ggplot() +
          geom_line(aes(x = !!sym(paste0(vars[1], "_grp")), 
                        y = 100 * estimate,
                        color = !!sym(vars[2])),
                    position = position_dodge(width = dodge_width),
                    linewidth = 0.2, linetype = 2) +
          geom_errorbar(aes(x = !!sym(paste0(vars[1], "_grp")),
                            ymin = 100 * `95% CI - lower`, 
                            ymax = 100 * `95% CI - upper`,
                            color = !!sym(vars[2])),
                        position = position_dodge(width = dodge_width),
                        linewidth = 0.2,
                        width = bar_width) +
          geom_point(aes(x = !!sym(paste0(vars[1], "_grp")),
                         y = 100 * estimate, 
                         color = !!sym(vars[2])), 
                     position = position_dodge(width = dodge_width),
                     size = 1.1) +
          geom_hline(yintercept = 100 * ate_all, linetype = 2, linewidth = 0.3) +
          xlab(glue("{vars[1]} ({if(vars[1] %in% c('Age', 'Years of education')) 'years' else if (vars[1] == 'C-reactive protein') 'mg/L' else if (vars[1] == 'Lactate dehydrogenase') 'U/L' else if (vars[1] == 'Thrombocytes') 'x 10\u2079/L' else if (vars[1] == 'Sodium') 'mmol/L' else 'days'})")) +
          ylab("Excess risk") +
          scale_x_continuous(
            breaks = cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table[[paste0(vars[1], "_grp")]],
            labels = str_extract(cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table[[vars[1]]], ".+?(?=\\s[a-zA-Z])"),
            guide = guide_axis(angle = 0)
          ) +
          scale_y_continuous(
            breaks = seq(-5, 20, 5),
            labels = paste0(seq(-5, 20, 5), "%")
          ) +
          coord_cartesian(xlim = if(vars[1] == 'Age') c(45, 80) else if (vars[1] == 'Sodium') c(131, 144) else if (vars[1] == 'Hemoglobin') c(7, 11) else if (vars[1] == 'Lactate dehydrogenase') c(130, 290) else if (vars[1] == "Thrombocytes") c(155, 390) else if (vars[1] == "C-reactive protein") c(0, 40) else if (vars[1] == "Years of education") c(4, 18) else c(-2, 13), ylim = c(-5, 20)) +
          scale_color_jama() +
          guides(color = guide_legend(nrow = guide_nrow, ncol = guide_ncol, byrow = TRUE)) +
          theme(
            legend.position = "top",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25)
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/subgroup_cate_", gsub(" ", "", vars[1]), "_", gsub(" ", "", vars[2]), "_plot_test", ".jpg"),
          plot = plot,
          device = "jpeg",
          width = 9,
          height = 7,
          unit = "cm",
          scale = 1,
          dpi = 240,
          quality = 90
        )
        plot_nl <- cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table |>
          ggplot() +
          geom_errorbar(aes(x = !!sym(paste0(vars[1], "_grp")),
                            ymin = 100 * `95% CI - lower`, 
                            ymax = 100 * `95% CI - upper`,
                            color = !!sym(vars[2])),
                        position = position_dodge(width = dodge_width),
                        linewidth = 0.7,
                        width = 0) +
          geom_point(aes(x = !!sym(paste0(vars[1], "_grp")),
                         y = 100 * estimate, 
                         color = !!sym(vars[2])), 
                     position = position_dodge(width = dodge_width),
                     size = 1.2) +
          geom_hline(yintercept = 100 * ate_all, linetype = 2, linewidth = 0.3) +
          xlab(glue("{vars[1]} ({if(vars[1] %in% c('Age', 'Years of education')) 'years' else if (vars[1] == 'C-reactive protein') 'mg/L' else if (vars[1] == 'Lactate dehydrogenase') 'U/L' else if (vars[1] == 'Thrombocytes') 'x 10\u2079/L' else if (vars[1] == 'Sodium') 'mmol/L' else 'days'})")) +
          ylab("Excess risk") +
          scale_x_continuous(
            breaks = cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table[[paste0(vars[1], "_grp")]],
            labels = str_extract(cfw_subgroup_cate_test_agg$twoway[[paste0(vars[1], "_", vars[2])]]$table[[vars[1]]], ".+?(?=\\s[a-zA-Z])"),
            guide = guide_axis(angle = 0)
          ) +
          scale_y_continuous(
            breaks = seq(-5, 20, 5),
            labels = paste0(seq(-5, 20, 5), "%")
          ) +
          coord_cartesian(xlim = if(vars[1] == 'Age') c(45, 80) else if (vars[1] == 'Sodium') c(131, 144) else if (vars[1] == 'Hemoglobin') c(7, 11) else if (vars[1] == 'Lactate dehydrogenase') c(130, 290) else if (vars[1] == "Thrombocytes") c(155, 390) else if (vars[1] == "C-reactive protein") c(0, 40) else if (vars[1] == "Years of education") c(4, 18) else c(-2, 13), ylim = c(-5, 20)) +
          scale_color_jama() +
          guides(color = guide_legend(nrow = guide_nrow, ncol = guide_ncol, byrow = TRUE)) +
          theme(
            legend.position = "top",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25)
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/subgroup_cate_nl_", gsub(" ", "", vars[1]), "_", gsub(" ", "", vars[2]), "_plot_test", ".jpg"),
          plot = plot_nl,
          device = "jpeg",
          width = 9,
          height = 7,
          unit = "cm",
          scale = 1,
          dpi = 240,
          quality = 90
        )
        return(list(line = plot, no_line = plot_nl))
      }
    ) |>
      list_transpose()
    
    subgroup_plot_comb <- cowplot::plot_grid(
      plotlist = plot_list$line,
      nrow = 3, 
      ncol = 2,
      labels = LETTERS[1:6],
      label_size = 10
    )
    subgroup_plot_nl_comb <- cowplot::plot_grid(
      plotlist = plot_list$no_line,
      nrow = 3, 
      ncol = 2,
      labels = LETTERS[1:6],
      label_size = 10
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot_test", ".jpg"),
      plot = subgroup_plot_comb,
      device = "jpeg",
      width = 18,
      height = 21,
      unit = "cm",
      scale = 1,
      dpi = 240,
      quality = 90
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot_nl_test", ".jpg"),
      plot = subgroup_plot_nl_comb,
      device = "jpeg",
      width = 18,
      height = 21,
      unit = "cm",
      scale = 1,
      dpi = 240,
      quality = 90
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot_test", ".pdf"),
      plot = subgroup_plot_comb,
      device = Cairo::CairoPDF,
      width = 18,
      height = 21,
      unit = "cm",
      scale = 1,
      family = "Arial Unicode MS"
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot_nl_test", ".pdf"),
      plot = subgroup_plot_nl_comb,
      device = Cairo::CairoPDF,
      width = 18,
      height = 21,
      unit = "cm",
      scale = 1,
      family = "Arial Unicode MS"
    )
    
    return(NULL)
  }
)

### Compare predictions from bfz and hctz models on combined thiazide validation cohort --------------------------- ----
mice <- length(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")) > 1
if (num_clusters == 1 || !mice) {
  plan(sequential)
} else {
  plan(multisession, workers = num_clusters)
}
cfw_cate_test_bfz_hctz_comp <- future_pmap(
  .l = list(
    i = seq_along(list.files(path, pattern = "^cfw_full_test_mod_\\d{1,}")),
    cohort = tih_tbi_1_test_cohort_imputed
  ),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", 
      "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
      "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
      "excess_risk_cutoff", "path", "result_path", "vi_models", "continuous", "discrete"
    ),
    seed = seed
  ),
  .f = \(i, cohort)  {
    options("parallelly.availableCores.methods" = "system")
    cfw_cate_test <- map(
      .x = seq_along(
        list.files(
          "data/mice/4mo_bfz_all/1year/", 
          pattern = "^cfw_full_mod_\\d{1,}"
        )
      ),
      .f = \(j) {
        workers <- min(floor(num_threads / num_clusters), length(vi_models) + 1)
        threads_per_worker <- floor(num_threads / num_clusters / workers)
        if(workers == 1) {
          plan(sequential)
        } else {
          plan(multisession, workers = workers)
        }
        out <- furrr::future_map(
          .x = c(0, vi_models),
          .options = furrr_options(
            packages = c(
              "dplyr", "stringr", "tibble", "grf", "glue", "stats", "future", "furrr"
            ),
            globals = c("i", "j", "path", "num_threads", "num_clusters", "threads_per_worker",
                        "continuous", "discrete", "cohort", "DiscreteCovariatesToOneHot"),
            seed = TRUE
          ),
          .f = \(vi_model) {
            if (vi_model == 0) {
              cf_bfz <- readRDS(paste0("data/mice/4mo_bfz_all/1year/cfw_full_mod_", j, ".rds"))
              cf_hctz <- readRDS(paste0("data/mice/4mo_hctz_all/1year/cfw_full_mod_", j, ".rds"))
            } else {
              cf_bfz <- readRDS(paste0("data/mice/4mo_bfz_all/1year/cfw_vi_mod_", vi_model, "_", j, ".rds"))
              cf_hctz <- readRDS(paste0("data/mice/4mo_hctz_all/1year/cfw_vi_mod_", vi_model, "_", j, ".rds"))
            }
            continuous_cf <- expr(matches("^X_(0[12789]|[124][0-9]|3[^2]|50|69|7[0-9]|8[0-2])$"))
            discrete_cf <- expr(matches("^X_(0[45]|32|5[1-4])$"))
            predictions_bfz <- predict(
              cf_bfz, 
              newdata = cohort |>
                select({{ continuous_cf }}) |>
                bind_cols(
                  cohort |>
                    select({{ discrete_cf }}) |>
                    DiscreteCovariatesToOneHot()
                ) |>
                select(all_of(names(cf_bfz$X.orig))),
              estimate.variance = TRUE, 
              num.threads = threads_per_worker
            )
            predictions_hctz <- predict(
              cf_hctz, 
              newdata = cohort |>
                select({{ continuous_cf }}) |>
                bind_cols(
                  cohort |>
                    select({{ discrete_cf }}) |>
                    DiscreteCovariatesToOneHot()
                ) |>
                select(all_of(names(cf_hctz$X.orig))),
              estimate.variance = TRUE, 
              num.threads = threads_per_worker
            )
            return(
              list(
                bfz = predictions_bfz,
                hctz = predictions_hctz
              )
            )
          }
        )
        plan(sequential)
        names(out) <- c(
          "full_mod",
          paste0("vi_mod_", vi_models)
        )
        return(out)
      }
    )
    
    # pool CATE estimates from each development model
    plan(multisession, workers = min(num_threads / num_clusters, length(vi_models) + 1))
    cfw_cate_test_pooled <- cfw_cate_test |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          model |>
            list_transpose(simplify = FALSE) |>
            map(
              \(drug) {
                drug |>
                  map(\(x) list_transpose(as.list(x), simplify = TRUE)) |>
                  list_transpose(simplify = FALSE) |>
                  map(\(x) list_transpose(as.list(x), simplify = TRUE)) |>
                  (\(x) {
                    n <- length(x)
                    map(
                      x,
                      \(x, n) {
                        if (is.list(x) && length(x[[1]]) > 1) {
                          pool <- pool.scalar(Q = x$predictions, U = x$variance.estimates, n = n, k = 1, rule = "rubin1987")
                          out <- tibble(
                            predictions = pool$qbar,
                            variance_estimates = pool$t
                          )
                        } else {
                          pool <- NULL
                          out <- tibble(
                            predictions = x$predictions,
                            variance_estimates = x$variance.estimates
                          )
                        }
                        return(list(table = out, pool_results = pool))
                      }, 
                      n = n
                    )
                  })() |>
                  list_transpose()
              }
            )
        }
      )
    plan(sequential)
    
    saveRDS(
      cfw_cate_test_pooled, 
      paste0(
        "results/data/mice/4mo_bfz_hctz_comp/1year/cfw_cate_test_bfz_hctz_comp_",
        i,
        ".rds"
      )
    )
    
    return(map(cfw_cate_test_pooled, \(x) map(x, \(y) y$table)))
  }
)

# aggregate cate estimates
cfw_cate_test_bfz_hctz_comp_agg <- cfw_cate_test_bfz_hctz_comp |>
  list_transpose(simplify = FALSE) |>
  future_map(
    .options = furrr_options(
      packages = c("mice", "dplyr", "purrr")
    ),
    .f = \(model) {
      model |>
        list_transpose(simplify = FALSE) |>
        map(
          \(drug) {
            drug |>
              map(\(x) list_transpose(as.list(x), simplify = TRUE)) |>
              list_transpose(simplify = FALSE) |>
              map(\(x) list_transpose(as.list(x), simplify = TRUE)) |>
              (\(x) {
                n <- length(x)
                map(
                  x,
                  \(x, n) {
                    if (is.list(x) && length(x[[1]]) > 1) {
                      pool <- pool.scalar(Q = x$predictions, U = x$variance_estimates, n = n, k = 1, rule = "rubin1987")
                      out <- tibble(
                        estimate = pool$qbar,
                        std_err = sqrt(pool$t),
                        lower = estimate - qnorm(0.975) * std_err,
                        upper = estimate + qnorm(0.975) * std_err
                      )
                    } else {
                      pool <- NULL
                      out <- tibble(
                        estimate = x$predictions,
                        std_err = sqrt(x$variance_estimates),
                        lower = estimate - qnorm(0.975) * std_err,
                        upper = estimate + qnorm(0.975) * std_err
                      )
                    }
                    return(list(table = out, pool_results = pool))
                  }, 
                  n = 1
                )
              })() |>
              list_transpose()
          }
        )
    }
  )
plan(sequential)

saveRDS(
  cfw_cate_test_bfz_hctz_comp_agg, 
  paste0("results/data/mice/4mo_bfz_hctz_comp/1year/cfw_cate_test_bfz_hctz_comp_agg.rds"),
)

# scatter plots with bfz vs. hctz predictions
imap(
  cfw_cate_test_bfz_hctz_comp_agg,
  \(mod_data, nm) {
    bfz_cutoffs <- as.numeric(quantile(mod_data$bfz$table$estimate, excess_risk_quantiles))
    hctz_cutoffs <- as.numeric(quantile(mod_data$hctz$table$estimate, excess_risk_quantiles))
    plot_data <- mod_data$bfz$table |>
      select("bfz" = estimate) |>
      bind_cols(
        select(mod_data$hctz$table, "hctz" = estimate)
      )
    plot <- plot_data |>
      ggplot(aes(x = bfz, y = hctz)) +
      geom_hex(bins = 100) +
      geom_vline(xintercept = bfz_cutoffs, linetype = 2) +
      geom_hline(yintercept = hctz_cutoffs, linetype = 2) +
      geom_label_repel(
        aes(x = x, y = y, label = label),
        data = tibble(
          x = bfz_cutoffs,
          y = range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))[2] - diff(range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))) * 0.05,
          label = sprintf("%.2f", excess_risk_quantiles)
        ),
        direction = "both",
        force = 2,
        force_pull = 0.25
      ) +
      geom_label_repel(
        aes(x = x, y = y, label = label),
        data = tibble(
          x = range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))[2] - diff(range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))) * 0.05,
          y = hctz_cutoffs,
          label = sprintf("%.2f", excess_risk_quantiles)
        ),
        direction = "both",
        force = 2,
        force_pull = 0.25
      ) +
      geom_label(
        aes(x = x, y = y, label = label),
        data = tibble(
          x = range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))[2] - diff(range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))) * 0.2,
          y = range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))[2] - diff(range(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))) * 0.05,
          label = paste("Correlation:", sprintf("%.2f", cor(mod_data$hctz$table$estimate, mod_data$bfz$table$estimate)))
        )
      ) +
      xlab("Excess risk from BFZ model") + ylab("Excess risk from HCTZ model") +
      coord_cartesian(
        xlim = c(min(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate)), max(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate))),
        ylim = c(min(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate)), max(c(mod_data$bfz$table$estimate, mod_data$hctz$table$estimate)))
      ) +
      theme(
        legend.position = "none"
      )
    ggplot2::ggsave(
      filename = paste0("results/figures/mice/4mo_bfz_hctz_comp/1year/cate_scatter_plot_bfz_hctz_comp_", nm, ".jpg"),
      plot = plot,
      device = "jpeg",
      width = 10,
      height = 7,
      unit = "cm",
      scale = 1.3,
      dpi = 240,
      quality = 90,
      create.dir = TRUE
    )
    return(NULL)
  }
)


### a priori subgroup AIPW estimates ------------------------------------------------------------------------------ ----
plan(multisession, workers = num_clusters)
# prepare plot data tables for each imputation
cfw_subgroup_cate_apriori <- future_pmap(
  .l = list(
    i = seq_along(
      list.files(
        "data/mice/4mo_thiazide_all/1year/", 
        pattern = "^cfw_full_test_mod_\\d{1,}"
      )
    )
  ),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
      "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
      "callr", "fuzzyjoin"
    ),
    globals = c(
      "DiscreteCovariatesToOneHot", "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results",
      "names_index_comb", "mice", "excess_risk_cutoff", "path", "result_path", "vi_models", "continuous", "discrete", 
      "plot_model", "plot_covs", "cutpoints", "cut_labels", "discrete_labels", "CausalForestCATEAllTable", "ATEAll", 
      "CausalForestATEAllSubgroupTable"
    ),
    seed = TRUE
  ),
  .f = \(i)  {
    # read in models with combined thiazides
    cfw_mod <- readRDS(paste0("data/mice/4mo_thiazide_all/1year/cfw_full_mod_", i, ".rds"))
    cfw_cate <- readRDS(paste0("results/data/mice/4mo_thiazide_all/1year/cfw_cate_", i, ".rds"))$full_mod
    cfw_test_mod <- readRDS(paste0("data/mice/4mo_thiazide_all/1year/cfw_full_test_mod_", i, ".rds"))
    cfw_cate_test <- readRDS(paste0("results/data/mice/4mo_thiazide_all/1year/cfw_cate_test_", i, ".rds"))$full_mod$table
    cate_tables <- pmap(
      .l = list(
        cov_list = list(
          "Age" = list(
            "Age" = list(
              "40-55 years" = c("-Inf", "55"),
              "55-70 years" = c("55", "70"),
              "\u2265 70 years" = c("70", "Inf")
            )
          ),
          "House hold income" = list(
            "House hold income" = list(
              "Q1" = "House hold income - Q1",
              "Q2" = "House hold income - Q2",
              "Q3" = "House hold income - Q3",
              "Q4" = "House hold income - Q4",
              "Q5" = "House hold income - Q5"
            )
          ),
          "Heart failure" = list(
            "Heart failure" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Ischemic heart disease" = list(
            "Ischemic heart disease" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Cerebrovascular disease" = list(
            "Cerebrovascular disease" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Any malignancy" = list(
            "Any malignancy" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Malignancy associated with hyponatraemia" = list(
            "Malignancy associated with hyponatraemia" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Other malignancy" = list(
            "Other malignancy" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Liver disease and peritonitis" = list(
            "Liver disease and peritonitis" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Pancreatitis" = list(
            "Pancreatitis" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Liver disease, perito- and pancreatitis" = list(
            "Liver disease and peritonitis" = list(
              "yes" = c("Liver disease and peritonitis", "Pancreatitis"),
              "no" = 0
            ),
            "Pancreatitis" = list(
              "yes" = c(0, 1),
              "no" = 0
            )
          ),
          "Chronic obstructive pulmonary disease" = list(
            "Chronic obstructive pulmonary disease" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Agents for obstructive pulmonary disease" = list(
            "Agents for obstructive pulmonary disease" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          # COPD registered either by use of agent for COPD or by diagnosis code - yes
          # COPD not registers by either agent or diagnosis - no
          "Agents for and chronic obstructive pulmonary disease" = list(
            "Chronic obstructive pulmonary disease" = list(
              "yes" = c("Chronic obstructive pulmonary disease", "Agents for obstructive pulmonary disease"),
              "no" = 0
            ),
            "Agents for obstructive pulmonary disease" = list(
              "yes" = c(0, 1),
              "no" = 0
            )
          ),
          "Diabetes" = list(
            "Diabetes" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Antidiabetics" = list(
            "Antidiabetic (not insulin)" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Insulin" = list(
            "Insulin" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Diabetes including insulin and antidiabetics" = list(
            "Diabetes" = list(
              "yes" = c("Diabetes", "Antidiabetic (not insulin)", "Insulin"),
              "no" = 0
            ),
            "Antidiabetic (not insulin)" = list(
              "yes" = c(0, 1),
              "no" = 0
            ),
            "Insulin" = list(
              "yes" = c(0, 1),
              "no" = 0
            )
          ),
          "Dehydration" = list(
            "Dehydration" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Frail general health" = list(
            "Frail general health" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Physical impairment" = list(
            "Physical impairment" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Mental impairment" = list(
            "Mental impairment" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Rehabilitation contacts" = list(
            "Rehabilitation contacts" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Podiatric contacts" = list(
            "Podiatric contacts" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Frailty conditions" = list(
            "Frail general health" = list(
              "yes" = c("Frail general health", "Physical impairment", "Mental impairment", "Rehabilitation contacts", "Podiatric contacts"),
              "no" = 0
            ),
            "Physical impairment" = list(
              "yes" = c(0, 1),
              "no" = 0
            ),
            "Mental impairment" = list(
              "yes" = c(0, 1),
              "no" = 0
            ),
            "Rehabilitation contacts" = list(
              "yes" = c(0, 1),
              "no" = 0
            ),
            "Podiatric contacts" = list(
              "yes" = c(0, 1),
              "no" = 0
            )
          ),
          "Alcohol abuse" = list(
            "Alcohol abuse" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "History of hyponatraemia" = list(
            "History of hyponatraemia" = list(
              "never" = "History of hyponatraemia - never",
              "<4 months" = "History of hyponatraemia - <4 months",
              "\u22654 months" = "History of hyponatraemia - \u22654 months"
            )
          ),
          "Opioids" = list(
            "Opioids" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Antidepressives" = list(
            "Antidepressives" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          "Antipsychotics" = list(
            "Antipsychotics" = list(
              "yes" = 1,
              "no" = 0
            )
          ),
          # c("< 135 mmol/L", "135-138 mmol/L", "138-141 mmol/L", "\u2265 141 mmol/L")
          "Sodium" = list(
            "Sodium" = list(
              "< 135 mmol/L" = c("-Inf", "135"),
              "135-138 mmol/L" = c("135", "138"),
              "138-141 mmol/L" = c("138", "141"),
              "\u2265 141 mmol/L" = c("141", "Inf")
            )
          ),
          "eGFR" = list(
            "eGFR" = list(
              "30-45 mL/min/1.73 m\u00B2" = c("30", "45"),
              "45-60 mL/min/1.73 m\u00B2" = c("45", "60"),
              "60-90 mL/min/1.73 m\u00B2" = c("60", "90"),
              "\u2265 90 mL/min/1.73 m\u00B2" = c("90", "Inf")
            )
          ),
          "Potassium" = list(
            "Potassium" = list(
              "< 3.5 mmol/L" = c("-Inf", "3.5"),
              "3.5-4.0 mmol/L" = c("3.5", "4"),
              "4.0-4.5 mmol/L" = c("4", "4.5"),
              "\u2265 4.5 mmol/L" = c("4.5", "Inf")
            )
          )
        )
      ),
      .f = \(cov_list) {
        CausalForestCATEAllTable(
          X.orig = bind_rows(
            cfw_mod$X.orig,
            cfw_test_mod$X.orig
          ) |>
            as_tibble() |>
            structure(
              names = slice(names_index_comb, (\(tbl) map_int(names(cfw_mod$X.orig), \(nm) which(tbl$short == nm)))(tbl = names_index_comb))$full
            ),
          Y.orig = c(cfw_mod$Y.orig, cfw_test_mod$Y.orig),
          Y.hat = c(cfw_mod$Y.hat, cfw_test_mod$Y.hat),
          W.orig = c(cfw_mod$W.orig, cfw_test_mod$W.orig),
          W.hat = c(cfw_mod$W.hat, cfw_test_mod$W.hat),
          tau.hat.pointwise = c(pull(filter(cfw_cate, observed), predictions), cfw_cate_test$predictions),
          sample.weights = c(cfw_mod$sample.weights, cfw_test_mod$sample.weights),
          cov_list = cov_list
        )
      }
    )
    return(cate_tables)
  }
)
plan(sequential)
saveRDS(
  cfw_subgroup_cate_apriori, 
  "results/data/mice/4mo_thiazide_all/1year/cfw_subgroup_cate_apriori.rds"
)

# Aggregate results
cfw_subgroup_cate_apriori_agg <- cfw_subgroup_cate_apriori |>
  list_transpose(simplify = FALSE) |>
  future_map(
    .options = furrr_options(
      packages = c("mice", "dplyr", "purrr")
    ),
    .f = \(model) {
      subgroup <- names(model[[1]])[1]
      if (length(model) > 1) {
        out <- model |>
          map(\(df) df |> arrange(subgroup, .locale = "en")) |>
          list_transpose(simplify = FALSE) |>
          (\(list) {
            list <- list[c(subgroup, "estimate", "std_err", "n")]
            i <- 1
            repeat {
              var_name <- c()
              for (j in seq_along(list[[subgroup]])) {
                var_name[j] <- list[[subgroup]][[j]][i]
              }
              if (length(unique(var_name)) > 1) {
                name_to_append <- sort(var_name)[1]
              } else if (!is.na(unique(var_name))){
                i <- i + 1
                next
              } else {
                break
              }
              for (j in seq_along(list[[subgroup]])) {
                if (
                  !is.na(list[[subgroup]][[j]][i]) && 
                  list[[subgroup]][[j]][i] == name_to_append
                ) {
                  n <- list$n[[j]][i]
                  break
                }
              }
              i <- i + 1
            }
            return(list)
          })() |>
          map(
            \(x) {
              out <- list()
              for(i in seq_along(x[[1]])) {
                out[[i]] <- map_vec(x, \(y) y[i])
              }
              return(out)
            }
          ) |>
          list_transpose(simplify = FALSE) |>
          map(
            \(x) {
              pool <- pool.scalar(Q = x$estimate, U = (x$std_err)^2, n = floor(mean(x$n)), k = 1, rule = "rubin1987")
              out <- tibble(
                name = x[[subgroup]][1],
                n = floor(mean(x$n)),
                estimate = pool$qbar,
                std_err = sqrt(pool$t),
                `95% CI - lower` = estimate - qnorm(0.975) * std_err,
                `95% CI - upper` = estimate + qnorm(0.975) * std_err
              )
              return(list(table = out, pool = pool))
            }
          ) |>
          list_transpose()
        names(out$pool) <- out$table$name
      } else {
        out <- list(
          pool = NULL,
          table = model[[1]] |>
            rename(
              "name" = subgroup
            )
        )
      }
      return(out)
    }
  )

writexl::write_xlsx( 
  map(cfw_subgroup_cate_apriori_agg, \(x) x$table),
  path = "results/tables/mice/4mo_thiazide_all/1year/cfw_subgroup_cate_apriori.xlsx"
)
