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
library(ragg)
library(cowplot)

## source helper functions ..................................................................................... .. ----
source("functions/CATEPlot.R", encoding = 'UTF-8')
source("functions/CATEPlots.R", encoding = 'UTF-8')
source("functions/CausalForestATESubgroupTable.R", encoding = 'UTF-8')
source("functions/CausalForestCATETable.R", encoding = 'UTF-8')
source("functions/CausalForestDynamicSubgroups.R", encoding = 'UTF-8')
source("functions/CForBenefit.R", encoding = 'UTF-8')
source("functions/CovariateBalance.R", encoding = 'UTF-8')
source("functions/DiscreteCovariatesToOneHot.R", encoding = 'UTF-8')
source("functions/futuremice.R", encoding = 'UTF-8')
source("functions/GRFAnalysisWrapper.R", encoding = 'UTF-8')
source("functions/MakeTIHCohort.R", encoding = 'UTF-8')
source("functions/MapCovariateBalance.R", encoding = 'UTF-8')
source("functions/mbcal.R", encoding = 'UTF-8')
source("functions/mem_avail.R", encoding = 'UTF-8')
source("functions/OneHotToFactor.R", encoding = 'UTF-8')
source("functions/OverlapPlot.R", encoding = 'UTF-8')
source("functions/predtools_calibration_plot.R", encoding = 'UTF-8')
source("functions/PrepareTIHCovariates.R", encoding = 'UTF-8')
source("functions/RATETest.R", encoding = 'UTF-8')
source("functions/RATETestcfwHelper.R", encoding = 'UTF-8')
source("functions/RATETestcfwHelperParallel.R", encoding = 'UTF-8')
source("functions/SubgroupATETable.R", encoding = 'UTF-8')
source("functions/TestCalibrationData.R", encoding = 'UTF-8')
source("functions/TIH_DynamicSubgroups.R", encoding = 'UTF-8')
source("functions/TimeDiff.R", encoding = 'UTF-8')
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
num_fold_par <- 2L
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
dynamic_subgroups_temp_path <- "PATH" 
## list of seeds used for reproducibility ......................................................................... ----
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

## load thiazide cohort ........................................................................................... ----
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
  full = c(# event indicator: (0: censored, 1: hyponatraemia)
    "Follow-up for <130mM event", "hyponatraemia (<130mM) event", 
    "Follow-up for <125mM event", "hyponatraemia (<125mM) event", 
    "hyponatraemia (<130 mM) event inside 90 days", "hyponatraemia (<130 mM) event or >90 days",
    "hyponatraemia (<125 mM) event inside 90 days", "hyponatraemia (<125 mM) event or >90 days",
    "hyponatraemia (<130 mM) event inside 120 days", "hyponatraemia (<130 mM) event or >120 days",
    "hyponatraemia (<125 mM) event inside 120 days", "hyponatraemia (<125 mM) event or >120 days",
    "BFZ initiation", "HCTZ initiation", "thiazide drug initiation", 
    "Age", "Sex", "Region of residence", 
    "House hold income", "Years of education",  
    "Calendar year of inclusion", "Heart failure", "Essential hypertension", 
    "Secondary hypertension", "Ischemic heart disease", "Cerebrovascular disease", 
    "Other CNS disorders", "Arrhythmias", "Any malignancy", 
    "Malignancy associated with hyponatraemia",  "Other malignancy", 
    "Renal disorders", "Liver disease and peritonitis", "Pancreatitis", 
    "Chronic obstructive pulmonary disease", "Diabetes", "Dehydration", 
    "Frail general health", "Physical impairment", "Mental impairment", 
    "Rehabilitation contacts", "Podiatric contacts", "Alcohol abuse", 
    "drug abuse", "HIV", "Anorexia and primary polydipsia",  
    "History of hyponatraemia", "Drugs for peptic ulcer and gastroesophageal reflux", 
    "Drugs used for obstipation or diarrhea", "Antidiabetic (not insulin)", 
    "Insulin", "Oral anticoagulants", "Aspirin", "ADPi", "Nitrates",
    "Beta-blockers", "Lipid lowering drugs", "Desmopressin", "Eltroxin", 
    "NSAIDs", "Antiepileptics", "Opioids", "Antidepressives", "Antipsychotics",
    "Agents for obstructive pulmonary disease",
    "No. of different prescription drugs used",
    "Days of hospitalization in the past year prior to index data",
    "No. of outpatient contacts in the past year prior to index data",
    "No. of primary care contacts in the last 2 years", "cSodium", "ceGRF",
    "cPotassium", "cHemoglobin", "cALAT", "cAlkaline phosphatase", 
    "cAlbumin", "cLactate dehydrogenase", "cCarbamide", "cThrombocytes",
    "cLeukocytes", "cC-reactive protein", "cCholesterol", "cLDL-C", "Sodium", 
    "eGFR", "Potassium", "Hemoglobin", "ALAT", "Alkaline phosphatase", 
    "Albumin", "Lactate dehydrogenase", "Carbamide", "Thrombocytes",
    "Leukocytes", "C-reactive protein", "Cholesterol", "LDL-C"
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
# indication codes 54 and 57 are used for hypertension:
# 54: for blood pressure (for blodtrykket)              
# 57: against hypertension (mod forhøjet blodtryk)
# cohort with all covariates included. Selection of different subsets of the covariates happens later.
tih_cohort_all <- study_cohort |>
  filter(grepl("5(4|7)$", indo)) |>
  MakeTIHCohort(h = horizons)

## Split into training and test data .............................................................................. ----
# split temporally. Training data 2014-2018 and test data 2019-2020
tih_cohort <- tih_cohort_all |>
  filter(X_06 %in% c("2014", "2015", "2016", "2017", "2018"))
tih_cohort$X_06 <- droplevels(tih_cohort$X_06)
attr(tih_cohort$X_06, "label") <- "Calendar year of inclusion"

# filter on complete cases and four least missing blood test results
tih_cohort_complete_case <- 
  tih_cohort |>
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
  filter(!is.na(rowSums(select(tih_cohort, "X_69", "X_70", "X_71", "X_73")))) |>
  replace_na(
    list(
      X_04 = "Q4",
      X_05 = "10-12"
    )
  )

## Construct a table 1 from development and validation cohorts .................................................... ----
table_one_vars <- tih_cohort |>
  select(starts_with("X")) |>
  select(!ends_with("tbi")) |>
  select(matches("[0-4]\\d|5[0-4]|69|7\\d|8[0-2]")) |>
  map_chr(\(x) attr(x, "label"))
table_one_vars_cc <- tih_cohort_complete_case |>
  select(starts_with("X")) |>
  select(!ends_with("tbi")) |>
  select(matches("[0-4]\\d|5[0-4]|69|7\\d|8[0-2]")) |>
  map_chr(\(x) attr(x, "label"))


table_one_dev <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_dev_cc <- CreateTableOne(
  vars = as.character(table_one_vars_cc), 
  strata = "W", 
  data = tih_cohort_complete_case |>
    rename(
      !!!setNames(as.list(names(table_one_vars_cc)), table_one_vars_cc)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_dev_bfz <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort |>
    filter(!is.na(W_bvc)) |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)
table_one_dev_hctz <- CreateTableOne(
  vars = as.character(table_one_vars), 
  strata = "W", 
  data = tih_cohort |>
    filter(!is.na(W_hvr)) |>
    rename(
      !!!setNames(as.list(names(table_one_vars)), table_one_vars)
    ) |>
    mutate(Sex = factor(Sex, levels = c(0, 1), labels = c("female", "male")))
)

table_one_dev_exp <- print(
  table_one_dev,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_dev_cc_exp <- print(
  table_one_dev_cc,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_dev_bfz_exp <- print(
  table_one_dev_bfz,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)
table_one_dev_hctz_exp <- print(
  table_one_dev_hctz,
  nonnormal = TRUE,
  missing = TRUE,
  printToggle = FALSE
)

for (col in seq_len(ncol(table_one_dev_exp))) {
  for (row in seq_len(nrow(table_one_dev_exp))) {
    cell <- table_one_dev_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_dev_exp[row, col] <- ifelse(as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, "<5 ( 0.0)", cell)
    }
  }
}
for (col in seq_len(ncol(table_one_dev_cc_exp))) {
  for (row in seq_len(nrow(table_one_dev_cc_exp))) {
    cell <- table_one_dev_cc_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_dev_cc_exp[row, col] <- ifelse(as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, "<5 ( 0.0)", cell)
    }
  }
}
for (col in seq_len(ncol(table_one_dev_bfz_exp))) {
  for (row in seq_len(nrow(table_one_dev_bfz_exp))) {
    cell <- table_one_dev_bfz_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_dev_bfz_exp[row, col] <- ifelse(as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, "<5 ( 0.0)", cell)
    }
  }
}
for (col in seq_len(ncol(table_one_dev_hctz_exp))) {
  for (row in seq_len(nrow(table_one_dev_hctz_exp))) {
    cell <- table_one_dev_hctz_exp[row, col] 
    if (cell != "" && 
        !str_detect(cell, "^\\s{1,}$") && 
        !str_detect(cell, "^(<|>|[a-z]|[A-Z])") &&
        str_detect(cell, "^\\s{0,}\\d{1,}\\s{1}\\(") &&
        str_extract(cell, "^\\s{0,}\\d{1,}") == str_extract(cell, "^\\s{0,}\\d{1,}\\.?")) {
      table_one_dev_hctz_exp[row, col] <- ifelse(as.numeric(str_extract(cell, "^\\s{0,}\\d{1,}")) < 5, "<5 ( 0.0)", cell)
    }
  }
}

write.csv2(table_one_dev_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_dev.csv"))
write.csv2(table_one_dev_cc_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_dev_cc.csv"))
write.csv2(table_one_dev_bfz_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_dev_bfz.csv"))
write.csv2(table_one_dev_hctz_exp, file = paste0(str_sub(tables_path, 1, -6), "table1_dev_hctz.csv"))

## proportion of thiazide vs. non-thiazide antihypertensive prescriptions by year ................................. ---- 
tih_cohort_all |>
  summarise(thiazide_prop = sum(W) / n(), .by = "X_06") |>
  arrange(X_06)

## count of censoring by thiazide group ........................................................................... ----
censtab <- left_join(
  tih_cohort,
  study_cohort |> 
    select("pnr", "dead_exit", "emigrated_exit", "hyp_stop_add_on_exit", "reg_midt_censor_date_exit"), 
  by = "pnr"
)
censtab |> filter(Y_130 < 120) |> count(W, Y_130_120)
censtab |> filter(is.na(Y_130_120) & Y_130 < 120) |> count(W, dead_exit)
censtab |> filter(is.na(Y_130_120) & Y_130 < 120) |> count(W, emigrated_exit)
censtab |> filter(is.na(Y_130_120) & Y_130 < 120) |> count(W, hyp_stop_add_on_exit)
censtab |> filter(is.na(Y_130_120) & Y_130 < 120) |> count(W, reg_midt_censor_date_exit)

## correlations between covariates and study exposure/outcome ..................................................... ----
# continuous covariates (pearson correlation)
cor_W_continuous <- cor(
  select(tih_cohort, c(17, seq(71, 110, 3))), select(tih_cohort, "W"), 
  use = "complete.obs", method = "pearson"
)
cor_Y_130_120_continuous <- cor(
  select(tih_cohort, c(17, seq(71, 110, 3))), select(tih_cohort, "Y_130_120"), 
  use = "complete.obs", method = "pearson"
)
cor_Y_125_90_continuous <- cor(
  select(tih_cohort, c(17, seq(71, 110, 3))), select(tih_cohort, "Y_125_90"), 
  use = "complete.obs", method = "pearson"
)
# discrete (rank correlation)
cor_W_discrete <- map_dbl(
  select(tih_cohort, c(18, 23:47, 49:66)),
  \(x) pcaPP::cor.fk(x, pull(tih_cohort, "W"))
) |> as.matrix()
cor_Y_130_120_discrete <- map_dbl(
  select(filter(tih_cohort, !is.na(tih_cohort$Y_130_120)), c(18, 23:47, 49:66)),
  \(x) pcaPP::cor.fk(x, pull(filter(tih_cohort, !is.na(tih_cohort$Y_130_120)), "Y_130_120"))
) |> as.matrix()
cor_Y_125_90_discrete <- map_dbl(
  select(filter(tih_cohort, !is.na(tih_cohort$Y_125_90)), c(18, 23:47, 49:66)),
  \(x) pcaPP::cor.fk(x, pull(filter(tih_cohort, !is.na(tih_cohort$Y_125_90)), "Y_125_90"))
) |> as.matrix()

# Imputation of missing blood test measurements ------------------------------------------------------------------- ----
## Imputation using multiple imputation ........................................................................... ----
# use multiple imputation to impute m complete data matrices, then train a causal forest model to each, 
# and finally aggregating the CATE estimates from each. 
if (impute_data) {
  plan(multisession, workers = 3)
  future_pmap(
    .l = list(
      cohort = list(
        tih_cohort, 
        tih_cohort |> filter(!is.na(W_bvc)),
        tih_cohort |> filter(!is.na(W_hvr))
      ),
      exposure = c("thiazide", "bfz", "hctz"),
      tbi = rep(years(1), 3),
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
            # These arguments are only relevant if polr fails and imputation falls back on nnet::multinom()
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
            donors = 5L, # size of donor pool. Use default 5L.
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
        paste0("data/cohorts/", exposure, "_tbi_", tbi@year, "_cohort_imputed_mice.rds")
      )
    } 
  )
  plan(sequential)
}

# save plot of trace lines from MICE
jpeg(
  filename = paste0(figures_path, "diagnostics/trace_plot.jpg"),
  width = 15,
  height = 21,
  units = "cm",
  res = 240,
  quality = 90
)
plot(tih_mice_imp, layout = c(4, 8))
dev.off()

## read in imputed data from multiple imputation .................................................................. ----
tih_mice_imp <- readRDS("data/cohorts/thiazide_tbi_1_cohort_imputed_mice.rds")
tih_tbi_1_cohort_imputed_mice <- list()
for (i in seq_len(tih_mice_imp$m)) {
  tih_tbi_1_cohort_imputed_mice[[i]] <- complete(tih_mice_imp, i)
}
bfz_mice_imp <- readRDS("data/cohorts/bfz_tbi_1_cohort_imputed_mice.rds")
bfz_tbi_1_cohort_imputed_mice <- list()
for (i in seq_len(bfz_mice_imp$m)) {
  bfz_tbi_1_cohort_imputed_mice[[i]] <- complete(bfz_mice_imp, i)
}
hctz_mice_imp <- readRDS("data/cohorts/hctz_tbi_1_cohort_imputed_mice.rds")
hctz_tbi_1_cohort_imputed_mice <- list()
for (i in seq_len(hctz_mice_imp$m)) {
  hctz_tbi_1_cohort_imputed_mice[[i]] <- complete(hctz_mice_imp, i)
}

tih_tbi_1_cohort_imputed <- tih_tbi_1_cohort_imputed_mice
bfz_tbi_1_cohort_imputed <- bfz_tbi_1_cohort_imputed_mice
hctz_tbi_1_cohort_imputed <- hctz_tbi_1_cohort_imputed_mice

# filter on complete cases and four least missing blood test results
tih_cohort_cc_comp <-  map(
  tih_tbi_1_cohort_imputed,
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

# train causal forest with censoring weights ---------------------------------------------------------------------- ----
# Training a causal forest model requires estimates of sample weights, expected outcomes and expected propensities to 
# treat. These estimates will be obtained from survival and regression forest models, using automatic parameter tuning 
# and all covariates. We will use a model using continuous blood test measurements with missing values imputed. 
# Parameter tuning will be used with all covariates, and if values better than defaults are found, we will tune all 
# following models as well. 

# To allow multiple imputation the future_pmap function will be used to allow the analysis pipeline to be run 
# multiple times on different datasets.
if(train_cf) {
  # set seed to ensure consistent indices for dynamic subgroups
  set.seed(seeds[4])
  plan(sequential)
  # map over different analyses (e.g. thiazide combined, BFZ, and HCTZ)
  future_pmap(
    .l = list(
      name = list(
        # tih = all thiazide drugs
        # bfz = only bfz's
        # hctz = only hctz's
        # 4mo = events within 4 months of index date
        # abt(n) = all blood tests n years back
        # nbt = no blood tests
        # kbt(n) kidney related blood tests n years back 
        # rds = removed Denmark specific covariates
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
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
        list(tih_cohort), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        bfz_tbi_1_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        bfz_tbi_1_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        hctz_tbi_1_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        hctz_tbi_1_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        list(tih_cohort_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
        tih_cohort_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
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
        paste0(training_path, "4mo_thiazide_complete_case/1year/"),
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
                   "furrr", "future.callr", "DescTools"),
      globals = c("num_threads", "num_clusters", "names_index", "names_index_comb", "GRFAnalysisWrapper", 
                  "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel", "VariableImportanceWrapper", "tolerence",
                  "save_sampleweight_forest", "save_exposure_forest", "save_outcome_forest", "n_folds", "n_rankings",
                  "dynamic_subgroups", "vcovHC", "vi_models", "TIH_DynamicSubgroups", "num_fold_par",
                  "dynamic_subgroups_temp_path", "cluster_name"),
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
          packages = c("grf", "rlang", "dplyr", "stringr", "tidyr", "purrr", "glue", 
                       "stats", "writexl", "future", "furrr", "future.callr", "DescTools"),
          globals = c("path", "num_threads", "num_clusters", "names_index", "names_index_comb", "horizon", 
                      "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel",
                      "VariableImportanceWrapper", "save_sampleweight_forest", "save_exposure_forest", 
                      "save_outcome_forest", "continuous", "discrete", "tune_parameters", "indices_dyn", 
                      "n_folds", "n_rankings", "name", "vcovHC", "vi_models", "TIH_DynamicSubgroups",
                      "dynamic_subgroups", "num_fold_par", "tolerence", "dynamic_subgroups_temp_path",
                      "thiazide_limit", "cluster_name"),
          seed = seed
        ),
        .f = \(cohort, data_index) {
          # create seeds for grf forest generation
          rng_seed <- sample(-999999999:999999999, 4)
          
          if (file.exists(paste0(path, "cfw_full_mod_", data_index, ".rds"))) {
            # loading existing trained model
            cfw_full <- readRDS(paste0(path, "cfw_full_mod_", data_index, ".rds"))
          } else {
            ### Survival forest model to predict sample weights ................................................... ----
            # Note: automatic tuning not implemented for survival forests
            cfw_sample_weight <-
              (function(data) {
                # set up variables based on the covariates used
                continuous <- continuous$sample_weight
                discrete <- discrete$sample_weight
                # calculate sample weights based on censoring process.
                # The censoring probabilities are P(C > T|X, W), therefore the variables used in the survival forest
                # include both X and W.
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
                # we use the survival probability of the censoring process at the "outcome time", which is the observed 
                # time when considering a hyponatraemia event and at the horizon when observing no hyponatraemia up to 
                # the horizon. If one makes it past the horizon, we know they survived up to this time.
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
              saveRDS(cfw_sample_weight, paste0(path, "cfw_sample_weight_", data_index, ".rds"))
            }
            
            sample_weights <- cfw_sample_weight$sample_weights
            
            # remove survival forest from memory
            rm("cfw_sample_weight")
            gc()
            
            # create data without censored observations (complete-case)
            data_obs <- filter(cohort, !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
            
            ### Regression forest models to predict expected outcome and propensity ................................. ----
            cfw_exp_out <- (function(data, sample_weights, tolerence) {
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
                tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
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
                tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
                seed = rng_seed[3],
                num.threads = floor(num_threads / num_clusters)
              )
              # expected outcome and treatment propensity in non-censored
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
                out$exp_out$W_hat[which(out$exp_out$W_hat) < tolerence] <- tolerence
                out$exp_out$W_hat[which(out$exp_out$W_hat) > 1 - tolerence] <- 1 - tolerence
                warning("Some units have propensity scores close to 0 or 1.",
                        "Scores are set to be at least", tolerence, " away from 0 and 1.")
              }
              return(out)
            })(data = data_obs, sample_weights = sample_weights, tolerence = tolerence)
            
            if (save_exposure_forest) {
              saveRDS(cfw_exp_out$exposure_forest, paste0(path, "cfw_exposure_", data_index, ".rds"))
            }
            if (save_outcome_forest) {
              saveRDS(cfw_exp_out$outcome_forest, paste0(path, "cfw_outcome_", data_index, ".rds"))
            }
            
            # save table with expected outcome and propensity scores in the full cohort
            saveRDS(cfw_exp_out$exp_out, paste0(path, "cfw_exp_out_", data_index, ".rds"))
            
            # extract expected outcome and propensity scores for those non-censored
            Y_hat <- filter(cfw_exp_out$exp_out, observed)$Y_hat
            W_hat <- filter(cfw_exp_out$exp_out, observed)$W_hat
            
            # remove outcome and exposure forests from memory
            rm("cfw_exp_out")
            gc()
            
            ### Causal forest models with all covariates .......................................................... ----
            cfw_full <-
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
                  tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
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
                data = data_obs,
                sample_weights = sample_weights,
                Y_hat = Y_hat,
                W_hat = W_hat
              )
            
            saveRDS(cfw_full, paste0(path, "cfw_full_mod_", data_index, ".rds"))
          }
          
          # calculate variable importance from model with all covariates
          cfw_full_vi <- VariableImportanceWrapper(cfw_full, names_index)
          saveRDS(cfw_full_vi, paste0(path, "cfw_full_vi_", data_index, ".rds"))
          writexl::write_xlsx(cfw_full_vi, paste0(path, "cfw_full_vi_", data_index, ".xlsx"))
          
          # remove causal forest from memory
          rm("cfw_full")
          
          # combine variable importance from different levels of categorical variables
          cfw_full_vi_comb <- (\(vi, nm_index) {
            for (nm in nm_index$full) {
              if (
                nrow(
                  filter(
                    cfw_full_vi,
                    str_detect(cfw_full_vi$variable_name, str_c("^", nm, " "))
                  )
                ) > 1
              ) {
                vi <- bind_rows(
                  vi,
                  cfw_full_vi |>
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
            vi = cfw_full_vi[, 1:2],
            nm_index = filter(names_index, !grepl("^(X_(5[5-9]|6[0-8])$|[^X])", short))
          )
          writexl::write_xlsx(cfw_full_vi_comb, paste0(path, "cfw_full_vi_comb_", data_index, ".xlsx"))
          
          # Compute tau_hat and rankings within each fold using remaining folds to train a causal forest
          if (dynamic_subgroups) {
            if (!file.exists(paste0(path, "cfw_full_mod_dyn_", data_index, ".rds"))) {
              TIH_DynamicSubgroups(
                temp_path = dynamic_subgroups_temp_path,
                remove_temp = "dynamic_subgroups",
                cf_append = "full",
                model_name = "cfw_full_mod",
                setup_name = name,
                cohort = cohort,
                cluster_name = cluster_name,
                n_folds = n_folds,
                n_rankings = n_rankings,
                continuous = continuous,
                discrete = discrete,
                tune_parameters = tune_parameters,
                data_index = data_index,
                num_threads = num_threads,
                num_clusters = num_clusters,
                num_fold_par = num_fold_par,
                horizon = horizon,
                thiazide_limit = thiazide_limit
              )
            }
          }
          
          return(invisible(NULL))
        }
      )
      
      # aggregate variable importance from each imputed full model:
      if (file.exists(paste0(path, "cfw_full_vi_agg.rds"))) {
        # load existing aggregated variable importance values
        cfw_full_vi_agg <- readRDS(paste0(path, "cfw_full_vi_agg.rds"))
        cfw_full_vi_comb <- readRDS(paste0(path, "cfw_full_vi_comb_agg.rds"))
      } else {
        # load variable importance values from each imputed model
        cfw_full_vi <- list()
        for (i in seq_along(list.files(path, pattern = "^cfw_full_vi_\\d{1,}.rds"))) {
          cfw_full_vi[[i]] <- readRDS(paste0(path, "cfw_full_vi_", i, ".rds"))
        }
        # pool variable importance estimates
        cfw_full_vi_agg <- cfw_full_vi |>
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
        saveRDS(cfw_full_vi_agg, paste0(path, "cfw_full_vi_agg.rds"))
        # combine variable importance from different levels of categorical variables
        cfw_full_vi_comb <- (\(vi, nm_index) {
          for (nm in names_index$full) {
            if (
              nrow(
                filter(
                  cfw_full_vi_agg$table,
                  str_detect(cfw_full_vi_agg$table$variable_name, str_c("^", nm, " "))
                )
              ) > 1
            ) {
              vi <- bind_rows(
                vi,
                cfw_full_vi_agg$table |>
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
          vi = cfw_full_vi_agg[["table"]][, 1:2],
          nm_index = filter(names_index, !grepl("^(X_(5[5-9]|6[0-8])$|[^X])", short))
        )
        saveRDS(cfw_full_vi_comb, paste0(path, "cfw_full_vi_comb_agg.rds"))
        writexl::write_xlsx(cfw_full_vi_comb, paste0(path, "cfw_full_vi_comb_agg.xlsx"))
      }
      
      future_pmap(
        .l = list(
          cohort = cohort,
          data_index = seq_along(cohort)
        ),
        .options = furrr_options(
          packages = c("grf", "rlang", "dplyr", "stringr", "tidyr", "purrr", "glue", 
                       "stats", "writexl", "future", "furrr", "future.callr", "DescTools"),
          globals = c("path", "num_threads", "num_clusters", "names_index", "names_index_comb", "horizon", 
                      "GRFAnalysisWrapper", "DiscreteCovariatesToOneHot", "RATETestcfwHelperParallel",
                      "VariableImportanceWrapper", "save_sampleweight_forest", "save_exposure_forest", 
                      "save_outcome_forest", "continuous", "discrete", "tune_parameters", "indices_dyn", 
                      "n_folds", "n_rankings", "name", "vcovHC", "vi_models", "TIH_DynamicSubgroups",
                      "dynamic_subgroups", "cfw_full_vi_comb", "num_fold_par", "tolerence",
                      "dynamic_subgroups_temp_path", "thiazide_limit", "cluster_name"),
          seed = seed
        ),
        .f = \(cohort, data_index) {
          ### read in data from model with all covariates ......................................................... ----
          # loading existing trained model
          cfw_full <- readRDS(paste0(path, "cfw_full_mod_", data_index, ".rds"))
          sample_weights <- cfw_full$sample.weights
          W_hat <- cfw_full$W.hat
          Y_hat <- cfw_full$Y.hat
          # create data without censored observations (complete-case)
          data_obs <- filter(cohort, !is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon))))
          # set flag for further tuning if successful for all covariates
          tune <- cfw_full$tuning.output$status != "default"
          tune_param <- if (tune)  tune_parameters$cf else "none"
          # remove causal forest
          rm("cfw_full")
          gc()
          
          ### Models with highest variable importance from full_cont model ........................................ ----
          # create seeds for grf forest generation
          rng_seed <- sample(-999999999:999999999, length(vi_models))
          
          # run causal forest analysis for n most important covariates
          for (n in vi_models) {
            (function(data, Y_hat, W_hat, sample_weights, continuous, discrete) {
              if (!file.exists(paste0(path, "cfw_vi_mod_", n, "_", data_index, ".rds"))) {
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
                  tune.num.trees = 200, tune.num.reps = 100, tune.num.draws = 1000,
                  seed = rng_seed[which(vi_models == n)],
                  num.threads = floor(num_threads / num_clusters)
                )
                if("X_03_SjÃ¦lland" %in% names(out$X.orig)) {
                  names(out$X.orig)[
                    which(names(out$X.orig) == "X_03_SjÃ¦lland")
                  ] <- "X_03_Sjælland"
                }
                saveRDS(out, paste0(path, "cfw_vi_mod_", n, "_", data_index, ".rds"))
                
                # remove causal forest from memory
                rm(list = "out")
              }
              
              # Compute tau_hat and rankings within each fold using remaining folds to train a causal forest
              if (dynamic_subgroups) {
                if (!file.exists(paste0(path, "cfw_vi_mod_", n, "_dyn_", data_index, ".rds"))) {
                  TIH_DynamicSubgroups(
                    temp_path = dynamic_subgroups_temp_path,
                    remove_temp = "dynamic_subgroups",
                    cf_append = paste0("vi_", n),
                    model_name = paste0("cfw_vi_mod_", n),
                    setup_name = name,
                    cohort = cohort,
                    cluster_name = cluster_name,
                    n_folds = n_folds,
                    n_rankings = n_rankings,
                    continuous = continuous,
                    discrete = discrete,
                    tune_parameters = list(cf = tune_param),
                    data_index = data_index,
                    num_threads = num_threads,
                    num_clusters = num_clusters,
                    num_fold_par = num_fold_par,
                    horizon = horizon,
                    thiazide_limit = thiazide_limit
                  )
                }
              }
              return(invisible(NULL))
            })(data = data_obs,
               Y_hat = Y_hat,
               W_hat = W_hat,
               sample_weights = sample_weights,
               continuous = continuous,
               discrete = discrete)
          }
          
          # empty return
          return(invisible(NULL))
        }
      )
    }
  )
}

# check forests -------------------------------------------------------------------------------------------  -----------
# each setup is analyzed in turn. For each setup, each imputation is analyzed in turn.
# A flag is created if mice was used to indicate use of multicore processing.

if (compute_results) {
  plan(sequential)
  future_pmap(
    .l = list(
      cohort = list(
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 3 months, outcome cutoff 125 mmol/L, includes all blood test results 1 year back
        list(tih_cohort), # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results excluded
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back
        tih_tbi_1_cohort_imputed, # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        bfz_tbi_1_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        bfz_tbi_1_cohort_imputed, # analysis bfz vs ccb, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        hctz_tbi_1_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, includes all blood test results 1 year back
        hctz_tbi_1_cohort_imputed, # analysis hctz vs ras, outcome within 4 months, outcome cutoff 130 mmol/L, blood test results except sodium, eGRF, and potassium excluded 1 year back, remove Denmark specific covariates
        list(tih_cohort_complete_case),  # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, only complete cases
        tih_cohort_cc_comp         # analysis thiazide combined, outcome within 4 months, outcome cutoff 130 mmol/L, includes 4 least missing variables, with MI for comparison with CC.
      ),
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
        paste0(training_path, "4mo_thiazide_complete_case/1year/"),
        paste0(training_path, "4mo_thiazide_cc_comp/1year/")
      ),
      result_path = list(
        paste0(results_path, "4mo_thiazide_all/1year/"),
        paste0(results_path, "3mo_tiazide_all_severe/1year/"),
        paste0(results_path, "4mo_thiazide_no_lab/"),
        paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
        paste0(results_path, "4mo_thiazide_international/1year/"),
        paste0(results_path, "4mo_bfz_all/1year/"),
        paste0(results_path, "4mo_bfz_international/1year/"),
        paste0(results_path, "4mo_hctz_all/1year/"),
        paste0(results_path, "4mo_hctz_international/1year/"),
        paste0(results_path, "4mo_thiazide_complete_case/1year/"),
        paste0(results_path, "4mo_thiazide_cc_comp/1year/")
      ),
      seed = seeds[16:26],
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
      plot_model = list(
        "vi_mod_4", # thiazide all
        "vi_mod_4", # tiazide all severe
        "vi_mod_4", # thiazide no lab
        "vi_mod_4", # thiazide only fluid lab results
        "vi_mod_4", # thiazide only fluid lab results and no Denmark specific covariates
        "vi_mod_4", # bfz all
        "vi_mod_4", # bfz only fluid lab results and no Denmark specific covariates
        "vi_mod_4", # hctz all
        "vi_mod_4", # hctz only fluid lab results and no Denmark specific covariates
        "vi_mod_4", # thiazide complete case
        "vi_mod_4", # thiazide complete case comparison
      ),
      plot_covs = list(
        c("X_01", "X_69", "X_72", "X_80"), # thiazide all
        c("X_01", "X_69", "X_72", "X_80"), # tiazide all severe
        c("X_01", "X_69", "X_72", "X_80"), # thiazide no lab
        c("X_01", "X_69", "X_72", "X_80"), # thiazide only fluid lab results
        c("X_01", "X_69", "X_72", "X_80"), # thiazide only fluid lab results and no Denmark specific covariates
        c("X_01", "X_69", "X_76", "X_80"), # bfz all
        c("X_01", "X_69", "X_72", "X_80"), # bfz only fluid lab results and no Denmark specific covariates
        c("X_01", "X_69", "X_72", "X_78"), # hctz all
        c("X_01", "X_69", "X_72", "X_80"), # hctz only fluid lab results and no Denmark specific covariates
        c("X_01", "X_05", "X_52", "X_69"), # thiazide complete case
        c("X_01", "X_05", "X_52", "X_69"), # thiazide complete case comparison
      ),
      cutpoints = list(
        list(
          Age = c(55, 70),
          Sodium = c(135, 138, 141),
          Hemoglobin = c(8, 9, 10),
          `C-reactive protein` = c(5, 20)
        ), # thiazide all
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
        ), # thiazide no lab
        list(
          Age = c(55, 70),
          Sodium = c(135, 138, 141),
          Hemoglobin = c(8, 9, 10),
          `C-reactive protein` = c(5, 20)
        ), # thiazide only fluid lab results
        list(
          Age = c(55, 70),
          Sodium = c(135, 138, 141),
          Hemoglobin = c(8, 9, 10),
          `C-reactive protein` = c(5, 20)
        ), # thiazide only fluid lab results and no Denmark specific covariates
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
        ), # thiazide complete case
        list(
          Age = c(55, 70),
          `Years of education` = list("<10", "10-12", c("13-15", ">15")),
          `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
          Sodium = c(135, 138, 141)
        ) # thiazide complete case comparison
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
        ), # thiazide all
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
        ), # thiazide no lab
        list(
          Age = NA,
          Sodium = NA,
          Hemoglobin = NA,
          `C-reactive protein` = NA
        ), # thiazide only fluid lab results
        list(
          Age = NA,
          Sodium = NA,
          Hemoglobin = NA,
          `C-reactive protein` = NA
        ), # thiazide only fluid lab results and no Denmark specific covariates
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
        ), # thiazide complete case
        list(
          Age = NA,
          `Years of education` = list("<10", "10-12", c("13-15", ">15")),
          `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
          Sodium = NA
        ) # thiazide complete case comparison
      )
    ) |> map(\(x) x[setup_index]),
    .options = furrr_options(
      packages = c(
        "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
        "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
        "callr", "fuzzyjoin"
      ),
      globals = c(
        "MapCovariateBalance", "CovariateBalance", "CForBenefit", "RATETestcfwHelperParallel", "vi_models",
        "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
        "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "names_index_comb", "mice", 
        "excess_risk_quantiles", "excess_risk_cutoff", "CausalForestCATETable", "CausalForestATESubgroupTable",
        "horizon", "thiazide_limit"
      ),
      seed = TRUE
    ),
    .f = \(cohort, path, result_path, seed, continuous, discrete, plot_model, plot_covs, cutpoints, cut_labels, discrete_labels) {
      mice <- length(list.files(path, pattern = "^cfw_full_mod_\\d{1,}")) > 1
      if (num_clusters == 1 || !mice) {
        plan(sequential)
      } else {
        plan(multisession, workers = num_clusters)
      }
      future_pmap(
        .l = list(
          i = seq_along(list.files(path, pattern = "^cfw_full_mod_\\d{1,}")),
          cohort = cohort
        ),
        .options = furrr_options(
          packages = c(
            "dplyr", "stringr", "purrr", "tibble", "tidyr", "ggsci", "ggplot2", "data.table", "cowplot", "gridExtra", 
            "patchwork", "grf", "glue", "stats", "future", "furrr", "rlang", "MatchIt", "Hmisc", "cli", "future.callr", 
            "callr", "fuzzyjoin"
          ),
          globals = c(
            "MapCovariateBalance", "CovariateBalance", "CForBenefit","RATETestcfwHelperParallel", "vi_models",
            "VariableImportanceWrapper", "DiscreteCovariatesToOneHot", "OverlapPlot", "CATEPlot", "SubgroupATETable", 
            "n_rankings", "num_threads", "num_clusters", "names_index", "compute_results", "excess_risk_quantiles",
            "names_index_comb", "mice", "excess_risk_cutoff", "path", "result_path", "continuous", "discrete", 
            "CausalForestCATETable", "CausalForestATESubgroupTable", "plot_model", "plot_covs", "cutpoints", 
            "cut_labels", "discrete_labels", "horizon", "thiazide_limit"
          ),
          seed = seed
        ),
        .f = \(i, cohort)  {
          options("parallelly.availableCores.methods" = "system")
          
          ### read in trained models .............................................................................. ----
          # with trees
          cfw <- list()
          cfw[["full_mod"]] <- readRDS(paste0(path, "cfw_full_mod_", i, ".rds"))
          for (j in vi_models) {
            cfw[[paste0("vi_mod_", j)]] <- readRDS(paste0(path, "cfw_vi_mod_", j, "_", i, ".rds"))
          }
          # without trees
          cfw_nt <- map(
            cfw,
            \(cf) cf[-c(4:11)] |> structure(class = c("causal_forest", "grf"))
          ) 
          
          # set up parallel workers over each model setup
          workers <- min(floor(num_threads / num_clusters), length(cfw_nt))
          threads_per_worker <- floor(num_threads / num_clusters / workers)
          
          ### covariate tables .................................................................................... ----
          cfw_covariates <- map(cfw_nt, \(cf) cf$X)
          saveRDS(cfw_covariates, paste0(result_path, "cfw_covariates_", i, ".rds"))
          rm(cfw_covariates)
          
          ### out-of-bag CATE estimates ........................................................................... ----
          cfw_cate <- map(
            .x = cfw,
            .f = \(cf) {
              pred <- predict(
                cf, 
                estimate.variance = TRUE, 
                num.threads = floor(num_threads / num_clusters)
              )
              continuous_cf <- continuous$cf
              discrete_cf <- discrete$cf
              pred_cens <- predict(
                cf,
                newdata = cohort |> 
                  filter(is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon)))) |>
                  select({{ continuous_cf }}) |>
                  bind_cols(
                    cohort |>
                      filter(is.na(!!sym(paste0("Y_", thiazide_limit, "_", horizon)))) |>
                      select({{ discrete_cf }}) |>
                      DiscreteCovariatesToOneHot()
                  ) |>
                  select(all_of(names(cf$X.orig))),
                estimate.variance = TRUE,
                num.threads = floor(num_threads / num_clusters)
              )
              return(
                tibble(
                  id = c(
                    which(!is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]])), 
                    which(is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]]))
                  ),
                  observed = !is.na(cohort[[paste0("Y_", thiazide_limit, "_", horizon)]]),
                  predictions = c(pred$predictions, pred_cens$predictions),
                  variance_estimates = c(pred$variance.estimates, pred_cens$variance.estimates)
                ) |>
                  arrange(id)
              )
            }
          )
          saveRDS(cfw_cate, paste0(result_path, "cfw_cate_", i, ".rds"))
          
          ### Predicted outcome and propensity to treat estimates ................................................. ----
          cfw_exp_out <- readRDS(paste0(path, "cfw_exp_out_", i, ".rds"))
          saveRDS(cfw_exp_out, paste0(result_path, "cfw_exp_out_", i, ".rds"))
          
          ### Expected outcomes in each treatment group ........................................................... ----
          cfw_outcome_trt_grp <- map(
            .x = cfw_cate,
            .f = \(cate, exp_out) {
              mu_hat_0 <- as.numeric(exp_out$Y_hat - exp_out$W_hat * cate$predictions)
              mu_hat_1 <- as.numeric(exp_out$Y_hat + (1 - exp_out$W_hat) * cate$predictions)
              
              return(
                tibble(
                  id = cate$id,
                  observed = cate$observed,
                  mu_hat_0 = mu_hat_0,
                  mu_hat_1 = mu_hat_1
                )
              )
            },
            exp_out = cfw_exp_out
          )
          saveRDS(cfw_outcome_trt_grp, paste0(result_path, "cfw_outcome_trt_grp_", i, ".rds"))
          rm(cfw_outcome_trt_grp)
          
          ### variable importance ................................................................................. ----
          cfw_vi <- map(cfw, \(cf) VariableImportanceWrapper(cf, names_index))
          saveRDS(cfw_vi, paste0(result_path, "cfw_vi_", i, ".rds"))
          rm(cfw_vi)
          rm(cfw)
          gc()
          
          ### test calibration .................................................................................... ----
          # Best linear fit regressing the residual outcome against the forest average 
          # treatment effect (times residual propensity) and the residual treatment effect
          # (times residual propensity) (not available for causal survival forest)
          cfw_cal <- map(cfw_nt, \(cf) test_calibration(cf, "HC0"))
          saveRDS(cfw_cal, paste0(result_path, "cfw_cal_", i, ".rds"))
          rm(cfw_cal)
          
          ### non-censoring ....................................................................................... ----
          # During the training process, causal_survival_forest() will print a warning
          # if any censoring probabilities are close to 1 for Y < h.
          
          ### overlap ............................................................................................. ----
          cfw_overlap <- (\(cf) {
            plot <- OverlapPlot(cf, draw = FALSE)
            rlang::env_bind(plot$plot_env, csf = NULL)
            return(plot)
          })(cf = cfw_nt$full_mod)
          saveRDS(cfw_overlap, paste0(result_path, "cfw_overlap_", i, ".rds"))
          rm(cfw_overlap)
          
          ### CATE distribution ................................................................................... ----
          cfw_cate_plot <- map(
            cfw_nt,
            \(cf) {
              plot <- CATEPlot(cf, draw = FALSE, xlim = c(-5, 20), ylim = c(0, 1))
              plot$plot_env |> rlang::env_bind(csf = NULL)
              return(plot)
            }
          )
          saveRDS(cfw_cate_plot, paste0(result_path, "cfw_cate_plot_", i, ".rds"))
          rm(cfw_cate_plot)
          
          ### Covariate balance ................................................................................... ----
          cfw_balance_data <- map(
            .x = cfw_nt,
            .f = \(cf) MapCovariateBalance(cf, names_index, everything(), plot_data_only = TRUE)
          )
          saveRDS(cfw_balance_data, paste0(result_path, "cfw_balance_data_", i, ".rds"))
          rm(cfw_balance_data)
          
          ### expected loss using R-loss criterion ................................................................ ----
          # An error measure based on the R-loss criterion (residual outcome minus residual treatment indicator times 
          # individual level treatment effect) is directly supported by the grf package on the training sample. 
          # The package also calculates a bias term and returns the debiased error (raw error - bias).
          cfw_rloss_error <- map(
            cfw_nt,
            \(cf) {
              list(
                rloss_error = predict(cf)$debiased.error,
                rloss_error_avg = weighted.mean(predict(cf)$debiased.error, cf$sample.weights)
              )
            }
          )
          saveRDS(cfw_rloss_error, paste0(result_path, "cfw_rloss_error_", i, ".rds"))
          rm(cfw_rloss_error)
          cfw_excess_error <- map(
            cfw_nt,
            \(cf) {
              list(
                excess_error = predict(cf)$excess.error,
                excess_error_avg = weighted.mean(predict(cf)$excess.error, cf$sample.weights)
              )
            }
          )
          saveRDS(cfw_excess_error, paste0(result_path, "cfw_excess_error_", i, ".rds"))
          rm(cfw_excess_error)
          
          cfw_r2 <- map(
            cfw_nt,
            \(cf) {
              ss_tot <- with(
                cf,
                weighted.mean(((Y.orig - Y.hat) - (W.orig - W.hat) * mean(predictions))^2, sample.weights)
              )
              ss_res <- with(
                cf, 
                weighted.mean(((Y.orig - Y.hat) - (W.orig - W.hat) * predictions)^2, sample.weights)
              )
              return(1 - ss_res / ss_tot)
            }
          )
          saveRDS(cfw_r2, paste0(result_path, "cfw_r2_", i, ".rds"))
          rm(cfw_r2)
          
          # Manually calculate raw R-loss errors. These are not equivalent to the debiased error, firstly because the 
          # grf package uses a bias corrected version of the R-loss estimating equation (using (centered_outcome - 
          # mean_of_centered_outcome)), second because the bias term calculated and corrected to get the debiased error 
          # is based on a jackknife approach, leaving out one tree at a time to estimate the excess error from
          # using finitely many trees.
          cfw_mse <- map(
            cfw_nt,
            \(cf) {
              pred <- with(cf, predictions * (W.orig - W.hat))
              obs <- with(cf, Y.orig - Y.hat)
              sq_err <- (pred - obs)^2
              return(weighted.mean(sq_err, cf$sample.weights))
            }
          )
          saveRDS(cfw_mse, paste0(result_path, "cfw_mse_", i, ".rds"))
          rm(cfw_mse)
          
          ### tests of heterogeneity (RATE-based and BLP-based) ................................................... ----
          # best linear projection for weighted causal forest
          cfw_blp <- map(
            cfw_nt,
            \(cf) {
              best_linear_projection(cf, cf$X.orig, vcov.type = "HC0") |>
                (\(blp) {
                  cbind(blp, p.adjust(blp[,"Pr(>|t|)"], "BH")) |>
                    structure(
                      dimnames = list(
                        rownames(blp),
                        c(colnames(blp), "pval_bh")
                      )
                    ) |>
                    as_tibble() |>
                    structure(
                      nobs = attr(blp, "nobs"),
                      df = attr(blp, "df")
                    ) |>
                    mutate(
                      name = c(
                        "(Intercept)",
                        fuzzyjoin::stringdist_left_join(
                          tibble(err = stringr::str_remove(rownames(blp), "TRUE")),
                          tibble(err = names(cf$X.orig), name = names(cf$X.orig)),
                          by = "err",
                          distance_col = "dist"
                        ) |>
                          group_by(err.x) |>
                          filter(dist == min(dist)) |>
                          ungroup() |>
                          pull(name)
                      )
                    )|>
                    left_join(names_index_comb, by = c("name" = "short")) |>
                    select(name, full, everything())
                })() |>
                arrange(pval_bh)
            }
          )
          saveRDS(cfw_blp, paste0(result_path, "cfw_blp_", i, ".rds"))
          rm(cfw_blp)
          
          # RATE test for censoring weighted causal forest
          # Run test for all one-hot encoded covariates. This means for factors we test
          # each level against all other levels. Also create tests of ordered factors
          # when appropriate.
          cfw_rate <- map(
            cfw_nt,
            \(cf) {
              RATETestcfwHelperParallel(
                cf,
                parallel = TRUE,
                n_cores = workers
              ) |>
                left_join(
                  names_index_comb,
                  by = c("covariate" = "short")
                ) |>
                mutate(full = ifelse(is.na(full), covariate, full)) |>
                select("covariate", "full", everything())
            }
          )
          saveRDS(cfw_rate, paste0(result_path, "cfw_rate_", i, ".rds"))
          rm(cfw_rate)
          
          ### C-for-benefit discrimination measure ................................................................ ----
          # Calculate C-for-benefit using the training data. Since the models are trained
          # on this data, we expect good discrimination performance.
          
          # Get estimated risk-difference from the causal forest models
          cfw_tau_hat_train <- map(cfw_nt, \(cf) as.numeric(cf$predictions))
          saveRDS(cfw_tau_hat_train, paste0(result_path, "cfw_tau_hat_train_", i, ".rds"))
          
          if (workers == 1) {
            plan(sequential)
          } else {
            plan(multisession, workers = workers)
          }
          
          # C-for-benefit using CATE matching and modified risk-difference using predicted risk under observed
          # treatment from each patient in a pair to obtain the RD. 
          cfw_cfb_cate <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("grf", "dplyr", "rlang", "glue", "stats", "MatchIt", 
                           "Hmisc", "cli", "tidyr", "stringr"),
              globals = c("CForBenefit"),
              seed = TRUE
            ),
            .f = \(cf) {
              CForBenefit(
                cf,
                match = "CATE",
                match_method = "nearest",
                match_distance = "mahalanobis",
                tau_hat_method = "risk_diff",
                CI = "simple",
                estimand = "ATT"
              )
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_cfb_cate, paste0(result_path, "cfw_cfb_cate_", i, ".rds"))
          rm(cfw_cfb_cate)
          
          # C-for-benefit using control-risk matching and predicted risk-difference of treated patient in each pair
          cfw_cfb_control_risk <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("grf", "dplyr", "rlang", "glue", "stats", "MatchIt", 
                           "Hmisc", "cli", "tidyr", "stringr"),
              globals = c("CForBenefit"),
              seed = TRUE
            ),
            .f = \(cf) {
              CForBenefit(
                cf,
                match = "control_risk",
                match_method = "nearest",
                match_distance = "mahalanobis",
                tau_hat_method = "tau_treated",
                CI = "simple",
                estimand = "ATT"
              )
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_cfb_control_risk, paste0(result_path, "cfw_cfb_control_risk_", i, ".rds"))
          rm(cfw_cfb_control_risk)
          
          # C-for-benefit using CATE + control risk full matching and modified risk-difference using predicted risk 
          # under observed treatment from each patient in a pair to obtain the RD. 
          cfw_cfb_combined <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("grf", "dplyr", "rlang", "glue", "stats", "MatchIt", 
                           "Hmisc", "cli", "tidyr", "stringr"),
              globals = c("CForBenefit"),
              seed = TRUE
            ),
            .f = \(cf) {
              CForBenefit(
                cf,
                match = "combined_outcome_risk",
                match_method = "nearest",
                match_distance = "mahalanobis",
                tau_hat_method = "risk_diff",
                CI = "simple",
                estimand = "ATE",
                replace = TRUE
              )
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_cfb_combined, paste0(result_path, "cfw_cfb_combined_", i, ".rds"))
          rm(cfw_cfb_combined)
          
          # model-based C-for-benefit
          cfw_mbcfb <- future_pmap(
            .l = list(
              cf = cfw_nt
            ),
            .options = furrr_options(
              packages = c("grf", "dplyr", "rlang", "glue", "stats", "MatchIt", 
                           "Hmisc", "cli", "tidyr", "stringr", "Rcpp"),
              globals = c("CForBenefit"),
              seed = TRUE
            ),
            .f = \(cf) {
              CForBenefit(
                cf,
                match = "none",
                CI = "simple"
              )
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_mbcfb, paste0(result_path, "cfw_mbcfb_", i, ".rds"))
          rm(cfw_mbcfb)
          rm(cfw_tau_hat_train)
          
          ### AIPW average treatment effects ...................................................................... ----
          # sex = 0 for females and sex = 1 for males
          cfw_ate <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("stringr", "purrr", "dplyr", "grf"),
              globals = c("SubgroupATETable", "names_index_comb"),
              seed = TRUE
            ),
            .f = \(model) {
              c(
                list(SubgroupATETable(NULL, NULL, model, level = 0.95)),
                map(
                  names(model$X.orig),
                  \(nm) {
                    if (str_detect(nm, "^X_\\d{2}$")) {
                      if (nm == "X_01") {
                        x <- c(map(seq(40, 70, 10), \(n) list(n, n + 10)), list(list(80, 110)))
                      } else if (str_detect(nm, "X_(69|7[0-9]|8[012])")) {
                        x <- map(1:4, \(n) as.list(quantile(model$X.orig[[nm]], seq(0, 1, 0.25))[c(n, n+1)]))
                      } else {
                        x <- c(0, 1)
                      }
                    } else {
                      x <- 1
                    }
                    map(x, \(x) SubgroupATETable(x, nm, model, names_index = names_index_comb, level = 0.95))
                  }
                )
              ) |>
                list_flatten() |>
                list_rbind()
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          
          saveRDS(cfw_ate, paste0(result_path, "cfw_ate_", i, ".rds"))
          
          ### excess risk strategies .............................................................................. ----
          # Look at the expected excess risk among those with a CATE above certain
          # cutoff (0%, 1%, 5%). Also look at population excess risk when applying 
          # treatment strategy where treatment is avoided based on these cut-offs.
          
          cfw_risk_strategy <- future_map2(
            .x = cfw_nt,
            .y = cfw_cate,
            .options = furrr_options(
              packages = c("stringr", "purrr", "dplyr", "grf"),
              globals = c("workers", "excess_risk_cutoff", "excess_risk_quantiles"),
              seed = TRUE
            ),
            .f = \(model, cate_all) {
              CFRiskStrategy(
                model, 
                cate_all$predictions,
                W,
                type = "train",
                level = 0.95,
                excess_risk_cutoff = excess_risk_cutoff,
                excess_risk_quantiles = excess_risk_quantiles
              )
            },
            W = cohort$W
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_risk_strategy, paste0(result_path, "cfw_risk_strategy_", i, ".rds"))
          rm(cfw_risk_strategy)
          
          ### data for excess risk reduction plot ................................................................. ----
          cfw_risk_reduction_plot <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("stringr", "purrr", "dplyr", "grf"),
              globals = c("workers", "excess_risk_cutoff", "excess_risk_quantiles"),
              seed = TRUE
            ),
            .f = \(model) {
              aipw <- map(
                seq(0, 1, 0.02),
                \(percentile) {
                  # determine the cutoff in excess risk among treated units to treat with new strategy
                  cutoff <- tibble(
                    predictions = model$predictions[,1],
                    W = model$W.orig,
                    sample_weight = model$sample.weights
                  ) |>
                    filter(W == 1) |>
                    arrange(predictions) |>
                    mutate(cum_weight = cumsum(sample_weight)) |>
                    (\(x) slice_head(x, n = sum(x$cum_weight <= percentile * sum(x$sample_weight))))() |>
                    pull(predictions) |>
                    max()
                  # ATE among samples changed from treated with new strategy
                  excess_risk <- tryCatch(
                    average_treatment_effect(
                      model, 
                      target.sample = "treated", 
                      subset = model$predictions[,1] > cutoff
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
                  # Difference in risk between old and new strategies 
                  risk_avoid <- tryCatch(
                    tibble(
                      percentile = percentile,
                      lower_0.025 = (1 - percentile) * excess_risk[["estimate"]] + qnorm(0.025) * (1 - percentile) * excess_risk[["std_err"]],
                      estimate = (1 - percentile) * excess_risk[["estimate"]],
                      upper_0.975 = (1 - percentile) * excess_risk[["estimate"]] + qnorm(0.975) * (1 - percentile) * excess_risk[["std_err"]],
                      std_err = (1 - percentile) * excess_risk[["std_err"]]
                    ),
                    error = \(e) {
                      tibble(
                        percentile = percentile,
                        lower_0.025 = NA_real_,
                        estimate = (1 - percentile) * excess_risk[["estimate"]],
                        upper_0.975 = NA_real_,
                        std_err = (1 - percentile) * excess_risk[["std_err"]]
                      )
                    }
                  )
                  
                  return(risk_avoid)
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
                    predictions = model$predictions[, 1],
                    W = model$W.orig,
                    sample_weight = model$sample.weights
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
                  ate_treated <- tibble(estimate = weighted.mean(pred$predictions, pred$sample_weight))
                  risk_avoid <- tryCatch(
                    tibble(
                      percentile = percentile,
                      estimate = ate_treated[["estimate"]] - ate_treated_cutoff[["estimate"]]
                    ),
                    error = \(e) {
                      tibble(
                        percentile = percentile,
                        estimate = ate_treated[["estimate"]]
                      )
                    }
                  )
                }
              ) |>
                list_rbind()
              ols <- add_row(
                ols, 
                percentile = 0,
                estimate =  mean(model$predictions[model$W.orig == 1]),
                .before = 2L
              ) |>
                slice(2:n())
              
              out <- list(aipw = aipw, ols = ols)
              return(out)
            }
          )
          if (workers > 1) {
            plan(sequential)
            plan(multisession, workers = workers)
          }
          saveRDS(cfw_risk_reduction_plot, paste0(result_path, "cfw_risk_reduction_plot_", i, ".rds"))
          rm(cfw_risk_reduction_plot)
          
          ## determine proportion of treated with predicted excess risk above zero
          cfw_risk_reduction_plot_zp <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("stringr", "purrr", "dplyr", "grf"),
              globals = c("workers", "excess_risk_cutoff", "excess_risk_quantiles"),
              seed = TRUE
            ),
            .f = \(model) {
              mean(model$predictions[,1] <= 0)
            }
          )
          saveRDS(cfw_risk_reduction_plot_zp, paste0(result_path, "cfw_risk_reduction_plot_zp_", i, ".rds"))
          rm(cfw_risk_reduction_plot_zp)
          
          ### data for bencalibr calibration plot ................................................................. ----
          cfw_bencalibr <- future_map(
            .x = cfw_nt,
            .options = furrr_options(
              packages = c("stringr", "purrr", "dplyr", "grf"),
              globals = c("workers"),
              seed = TRUE
            ),
            .f = \(model) {
              out <- tibble(
                y.observed = model$Y.orig,
                treat = model$W.orig,
                predicted.treat.0 = model$Y.hat - model$W.hat * model$predictions[, 1],
                predicted.treat.1 = model$Y.hat + (1 - model$W.hat) * model$predictions[, 1]
              )
              return(out)
            }
          )
          saveRDS(cfw_bencalibr, paste0(result_path, "cfw_bencalibr_", i, ".rds"))
          rm(cfw_bencalibr)
          
          ### data for dynamic subgroups .......................................................................... ----
          cfw_dyn <- map(names(cfw_nt), \(nm) readRDS(paste0(path, "cfw_", nm, "_dyn_", i, ".rds"))) |>
            structure(names = names(cfw_nt))
          saveRDS(cfw_dyn, paste0(result_path, "cfw_dyn_", i, ".rds"))
          rm(cfw_dyn)
          
          ### Excess risk distribution plots by covariate combinations ............................................ ----
          # covariate data for all (censored and uncensored) observations
          continuous_cf <- continuous$cf
          discrete_cf <- discrete$cf
          cov_data <- cohort |>
            select({{ continuous_cf }}) |>
            bind_cols(
              cohort |>
                select({{ discrete_cf }}) |>
                DiscreteCovariatesToOneHot()
            ) |>
            select(all_of(names(cfw_nt[[plot_model]]$X.orig))) |>
            as_tibble()
          
          # plot data 
          plot_data <- cov_data |>
            bind_cols(cfw_cate[[plot_model]]) |>
            mutate(
              ate_diff_signif = 
                (cfw_ate[[plot_model]] |> filter(subgroup == "Full population") |> pull(estimate) < 
                   predictions - qnorm(0.975) * sqrt(variance_estimates)) |
                (cfw_ate[[plot_model]] |> filter(subgroup == "Full population") |> pull(estimate) >
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
          saveRDS(plot_data, paste0(result_path, "cfw_cate_histogram_data_", plot_model, "_",i, ".rds"))
          
          ### Average risk of thiazide-induced hyponatraemia in subgroups by covariate combinations ............... ----
          cfw_subgroup_cate <- list()
          cfw_subgroup_cate$oneway <- map(
            as.list(names(cutpoints)),
            \(var) {
              cate_table <- CausalForestCATETable(
                (\(mod) {
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
                  return(mod)
                })(mod = cfw_nt[[plot_model]]),
                if (sum(grepl(plot_covs[names(cutpoints) == var], names(cfw_nt[[plot_model]]$X.orig))) == 1) {
                  list2(
                    "{var[1]}" := 
                      structure(
                        map2(c(-Inf, cutpoints[[var[1]]]), c(cutpoints[[var[1]]], Inf), \(x, y) as.character(c(x, y))), 
                        names = cut_labels[[var[1]]]
                      )
                  )
                } else {
                  list2(
                    "{var[1]}" := 
                      structure(
                        map(discrete_labels[[var]], \(x) paste0(var, "_", x)), 
                        names = cut_labels[[var]]
                      )
                  )
                }
              ) |>
                select(1:4) |>
                left_join(group_median[[var]], by = var)
            }
          ) |>
            structure(names = names(cutpoints))
          
          cfw_subgroup_cate$twoway <-
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
                cate_table <- CausalForestCATETable(
                  (\(mod) {
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
                    return(mod)
                  })(mod = cfw_nt[[plot_model]]),
                  if (
                    sum(grepl(plot_covs[names(cutpoints) == vars[1]], names(cfw_nt[[plot_model]]$X.orig))) == 1 &
                    sum(grepl(plot_covs[names(cutpoints) == vars[2]], names(cfw_nt[[plot_model]]$X.orig))) == 1
                  ) {
                    list2(
                      "{vars[1]}" := 
                        rep(
                          structure(
                            map2(
                              c(-Inf, cutpoints[[vars[1]]]), 
                              c(cutpoints[[vars[1]]], Inf), 
                              \(x, y) as.character(c(x, y))
                            ), 
                            names = cut_labels[[vars[1]]]
                          ), 
                          times = length(cutpoints[[vars[2]]]) + 1
                        ),
                      "{vars[2]}" := 
                        rep(
                          structure(
                            map2(
                              c(-Inf, cutpoints[[vars[2]]]), 
                              c(cutpoints[[vars[2]]], Inf), 
                              \(x, y) as.character(c(x, y))
                            ), 
                            names = cut_labels[[vars[2]]]
                          ), 
                          each = length(cutpoints[[vars[1]]]) + 1
                        )
                    )
                  } else if (
                    sum(grepl(plot_covs[names(cutpoints) == vars[1]], names(cfw_nt[[plot_model]]$X.orig))) == 1
                  ) {
                    list2(
                      "{vars[1]}" := 
                        rep(
                          structure(
                            map2(
                              c(-Inf, cutpoints[[vars[1]]]), 
                              c(cutpoints[[vars[1]]], Inf), 
                              \(x, y) as.character(c(x, y))
                            ), 
                            names = cut_labels[[vars[1]]]), 
                          times = length(cutpoints[[vars[2]]])
                        ),
                      "{vars[2]}" := 
                        rep(
                          structure(
                            map(
                              discrete_labels[[vars[2]]], 
                              \(x) paste0(vars[2], "_", x)
                            ), 
                            names = cut_labels[[vars[2]]]
                          ), 
                          each = length(cutpoints[[vars[1]]]) + 1
                        )
                    )
                  } else if (
                    sum(grepl(plot_covs[names(cutpoints) == vars[2]], names(cfw_nt[[plot_model]]$X.orig))) == 1
                  ) {
                    list2(
                      "{vars[1]}" := 
                        rep(
                          structure(
                            map(
                              discrete_labels[[vars[1]]], 
                              \(x) paste0(vars[1], "_", x)
                            ), 
                            names = cut_labels[[vars[1]]]
                          ), 
                          times = length(cutpoints[[vars[2]]]) + 1
                        ),
                      "{vars[2]}" := 
                        rep(
                          structure(
                            map2(
                              c(-Inf, cutpoints[[vars[2]]]), 
                              c(cutpoints[[vars[2]]], Inf), 
                              \(x, y) as.character(c(x, y))
                            ), 
                            names = cut_labels[[vars[2]]]
                          ), 
                          each = length(cutpoints[[vars[1]]])
                        )
                    )
                  } else {
                    list2(
                      "{vars[1]}" := 
                        rep(
                          structure(
                            map(
                              discrete_labels[[vars[1]]], 
                              \(x) paste0(vars[1], "_", x)
                            ), 
                            names = cut_labels[[vars[1]]]
                          ), 
                          times = length(cutpoints[[vars[2]]])
                        ),
                      "{vars[2]}" := 
                        rep(
                          structure(
                            map(
                              discrete_labels[[vars[2]]], 
                              \(x) paste0(vars[2], "_", x)
                            ), 
                            names = cut_labels[[vars[2]]]
                          ), 
                          each = length(cutpoints[[vars[1]]])
                        )
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
          
          saveRDS(cfw_subgroup_cate, paste0(result_path, "cfw_subgroup_cate_", plot_model, "_", i, ".rds"))
          rm(cfw_subgroup_cate)
          rm(cfw_ate)
          rm(cfw_cate)
        }
      )
      plan(sequential)
      gc()
    }
  )
}

## aggregate results from multiple imputation ----------------------------------------------------------------- -------- 
# map over each analysis setup
plan(sequential)
future_pmap(
  .l = list(
    table_path = list(
      paste0(tables_path, "4mo_thiazide_all/1year/"),
      paste0(tables_path, "4mo_thiazide_no_lab/"),
      paste0(tables_path, "4mo_thiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_thiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_thiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_thiazide_cc_comp/1year/")
    ),
    figure_path = list(
      paste0(figures_path, "4mo_thiazide_all/1year/"),
      paste0(figures_path, "4mo_thiazide_no_lab/"),
      paste0(figures_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_thiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_thiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_thiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_thiazide_all/1year/"),
      paste0(results_path, "4mo_thiazide_no_lab/"),
      paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_thiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_thiazide_complete_case/1year/"),
      paste0(results_path, "4mo_thiazide_cc_comp/1year/")
    ),
    W = list(
      tih_tbi_1_cohort_imputed[[1]]$W,
      tih_cohort$W,
      tih_tbi_1_cohort_imputed[[1]]$W, 
      tih_tbi_1_cohort_imputed[[1]]$W,
      bfz_tbi_1_cohort_imputed[[1]]$W_bvc,
      bfz_tbi_1_cohort_imputed[[1]]$W_bvc,
      hctz_tbi_1_cohort_imputed[[1]]$W_hvr,
      hctz_tbi_1_cohort_imputed[[1]]$W_hvr,
      tih_cohort_complete_case$W,
      tih_cohort_cc_comp[[1]]$W,
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "stringr", "purrr", "tibble", "tidyr", "glue", "rlang",
      "mice",
      "ggsci", "ggplot2", 
      "data.table", "cowplot", "gridExtra", "patchwork", "grf", 
      "stats", "future", "furrr",  "MatchIt", "Hmisc", "cli",
      "future.callr", "callr", "fuzzyjoin"
    ),
    globals = c("MapCovariateBalance", "CovariateBalance", "CForBenefit",
                "RATETestcfwHelperParallel", "VariableImportanceWrapper", 
                "DiscreteCovariatesToOneHot",
                "OverlapPlot", "CATEPlot", "SubgroupATETable", "n_rankings",
                "num_threads", "num_clusters", "names_index", "compute_results",
                "names_index_comb", "mice"),
    seed = TRUE
  ),
  .f = \(table_path, figure_path, result_path, W) {
    # read in results ............................................................................................. ----
    cfw_covariates <- list()
    cfw_cate <- list()
    cfw_exp_out <- list()
    cfw_outcome_trt_grp <- list()
    cfw_vi <- list()
    cfw_cal <- list()
    cfw_overlap <- list()
    cfw_cate_plot <- list()
    cfw_balance_data <- list()
    cfw_rloss_error <- list()
    cfw_mse <- list()
    cfw_blp <- list()
    cfw_rate <- list()
    cfw_tau_hat_train <- list()
    cfw_cfb_cate <- list()
    cfw_cfb_control_risk <- list()
    cfw_cfb_combined <- list()
    cfw_mbcfb <- list()
    cfw_ate <- list()
    cfw_risk_strategy <- list()
    cfw_risk_reduction_plot <- list()
    cfw_risk_reduction_plot_zp <- list()
    cfw_bencalibr <- list()
    cfw_dyn <- list()
    cfw_cate_histogram_data <- list()
    cfw_subgroup_cate <- list()
    for (i in seq_along(list.files(result_path, pattern = "^cfw_vi_\\d{1,}"))) {
      cfw_covariates[[i]] <- readRDS(paste0(result_path, "cfw_covariates_", i, ".rds"))
      cfw_cate[[i]] <- readRDS(paste0(result_path, "cfw_cate_", i, ".rds"))
      cfw_exp_out[[i]] <- readRDS(paste0(result_path, "cfw_exp_out_", i, ".rds"))
      cfw_outcome_trt_grp[[i]] <- readRDS(paste0(result_path, "cfw_outcome_trt_grp_", i, ".rds"))
      cfw_vi[[i]] <- readRDS(paste0(result_path, "cfw_vi_", i, ".rds"))
      cfw_cal[[i]] <- readRDS(paste0(result_path, "cfw_cal_", i, ".rds"))
      cfw_overlap[[i]] <- readRDS(paste0(result_path, "cfw_overlap_", i, ".rds"))
      cfw_cate_plot[[i]] <- readRDS(paste0(result_path, "cfw_cate_plot_", i, ".rds"))
      cfw_balance_data[[i]] <- readRDS(paste0(result_path, "cfw_balance_data_", i, ".rds"))
      cfw_rloss_error[[i]] <- readRDS(paste0(result_path, "cfw_rloss_error_", i, ".rds"))
      cfw_mse[[i]] <- readRDS(paste0(result_path, "cfw_mse_", i, ".rds"))
      cfw_blp[[i]] <- readRDS(paste0(result_path, "cfw_blp_", i, ".rds"))
      cfw_rate[[i]] <- readRDS(paste0(result_path, "cfw_rate_", i, ".rds"))
      cfw_tau_hat_train[[i]] <- readRDS(paste0(result_path, "cfw_tau_hat_train_", i, ".rds"))
      cfw_cfb_cate[[i]] <- readRDS(paste0(result_path, "cfw_cfb_cate_", i, ".rds"))
      cfw_cfb_control_risk[[i]] <- readRDS(paste0(result_path, "cfw_cfb_control_risk_", i, ".rds"))
      cfw_cfb_combined[[i]] <- readRDS(paste0(result_path, "cfw_cfb_combined_", i, ".rds"))
      cfw_mbcfb[[i]] <- readRDS(paste0(result_path, "cfw_mbcfb_", i, ".rds"))
      cfw_ate[[i]] <- readRDS(paste0(result_path, "cfw_ate_", i, ".rds"))
      cfw_risk_strategy[[i]] <- readRDS(paste0(result_path, "cfw_risk_strategy_", i, ".rds"))
      cfw_risk_reduction_plot[[i]] <- readRDS(paste0(result_path, "cfw_risk_reduction_plot_", i, ".rds"))
      cfw_risk_reduction_plot_zp[[i]] <- readRDS(paste0(result_path, "cfw_risk_reduction_plot_zp_", i, ".rds"))
      cfw_bencalibr[[i]] <- readRDS(paste0(result_path, "cfw_bencalibr_", i, ".rds"))
      cfw_dyn[[i]] <- readRDS(paste0(result_path, "cfw_dyn_", i, ".rds"))
      cfw_cate_histogram_data[[i]] <- readRDS(paste0(result_path, "cfw_cate_histogram_data_vi_mod_4_", i, ".rds"))
      cfw_subgroup_cate[[i]] <- readRDS(paste0(result_path, "cfw_subgroup_cate_vi_mod_4_", i, ".rds"))
    }
    
    # aggregate results ........................................................................................... ----
    plan(multisession, workers = min(num_clusters, length(cfw_covariates[[1]])))
    cfw_agg <- list()
    ### pool covariate tables ----
    cfw_agg$covariates <- cfw_covariates |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          # add NA column when used in some but not all imputation models
          cov_nm <- map(model, names)
          for (i in seq_along(cov_nm)) {
            for (j in seq_along(cov_nm)[-i]) {
              vars <- setdiff(cov_nm[[j]], names(model[[i]]))
              if (length(vars) > 0) {
                model[[i]] <- cbind(
                  model[[i]],
                  select(model[[j]], !!vars)|> 
                    mutate(across(everything(), ~ NA))
                )
              }
            }
          }
          # replace values with NA if they are imputed
          out <- model |> 
            list_transpose(simplify = FALSE) |>
            map(
              \(cov) {
                # cov |> 
                #   list_transpose() |>
                #   map_dbl(
                #     \(x) ifelse(length(unique(x)) == 1, x[1], NA_real_)
                #     )
                cov |>
                  (\(x) structure(x, names = seq_along(x)))() |>
                  as_tibble() |> 
                  (\(x) {
                    rowmeans <- rowMeans(x)
                    ifelse(
                      !is.na(rowmeans) & rowmeans == x[[1]], 
                      x[[1]], 
                      NA_real_
                    )
                  })()
              }
            ) |>
            as_tibble()
          # NOTE: the above map is slow if checking equality. The current use of rowMeans is faster, but 
          # may not be correct in all cases.
          
          return(out)
        }
      )
    
    ### pool CATE estimates ----
    cfw_agg$cate <- cfw_cate |>
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
                    pool <- pool.scalar(Q = x$predictions, U = x$variance_estimates, n = n, k = 1, rule = "rubin1987")
                    out <- tibble(
                      id = x$id[1],
                      observed = x$observed[1],
                      estimate = pool$qbar,
                      std_err = sqrt(pool$t),
                      lower = estimate - qnorm(0.975) * std_err,
                      upper = estimate + qnorm(0.975) * std_err
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      id = x$id,
                      observed = x$observed,
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
    
    ### pool predicted outcome and propensity scores ----
    cfw_agg$exp_out <- cfw_exp_out |>
      list_transpose(simplify = FALSE) |>
      map(\(x) list_transpose(x, simplify = TRUE)) |>
      list_transpose(simplify = FALSE) |>
      (\(x) {
        future_map(
          .x = x,
          .options = furrr_options(
            packages = c("mice", "dplyr", "purrr")
          ),
          .f = \(x) {
            if (is.list(x) && length(x[[1]]) > 1) {
              pool_yhat <- pool.scalar(Q = x$Y_hat, U = 0, rule = "rubin1987")
              pool_what <- pool.scalar(Q = x$W_hat, U = 0, rule = "rubin1987")
              out <- tibble(
                id = x$id[1],
                observed = x$observed[1],
                Y_hat = pool_yhat$qbar,
                W_hat = pool_what$qbar
              )
            } else {
              pool_yhat <- NULL
              pool_what <- NULL
              out <- tibble(
                id = x$id,
                observed = x$observed,
                Y_hat = x$Y_hat,
                W_hat = x$W_hat
              )
            }
            return(list(table = out, pool_results = list(Y_hat = pool_yhat, W_hat = pool_what)))
          }
        )
      })() |>
      list_transpose()
    
    ### pool estimated outcomes from each treatment group ----
    cfw_agg$outcome_trt_grp <- cfw_outcome_trt_grp |>
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
                      mu_hat_0 = pool.scalar(Q = x$mu_hat_0, U = 0, n = n, k = 1, rule = "rubin1987"),
                      mu_hat_1 = pool.scalar(Q = x$mu_hat_1, U = 0, n = n, k = 1, rule = "rubin1987")
                    )
                    out <- tibble(
                      id = x$id[1],
                      observed = x$observed[1],
                      mu_hat_0 = pool$mu_hat_0$qbar,
                      mu_hat_1 = pool$mu_hat_1$qbar
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      id = x$id,
                      observed = x$observed,
                      mu_hat_0 = x$mu_hat_0,
                      mu_hat_1 = x$mu_hat_1
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
    
    ### pool variable importance estimates ----
    cfw_agg$variable_importance <- cfw_vi |>
      list_transpose(simplify = FALSE) |>
      (\(x) {
        future_map(
          .x = x,
          .options = furrr_options(
            packages = c("mice", "dplyr", "purrr")
          ),
          .f = \(model) {
            model |>
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
              )
          }
        )
      })() |> 
      map(
        \(model) {
          tbl <- map(model,\(x) x[[1]]) |> 
            list_rbind() |>
            arrange(desc(variable_importance_num))
          pool <- map(model, \(x) x[[2]])
          names(pool) <- map_chr(model, \(x) x[[1]]$variable_name)
          return(list(table = tbl, pool_results = pool))
        }
      )
    
    # aggregate variable importance by covariates
    cfw_agg$variable_importance_comb <-
      map(
        cfw_agg$variable_importance,
        \(vi) {
          (\(vi, nm_index) {
            vi_out <- vi[["table"]][, 1:2]
            for (nm in names_index$full) {
              if (
                nrow(
                  filter(
                    vi$table,
                    str_detect(vi$table$variable_name, str_c("^", nm, " "))
                  )
                ) > 1
              ) {
                vi_out <- bind_rows(
                  vi_out,
                  vi$table |>
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
              vi_out |>
                mutate(
                  variable_name = trimws(variable_name),
                  variable_importance = sprintf("%.5f", variable_importance_num)
                ) |>
                arrange(desc(variable_importance_num)) |>
                left_join(
                  select(nm_index, "short", "full"),
                  by = c("variable_name" = "full")
                )
            )
          })(
            vi = vi,
            nm_index = filter(names_index, !grepl("^(X_(5[5-9]|6[0-8])$|[^X])", short))
          )
        }
      )
    
    ### pool test calibration results ----
    cfw_agg$calibration <- cfw_cal |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          map(
            model, 
            \(iter) {
              dplyr::tibble(
                mean.forest.prediction = iter[1, 1:2], 
                differential.forest.prediction = iter[2, 1:2]
              )
            }
          ) |>
            purrr::list_transpose(simplify = FALSE) |>
            (\(x) {
              n <- attr(cfw_cal[[1]]$full, "nobs")
              df <- attr(cfw_cal[[1]]$full, "df")
              if (length(x$mean.forest.prediction) > 1) {
                pool <- map(
                  x,
                  \(x, n) {
                    x <- list_transpose(x)
                    pool.scalar(Q = x$Estimate, U = x$`Std. Error`, n = n, k = 2, rule = "rubin1987")
                  }, 
                  n = n
                )
                table <- matrix(
                  data = c(
                    pool$mean.forest.prediction$qbar, 
                    pool$differential.forest.prediction$qbar,
                    pool$mean.forest.prediction$t,
                    pool$differential.forest.prediction$t,
                    pool$mean.forest.prediction$qbar / pool$mean.forest.prediction$t,
                    pool$differential.forest.prediction$qbar / pool$differential.forest.prediction$t,
                    pt(abs(pool$mean.forest.prediction$qbar / pool$mean.forest.prediction$t), df = pool$mean.forest.prediction$df, lower.tail = FALSE),
                    pt(abs(pool$differential.forest.prediction$qbar / pool$differential.forest.prediction$t), df = pool$differential.forest.prediction$df, lower.tail = FALSE)
                  ),
                  nrow = 2,
                  ncol = 4
                ) |>
                  structure(
                    class = "coeftest",
                    dimnames = list(
                      c("mean.forest.prediction", "differential.forest.prediction"),
                      c("Estimate", "Std. Error", "t value", "Pr(>t)")
                    ),
                    method = "Best linear fit using forest predictions (on held-out data)\nas well as the mean forest prediction as regressors, along\nwith one-sided heteroskedasticity-robust (HC3) SEs",
                    df = c(pool$mean.forest.prediction$df, pool$differential.forest.prediction$df),
                    nobs = n
                  )
              } else {
                pool <- NULL
                table <- matrix(
                  data = c(
                    x$mean.forest.prediction[[1]][1], 
                    x$differential.forest.prediction[[1]][1],
                    x$mean.forest.prediction[[1]][2],
                    x$differential.forest.prediction[[1]][2],
                    x$mean.forest.prediction[[1]][1] / x$mean.forest.prediction[[1]][2],
                    x$differential.forest.prediction[[1]][1] / x$differential.forest.prediction[[1]][2],
                    pt(abs(x$mean.forest.prediction[[1]][1] / x$mean.forest.prediction[[1]][2]), df = n, lower.tail = FALSE),
                    pt(abs(x$differential.forest.prediction[[1]][1] / x$differential.forest.prediction[[1]][2]), df = df, lower.tail = FALSE)
                  ),
                  nrow = 2,
                  ncol = 4
                ) |>
                  structure(
                    class = "coeftest",
                    dimnames = list(
                      c("mean.forest.prediction", "differential.forest.prediction"),
                      c("Estimate", "Std. Error", "t value", "Pr(>t)")
                    ),
                    method = "Best linear fit using forest predictions (on held-out data)\nas well as the mean forest prediction as regressors, along\nwith one-sided heteroskedasticity-robust (HC3) SEs",
                    df = df,
                    nobs = n
                  )
              }
              return(c(list(table = table), pool))
            })()
        }
      )
    
    ### overlap ----
    cfw_overlap <- tibble(W_hat = cfw_agg$exp_out$table$W_hat) |>
      ggplot() + 
      geom_density(
        aes(x = W_hat, y = after_stat(density)),
        color = "white",
        fill = "darkgrey"
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
      labs(title = "Estimates propensity score for exposure") + 
      xlab("Propensity score") +
      ylab("")
    cfw_overlap$plot_env <- new_environment(data = list(W_hat = cfw_agg$exp_out$table$W_hat))
    cfw_agg$overlap <- cfw_overlap
    
    ### pool balance plot data ----
    cfw_agg$balance_data <- cfw_balance_data |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          model |>
            list_transpose(simplify = FALSE) |>
            imap(
              \(x, nm) {
                if (nm %in% c("grepl_string", "names_string")) return(x[[1]])
                x <- list_transpose(x)
                x$love_data <- x$love_data |>
                  list_transpose() |>
                  (\(x) {
                    out <- list()
                    out$covariate_name <- x$covariate_name[[1]]
                    out$Cohort <- x$Cohort[[1]]
                    out$value <- list_transpose(x$value) |>
                      map_dbl(
                        \(val) {
                          pool <- pool.scalar(Q = val, U = 0, n = length(val), k = 1, rule = "rubin1987")
                          return(pool$qbar)
                        }
                      )
                    return(as_tibble(out))
                  })()
                x$cd_data <- x$cd_data |>
                  list_transpose() |>
                  imap(
                    \(x, nm) {
                      if (nm == "IPW") {
                        return(
                          list_transpose(x) |>
                            map_dbl(
                              \(val) {
                                pool <- pool.scalar(Q = val, U = 0, n = length(val), k = 1, rule = "rubin1987")
                                return(pool$qbar)
                              }
                            )
                        )
                      }
                      if (nm == "W") return(x[[1]])
                      return(x[[1]])
                    }
                  ) |>
                  as_tibble()
                x$cd_type <- dplyr::coalesce(!!!list_transpose(as.list(x$cd_type)))
                x$X_type <- dplyr::coalesce(!!!list_transpose(as.list(x$X_type)))
                x$X_orig_old <- x$X_orig_old[[1]]
                
                return(x)
              }
            )
        }
      )
    
    ### pool r-loss error results ----
    cfw_agg$rloss_error <- cfw_rloss_error |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          n <- length(model[[1]])
          if (length(model) > 1) {
            pool <- map(model, mean) |>
              list_transpose(simplify = TRUE) |>
              (\(x) pool.scalar(Q = x[[1]], U = 0, n = n, k = 1, rule = "rubin1987"))()
            table <- tibble(rloss_error = pool$qbar)
          } else {
            pool <- NULL
            table <- tibble(rloss_error = map(model, mean)[[1]])
          }
          return(list(table = table, pool = pool))
        }
      )
    
    ### pool mse results ----
    cfw_agg$mean_square_error <- cfw_mse |>
      list_transpose(simplify = TRUE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            pool <- pool.scalar(Q = model, U = 0, n = 1, k = 1, rule = "rubin1987")
            table <- tibble(mean_square_error = pool$qbar)
          } else {
            pool <- NULL
            table <- tibble(mean_square_error = model)
          }
          return(list(table = table, pool = pool))
        }
      )
    
    ### pool tests of heterogeneity ----
    # best linear projection
    cfw_agg$best_linear_projection <- cfw_blp |>
      purrr::list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            n <- attr(model[[1]], "nobs")
            out <- model |>
              purrr::map(\(df) df |> arrange(name, .locale = "en")) |>
              purrr::list_transpose(simplify = FALSE) |>
              (\(list) {
                list <- list[c("name", "Estimate", "Std. Error")]
                i <- 1
                repeat {
                  var_name <- c()
                  for (j in seq_along(list$name)) {
                    var_name[j] <- list$name[[j]][i]
                  }
                  if (length(unique(var_name)) > 1) {
                    name_to_append <- sort(var_name)[1]
                  } else if (!is.na(unique(var_name))){
                    i <- i + 1
                    next
                  } else {
                    break
                  }
                  for (j in seq_along(list$name)) {
                    if (is.na(list$name[[j]][i]) || list$name[[j]][i] != name_to_append) {
                      list$name[[j]] <- append(list$name[[j]], name_to_append, i - 1)
                      list$Estimate[[j]] <- append(list$Estimate[[j]], 0, i - 1)
                      list$`Std. Error`[[j]] <- append(list$`Std. Error`[[j]], 0, i - 1)
                    }
                  }
                  i <- i + 1
                }
                return(list)
              })() |>
              purrr::map(\(x) list_transpose(x, simplify = TRUE)) |>
              purrr::list_transpose(simplify = FALSE) |>
              purrr::map(
                \(x) {
                  pool <- pool.scalar(Q = x$Estimate, U = (x$`Std. Error`)^2, n = n, k = 1, rule = "rubin1987")
                  out <- dplyr::tibble(
                    name = x$name[1],
                    Estimate = pool$qbar,
                    std_err = sqrt(pool$t),
                    t_value = Estimate / std_err,
                    pval = 2 * pt(abs(t_value), pool$df, lower.tail = FALSE)
                  )
                  return(list(table = out, pool = pool))
                }
              ) |>
              list_transpose() |>
              (\(x) {
                x$table <- dplyr::mutate(x$table, pval_bh = p.adjust(pval, "BH"))
                return(x)
              })()
            names(out$pool) <- out$table$name
            out$table <- out$table |>
              dplyr::left_join(names_index_comb, by = c("name" = "short")) |>
              dplyr::select(name, full, everything()) |>
              dplyr::arrange(pval_bh)
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                rename(
                  "std_err" = "Std. Error",
                  "t_value" = "t value",
                  "pval" = "Pr(>|t|)"
                )
            )
          }
          return(out)
        } 
      )
    
    ### pool rank average treatment effect ----
    cfw_agg$rate <- cfw_rate |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              map(\(df) arrange(df, covariate, .locale = "en")) |>
              list_transpose(simplify = FALSE) |>
              (\(list) {
                list <- list[c("covariate", "estimate", "std.err")]
                i <- 1
                repeat {
                  var_name <- c()
                  for (j in seq_along(list$covariate)) {
                    var_name[j] <- list$covariate[[j]][i]
                  }
                  if (length(unique(var_name)) > 1) {
                    name_to_append <- sort(var_name)[1]
                  } else if (!is.na(unique(var_name))){
                    i <- i + 1
                    next
                  } else {
                    break
                  }
                  for (j in seq_along(list$covariate)) {
                    if (is.na(list$covariate[[j]][i]) || list$covariate[[j]][i] != name_to_append) {
                      list$covariate[[j]] <- append(list$covariate[[j]], name_to_append, i - 1)
                      list$estimate[[j]] <- append(list$estimate[[j]], 0, i - 1)
                      list$std.err[[j]] <- append(list$std.err[[j]], 0, i - 1)
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
                  pool <- pool.scalar(Q = x$estimate, U = (x$std.err)^2, n = 1, k = 1, rule = "rubin1987")
                  out <- tibble(
                    name = x$covariate[1],
                    estimate = pool$qbar,
                    std_err = sqrt(pool$t),
                    t_value = estimate / std_err,
                    pval = 2 * pnorm(abs(t_value), lower.tail = FALSE)
                  )
                  return(return(list(table = out, pool = pool)))
                }
              ) |>
              list_transpose() |>
              (\(x) {
                x$table <- mutate(x$table, pval_bh = p.adjust(pval, "BH"))
                return(x)
              })()
            names(out$pool) <- out$table$name
            out$table <- out$table |>
              left_join(names_index_comb, by = c("name" = "short")) |>
              mutate(full = ifelse(is.na(full), name, full)) |>
              select(name, full, everything()) |>
              arrange(pval_bh)
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                reframe(
                  name = covariate,
                  full = full,
                  estimate = estimate,
                  std_err = std.err,
                  t_value = estimate / std_err,
                  pval = pval,
                  pval_bh = pval_bh
                )
            )
          }
          return(out)
        }
      )
    saveRDS(cfw_agg, paste0(result_path, "cfw_agg.rds"))
    
    ### pool C-for-benefit ----
    cfw_agg$c_for_benefit_cate <- cfw_cfb_cate |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose() |>
              (\(x) x[c("c_for_benefit", "c_for_benefit_se")])() |>
              (\(x) {
                pool <- pool.scalar(Q = x$c_for_benefit, U = (x$c_for_benefit_se)^2, n = 1, k = 1, rule = "rubin1987")
                out <- tibble(
                  c_for_benefit = pool$qbar,
                  c_for_benefit_se = sqrt(pool$t),
                  lower_CI = c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se,
                  upper_CI = c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
                )
                return(return(list(table = out, pool = pool)))
              })()
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                (\(x) {
                  tibble(
                    "c_for_benefit" = x$c_for_benefit,
                    "c_for_benefit_se" = x$c_for_benefit_se, 
                    "lower_ci" = x$lower_ci, 
                    "upper_ci" = x$upper_ci
                  )
                })()
            )
          }
          return(out)
        }
      )
    
    cfw_agg$c_for_benefit_control_risk <- cfw_cfb_control_risk |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose() |>
              (\(x) x[c("c_for_benefit", "c_for_benefit_se")])() |>
              (\(x) {
                pool <- pool.scalar(Q = x$c_for_benefit, U = (x$c_for_benefit_se)^2, n = 1, k = 1, rule = "rubin1987")
                out <- tibble(
                  c_for_benefit = pool$qbar,
                  c_for_benefit_se = sqrt(pool$t),
                  lower_CI = c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se,
                  upper_CI = c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
                )
                return(return(list(table = out, pool = pool)))
              })()
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                (\(x) {
                  tibble(
                    "c_for_benefit" = x$c_for_benefit,
                    "c_for_benefit_se" = x$c_for_benefit_se, 
                    "lower_ci" = x$lower_ci, 
                    "upper_ci" = x$upper_ci
                  )
                })()
            )
          }
          return(out)
        }
      )
    
    cfw_agg$c_for_benefit_combined <- cfw_cfb_combined |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose() |>
              (\(x) x[c("c_for_benefit", "c_for_benefit_se")])() |>
              (\(x) {
                pool <- pool.scalar(Q = x$c_for_benefit, U = (x$c_for_benefit_se)^2, n = 1, k = 1, rule = "rubin1987")
                out <- tibble(
                  c_for_benefit = pool$qbar,
                  c_for_benefit_se = sqrt(pool$t),
                  lower_CI = c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se,
                  upper_CI = c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
                )
                return(return(list(table = out, pool = pool)))
              })()
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                (\(x) {
                  tibble(
                    "c_for_benefit" = x$c_for_benefit,
                    "c_for_benefit_se" = x$c_for_benefit_se, 
                    "lower_ci" = x$lower_ci, 
                    "upper_ci" = x$upper_ci
                  )
                })()
            )
          }
          return(out)
        }
      )
    
    cfw_agg$model_based_c_for_benefit <- cfw_mbcfb |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          if (length(model) > 1) {
            out <- model |>
              list_transpose() |>
              (\(x) x[c("c_for_benefit", "c_for_benefit_se", "lower_CI", "upper_CI")])() |>
              (\(x) {
                pool <- pool.scalar(Q = x$c_for_benefit, U = (x$c_for_benefit_se)^2, n = 1, k = 1, rule = "rubin1987")
                out <- tibble(
                  c_for_benefit = pool$qbar,
                  c_for_benefit_se = sqrt(pool$t),
                  lower_CI = c_for_benefit + qnorm(0.5 - 0.95 / 2) * c_for_benefit_se,
                  upper_CI = c_for_benefit + qnorm(0.5 + 0.95 / 2) * c_for_benefit_se
                )
                return(return(list(table = out, pool = pool)))
              })()
          } else {
            out <- list(
              pool = NULL,
              table = model[[1]] |>
                (\(x) {
                  tibble(
                    "c_for_benefit" = x$c_for_benefit,
                    "c_for_benefit_se" = x$c_for_benefit_se, 
                    "lower_ci" = x$lower_ci, 
                    "upper_ci" = x$upper_ci
                  )
                })()
            )
          }
          return(out)
        }
      )
    
    ### pool ATE results ----
    # change names of blood test quartiles so they are consistent
    cfw_agg$ate <- cfw_ate |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          nm_qrt <- names_index_comb |>
            filter(grepl("^X_(69|7[0-9]|8[0-2])", short)) |>
            pull(full)
          
          model <- map(
            model,
            \(x) {
              x|>
                mutate(
                  cov = str_extract(subgroup, paste0(nm_qrt, collapse = "|")),
                  quartile_low = ifelse(
                    is.na(cov), 
                    NA_character_, 
                    str_extract(
                      str_remove(
                        str_remove(
                          subgroup, 
                          cov
                        ), 
                        "-\\s"
                      ), 
                      "-?\\d{1,}\\.?\\d?"
                    )
                  ),
                  quartile_high = ifelse(
                    is.na(cov), 
                    NA_character_, 
                    str_remove(
                      str_extract(
                        str_remove(
                          str_remove(
                            subgroup,
                            cov
                          ), 
                          "-\\s"
                        ), 
                        "-{1,2}\\d{1,}\\.?\\d?$"
                      ), 
                      "^-"
                    )
                  ) 
                ) |>
                mutate(
                  cov_grp = ifelse(
                    is.na(cov), 
                    NA_character_, 
                    paste0(
                      cov, 
                      " - Q", 
                      which(
                        sort.int(as.numeric(quartile_low)) == 
                          as.numeric(quartile_low)
                      )
                    )
                  ),
                  .by = "cov"
                ) |>
                mutate(
                  subgroup = ifelse(is.na(cov), subgroup, cov_grp)
                ) |>
                select(-c("cov", cov_grp))
            }
          )
          
          if (length(model) > 1) {
            out <- model |>
              map(\(df) df |> arrange(subgroup, .locale = "en")) |>
              list_transpose(simplify = FALSE) |>
              (\(list) {
                list <- list[c("subgroup", "estimate", "std_err", "n", "quartile_low", "quartile_high")]
                i <- 1
                repeat {
                  var_name <- c()
                  for (j in seq_along(list$subgroup)) {
                    var_name[j] <- list$subgroup[[j]][i]
                  }
                  if (length(unique(var_name)) > 1) {
                    name_to_append <- sort(var_name)[1]
                  } else if (!is.na(unique(var_name))){
                    i <- i + 1
                    next
                  } else {
                    break
                  }
                  for (j in seq_along(list$subgroup)) {
                    if (
                      !is.na(list$subgroup[[j]][i]) && 
                      list$subgroup[[j]][i] == name_to_append
                    ) {
                      n <- list$n[[j]][i] # number of samples in group to append
                      break
                    }
                  }
                  for (j in seq_along(list$subgroup)) {
                    if (is.na(list$subgroup[[j]][i]) || list$subgroup[[j]][i] != name_to_append) {
                      list$subgroup[[j]] <- append(list$subgroup[[j]], name_to_append, i - 1)
                      # when not in the model, assume estimate is the average of full population
                      list$estimate[[j]] <- append(
                        list$estimate[[j]], 
                        list$estimate[[j]][list$subgroup[[j]] == "Full population"], 
                        i - 1
                      )
                      list$std_err[[j]] <- append(
                        list$std_err[[j]], 
                        list$std_err[[j]][list$subgroup[[j]] == "Full population"], 
                        i - 1
                      )
                      list$n[[j]] <- append(
                        list$n[[j]],
                        n,
                        i - 1
                      )
                      list$quartile_low[[j]] <- append(
                        list$quartile_low[[j]],
                        NA_character_,
                        i - 1
                      )
                      list$quartile_high[[j]] <- append(
                        list$quartile_high[[j]],
                        NA_character_,
                        i - 1
                      )
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
                    name = x$subgroup[1],
                    # When subgroup along an imputed variable, n may differ for each imputation. 
                    # How to define n? Currently use the mean. 
                    n = floor(mean(x$n)),
                    estimate = pool$qbar,
                    std_err = sqrt(pool$t),
                    `95% CI - lower` = estimate - qnorm(0.975) * std_err,
                    `95% CI - upper` = estimate + qnorm(0.975) * std_err,
                    quartile_low = list(x$quartile_low),
                    quartile_high = list(x$quartile_high)
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
                  "name" = "subgroup"
                )
            )
          }
          return(out)
        }
      )
    
    ### pool excess risk strategy results ----
    cfw_agg$risk_strategy <- cfw_risk_strategy |>
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
                    pool <- pool.scalar(Q = x$estimate, U = (x$std_err)^2, n = 1, k = 1, rule = "rubin1987")
                    out <- tibble(
                      "name" = x$name[1],
                      "lower_0.025" = pool$qbar + qnorm(0.5 - 0.95 / 2) * sqrt(pool$t),
                      "estimate" = pool$qbar,
                      "upper_0.975" = pool$qbar + qnorm(0.5 + 0.95 / 2) * sqrt(pool$t),
                      "std_err" = sqrt(pool$t)
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
    
    ### pool data for excess risk reduction plots ----
    cfw_agg$cfw_risk_reduction_plot <- cfw_risk_reduction_plot |>
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
                            "upper_0.975" = pool.scalar(Q = x[["upper_0.975"]], U = 0, n = 1, k = 1, rule = "rubin1987")
                          )
                          out <- tibble(
                            "percentile" = x$percentile[1],
                            "lower_0.025" = pool[["lower_0.025"]]$qbar,
                            "estimate" = pool$estimate$qbar,
                            "upper_0.975" = pool[["upper_0.975"]]$qbar,
                            "std_err" = NA
                          )
                          return(list(table = out, pool = pool))
                        } else {
                          pool <- pool.scalar(Q = x$estimate, U = (x$std_err)^2, n = 1, k = 1, rule = "rubin1987")
                          out <- tibble(
                            "percentile" = x$percentile[1],
                            "lower_0.025" = pool$qbar + qnorm(0.5 - 0.95 / 2) * sqrt(pool$t),
                            "estimate" = pool$qbar,
                            "upper_0.975" = pool$qbar + qnorm(0.5 + 0.95 / 2) * sqrt(pool$t),
                            "std_err" = sqrt(pool$t)
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
                        "std_err" = x$std_err
                      )
                    })()
                )
              }
            )
          }
          return(out)
        }
      )
    cfw_agg$cfw_risk_reduction_plot_zp <- cfw_risk_reduction_plot_zp |>
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
    
    ### pool data for bencalibr calibration plots ----
    cfw_agg$bencalibr <- cfw_bencalibr |>
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
                      y.observed = pool.scalar(Q = x$y.observed, U = 0, n = n, k = 1, rule = "rubin1987"),
                      treat = pool.scalar(Q = x$treat, U = 0, n = n, k = 1, rule = "rubin1987"),
                      predicted.treat.0 = pool.scalar(Q = x$predicted.treat.0, U = 0, n = n, k = 1, rule = "rubin1987"),
                      predicted.treat.1 = pool.scalar(Q = x$predicted.treat.1, U = 0, n = n, k = 1, rule = "rubin1987")
                    )
                    out <- tibble(
                      y.observed = pool$y.observed$qbar,
                      treat = pool$treat$qbar,
                      predicted.treat.0 = pool$predicted.treat.0$qbar,
                      predicted.treat.1 = pool$predicted.treat.1$qbar
                    )
                  } else {
                    pool <- NULL
                    out <- tibble(
                      y.observed = x$y.observed,
                      treat = x$treat,
                      predicted.treat.0 = x$predicted.treat.0,
                      predicted.treat.1 = x$predicted.treat.1
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
    
    ### pool data for dynamic subgroups ----
    cfw_agg$dyn <- cfw_dyn |>
      list_transpose(simplify = FALSE) |>
      future_map(
        .options = furrr_options(
          packages = c("mice", "dplyr", "purrr")
        ),
        .f = \(model) {
          out <- model |>
            list_transpose(simplify = FALSE) |>
            (\(list) {
              rank_tbl_pool <- list$rank_tbl |>
                list_transpose(simplify = FALSE) |>
                (\(x) x[c("id", "sample_weights", "Y_hat", "W_hat", "mu_hat_0", "mu_hat_1", "tau_hat", "aipw_scores", "fold")])() |> # , "tau_hat_var"
                map(\(x) list_transpose(x, simplify = TRUE)) |>
                list_transpose(simplify = FALSE) |>
                (\(x) {
                  n <- length(x)
                  map(
                    x,
                    \(x, n) {
                      nm <- names(x)
                      tmp <- tryCatch(as.double(x), error = \(e) {})
                      if (!is.null(tmp)) {x <- tmp; names(x) <- nm}
                      if (is.list(x)) {
                        pool <- list(
                          "sample_weights" = pool.scalar(Q = x[["sample_weights"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "Y_hat" = pool.scalar(Q = x[["Y_hat"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "W_hat" = pool.scalar(Q = x[["W_hat"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "mu_hat_0" = pool.scalar(Q = x[["mu_hat_0"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "mu_hat_1" = pool.scalar(Q = x[["mu_hat_1"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "tau_hat" = pool.scalar(Q = x[["tau_hat"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "aipw_scores" = pool.scalar(Q = x[["aipw_scores"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                          "fold" = pool.scalar(Q = x[["fold"]], U = 0, n = 1, k = 1, rule = "rubin1987")
                        )
                        tbl <- tibble(
                          "id" = x[["id"]][[1]],
                          "sample_weights" = pool$sample_weights$qbar,
                          "Y_hat" = pool$Y_hat$qbar,
                          "W_hat" = pool$W_hat$qbar,
                          "mu_hat_0" = pool$mu_hat_0$qbar,
                          "mu_hat_1" = pool$mu_hat_1$qbar,
                          "tau_hat" = pool$tau_hat$qbar,
                          "aipw_scores" = pool$aipw_scores$qbar,
                          "fold" = pool$fold$qbar
                        )
                      } else {
                        pool <- NULL
                        tbl <- tibble(
                          "id" = x["id"],
                          "sample_weights" = x["sample_weights"],
                          "Y_hat" = x["Y_hat"],
                          "W_hat" = x["W_hat"],
                          "mu_hat_0" = x["mu_hat_0"],
                          "mu_hat_1" = x["mu_hat_1"],
                          "tau_hat" = x["tau_hat"],
                          "aipw_scores" = x["aipw_scores"],
                          "fold" = x["fold"]
                        )
                      }
                      return(list(tbl = tbl, pool = pool))
                    },
                    n = n)
                })() |>
                list_transpose()
              
              out_heatmap_tbl <- list()
              out_heatmap_pool <- list()
              for (i in seq_along(n_rankings)) {
                rank_tbl_pool$tbl <- rank_tbl_pool$tbl |>
                  mutate(
                    "ranking_{n_rankings[i]}" := dplyr::tibble(
                      id = seq_along(tau_hat),
                      tau_hat = tau_hat
                    ) |>
                      dplyr::arrange(.data$tau_hat) |>
                      dplyr::mutate(
                        id_2 = seq_along(.data$tau_hat),
                        rank = cut(
                          x = .data$id_2,
                          breaks = c(
                            seq(0, length(tau_hat) %% n_rankings[i], by = 1) *
                              (length(tau_hat) %/% n_rankings[i] + 1),
                            seq(
                              length(tau_hat) %% n_rankings[i] + 1,
                              n_rankings[i],
                              length.out = n_rankings[i] - length(tau_hat) %% n_rankings[i]
                            ) *
                              length(tau_hat) %/% n_rankings[i] + length(tau_hat) %% n_rankings[i]
                          ),
                          include.lowest = TRUE,
                          labels = seq_len(n_rankings[i])
                        )
                      ) |>
                      dplyr::arrange(.data$id) |>
                      dplyr::pull(.data$rank)
                  )
                heatmap_data_pool <- map(
                  list[[glue::glue("heatmap_data_{n_rankings[i]}")]],
                  \(tbl) {
                    map(
                      tbl,
                      \(var) {
                        attributes(var) <- NULL
                        return(var)
                      }
                    ) |>
                      as_tibble()
                  }
                ) |>
                  list_transpose(simplify = FALSE) |>
                  (\(x) x[c("covariate", "avg", "stderr", "scaling", "variation")])() |>
                  map(\(x) list_transpose(x, simplify = TRUE)) |>
                  list_transpose(simplify = FALSE) |>
                  (\(x) {
                    n <- length(x)
                    map(
                      x,
                      \(x, n) {
                        if (is.list(x) && length(x[["covariate"]]) > 1) {
                          pool <- list(
                            "avg" = pool.scalar(Q = x[["avg"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "stderr" = pool.scalar(Q = x[["stderr"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "scaling" = pool.scalar(Q = x[["scaling"]], U = 0, n = 1, k = 1, rule = "rubin1987"),
                            "variation" = pool.scalar(Q = x[["variation"]], U = 0, n = 1, k = 1, rule = "rubin1987")
                          )
                          out <- tibble(
                            "covariate" = x[["covariate"]][[1]],
                            "avg" = pool$avg$qbar,
                            "stderr" = pool$stderr$qbar,
                            "scaling" = pool$scaling$qbar,
                            "variation" = pool$variation$qbar
                          )
                        } else {
                          pool <- NULL
                          out <- tibble(
                            "covariate" = x[["covariate"]],
                            "avg" = x[["avg"]],
                            "stderr" = x[["stderr"]],
                            "scaling" = x[["scaling"]],
                            "variation" = x[["variation"]]
                          )
                        }
                        return(list(tbl = out, pool = pool))
                      },
                      n = n)
                  })()|>
                  list_transpose()
                heatmap_data_pool$tbl <- heatmap_data_pool$tbl |>
                  mutate(
                    labels = paste0(
                      sprintf("%.3f", avg),
                      "\n",
                      "(",
                      sprintf("%.3f", stderr),
                      ")"
                    )
                  )
                heatmap_data_pool$tbl <- heatmap_data_pool$tbl |>
                  mutate(
                    "ranking_{n_rankings[i]}" := rep(paste0("Q", seq_len(n_rankings[i])), n() / n_rankings[i])
                  )
                
                out_heatmap_tbl <- c(
                  out_heatmap_tbl,
                  list2("heatmap_data_{n_rankings[i]}" := heatmap_data_pool$tbl)
                )
                out_heatmap_pool <- c(
                  out_heatmap_pool,
                  list2("heatmap_data_{n_rankings[i]}" := heatmap_data_pool$pool)
                )
              }
              
              return(
                c(
                  list(rank_tbl = rank_tbl_pool$tbl), 
                  out_heatmap_tbl, 
                  list(
                    pool = c(
                      list(rank_tbl = rank_tbl_pool$pool),
                      out_heatmap_pool
                    )
                  )
                )
              )
            })()
          return(out)
        }
      )
    
    # pool data for excess risk distribution plots by covariate combinations ----
    cfw_agg$cate_histogram_data <- list_rbind(cfw_cate_histogram_data) |>
      select(-starts_with("X"))
    
    # pool data for average risk of thiazide-induced hyponatraemia in subgroups by covariate combinations ----
    cfw_agg$subgroup_cate <- cfw_subgroup_cate |>
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
                    for (j in seq_along(x)) {
                      out[[i]][j] <- x[[j]][i]
                    }
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
                      out <- mutate(out, "{var1}" := str_extract(x[[1]][1], ".+?(?=_)"))
                      out <- mutate(out, "{var2}" := str_extract(x[[1]][1], "(?<=_).+"))
                      out <- mutate(out, "{var1}_grp" := x[[paste0(var1, "_grp")]][1])
                      out <- mutate(out, "{var2}_grp" := x[[paste0(var2, "_grp")]][1])
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
    
    # close multisession connections
    plan(sequential)
    
    # save aggregated results
    saveRDS(cfw_agg, paste0(result_path, "cfw_agg.rds"))
    
    # remove un-pooled data
    rm(list = c(
      "cfw_covariates", "cfw_cate", "cfw_outcome_trt_grp", "cfw_vi", "cfw_cal", "cfw_overlap",
      "cfw_cate_plot", "cfw_balance_data", "cfw_rloss_error", "cfw_mse", "cfw_blp", "cfw_rate",
      "cfw_tau_hat_train", "cfw_cfb_cate", "cfw_cfb_control_risk", "cfw_mbcfb","cfw_ate",
      "cfw_risk_strategy", "cfw_bencalibr", "cfw_dyn"
    ))
    
    
    
    # save tables with results .................................................................................... ----
    # variable importance
    cfw_agg$variable_importance |>
      map(\(x) select(x$table, variable_name, variable_importance, variable_importance_num)) |>
      dfs2xlsx(paste0(table_path, "cfw_vi.xlsx"), rowNames = FALSE)
    cfw_agg$variable_importance_comb |>
      map(\(x) select(x, variable_name, variable_importance, variable_importance_num)) |>
      dfs2xlsx(paste0(table_path, "cfw_vi_comb.xlsx"), rowNames = FALSE)
    # calibration from R-learner
    cfw_agg$calibration|>
      map(\(x) {
        tbl <- x$table
        tibble(
          coefficient = rownames(tbl),
          estimate = tbl[, 1],
          std_err = tbl[, 2],
          t_value = tbl[, 3],
          p_value = tbl[, 4],
          df = attr(tbl, "df"),
          nobs = rep(attr(tbl, "nobs"), 2)
        )
      }) |>
      dfs2xlsx(paste0(table_path, "cfw_cal.xlsx"), rowNames = FALSE)
    # rloss-error
    cfw_agg$rloss_error |>
      map(
        \(x) x$table
      ) |>
      list_rbind() |>
      mutate(model = names(cfw_agg$rloss_error), .before = 1) |>
      mutate(relative_error = rloss_error / min(rloss_error)) |>
      (
        \(x) {
          dfs2xlsx(
            list(`rloss error` = x), 
            paste0(table_path, "cfw_rloss_error.xlsx"), 
            rowNames = FALSE
          )
        }
      )()
    # mean squared error
    cfw_agg$mean_square_error |>
      map(
        \(x) x$table
      ) |>
      list_rbind() |>
      mutate(model = names(cfw_agg$mean_square_error), .before = 1) |>
      mutate(relative_error = mean_square_error / min(mean_square_error)) |>
      (
        \(x) {
          dfs2xlsx(
            list(`mean square error` = x), 
            paste0(table_path, "cfw_mse.xlsx"), 
            rowNames = FALSE
          )
        }
      )()
    # best linear projection of CATEs
    cfw_agg$best_linear_projection |>
      map(\(x) x$table) |>
      dfs2xlsx(paste0(table_path, "cfw_blp.xlsx"), rowNames = FALSE)
    # RATE tests of presence of heterogeneity
    cfw_agg$rate |> 
      map(\(x) x$table) |>
      dfs2xlsx(paste0(table_path, "cfw_rate.xlsx"), rowNames = FALSE)
    # C-for-benefit
    list(
      `CATE matched` = cfw_agg$c_for_benefit_cate |>
        map(\(x) x$table) |>
        list_rbind() |>
        mutate(model = names(cfw_agg$c_for_benefit_cate), .before = 1),
      `Control risk matched` = cfw_agg$c_for_benefit_control_risk |>
        map(\(x) x$table) |>
        list_rbind() |>
        mutate(model = names(cfw_agg$c_for_benefit_control_risk), .before = 1),
      `Combined` = cfw_agg$c_for_benefit_combined |>
        map(\(x) x$table) |>
        list_rbind() |>
        mutate(model = names(cfw_agg$c_for_benefit_combined), .before = 1),
      `Model based` = cfw_agg$model_based_c_for_benefit |>
        map(\(x) x$table) |>
        list_rbind() |>
        mutate(model = names(cfw_agg$model_based_c_for_benefit), .before = 1)
    ) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_c_for_benefit.xlsx"
        ), 
        rowNames = FALSE
      )
    # ATE in subgroups
    cfw_agg$ate |> 
      map(\(x) {
        x$table |>
          unnest(cols = c(quartile_low, quartile_high), keep_empty = TRUE) |>
          mutate(name_quartile = paste0(row_number()), .by = name) |>
          pivot_wider(names_from = name_quartile, values_from = c(quartile_low, quartile_high), names_glue = "{.value}_{name_quartile}")
      }) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_ate.xlsx"
        ), 
        rowNames = FALSE
      )
    # Risk strategy table
    cfw_agg$risk_strategy |> 
      map(\(x) {
        x$table
      }) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_strategy.xlsx"
        ), 
        rowNames = FALSE
      )
    # Excess risk reduction
    list_transpose(cfw_agg$cfw_risk_reduction_plot)$aipw |>
      map(\(x) mutate(x$table, percentile = 1 - percentile)) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_reduction_aipw.xlsx"
        ), 
        rowNames = FALSE
      )
    
    list_transpose(cfw_agg$cfw_risk_reduction_plot)$ols |>
      map(\(x) mutate(x$table, percentile = 1 - percentile)) |>
      dfs2xlsx(
        paste0(
          table_path, 
          "cfw_risk_reduction_ols.xlsx"
        ), 
        rowNames = FALSE
      )
    
    # save figures with results ................................................................................... ----
    # Propensity score distribution to access overlap
    ggsave(
      paste0(figure_path, "overlap_plot.jpg"),
      cfw_agg$overlap +
        xlab("Probability of initiating thiazide vs. non-thiazide treatment") +
        ylab("Density") +
        ggtitle(""),
      "jpeg", NULL, 1.3, 10, 5, "cm", 240,
      quality = 90, create.dir = TRUE
    )
    ggsave(
      paste0(figure_path, "overlap_plot.pdf"),
      cfw_agg$overlap +
        xlab("Probability of initiating thiazide vs. non-thiazide treatment") +
        ylab("Density") +
        ggtitle(""),
      "pdf", NULL, 1.3, 10, 5, "cm", 240, create.dir = TRUE
    )
    # Covariate balance (Love plots)
    imap(
      cfw_agg$balance_data,
      \(x, nm) {
        plots <- MapCovariateBalance(
          data = x,
          plots = "Love",
          nm_ind = names_index,
          plots_only = TRUE,
          balance_table = FALSE
        )
        love_plot_comb <- plots[[1]]$love
        love_plot_comb$data <- reduce(
          plots,
          \(x1, x2) {
            list(
              love = list(
                data = bind_rows(
                  x1$love$data,
                  x2$love$data
                )
              )
            )
          }
        )[["love"]][["data"]] |>
          mutate(cov_val = sum((Cohort == "Before adjustment") * value), .by = "covariate_name") |>
          arrange(cov_val, value) |>
          select(-c("cov_val")) |>
          mutate(covariate_name = factor(covariate_name, levels = unique(covariate_name))) |>
          filter(str_detect(covariate_name, " - ", negate = TRUE))
        nrow_data <- nrow(love_plot_comb$data)
        ggsave(paste0(figure_path, nm, "/love_plot.jpg"),
               love_plot_comb, "jpeg", NULL, 1.8, 16, 4 + ceiling(nrow_data / 9), "cm", 240,
               quality = 90, create.dir = TRUE)
        ggsave(paste0(figure_path, nm, "/love_plot.pdf"),
               love_plot_comb, "pdf", NULL, 1.8, 16, 4 + ceiling(nrow_data / 9), "cm", 240,
               create.dir = TRUE)
        return(NULL)
      }
    )
    
    # Distribution of CATEs
    imap(
      cfw_agg$cate,
      \(x, nm) {
        cate_table <- x$table
        cate_table <- cate_table |> 
          arrange(estimate) |>
          mutate(
            quantiles = cut(
              estimate, 
              breaks = quantile(
                estimate, 
                c(0, 0.5, 0.9, 0.95, 1)
              ),
              labels = c("0%-50%", "50%-90%", "90%-95%", "95%-100%"),
              include.lowest = TRUE
            ),
            order = seq_len(n())
          )
        plot_a <- ggplot(cate_table) +
          geom_histogram(aes(x = estimate, fill = quantiles), binwidth = 0.002) + 
          scale_x_continuous(breaks = seq(-0.025, 0.15, 0.025),
                             labels = c("-2.5%", "0%", "2.5%", "5%", "7.5%", "10%", "12.5%", "15%")) +
          scale_fill_jama() +
          coord_cartesian(xlim = c(-0.03, 0.16)) +
          xlab("Excess risk") + 
          ylab("Number of thiazide users") +
          theme(legend.title = element_blank())
        plot_a_trt <- ggplot(cate_table |> filter(W == 1)) +
          geom_histogram(aes(x = estimate, fill = quantiles), binwidth = 0.002) + 
          scale_x_continuous(breaks = seq(-0.025, 0.15, 0.025),
                             labels = c("-2.5%", "0%", "2.5%", "5%", "7.5%", "10%", "12.5%", "15%")) +
          scale_fill_jama() +
          coord_cartesian(xlim = c(-0.03, 0.16)) +
          xlab("Excess risk") + 
          ylab("Number of thiazide users") +
          theme(legend.title = element_blank())
        plot_b <- ggplot(cate_table) +
          geom_line(aes(x = order, y = estimate, color = quantiles), linewidth = 1.1) + 
          scale_y_continuous(breaks = seq(-0.025, 0.15, 0.025),
                             labels = c("-2.5%", "0%", "2.5%", "5%", "7.5%", "10%", "12.5%", "15%")) +
          scale_x_continuous(breaks = round(nrow(cate_table) * c(0, 0.25, 0.5, 0.75, 1))) +
          scale_color_jama() +
          coord_cartesian(ylim = c(-0.03, 0.16)) +
          xlab("Thiazide users ranked by excess risk") + 
          ylab("Excess risk") +
          theme(legend.title = element_blank())
        plot_b_trt <- ggplot(cate_table |> filter(W == 1) |> mutate(order = seq_len(n()))) +
          geom_line(aes(x = order, y = estimate, color = quantiles), linewidth = 1.1) + 
          scale_y_continuous(breaks = seq(-0.025, 0.15, 0.025),
                             labels = c("-2.5%", "0%", "2.5%", "5%", "7.5%", "10%", "12.5%", "15%")) +
          scale_x_continuous(breaks = round(nrow(cate_table |> filter(W == 1)) * c(0, 0.25, 0.5, 0.75, 1))) +
          scale_color_jama() +
          coord_cartesian(ylim = c(-0.03, 0.16)) +
          xlab("Thiazide users ranked by excess risk") + 
          ylab("Excess risk") +
          theme(legend.title = element_blank())
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/cate_density_plot", ".jpg"),
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
          filename = paste0(figure_path, nm, "/cate_density_plot_treated", ".jpg"),
          plot = plot_a_trt,
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
          filename = paste0(figure_path, nm, "/cate_rank_plot", ".jpg"),
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
          filename = paste0(figure_path, nm, "/cate_rank_plot_treated", ".jpg"),
          plot = plot_b_trt,
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
    
    # excess risk reduction plots
    pmap(
      list(
        cfw_agg$cfw_risk_reduction_plot,
        cfw_agg$cfw_risk_reduction_plot_zp,
        names(cfw_agg$cfw_risk_reduction_plot)
      ),
      \(x, zp, nm) {
        plot_aipw <- x$aipw$table |>
          slice(2:n()) |>
          ggplot() + 
          geom_linerange(
            aes(ymin = ymin, ymax = ymax, x = x),
            data = tibble(
              ymin = -Inf,
              ymax = x$aipw$table |> filter(percentile == (1 - 0.1)) |> pull(estimate),
              x = 100 * 0.1
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.7
          ) +
          geom_linerange(
            aes(xmin = xmin, xmax = xmax, y = y),
            data = tibble(
              xmin = -Inf,
              xmax = 100 * 0.1,
              y = x$aipw$table |> filter(percentile == (1 - 0.1)) |> pull(estimate)
            ),
            color = "red3",
            alpha = 0.8,
            linewidth = 0.7
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate), linewidth = 0.9) +
          geom_hline(yintercept = x$aipw$table$estimate[1], linetype = 2, linewidth = 0.7) +
          geom_vline(xintercept = 100 * (1 - zp$quantile), linetype = 2, linewidth = 0.7) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Percent not treated with thiazide") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(x$aipw$table$estimate))) +
          scale_x_continuous(breaks = seq(0, 100, 20))
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
            linewidth = 0.7
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
            linewidth = 0.7
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate), linewidth = 0.9) +
          geom_hline(yintercept = x$ols$table$estimate[1], linetype = 2, linewidth = 0.7) +
          geom_vline(xintercept = 100 * (1 - zp$quantile), linetype = 2, linewidth = 0.7) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Percent not treated with thiazide") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(x$ols$table$estimate))) +
          scale_x_continuous(breaks = seq(0, 100, 20))
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
            linewidth = 0.7,
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
            linewidth = 0.7,
            alpha = 0.8
          ) +
          geom_line(aes(x = 100 * (1 - percentile), y = estimate, color = est), linewidth = 0.9) +
          geom_hline(yintercept = x$aipw$table$estimate[1], linetype = 2, linewidth = 0.7) +
          geom_vline(xintercept = 100 * (1 - zp$quantile), linetype = 2, linewidth = 0.7) +
          ylab("Population reduction in hyponatraemia risk") + xlab("Percent not treated with thiazide") +
          scale_color_jama() +
          labs(color = "Estimator") +
          coord_cartesian(xlim = c(0, 100), ylim = c(0, max(c(x$aipw$table$estimate, x$ols$table$estimate)))) +
          scale_x_continuous(breaks = seq(0, 100, 20))
        ggplot2::ggsave(
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_aipw", ".jpg"),
          plot = plot_aipw,
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
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_ols", ".jpg"),
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
          filename = paste0(figure_path, nm, "/excess_risk_reduction_plot_comb", ".jpg"),
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
    
    return(NULL)
  }
)

## forest results ------------------------------------------------------ --------
### dynamic subgroups ------------------------------------------------------------
plan(sequential)
future_pmap(
  .l = list(
    W = list(
      tih_tbi_1_cohort_imputed[[1]]$W,
      tih_tbi_1_cohort_imputed[[1]]$W,
      tih_cohort$W,
      tih_tbi_1_cohort_imputed[[1]]$W, 
      tih_tbi_1_cohort_imputed[[1]]$W,
      bfz_tbi_1_cohort_imputed[[1]]$W_bvc,
      bfz_tbi_1_cohort_imputed[[1]]$W_bvc,
      hctz_tbi_1_cohort_imputed[[1]]$W_hvr,
      hctz_tbi_1_cohort_imputed[[1]]$W_hvr,
      tih_cohort_complete_case$W,
      tih_cohort_cc_comp[[1]]$W
    ),
    Y = list(
      tih_tbi_1_cohort_imputed[[1]]$Y_130_120,
      tih_tbi_1_cohort_imputed[[1]]$Y_125_90,
      tih_cohort$Y_130_120,
      tih_tbi_1_cohort_imputed[[1]]$Y_130_120, 
      tih_tbi_1_cohort_imputed[[1]]$Y_130_120,
      bfz_tbi_1_cohort_imputed[[1]]$Y_130_120,
      bfz_tbi_1_cohort_imputed[[1]]$Y_130_120,
      hctz_tbi_1_cohort_imputed[[1]]$Y_130_120,
      hctz_tbi_1_cohort_imputed[[1]]$Y_130_120,
      tih_cohort_complete_case$Y_130_120,
      tih_cohort_cc_comp[[1]]$Y_130_120
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
    )
  ) |> map(\(x) x[setup_index]),
  .options = furrr_options(
    packages = c(
      "dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid",
      "future", "furrr", "future.callr", "callr", "MatchIt"
    ),
    globals = c("n_rankings", "vcovHC", "horizon", "thiazide_limit"),
    seed = TRUE
  ),
  .f = \(W, Y, path, result_path, table_path, figure_path) {
    cfw_agg <- readRDS(paste0(result_path, "cfw_agg.rds"))
    plan(multisession, workers = min(num_clusters, length(cfw_agg$dyn)))
    cfw_dsg <- future_pmap(
      .l = list(
        agg = cfw_agg$dyn,
        nm = names(cfw_agg$dyn)
      ),
      .options = furrr_options(
        packages = c("dplyr", "purrr", "grf", "ggplot2", "ggsci", "ggpubr", "grid", "MatchIt"),
        globals = c("W", "Y", "n_rankings", "vcovHC", "horizon", "thiazide_limit"),
        seed = TRUE
      ),
      .f = \(agg, nm) {
        # add treatment indicator to dynamic subgroup data ----
        data_tbl <- agg$rank_tbl |>
          mutate(
            W = W[.data$id],
            Y = Y[.data$id]
          )
        out <- list()
        for (i in n_rankings) {
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
                ggplot2::theme_bw() +
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
                  W ~ tau_hat,
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
              arrange(!!sym(ranking), subclass),
            control_risk = map(
              levels(data_tbl[[ranking]]),
              \(lvl) {
                matched <- MatchIt::matchit(
                  W ~ mu_hat_0,
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
              mutate(Y_adj = Y - mu_hat_0) |>
              arrange(!!sym(ranking), subclass),
            combined = map(
              levels(data_tbl[[ranking]]),
              \(lvl) {
                matched_trt <- MatchIt::matchit(
                  W ~ mu_hat_0,
                  data = filter(data_tbl, (!!sym(ranking) == lvl) | W == 0),
                  method = "nearest",
                  distance = "mahalanobis",
                  estimand = "ATT",
                  replace = TRUE
                )
                matched_ctr <- MatchIt::matchit(
                  W ~ mu_hat_1,
                  data = filter(data_tbl, (!!sym(ranking) == lvl) | W == 1),
                  method = "nearest",
                  distance = "mahalanobis",
                  estimand = "ATC",
                  replace = TRUE
                )
                matched_patients_trt <- filter(data_tbl, (!!sym(ranking) == lvl) | W == 0) |>
                  filter(matched_trt$treat == 1) |>
                  mutate(subclass = seq_len(n())) |>
                  bind_rows(
                    filter(data_tbl, (!!sym(ranking) == lvl) | W == 0) |>
                      slice(as.integer(matched_trt$match.matrix[, 1])) |>
                      mutate(subclass = seq_len(n()))
                  ) |>
                  arrange(subclass, W) |>
                  mutate(
                    tau_hat = sum((W == 1) * tau_hat),
                    aipw_scores = sum((W == 1) * aipw_scores),
                    weights = sum((W == 1) * sample_weights),
                    .by = subclass
                  ) |>
                  mutate(Y_adj = Y - mu_hat_0)
                matched_patients_ctr <- filter(data_tbl, (!!sym(ranking) == lvl) | W == 1) |>
                  filter(matched_ctr$treat == 0) |>
                  mutate(subclass = seq_len(n()) + sum(filter(data_tbl, !!sym(ranking) == lvl)$W == 1)) |>
                  bind_rows(
                    filter(data_tbl, (!!sym(ranking) == lvl) | W == 1) |>
                      slice(as.integer(matched_ctr$match.matrix[, 1])) |>
                      mutate(subclass = seq_len(n()) + sum(filter(data_tbl, !!sym(ranking) == lvl)$W == 1))
                  ) |>
                  arrange(subclass, W) |>
                  mutate(
                    tau_hat = sum((W == 0) * tau_hat),
                    aipw_scores = sum((W == 0) * aipw_scores),
                    weights = sum((W == 0) * sample_weights),
                    .by = subclass
                  ) |>
                  mutate(Y_adj = Y - mu_hat_1)
                matched_patients <- bind_rows(
                  matched_patients_trt,
                  matched_patients_ctr
                )
                return(matched_patients)
              }
            ) |>
              list_rbind() |>
              arrange(!!sym(ranking), subclass)
          )
          
          # Mean bias
          # expected difference between observed effect and predicted effect in matched sample
          mean_bias <- map(
            data_tbl_matched,
            \(tbl) {
              tbl |> 
                summarise(
                  observed_treatment_effect = 
                    sum((W == 1) * weights * Y) / sum((W == 1) * weights) - 
                    sum((W == 0) * weights * Y) / sum((W == 0) * weights),
                  expected_treatment_effect = sum(
                    (W == 1) * weights * (mu_hat_1) / sum((W == 1) * weights) -
                      (W == 0) * weights * (mu_hat_0) / sum((W == 0) * weights)
                  ),
                  mean_bias = observed_treatment_effect - expected_treatment_effect
                )
            }
          )
          
          # calculate adjusted expected effect (expected outcome under treatment of treated minus expected outcome under control for non-treated)
          # Partition AIPW score into components for treated and untreated individual in each pair:
          # Treated: mu_1(X) + e(X)^-1 * (Y - mu_1(X)) 
          # Untreated: mu_0(X) + (1 - e(X))^-1 * (Y - mu_0(X)) 
          # expected_treatment_effect <- data_tbl_matched |>
          #   summarise(
          #     expected_treatment_effect = sum(
          #       W *       (mu_hat_1 ) - #+ W_hat^(-1) * (Y - mu_hat_1)
          #       (1 - W) * (mu_hat_0 )   #+ (1 - W_hat)^(-1) * (Y - mu_hat_0)
          #       ),
          #     .by = c("ranking", "subclass")
          #   ) |>
          #   summarise(
          #     expected_treatment_effect_se = sd(expected_treatment_effect),
          #     expected_treatment_effect = mean(expected_treatment_effect),
          #     .by = "ranking"
          #   )
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
                  sum((W == 1) * weights * mu_hat_1 / sum((W == 1) * weights)) -
                  sum((W == 0) * weights * mu_hat_0 / sum((W == 0) * weights)),
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
                  sum((W == 1) * weights * mu_hat_1 / sum((W == 1) * weights)) -
                  sum((W == 0) * weights * mu_hat_0 / sum((W == 0) * weights)),
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
          # observed_treatment_effect <- data_tbl_matched |> 
          #   summarise(
          #     observed_treatment_effect = sum(W * Y - (1 - W) * Y),
          #     .by = c("ranking", "subclass")
          #   ) |>
          #   summarise(
          #     observed_treatment_effect_se = sd(observed_treatment_effect),
          #     observed_treatment_effect = mean(observed_treatment_effect),
          #     .by = c("ranking")
          #   )
          # 
          # CATE matched:
          # The observed effect is calculated as the mean difference between the two treatment 
          # arms in the matched sample. The standard error is based on Y having a bernoulli distribution,
          # such that sum((W==1) * Y) is binomially distributed with estimated parameter 
          # sum((W == 1) * Y_obs) / sum(W == 1).
          # Then the variance is sum(W == 1) * sum((W == 1) * Y) / sum(W == 1) * (sum((W == 1) * (1 - Y)) / sum(W == 1))
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
                  sum((W == 1) * weights * Y) / sum((W == 1) * weights) - 
                  sum((W == 0) * weights * Y) / sum((W == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W == 1) * weights * (Y == 1)) * sum((W == 1) * weights * (Y == 0)) / sum((W == 1) * weights)^3 +
                    sum((W == 0) * weights * (Y == 1)) * sum((W == 0) * weights * (Y == 0)) / sum((W == 0) * weights)^3
                ),
                .by = all_of(ranking)
              ),
            control_risk = data_tbl_matched$control_risk |>
              summarise(
                observed_treatment_effect = 
                  sum((W == 1) * weights * Y_adj) / sum((W == 1) * weights) - 
                  sum((W == 0) * weights * Y_adj) / sum((W == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W == 1) * weights * (Y == 1)) * sum((W == 1) * weights * (Y == 0)) / sum((W == 1) * weights)^3 +
                    sum((W == 0) * weights * (Y == 1)) * sum((W == 0) * weights * (Y == 0)) / sum((W == 0) * weights)^3
                ),
                .by = all_of(ranking)
              ),
            combined = data_tbl_matched$combined |>
              summarise(
                observed_treatment_effect = 
                  sum((W == 1) * weights * Y_adj) / sum((W == 1) * weights) - 
                  sum((W == 0) * weights * Y_adj) / sum((W == 0) * weights),
                observed_treatment_effect_se = sqrt(
                  sum((W == 1) * weights * (Y == 1)) * sum((W == 1) * weights * (Y == 0)) / sum((W == 1) * weights)^3 +
                    sum((W == 0) * weights * (Y == 1)) * sum((W == 0) * weights * (Y == 0)) / sum((W == 0) * weights)^3
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
                W ~ mu_hat_0,
                data = data_tbl,
                method = "nearest",
                distance = "mahalanobis",
                estimand = "ATT",
                replace = TRUE
              )
              matched_ctr <- MatchIt::matchit(
                W ~ mu_hat_1,
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
                arrange(subclass, W) |>
                mutate(
                  tau_hat = sum((W == 1) * tau_hat),
                  aipw_scores = sum((W == 1) * aipw_scores),
                  weights = sum((W == 1) * sample_weights),
                  .by = subclass
                ) |>
                mutate(Y_adj = Y - mu_hat_0)
              matched_patients_ctr <- data_tbl |>
                filter(matched_ctr$treat == 0) |>
                mutate(subclass = seq_len(n()) + sum(data_tbl$W == 1)) |>
                bind_rows(
                  data_tbl |>
                    slice(as.integer(matched_ctr$match.matrix[, 1])) |>
                    mutate(subclass = seq_len(n()) + sum(data_tbl$W == 1))
                ) |>
                arrange(subclass, W) |>
                mutate(
                  tau_hat = sum((W == 0) * tau_hat),
                  aipw_scores = sum((W == 0) * aipw_scores),
                  weights = sum((W == 0) * sample_weights),
                  .by = subclass
                ) |>
                mutate(Y_adj = Y - mu_hat_1)
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
              observed_treatment_effect = sum((W == 1) * Y_adj - (W == 0) * Y_adj),
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
                geom_hex(binwidth = c(0.005, 0.04)) +
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
                  breaks = seq(-0.05, 0.2, 0.05),
                  labels = paste0(seq(-5, 20, 5), "%")
                ) +
                scale_y_continuous(
                  breaks = seq(-1, 1, 0.2),
                  labels = paste0(seq(-100, 100, 20), "%")
                ) +
                xlab("Predicted excess risk") + ylab("Observed excess risk") +
                coord_cartesian(xlim = c(-0.05, 0.2), ylim = c(-1, 1))
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
        cfw_dsg, # map over model, e.g. full, vi_1, etc.
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
      paste0(result_path, "cfw_dsg.rds")
    )
    
    # save calibration plots using dynamic subgroups
    cal_plot_combined_aipw_comb <- imap(
      cfw_dsg,
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
                    paste0("\u03B2\u2080 = ", sprintf('%.3f', coef_comb_ind$est[1]), "\n\u03B2\u2081 = ", sprintf('%.3f', coef_comb_ind$est[2]))
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
                    paste0("\u03B2\u2080 = ", sprintf('%.3f', coef_comb_ind$est[1]), "\n\u03B2\u2081 = ", sprintf('%.3f', coef_comb_ind$est[2]))
                  )
                )
              )
            ggplot2::ggsave(
              filename = paste0(figure_path, nm, "/cal_ate_plot_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/ols/cal_plot_cate_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/aipw/cal_plot_cate_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/ols/cal_plot_control_risk_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/aipw/cal_plot_control_risk_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/ols/cal_plot_combined_", rank, ".jpg"),
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
              filename = paste0(figure_path, nm, "/aipw/cal_plot_combined_", rank, ".jpg"),
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
    
    # save combined figure with calibration plots from all seven models
    dsg_plot_comb_5 <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb, 
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
    dsg_plot_comb_10 <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb, 
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
    dsg_plot_comb_20 <- cowplot::plot_grid(
      plotlist = imap(
        cal_plot_combined_aipw_comb, 
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
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_5.jpg"),
      plot = dsg_plot_comb_5,
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
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_10.jpg"),
      plot = dsg_plot_comb_10,
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
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_10.pdf"),
      plot = dsg_plot_comb_10,
      device = Cairo::CairoPDF,
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      family = "Arial Unicode MS",
      create.dir = TRUE
    )
    ggplot2::ggsave(
      filename = paste0(figure_path, "combined/cal_plot_aipw_combined_20.jpg"),
      plot = dsg_plot_comb_20,
      device = "jpeg",
      width = 12,
      height = 24,
      unit = "cm",
      scale = 1.7,
      dpi = 240,
      quality = 90,
      create.dir = TRUE
    )
    
    # save tables with coefficients from calibration plots
    cfw_dsg |> 
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
              nrow(cfw_dsg$full_mod$individual$calibration_data$combined)
            )
          ) |>
          arrange(num_ranks, exp_type, model)
      }) |>
      dfs2xlsx(paste0(table_path, "cfw_dyn_cal_coef.xlsx"), rowNames = FALSE)
    
    return(invisible(NULL))
  }
)

### gold standard comparison sensitivity analysis ----------------------------------------------------------------- ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_thiazide_all/1year/"),
      paste0(training_path, "4mo_thiazide_no_lab/"),
      paste0(training_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_thiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_thiazide_complete_case/1year/"),
      paste0(training_path, "4mo_thiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_thiazide_all/1year/"),
      paste0(results_path, "4mo_thiazide_no_lab/"),
      paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_thiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_thiazide_complete_case/1year/"),
      paste0(results_path, "4mo_thiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_thiazide_all/1year/"),
      paste0(tables_path, "4mo_thiazide_no_lab/"),
      paste0(tables_path, "4mo_thiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_thiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_thiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_thiazide_cc_comp/1year/")
    ),
    figure_path = c(
      paste0(figures_path, "4mo_thiazide_all/1year/"),
      paste0(figures_path, "4mo_thiazide_no_lab/"),
      paste0(figures_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_thiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_thiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_thiazide_cc_comp/1year/")
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
    cfw_agg <- readRDS(paste0(result_path, "cfw_agg.rds"))
    plan(multisession, workers = min(num_clusters, length(cfw_agg$cate)))
    cfw_gs_sens <- future_pmap(
      .l = list(
        gold_std = cfw_agg$cate, # gold standard model
        comp_list = map(seq_along(cfw_agg$cate), \(i) cfw_agg$cate[-i]), # list of comparison models
        gs_nm = names(cfw_agg$cate) # name of gold standard model
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
              filename = paste0(figure_path, "/", gs_nm, "/gold_standard/", comp_nm, "/cate_scatter_plot", ".jpg"),
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
                # What if point probability at the cutoff?
                comp_idx <- which(comp$table$estimate > comp_cutoff)
                comp_eidx <- which(comp$table$estimate == comp_cutoff)
                comp_nidx <- which(comp$table$estimate < comp_cutoff)
                
                eidx_split <- if (quantile > 0) {
                  ceiling(length(comp$table$estimate) * (1 - quantile)) - length(comp_idx)
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
    saveRDS(cfw_gs_sens, paste0(result_path, "cfw_gs_sens.rds"))
    
    # save tables with gold standard model comparison results
    cfw_gs_sens |>
      dfs2xlsx(paste0(table_path, "cfw_gs_sens.xlsx"), rowNames = FALSE)
    
    return(invisible(NULL))
  }
)

### Excess risk distribution plots by covariate combinations ------------------------------------------------------ ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_thiazide_all/1year/"),
      paste0(training_path, "4mo_thiazide_no_lab/"),
      paste0(training_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_thiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_thiazide_complete_case/1year/"),
      paste0(training_path, "4mo_thiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_thiazide_all/1year/"),
      paste0(results_path, "4mo_thiazide_no_lab/"),
      paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_thiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_thiazide_complete_case/1year/"),
      paste0(results_path, "4mo_thiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_thiazide_all/1year/"),
      paste0(tables_path, "4mo_thiazide_no_lab/"),
      paste0(tables_path, "4mo_thiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_thiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_thiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_thiazide_cc_comp/1year/")
    ),
    figure_path = c(
      paste0(figures_path, "4mo_thiazide_all/1year/"),
      paste0(figures_path, "4mo_thiazide_no_lab/"),
      paste0(figures_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_thiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_thiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_thiazide_cc_comp/1year/")
    ),
    cutpoints = list(
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide no lab
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide only fluid lab results
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide only fluid lab results and no Denmark specific covariates
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
      ), # thiazide complete case
      list(
        Age = c(55, 70),
        `Years of education` = list("<10", "10-12", c("13-15", ">15")),
        `Days of hospitalization in the past year prior to index data` = list("0", "1-7", c("8-14", "\u226515")),
        Sodium = c(135, 138, 141)
      ) # thiazide complete case comparison
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
  .f = \(path, result_path, table_path, figure_path, cutpoints) {
    cfw_agg <- readRDS(paste0(result_path, "cfw_agg.rds"))
    # obtain data for plot labels
    ate_all <- cfw_agg$ate[[plot_model]]$table |> filter(name == "Full population") |> pull(estimate)
    label_data_oneway <- map(
      as.list(names(cutpoints)),
      \(var) {
        cfw_agg$cate_histogram_data |>
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
        cfw_agg$cate_histogram_data |>
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
        round(100 * max(cfw_agg$cate_histogram_data$predictions) * 5) / 5 - 
          round(100 * min(cfw_agg$cate_histogram_data$predictions) * 5) / 5
      ) * 5
    ) + 
      1 +
      if (sign(min(label_data_oneway[[1]]$x_l) - 100 * min(cfw_agg$cate_histogram_data$predictions)) == -1) {
        round(
          (
            round(100 * min(cfw_agg$cate_histogram_data$predictions) * 5) / 5 - 
              round(min(label_data_oneway[[1]]$x_l) * 5) / 5
          ) * 5
        )
      } else {
        0
      }
    map(
      as.list(names(cutpoints)),
      \(var) {
        plot <- cfw_agg$cate_histogram_data |>
          arrange(!!sym(var)) |>
          ggplot() +
          geom_histogram(aes(x = 100 * predictions, 
                             y = after_stat(count) *
                               5 / 
                               rep(
                                 label_data_oneway[[var[1]]]$count,
                                 each = 2 * n_bins
                               ), 
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
          filename = paste0(figure_path, "vi_mod_4/density_", gsub(" ", "", var[1]), "_plot", ".jpg"),
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
        plot <- cfw_agg$cate_histogram_data |>
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
                                   if (summarise(cfw_agg$cate_histogram_data, sum(!!sym(vars[1]) == levels(!!sym(vars[1]))[i] & !!sym(vars[2]) == levels(!!sym(vars[2]))[j] & ate_diff_signif))[[1]] == 0) {
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
          filename = paste0(figure_path, "vi_mod_4/density_", gsub(" ", "", vars[1]), "_", gsub(" ", "", vars[2]), "_plot", ".jpg"),
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
    
    return(invisible(NULL))
  }
)

### plots with average risk of thiazide-induced hyponatraemia in subgroups by covariate combinations ---------------- ----
plan(sequential)
future_pmap(
  .l = list(
    path = list(
      paste0(training_path, "4mo_thiazide_all/1year/"),
      paste0(training_path, "4mo_thiazide_no_lab/"),
      paste0(training_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(training_path, "4mo_thiazide_international/1year/"),
      paste0(training_path, "4mo_bfz_all/1year/"),
      paste0(training_path, "4mo_bfz_international/1year/"),
      paste0(training_path, "4mo_hctz_all/1year/"),
      paste0(training_path, "4mo_hctz_international/1year/"),
      paste0(training_path, "4mo_thiazide_complete_case/1year/"),
      paste0(training_path, "4mo_thiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_thiazide_all/1year/"),
      paste0(results_path, "4mo_thiazide_no_lab/"),
      paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_thiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_thiazide_complete_case/1year/"),
      paste0(results_path, "4mo_thiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_thiazide_all/1year/"),
      paste0(tables_path, "4mo_thiazide_no_lab/"),
      paste0(tables_path, "4mo_thiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_thiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_thiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_thiazide_cc_comp/1year/")
    ),
    figure_path = c(
      paste0(figures_path, "4mo_thiazide_all/1year/"),
      paste0(figures_path, "4mo_thiazide_no_lab/"),
      paste0(figures_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_thiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_thiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_thiazide_cc_comp/1year/")
    ),
    cutpoints = list(
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide all
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide no lab
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide only fluid lab results
      list(
        Age = c(55, 70),
        Sodium = c(135, 138, 141),
        Hemoglobin = c(8, 9, 10),
        `C-reactive protein` = c(5, 20)
      ), # thiazide only fluid lab results and no Denmark specific covariates
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
      ), # thiazide complete case
      list(
        Age = c(55, 70),
        `Years of education` = c(1.5, 2.5),
        `Days of hospitalization in the past year prior to index data` = c(1.5, 2.5),
        Sodium = c(135, 138, 141)
      ) # thiazide complete case comparison
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
  .f = \(path, result_path, table_path, figure_path, cutpoints) {
    cfw_agg <- readRDS(paste0(result_path, "cfw_agg.rds"))
    ate_all <- cfw_agg$ate$vi_mod_4$table |> filter(name == "Full population") |> pull(estimate)
    map(
      as.list(names(cutpoints)),
      \(var) {
        plot <- cfw_agg$subgroup_cate$oneway[[var]]$table |>
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
            breaks = cfw_agg$subgroup_cate$oneway[[var]]$table[[paste0(var, "_grp")]],
            labels = str_extract(cfw_agg$subgroup_cate$oneway[[var]]$table[[var]], ".+?(?=\\s[a-zA-Z])"),
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
          filename = paste0(figure_path, "vi_mod_4/subgroup_cate_", gsub(" ", "", var[1]), "_plot", ".jpg"),
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
        plot <- cfw_agg$subgroup_cate$twoway[[paste0(vars[1], "_", vars[2])]]$table |>
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
                     size = 1) +
          geom_hline(yintercept = 100 * ate_all, linetype = 2, linewidth = 0.3) +
          xlab(glue("{vars[1]} ({if(vars[1] %in% c('Age', 'Years of education')) 'years' else if (vars[1] == 'C-reactive protein') 'mg/L' else if (vars[1] == 'Lactate dehydrogenase') 'U/L' else if (vars[1] == 'Thrombocytes') 'x 10\u2079/L' else if (vars[1] == 'Sodium') 'mmol/L' else 'days'})")) +
          ylab("Excess risk") +
          scale_x_continuous(
            breaks = cfw_agg$subgroup_cate$twoway[[paste0(vars[1], "_", vars[2])]]$table[[paste0(vars[1], "_grp")]],
            labels = str_extract(cfw_agg$subgroup_cate$twoway[[paste0(vars[1], "_", vars[2])]]$table[[vars[1]]], ".+?(?=\\s[a-zA-Z])"),
            guide = guide_axis(angle = 0)
          ) +
          scale_y_continuous(
            breaks = seq(-5, 20, 5),
            labels = paste0(seq(-5, 20, 5), "%")
          ) +
          coord_cartesian(xlim = if(vars[1] == 'Age') c(45, 80) else if (vars[1] == 'Sodium') c(131, 144) else if (vars[1] == 'Hemoglobin') c(7, 11) else if (vars[1] == 'Lactate dehydrogenase') c(130, 290) else if (vars[1] == "Thrombocytes") c(155, 390) else if (vars[1] == "C-reactive protein") c(0, 40) else if (vars[1] == "Years of education") c(4, 18) else c(-2, 13), ylim = c(-5, 15)) +
          scale_color_jama() +
          guides(color = guide_legend(nrow = guide_nrow, ncol = guide_ncol, byrow = TRUE)) +
          theme(
            legend.position = "top",
            text = element_text(family = "sans", size = 7.5),
            axis.line = element_line(linewidth = 0.25),
            axis.ticks = element_line(linewidth = 0.25)
          )
        
        ggplot2::ggsave(
          filename = paste0(figure_path, "vi_mod_4/subgroup_cate_", gsub(" ", "", vars[1]), "_", gsub(" ", "", vars[2]), "_plot", ".jpg"),
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
    ) |>
      list_transpose()
    
    subgroup_plot_comb <- cowplot::plot_grid(
      plotlist = plot_list$line,
      nrow = 3, 
      ncol = 2,
      labels = LETTERS[1:6],
      label_size = 10
    )
    
    ggplot2::ggsave(
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot", ".jpg"),
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
      filename = paste0(figure_path, "vi_mod_4/subgroup_cate_comb_plot", ".pdf"),
      plot = subgroup_plot_comb,
      device = Cairo::CairoPDF,
      width = 18,
      height = 21,
      unit = "cm",
      scale = 1,
      family = "Arial Unicode MS"
    )
    
    return(invisible(NULL))
  }
)

### association between sex and hyponatraemia ---------------------------------------------------------------------------
plan(sequential)
future_pmap(
  .l = list(
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
      paste0(training_path, "4mo_thiazide_complete_case/1year/"),
      paste0(training_path, "4mo_thiazide_cc_comp/1year/")
    ),
    result_path = list(
      paste0(results_path, "4mo_thiazide_all/1year/"),
      paste0(results_path, "3mo_tiazide_all_severe/1year/"),
      paste0(results_path, "4mo_thiazide_no_lab/"),
      paste0(results_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(results_path, "4mo_thiazide_international/1year/"),
      paste0(results_path, "4mo_bfz_all/1year/"),
      paste0(results_path, "4mo_bfz_international/1year/"),
      paste0(results_path, "4mo_hctz_all/1year/"),
      paste0(results_path, "4mo_hctz_international/1year/"),
      paste0(results_path, "4mo_thiazide_complete_case/1year/"),
      paste0(results_path, "4mo_thiazide_cc_comp/1year/")
    ),
    table_path = list(
      paste0(tables_path, "4mo_thiazide_all/1year/"),
      paste0(tables_path, "3mo_tiazide_all_severe/1year/"),
      paste0(tables_path, "4mo_thiazide_no_lab/"),
      paste0(tables_path, "4mo_thiazide_lab_fluid/1year"),
      paste0(tables_path, "4mo_thiazide_international/1year/"),
      paste0(tables_path, "4mo_bfz_all/1year/"),
      paste0(tables_path, "4mo_bfz_international/1year/"),
      paste0(tables_path, "4mo_hctz_all/1year/"),
      paste0(tables_path, "4mo_hctz_international/1year/"),
      paste0(tables_path, "4mo_thiazide_complete_case/1year/"),
      paste0(tables_path, "4mo_thiazide_cc_comp/1year/")
    ),
    figure_path = c(
      paste0(figures_path, "4mo_thiazide_all/1year/"),
      paste0(figures_path, "3mo_tiazide_all_severe/1year/"),
      paste0(figures_path, "4mo_thiazide_no_lab/"),
      paste0(figures_path, "4mo_thiazide_lab_fluid/1year/"),
      paste0(figures_path, "4mo_thiazide_international/1year/"),
      paste0(figures_path, "4mo_bfz_all/1year/"),
      paste0(figures_path, "4mo_bfz_international/1year/"),
      paste0(figures_path, "4mo_hctz_all/1year/"),
      paste0(figures_path, "4mo_hctz_international/1year/"),
      paste0(figures_path, "4mo_thiazide_complete_case/1year/"),
      paste0(figures_path, "4mo_thiazide_cc_comp/1year/")
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
    cfw_agg <- readRDS(paste0(result_path, "cfw_agg.rds"))
    # plot distibution of CATEs by sex
    cfw_cate_by_sex_plot <- tibble(
      cate = cfw_agg$cate$full_mod$table$estimate,
      sex = factor(tih_cohort$X_02, levels = c(0, 1), labels = c("female", "male"))
    ) |>
      ggplot(aes(x = cate, y = after_stat(density), color = sex)) +
      geom_freqpoly(binwidth = 0.002, alpha = 0.5)
    # The distribution for females has a lower peek but a heavier tail compared to males
    
    # ATE by sex
    cfw_ate_by_sex <- cfw_agg$ate$full_mod$table |>
      filter(str_detect(name, "^Sex")) |>
      select(c(1:3, 5:6))
    # AIPW ATE estimate is 0.17% higher for females than males, 
    # but the confidence intervals have a large degree of overlap.
    
    # RATE for sex
    cfw_rate_sex <- cfw_agg$rate$full_mod$table |>
      filter(str_detect(name, "^X_02"))
    # The p-value for the RATE test of heterogeneity for sex is 0.41, indicating that we cannot dismiss the 
    # hypothesis of no heterogeneity.
  }
)
