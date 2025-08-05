# ===   0. Setup   ================================================================================================ ====
# load packages
library(tidyverse)
library(dtplyr)
library(lubridate)
library(data.table)
library(future)
library(furrr)
library(zoo)

# Define data paths
data_path <- "PATH"
cohort_path <- "data/cohort_raw_data/"

# === 1. Data read-in ============================================================================================= ====
# replace PNR with PNR_active in read-in registers (PNR_active contains newest known PNR for a person and should
# uniquely identify a person, where PNR can change over time)
t_person <- readRDS(paste0(data_path, "CPR_t_person.rds"))
pnr_active <- select(t_person, "PNR", "PNR_active")
saveRDS(pnr_active, paste0(cohort_path, "pnr_active.rds"))

## Source population ---------------------------------------------------------------------------------------------- ----
# Age >40, period Jan 2014 - Dec 2020
sourcepop <- t_person |>
  select(
    "pnr" = "PNR",
    "pnr_active" = "PNR_active",
    "sex" = "C_KON",
    "dob" = "D_FODDATO",
    "stat" = "C_STATUS",
    "stdt" = "D_STATUS_HEN_START"
  ) |>
  # status codes: 01: active with danish address, 70-90: Inactive; disappeared, emigrated or dead
  filter(stat %in% c("01", "70", "80", "90")) |> 
  mutate(
    start = pmax(dob + years(40), dmy("01/01/2014")),
    stop = case_when(
      stat == "01" ~ dmy("31/12/2020"),
      TRUE ~ pmin((stdt + days(1)), dmy("31/12/2020"))
    )
  ) |>
  filter(start < stop)

saveRDS(sourcepop, paste0(cohort_path, "sourcepop.rds"))

## Basic cohort --------------------------------------------------------------------------------------------------- ----
# load in all LMDB data
lmdb_1995_2020 <- readRDS(paste0(data_path, "LMDB_1995_2020.rds"))

# select data within 1 year of study period (01-01-2013), inner join with source population (both pnr and follow-up 
# period), add drug categories about hypertension. 
all_lmdb <- lmdb_1995_2020 |>
  filter(eksd > dmy("01/01/2013")) |>
  left_join(pnr_active, by = c("PNR")) |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything()) |>
  inner_join(sourcepop |> select("pnr", "start", "stop"), by = c("PNR" = "pnr")) |>
  filter(start - years(1) <= eksd & eksd < stop) |>
  rename("pnr" = "PNR") |>
  select("pnr", "eksd", "atc", "apk", "packsize", "strnum", "strunit", "volume", "indo") |>
  map2(
    c(
      "pnr", "Ekspeditionsdato", "ATC-kode 5. niveau", "Antal pakninger", "Pakningsstørrelse", "Styrke, numerisk",
      "Styrke, enhed", "Måltype", "Indikation for ordination"
    ),
    \(x, label) {
      x <- structure(
        x, 
        label = label
      )
    }
  ) |>
  as_tibble() |>
  mutate(
    drugcat = case_when(
      stringi::stri_sub(atc, 1, 7) == "C03AB01" ~ "thz",
      stringi::stri_sub(atc, 1, 5) == "C09BA" ~ "thz",
      stringi::stri_sub(atc, 1, 5) == "C09DA" ~ "thz",
      stringi::stri_sub(atc, 1, 4) == "C09A" ~ "nthz",
      stringi::stri_sub(atc, 1, 4) == "C09C" ~ "nthz",
      stringi::stri_sub(atc, 1, 5) == "C08CA" ~ "nthz",
      TRUE ~ "other"
    ),
    drugcat2 = case_when(
      stringi::stri_sub(atc, 1, 7) == "C03AB01" ~ "bfz",
      stringi::stri_sub(atc, 1, 5) == "C09BA" ~ "hctz",
      stringi::stri_sub(atc, 1, 5) == "C09DA" ~ "hctz",
      stringi::stri_sub(atc, 1, 4) == "C09A" ~ "ras",
      stringi::stri_sub(atc, 1, 4) == "C09C" ~ "ras",
      stringi::stri_sub(atc, 1, 5) == "C08CA" ~ "ccb",
      TRUE ~ "other"
    ),
    drugcat3 = case_when(
      stringi::stri_sub(atc, 1, 3) == "C02" ~ "hyp",
      stringi::stri_sub(atc, 1, 3) == "C03" ~ "hyp",
      stringi::stri_sub(atc, 1, 3) == "C04" ~ "hyp",
      stringi::stri_sub(atc, 1, 3) == "C08" ~ "hyp",
      stringi::stri_sub(atc, 1, 3) == "C09" ~ "hyp",
      TRUE ~ "other"
    ),
    drugcat4 = case_when(
      stringi::stri_sub(atc, 1, 4) == "C09A" ~ "ace",
      stringi::stri_sub(atc, 1, 4) == "C09C" ~ "arb",
      TRUE ~ "other"
    ),
    days_supply = apk * packsize
  )

saveRDS(all_lmdb, paste0(cohort_path, "all_lmdb.rds"))

rm(lmdb_1995_2020)

# remove drugs not related to hypertension and sum any prescriptions on a specific drug given on the same day
hyp_lmdb <- data.table::as.data.table(all_lmdb)[
  drugcat3 != "other",
  list(
    tmp_duration = sum(days_supply),
    tmp_apk = sum(apk),
    no_of_prescriptions = .N
  ),
  by = list(pnr, eksd, atc, strnum, strunit, indo, drugcat, drugcat2, drugcat3, drugcat4)
] |>
  as_tibble() |>
  # allow a 50% grace on duration of a prescription, remove any duplicate rows
  mutate(
    duration = 1.5 * tmp_duration,
    apk = tmp_apk
  ) |>
  select(
    "pnr", "eksd", "atc", "drugcat", "drugcat2", "drugcat3", "drugcat4", "apk", "strnum", "strunit", "indo", "duration"
  ) |>
  distinct()

# filter on first prescription (if tied on date, select first alphabetically among drugcat2)
first_prescription <- inner_join(
  select(hyp_lmdb, "pnr", "eksd", "atc", "drugcat", "drugcat2", "drugcat3", "drugcat4", "indo", "strnum", "strunit"),
  select(sourcepop, "pnr", "start", "stop"),
  by = "pnr"
) |>
  filter(start <= eksd & eksd < stop & drugcat != "other") |>
  slice_min(tibble(eksd, drugcat2), by = "pnr", with_ties = FALSE) |>
  select(
    "pnr", "start", "stop", "eksd", "atc", "drugcat", "drugcat2", "drugcat3", "drugcat4", "indo", "strnum", "strunit"
  )

saveRDS(first_prescription, paste0(cohort_path, "first_prescription.rds"))

# join first prescription to source population  
pop_with_hyp_date <- left_join(
  sourcepop,
  first_prescription |> select(-c("start", "stop")) |> rename("hyp_date" = "eksd"),
  by = "pnr"
)

rm(sourcepop, first_prescription)

# create table with periods of different exposure status
tmp_episodes <- bind_rows(
  data.table::as.data.table(pop_with_hyp_date)[
    is.na(hyp_date), 
    list(`in` = start, out = stop, indexdrug = "none"),
    by = pnr
  ] |> 
    as_tibble() |>
    inner_join(pop_with_hyp_date, by = "pnr"),
  data.table::as.data.table(pop_with_hyp_date)[!is.na(hyp_date), list(`in` = c(start, hyp_date)), by = pnr] |>
    as_tibble() |>
    bind_cols(
      data.table::as.data.table(pop_with_hyp_date)[!is.na(hyp_date), list(out = c(hyp_date, stop)), by = pnr] |>
        as_tibble() |>
        select("out")
    ) |>
    bind_cols(
      data.table::as.data.table(pop_with_hyp_date)[!is.na(hyp_date), list(indexdrug = c("none", drugcat2)), by = pnr] |>
        as_tibble() |>
        select("indexdrug")
    ) |>
    inner_join(pop_with_hyp_date, by = "pnr")
)

# split non-exposed period into two year courses
episodes <- bind_rows(
  tmp_episodes |>
    filter(indexdrug != "none" | (out - `in`) < years(2)) |>
    mutate(course = 1),
  data.table::as.data.table(tmp_episodes)[
    (out - `in`) >= years(2) & indexdrug == "none", 
    list(`in` = c(`in`, NA_Date_)), 
    by = list(pnr, indexdrug)
  ] |>
    as_tibble() |>
    bind_cols(
      data.table::as.data.table(tmp_episodes)[
        (out - `in`) >= years(2) & indexdrug == "none",
        list(out = c(NA_Date_, out)), 
        by = list(pnr, indexdrug)
      ] |>
        as_tibble() |>
        select("out")
    ) |>
    bind_cols(
      data.table::as.data.table(tmp_episodes)[
        (out - `in`) >= years(2) & indexdrug == "none", 
        list(course = c(1, 2)), 
        by = list(pnr, indexdrug)
      ] |>
        as_tibble() |>
        select("course")
    ) |>
    inner_join(tmp_episodes |> select(-c("in", "out")), by = c("pnr", "indexdrug"))
) |> mutate(
  `in` = fcase(
    is.na(`in`), lag(`in`) + years(2),
    !is.na(`in`), `in`
  ),
  out = fcase(
    is.na(`out`), `in` + years(2) - days(1),
    !is.na(`out`), `out`
  )
)

while(episodes |> filter((out - `in`) >= years(2) & indexdrug == "none") |> nrow() > 0) {
  episodes <- bind_rows(
    episodes |>
      filter(indexdrug != "none" | (out - `in`) < years(2)),
    data.table::as.data.table(episodes)[
      (out - `in`) >= years(2) & indexdrug == "none", 
      list(`in` = c(`in`, NA_Date_)), 
      by = list(pnr, indexdrug)
    ] |>
      as_tibble() |>
      bind_cols(
        data.table::as.data.table(episodes)[
          (out - `in`) >= years(2) & indexdrug == "none",
          list(out = c(NA_Date_, out)), 
          by = list(pnr, indexdrug)
        ] |>
          as_tibble() |>
          select("out")
      ) |>
      bind_cols(
        data.table::as.data.table(episodes)[
          (out - `in`) >= years(2) & indexdrug == "none", 
          list(course = c(course, course + 1)), 
          by = list(pnr, indexdrug)
        ] |>
          as_tibble() |>
          select("course")
      ) |>
      inner_join(tmp_episodes |> select(-c("in", "out")), by = c("pnr", "indexdrug"))
  ) |>
    mutate(
      `in` = fcase(
        is.na(`in`), lag(`in`) + years(2),
        !is.na(`in`), `in`
      ),
      out = fcase(
        is.na(`out`), `in` + years(2) - days(1),
        !is.na(`out`), `out`
      )
    )
}

episodes <- arrange(episodes, "pnr")

rm(pop_with_hyp_date)

# table with non-exposed time in follow-up period split into two year courses, containing both users and non-users
all_users_and_non_users <- episodes |>
  filter(out > `in`) |>
  mutate(
    count = seq_len(n()),
    .by = "pnr"
  ) |>
  mutate(
    id = paste0(pnr, "_", count)
  ) |>
  select(
    "pnr", "sex", "dob", "stat", "stdt" ,"start", "stop", "hyp_date", "atc", "drugcat", 
    "drugcat2", "drugcat3", "drugcat4", "indo", "strnum", "strunit", "in", "out",
    "indexdrug", "course", "id"
  )

saveRDS(all_users_and_non_users, paste0(cohort_path, "all_users_and_non_users.rds"))

rm(episodes)

# identify new users
new_users <- all_users_and_non_users |>
  filter(indexdrug != "none") |>
  left_join(
    hyp_lmdb |> select("pnr", "eksd", "duration", "drugcat2_hyp" = "drugcat2", "drugcat3_hyp" = "drugcat3"),
    by = "pnr",
    relationship = "many-to-many"
  ) |>
  summarise(
    across(everything(), \(x) x[1]),
    # use of antihypertensive drug within 1 year of first prescription with study drug in study period
    x_hyp = sum((eksd < hyp_date) & ((eksd + duration) > (hyp_date - years(1))) & drugcat3_hyp == "hyp") > 0,
    x_multiple = sum(drugcat2_hyp %in% c("bfz", "hctz", "ras", "ccb") & drugcat2 != drugcat2_hyp & hyp_date == eksd) > 0,
    .by = "id"
  ) |>
  select(-c("eksd", "duration", "drugcat2_hyp", "drugcat3_hyp", "course")) |>
  # filter on no antihypertensive drugs in 1 year before indexdate
  filter(x_hyp == 0) |>
  select(
    "id", "pnr", "sex", "dob", "start", "stop", "stat", "stdt", "in", "out", 
    "indexdrug", "x_multiple", "indexdate" = "hyp_date", "atc", "drugcat", 
    "drugcat2", "drugcat3", "drugcat4", "indo", "strnum", "strunit"
  )

tmp_basic_cohort <- arrange(new_users, id, pnr)

saveRDS(tmp_basic_cohort, paste0(cohort_path, "tmp_basic_cohort.rds"))

rm(all_users_and_non_users, new_users)

### exclusion criteria -------------------------------------------------------------------------------------------- ----
# lpr diagnosis codes
# read in LPR2 information on contacts and diagnoses
t_diag <- bind_rows(
  readRDS(paste0(data_path, "LPR_t_diag_1977_2018.RDS")),
  readRDS(paste0(data_path, "LPR_t_diag_2019.RDS")) |>
    mutate(C_DIAGMOD = NA_character_) |>
    select("RECNUM", "C_DIAG", "C_DIAGMOD", "C_DIAGTYPE", "C_TILDIAG")
) |>
  filter(# Exlude non-confirmed diagnoses, reference diagnoses, temporary diagnoses, and add-on diagnoses
    (C_DIAGMOD != "2" | is.na(C_DIAGMOD)) &
      !(C_DIAGTYPE %in% c("+", "H", "M"))
  ) 
t_adm <- bind_rows(
  readRDS(paste0(data_path, "LPR_t_adm_1977_2018.RDS")) |>
    select("PNR", "RECNUM", "D_INDDTO", "D_UDDTO", "C_PATTYPE", "C_INDM"),
  readRDS(paste0(data_path, "LPR_t_adm_2019.RDS")) |>
    select("PNR", "RECNUM", "D_INDDTO", "D_UDDTO", "C_PATTYPE", "C_INDM")
)

# read in LPR3 information on contacts and diagnoses
t_diag_lpr3 <- readRDS(paste0(data_path, "LPR3_F_diagnoser.RDS")) |>
  filter(senere_afkraeftet == "Nej" & !(diagnosekode %in% c("+", "H", "M"))) |>
  select("dw_ek_kontakt", "C_DIAG" = "diagnosekode", "C_DIAGTYPE" = "diagnosetype")
t_kontakt_lpr3 <- readRDS(paste0(data_path, "LPR3_F_kontakter.RDS")) 
t_forloeb_lpr3 <- readRDS(paste0(data_path, "LPR3_F_forloeb.RDS")) 

# To determine if a diagnosis is present in a period before index date, 
# only diagnosis code, start- and end time of diagnosis are needed.
# For use with determining inpatient and outpatient hospitalization, 
# a different approach will be needed for LPR3.
all_diag_lpr3 <- t_diag_lpr3 |>
  inner_join(
    t_kontakt_lpr3 |>
      select("PNR", "dw_ek_kontakt", "dato_start", "dato_slut"),
    by = "dw_ek_kontakt"
  ) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  mutate(
    PNR = PNR,
    C_DIAG = diagnosekode,
    C_DIAGTYPE = diagnosetype,
    D_INDDTO = dato_start,
    D_UDDTO = dato_slut,
    C_PATTYPE = "NA",
    C_INDM = "NA",
    .keep = "none"
  )

all_diag <- t_diag |>
  left_join(
    t_adm,
    by = "RECNUM"
  ) |>
  select("PNR", "C_DIAG", "C_DIAGTYPE", "D_INDDTO", "D_UDDTO", "C_PATTYPE", "C_INDM") |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything()) |>
  bind_rows(all_diag_lpr3)

saveRDS(all_diag, paste0(cohort_path, "all_diag.rds"))

rm(all_diag_lpr3)

# operation, treatment, and examination codes
t_sksopr <- bind_rows(
  readRDS(paste0(data_path, "LPR_t_sksopr_1996_2018.RDS")),
  readRDS(paste0(data_path, "LPR_t_sksopr_2019.RDS"))
)
t_sksube <- bind_rows(
  readRDS(paste0(data_path, "LPR_t_sksube_1999_2018.RDS")),
  readRDS(paste0(data_path, "LPR_t_sksube_2019.RDS"))
)
all_opr <- bind_rows(t_sksopr, t_sksube)
x_opr <- all_opr |>
  left_join(t_adm, by = "RECNUM") |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", all_of(names(all_opr))) |>
  filter(# With current use, only need codes starting with BJF
    str_detect(C_OPR, "^(BJF[DZ]|K[ABEFGHIJKLMNOPQ])") |
      str_detect(C_TILOPR, "^(BJF[DZ]|K[ABEFGHIJKLMNOPQ])")
  )

# both procedure and add-on codes in same column of LPR3. 
t_proc_kirurgi_lpr3 <- readRDS(paste0(data_path, "LPR3_F_procedurer_kirurgi.RDS"))
t_proc_andre_lpr3 <- readRDS(paste0(data_path, "LPR3_F_procedurer_andre.RDS"))

x_opr_lpr3 <- bind_rows(
  t_proc_kirurgi_lpr3 |>
    select("dw_ek_forloeb", "procedurekode", "dato_start", "tidspunkt_start"),
  t_proc_andre_lpr3 |>
    select("dw_ek_forloeb", "procedurekode", "dato_start", "tidspunkt_start")
) |> 
  left_join(t_forloeb_lpr3 |> select("PNR", "dw_ek_forloeb"), by = "dw_ek_forloeb") |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything()) |>
  filter( # With current use, only need codes starting with BJF
    str_detect(procedurekode, "^(BJF[DZ]|K[ABEFGHIJKLMNOPQ])")
  ) |>
  mutate(
    PNR = PNR,
    C_OPR = procedurekode,
    C_TILOPR = "",
    D_ODTO = dato_start,
    .keep = "none"
  )

# combine LPR2 and LPR3
x_opr <- bind_rows(
  x_opr |> select("PNR", "C_OPR", "C_TILOPR", "D_ODTO"),
  x_opr_lpr3
)

rm(t_sksopr, t_sksube, t_proc_kirurgi_lpr3, t_proc_andre_lpr3, all_opr, x_opr_lpr3)
saveRDS(x_opr, paste0(cohort_path, "x_opr.rds"))

# Users with multiple prescriptions on the same day
x_multiple <- tmp_basic_cohort |>
  select(id, x_multiple)

saveRDS(x_multiple, paste0(cohort_path, "x_multiple.rds"))

# Users not living in Region Midt
komreg_v4 <- read.delim(paste0(cohort_path, "komreg_v4.txt", sep = "=")) |>
  mutate(
    C_KOM = str_remove_all(C_KOM, "[^\\w]"),
    region = str_remove_all(region, "[^\\w]")
  )
address <- bind_rows(
  readRDS(paste0(data_path, "CPR_t_adresse.rds")),
  readRDS(paste0(data_path, "CPR_t_adresse_hist.rds")),
  readRDS(paste0(data_path, "CPR_t_arkiv_adresse_hist.rds"))
) |>
  left_join(komreg_v4, by = "C_KOM") |>
  filter(str_detect(region, "^8[1-5]")) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything())

cohort_address <- tmp_basic_cohort |> # Last known address as index date
  select("id", "pnr", "indexdate") |>
  left_join(address, by = c("pnr" = "PNR")) |>
  filter(
    (D_TILFLYT_DATO <= indexdate | is.na(D_TILFLYT_DATO)) & 
      (indexdate <= pmin(D_FRAFLYT_DATO, dmy("31122022")) | is.na(D_FRAFLYT_DATO))
  ) |>
  slice_max(D_TILFLYT_DATO, n = 1, by = "id") |>
  mutate(obs_udlandsk = is.na(region))

tmp_x_region_midt <- cohort_address |> # check still living in DK at index date
  mutate(x_region_midt = str_detect(region, "^82")) |>
  select("id", "pnr", "x_region_midt")

tmp_x_region_midt_2 <- tmp_basic_cohort |> # only check address changes before index date
  select("id", "pnr", "indexdate") |>
  left_join(address, by = c("pnr" = "PNR")) |>
  filter(D_TILFLYT_DATO <= indexdate | is.na(D_TILFLYT_DATO)) |>
  slice_max(D_TILFLYT_DATO, n = 1, by = "id", with_ties = FALSE) |>
  mutate(x_region_midt_2 = str_detect(region, "^82"))

tmp_x_region_midt_3 <- tmp_basic_cohort |>
  select("id") |>
  left_join(tmp_x_region_midt |> select("id", "x_region_midt"), by = "id") |>
  left_join(tmp_x_region_midt_2 |> select("id", "x_region_midt_2"), by = "id")

x_region_midt <- tmp_x_region_midt_3 |> # go with only checking address before index date (not if moved away)
  mutate(x_region_midt = ifelse(x_region_midt_2, x_region_midt_2, x_region_midt))

saveRDS(x_region_midt, paste0(cohort_path, "x_region_midt.rds"))
rm(komreg_v4, cohort_address, tmp_x_region_midt_2, tmp_x_region_midt_3)

# users with no or less than 2 years of known address in Denmark
no_address <- address |> 
  inner_join(tmp_basic_cohort |> select("id", "pnr", "indexdate"), by = c("PNR" = "pnr")) |>
  filter(# address periods starting or ending within 2 years of indexdate
    (D_TILFLYT_DATO <= indexdate | is.na(D_TILFLYT_DATO)) & 
      ((indexdate - years(2)) < pmin(D_FRAFLYT_DATO, dmy("31122022")) | is.na(D_FRAFLYT_DATO))
  ) |>
  select("id", "pnr" = "PNR", "D_TILFLYT_DATO", "D_FRAFLYT_DATO", "region", "indexdate") |> 
  arrange(id, D_TILFLYT_DATO) |>
  mutate(
    fraflyt_lag = lag(D_FRAFLYT_DATO),
    .by = "id"
  ) |>
  summarise(# users with missing address information within 2 years of index date
    pnr = first(pnr),
    x_no_address = isTRUE(any(fraflyt_lag < D_TILFLYT_DATO)),
    .by = "id"
  )

x_no_address <- tmp_basic_cohort |>
  select("id", "pnr") |>
  left_join(no_address |> select(-c("pnr")), by = "id") |>
  mutate(
    x_no_address = ifelse(is.na(x_no_address), TRUE, x_no_address)
  )

saveRDS(x_no_address, paste0(cohort_path, "x_no_address.rds"))
rm(no_address)

# Diseases: Endstage illness, Adrenal, KAS, CKD
tmp_x_diseases <- tmp_basic_cohort |>
  left_join(
    bind_rows(
      all_diag |> select("PNR", "C_DIAG", "date" = "D_INDDTO") |> mutate(C_OPR = NA_character_, C_TILOPR = NA_character_),
      x_opr |> select("PNR", "C_OPR", "C_TILOPR", "date" = "D_ODTO") |> mutate(C_DIAG = NA_character_)
    ), 
    by = c("pnr" = "PNR")
  ) |>
  filter(indexdate - years(5) <= date & date < indexdate) |>
  summarise(
    pnr = first(pnr),
    x01 = sum(
      str_detect(C_DIAG, "^D(E4[123]|R(34|64|40[23]|029)|Z(49|991)|F051|I702A)"),
      na.rm = TRUE
    ) > 0,
    x02 = sum(
      str_detect(C_DIAG, "^DE27"),
      na.rm = TRUE
    ) > 0,
    x03 = sum(
      str_detect(C_DIAG, "^D(I701|N280|Q271)"),
      na.rm = TRUE
    ) > 0,
    x04 = sum(
      str_detect(C_DIAG, "^D(N18[45]|Z992)") |
        str_detect(C_OPR, "^BJF[DZ]") |
        str_detect(C_TILOPR, "^BJF[DZ]"),
      na.rm = TRUE
    ) > 0,
    .by = "id"
  ) |> 
  right_join(tmp_basic_cohort |> select("id", "pnr"), by = c("id", "pnr")) |>
  mutate(across(x01:x04, \(x) ifelse(is.na(x), FALSE, x)))
attr(tmp_x_diseases$x01, "label") <- "End-stage illness"
attr(tmp_x_diseases$x02, "label") <- "Adrenal insufficiency disorders"
attr(tmp_x_diseases$x03, "label") <- "Kidney artery stenosis"
attr(tmp_x_diseases$x04, "label") <- "Severe CKD"

# GFR
future::plan("multisession", workers = 10)
gfr <- furrr::future_map(
  .x = 1:10,
  .options = furrr_options(
    globals = c("data_path"),
    packages = c("glue", "dplyr")
  ),
  .f = \(i) {
    data <- readRDS(glue("{data_path}LAB_dm_forsker_{i}.rds")) |>
      filter(ANALYSISCODE %in% c("DNK35131", "DNK35302"))
    return(data)
  }
) |>
  list_rbind()|>
  mutate(
    value_num = as.numeric(VALUE),
    sampling_year = year(SAMPLINGDATE)
  ) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything())
future::plan("sequential")

gfr_pop <- tmp_basic_cohort |>
  select("id", "pnr", "indexdate") |>
  inner_join(gfr |> select("PNR", "SAMPLINGDATE", "ANALYSISCODE", "value_num"), by = c("pnr" = "PNR")) |>
  filter(SAMPLINGDATE <= indexdate) |>
  mutate(
    gfr_level = case_when(
      value_num >= 90 ~ "Stage 1",
      value_num < 90 & value_num >= 60 ~ "Stage 2",
      value_num < 60 & value_num >= 30 ~ "Stage 3",
      value_num < 30 ~ "Stage 4",
      TRUE ~ NA_character_
    )
  )

x_gfr <- tmp_basic_cohort |>
  select("id") |>
  left_join(gfr_pop, by = "id") |>
  summarise(
    x_gfr = sum(
      (SAMPLINGDATE <= indexdate) &
        str_detect(gfr_level, "Stage 4")
    ) > 0,
    .by = "id"
  ) |>
  replace_na(list(x_gfr = FALSE))

tmp_x_diseases_2 <- tmp_x_diseases |>
  left_join(x_gfr, by = "id")

x_diseases <- tmp_x_diseases_2 |>
  mutate(x04 = ifelse(x_gfr, TRUE, x04))
attr(x_diseases$x04, "label") <- "Severe CKD"

saveRDS(x_diseases, paste0(cohort_path, "x_diseases.rds"))
rm(tmp_x_diseases, tmp_x_diseases_2, gfr, x_gfr)

# Recent hyponatremia diagnoses
x_hypo_na_diag <- tmp_basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(all_diag, by = c("pnr" = "PNR")) |>
  mutate(
    hypo_na_diag_mt4 = indexdate - years(5) < D_INDDTO & D_INDDTO < indexdate - 120 & C_DIAG %in% c("DE871A", "DE222A"),
    hypo_na_diag_lt4 = indexdate - 120 <= D_INDDTO & D_INDDTO < indexdate & C_DIAG %in% c("DE871A", "DE222A")
  ) |>
  summarise(
    hypo_na_diag_mt4 = any(hypo_na_diag_mt4),
    hypo_na_diag_lt4 = any(hypo_na_diag_lt4),
    .by = "id"
  )

saveRDS(x_hypo_na_diag, paste0(cohort_path, "x_hypo_na_diag.rds"))

# measured low sodium
future::plan("multisession", workers = 10)
sodium <- furrr::future_map(
  .x = 1:10,
  .options = furrr_options(
    globals = c("data_path"),
    packages = c("glue", "dplyr")
  ),
  .f = \(i) {
    data <- readRDS(glue("{data_path}LAB_dm_forsker_{i}.rds")) |>
      filter(ANALYSISCODE == "NPU03429")
    return(data)
  }
) |>
  list_rbind()|>
  mutate(
    value_num = as.numeric(VALUE),
    sampling_year = year(SAMPLINGDATE)
  ) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything())
future::plan("sequential")

sodium_pop <- tmp_basic_cohort |>
  select("id", "pnr", "indexdate") |>
  inner_join(sodium |> select("PNR", "SAMPLINGDATE", "ANALYSISCODE", "value_num"), by = c("pnr" = "PNR")) |>
  filter(indexdate - years(5) <= SAMPLINGDATE & SAMPLINGDATE <= indexdate) |>
  mutate(
    hypo_na = case_when(
      value_num < 130 ~ "Low sodium",
      TRUE ~ NA_character_
    )
  )

x_sodium <- tmp_basic_cohort |>
  select("id") |>
  left_join(sodium_pop, by = "id") |>
  mutate(
    na_mt4 = indexdate - years(5) < SAMPLINGDATE & SAMPLINGDATE < indexdate - 120 & hypo_na == "Low sodium",
    na_lt4 = indexdate - 120 <= SAMPLINGDATE & SAMPLINGDATE < indexdate & hypo_na == "Low sodium"
  ) |>
  summarise(
    na_mt4 = any(ifelse(is.na(na_mt4), FALSE, na_mt4)),
    na_lt4 = any(ifelse(is.na(na_lt4), FALSE, na_lt4)),
    .by = "id"
  )

saveRDS(x_sodium, paste0(cohort_path, "x_sodium.rds"))
rm(sodium, sodium_pop)

x_hyponatremia <- x_hypo_na_diag |>
  left_join(x_sodium, by = "id") |>
  mutate(
    d26 = case_when(
      hypo_na_diag_lt4 | na_lt4 ~ "lt_4mths",
      hypo_na_diag_mt4 | na_mt4 ~ "mt_4mths",
      TRUE ~ "never"
    )
  )

saveRDS(x_hyponatremia, paste0(cohort_path, "x_hyponatremia.rds"))
rm(x_hypo_na_diag, x_sodium)

# Combining x criteria
tmp_basic_cohort_1 <- tmp_basic_cohort |>
  select(-c("x_multiple")) |> # x_multiple already in tmp_basic_cohort, so can skip this and next line
  left_join(x_multiple, by = "id") |>
  left_join(x_region_midt |> select("id", "x_region_midt"), by = "id") |>
  left_join(x_no_address |> select("id", "x_no_address"), by = "id") |>
  left_join(x_diseases |> select("id", "x01", "x02", "x03", "x04"), by = "id")

saveRDS(tmp_basic_cohort_1, paste0(cohort_path, "tmp_basic_cohort_1.rds"))
rm(x_region_midt, x_no_address, x_diseases)

tmp_basic_cohort_2 <- tmp_basic_cohort_1 |>
  select(
    "id", "pnr", "sex", "dob", "stat", "stdt", "start", "stop", "entry" = "in", "exit" = "out", "indexdate", "indexdrug"
  )

# Exiting based on moving to Central Region Denmark
cohort_address_cens <- tmp_basic_cohort_2 |>
  select("id", "pnr", "entry", "exit") |>
  left_join(address, by = c("pnr" = "PNR")) |>
  filter(entry < D_TILFLYT_DATO & D_TILFLYT_DATO <= exit) |>
  filter(str_detect(region, "^82")) |>
  slice_min(D_TILFLYT_DATO, by = "id") |>
  select("id", "pnr", "C_KOM", "D_TILFLYT_DATO", "D_FRAFLYT_DATO", "region")

# basic cohort
basic_cohort <- tmp_basic_cohort_2 |>
  left_join(cohort_address_cens |> select("id", "reg_midt_censor_date" = "D_TILFLYT_DATO"), by = "id") |>
  relocate("reg_midt_censor_date", .before = "indexdate") |>
  arrange(id, indexdate)

saveRDS(basic_cohort, paste0(cohort_path, "basic_cohort.rds"))
rm(tmp_basic_cohort_2, cohort_address_cens)

## Episodes ------------------------------------------------------------------------------------------------------- ----
new_user_only <- basic_cohort

# date of prescription fulfillment and expected date the drugs last until.
hyp_exposure <- new_user_only |>
  select("id", "pnr", "indexdate", "indexdrug", "stop") |>
  inner_join(hyp_lmdb, by = "pnr") |>
  filter(indexdate <= eksd & eksd < stop & indexdrug == drugcat2) |>
  mutate(
    x1 = eksd,
    x2 = eksd + duration
  ) |>
  select("id", "pnr", "drugcat2", "x1", "x2")

# Data with an episode counter for use of drug, and a count of fulfilled prescriptions in each episode
allexp <- hyp_exposure |>
  arrange(id, x1, x2) |>
  lazy_dt() |>
  mutate(
    new_episode = is.na(lag(x1)) | x1 > lag(x2),
    episode = cumsum(new_episode),
    .by = "id"
  ) |>
  summarise(
    hyp_start = first(x1),
    hyp_stop = last(x2),
    count = n(),
    .by = c("id", "episode")
  ) |>
  as_tibble()
attr(allexp$hyp_stop, "label") <- NULL

tmp_episodes_1 <- new_user_only |>
  select(-c("reg_midt_censor_date")) |>
  left_join(allexp, by = "id") |>
  filter(indexdate == hyp_start) |>
  select(-c("episode", "count"))

# data on hypertensive add-on drugs given during the main exposure period
episodes_add_on_any_hyp <- tmp_episodes_1 |>
  select("id", "pnr", "indexdrug", "hyp_start", "exit") |>
  left_join(hyp_lmdb |> select("pnr", "eksd", "drugcat2", "drugcat3"), by = "pnr") |>
  filter(indexdrug != drugcat2 & hyp_start < eksd & eksd < exit & drugcat3 == "hyp") |>
  summarise(
    add_on_date_any_hyp = min(eksd),
    .by = "id"
  ) |>
  right_join(tmp_episodes_1, by = "id") |>
  arrange(id) |>
  relocate("add_on_date_any_hyp", .after = "hyp_stop") |>
  mutate(
    add_on_any_hyp = !is.na(add_on_date_any_hyp)
  )
attr(episodes_add_on_any_hyp$add_on_date_any_hyp, "label") <- "Date of receiving any hyp add-on during FU"
attr(episodes_add_on_any_hyp$add_on_any_hyp, "label") <- "Receiving any hyp add-on during FU"

episodes_add_on_bfzccb <- tmp_episodes_1 |>
  select("id", "pnr", "indexdrug", "hyp_start", "exit") |>
  left_join(hyp_lmdb |> select("pnr", "eksd", "drugcat2", "drugcat3"), by = "pnr") |>
  filter(indexdrug != drugcat2 & hyp_start < eksd & eksd < exit & drugcat2 %in% c("bfz", "ccb")) |>
  summarise(
    add_on_date_drugcat2 = min(eksd),
    .by = "id"
  ) |>
  right_join(tmp_episodes_1, by = "id") |>
  filter(indexdrug %in% c("bfz", "ccb")) |>
  arrange(id) |>
  relocate("add_on_date_drugcat2", .after = "hyp_stop") |>
  mutate(
    add_on_drugcat2 = !is.na(add_on_date_drugcat2)
  )
attr(episodes_add_on_bfzccb$add_on_date_drugcat2, "label") <- "Date of receiving comparator during FU"
attr(episodes_add_on_bfzccb$add_on_drugcat2, "label") <- "Receiving comparator (B/C,T/R) during FU"

episodes_add_on_hctzras <- tmp_episodes_1 |>
  select("id", "pnr", "indexdrug", "hyp_start", "exit") |>
  left_join(hyp_lmdb |> select("pnr", "eksd", "drugcat2", "drugcat3"), by = "pnr") |>
  filter(indexdrug != drugcat2 & hyp_start < eksd & eksd < exit & drugcat2 %in% c("hctz", "ras")) |>
  summarise(
    add_on_date_drugcat2 = min(eksd),
    .by = "id"
  ) |>
  right_join(tmp_episodes_1, by = "id") |>
  filter(indexdrug %in% c("hctz", "ras")) |>
  arrange(id) |>
  relocate("add_on_date_drugcat2", .after = "hyp_stop") |>
  mutate(
    add_on_drugcat2 = !is.na(add_on_date_drugcat2)
  )
attr(episodes_add_on_hctzras$add_on_date_drugcat2, "label") <- "Date of receiving comparator during FU"
attr(episodes_add_on_hctzras$add_on_drugcat2, "label") <- "Receiving comparator (B/C,T/R) during FU"

episodes_add_on_drugcat2_any <- tmp_episodes_1 |>
  select("id", "pnr", "indexdrug", "hyp_start", "exit") |>
  left_join(hyp_lmdb |> select("pnr", "eksd", "drugcat2", "drugcat3"), by = "pnr") |>
  filter(indexdrug != drugcat2 & hyp_start < eksd & eksd < exit & drugcat2 %in% c("bfz", "ccb", "hctz", "ras")) |>
  summarise(
    add_on_date_drugcat2_any = min(eksd),
    .by = "id"
  ) |>
  right_join(tmp_episodes_1, by = "id") |>
  filter(indexdrug %in% c("bfz", "ccb", "hctz", "ras")) |>
  arrange(id) |>
  relocate("add_on_date_drugcat2_any", .after = "hyp_stop") |>
  mutate(
    add_on_drugcat2_any = !is.na(add_on_date_drugcat2_any)
  )
attr(episodes_add_on_drugcat2_any$add_on_date_drugcat2_any, "label") <- "Date of receiving any other study drug during FU"
attr(episodes_add_on_drugcat2_any$add_on_drugcat2_any, "label") <- "Receiving any other study drug during FU"

add_on_hyp_start_stop <- episodes_add_on_any_hyp |>
  select("id", "indexdrug", "hyp_start", "hyp_stop", "add_on_date_any_hyp", "add_on_any_hyp") |>
  left_join(
    episodes_add_on_bfzccb |> 
      select(
        "id", 
        "add_on_date_drugcat2_bfzccb" = "add_on_date_drugcat2",
        "add_on_drugcat2_bfzccb" = "add_on_drugcat2"
      ), 
    by = "id"
  ) |>
  left_join(
    episodes_add_on_hctzras |> 
      select(
        "id", 
        "add_on_date_drugcat2_hctzras" = "add_on_date_drugcat2",
        "add_on_drugcat2_hctzras" = "add_on_drugcat2"
      ), 
    by = "id"
  ) |>
  left_join(
    episodes_add_on_drugcat2_any |> 
      select("id", "add_on_date_drugcat2_any", "add_on_drugcat2_any"), 
    by = "id"
  ) |>
  mutate(
    add_on_date_drugcat2 = as.Date(ifelse(
      indexdrug %in% c("bfz", "ccb"),
      add_on_date_drugcat2_bfzccb,
      add_on_date_drugcat2_hctzras
    )),
    add_on_drugcat2 = ifelse(
      indexdrug %in% c("bfz", "ccb"),
      add_on_drugcat2_bfzccb,
      add_on_drugcat2_hctzras
    )
  ) |>
  select("id", "indexdrug", "hyp_start", "hyp_stop", "add_on_date_any_hyp", "add_on_any_hyp",
         "add_on_date_drugcat2", "add_on_drugcat2", "add_on_date_drugcat2_any", "add_on_drugcat2_any")
attr(add_on_hyp_start_stop$add_on_date_drugcat2, "label") <- "Date of receiving comparator during FU"
attr(add_on_hyp_start_stop$add_on_drugcat2, "label") <- "Receiving comparator (B/C,T/R) during FU"

saveRDS(add_on_hyp_start_stop, paste0(cohort_path, "add_on_hyp_start_stop.rds"))

rm(
  hyp_lmdb, new_user_only, hyp_exposure, allexp, tmp_episodes_1, 
  episodes_add_on_any_hyp, episodes_add_on_bfzccb, episodes_add_on_hctzras, episodes_add_on_drugcat2_any
)

## Covariates ----------------------------------------------------------------------------------------------------- ----
sssy <- readRDS(paste0(data_path, "SSSY2005_2020.rds")) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  select("PNR", everything())

## pediatric information
podiatry_cov <- sssy |>
  filter(SPEC2 %in% c("54", "59", "60")) |>
  mutate(
    d_date = as_date(HONUGE, format = "%y%V") # format date from year week format: yyww
  )

## table with covariates related to diagnoses
covs_diag <- basic_cohort |>
  inner_join(all_diag |> select("PNR", "C_DIAG", "D_INDDTO"), by = c("pnr" = "PNR")) |>
  filter(indexdate - years(5) <= D_INDDTO & D_INDDTO < indexdate) |>
  lazy_dt() |>
  summarise(
    d01 = any(str_detect(C_DIAG, "^D(I(110|13[02]|42[06789]|50[019])|J819)")),
    d02 = any(str_detect(C_DIAG, "^DI10")),
    d03 = any(str_detect(C_DIAG, "^DI15")),
    d04 = any(str_detect(C_DIAG, "^DI2[0-5]")),
    d05 = any(str_detect(C_DIAG, "^D(I6([134]|76)|G45)")),
    d06 = any(str_detect(C_DIAG, "^D(E23|G35|G91|Q03|G610|G93[56])")),
    d07 = any(str_detect(C_DIAG, "^DI(48|49)|DZ950")),
    d08 = any(str_detect(C_DIAG, "^DC(?!44)")),
    d09 = any(str_detect(C_DIAG, "^DC([01236]|4[0156789]|5[1-9]|7[012]|8[0-6]|884)")),
    d10 = any(str_detect(C_DIAG, "^(DC([012356]|50|7[3-9]|8[09]|9[0-7])|DC88(?!4))")),
    d12 = any(str_detect(C_DIAG, "^D(K65|K7[0-7]|I982|Z944|D684C)")),
    d13 = any(str_detect(C_DIAG, "^DK8(5|60|61)")),
    d14 = any(str_detect(C_DIAG, "^DJ4[456]")),
    d15 = any(str_detect(C_DIAG, "^DE1[0-4]")),
    d16 = any(str_detect(C_DIAG, "^DE86")),
    d17 = any(str_detect(C_DIAG, "^D(H54|L89|L97|R(32|5[34]|63[0346])|Z(74|96))")),
    d18 = any(str_detect(C_DIAG, "^D(R(26|296|67[1-9])|S(12|52|22[01]|320|42[2-4]|72[0-2])|T08|M625|Z993)")),
    d19 = any(str_detect(C_DIAG, "^D(F0[0-3]|G30)")),
    d20 = any(str_detect(C_DIAG, "^DZ50")),
    d23 = any(str_detect(C_DIAG, "^D(F[12345689]|T4[03]|R78[1-5])")),
    d24 = any(str_detect(C_DIAG, "^DB2[0-4]")),
    d25 = any(str_detect(C_DIAG, "^DR63[01]")),
    .by = "id"
  ) |>
  as_tibble() |>
  right_join(basic_cohort, by = "id") |>
  mutate(across(starts_with("d"), \(x) replace_na(x, FALSE))) |>
  arrange(id) |>
  select("id", starts_with("d"), -c("dob")) |>
  inner_join(
    basic_cohort |>
      inner_join(all_diag |> select("PNR", "C_DIAG", "D_INDDTO"), by = c("pnr" = "PNR")) |>
      filter(indexdate - years(5) <= D_INDDTO     & D_INDDTO     < indexdate) |>
      inner_join(gfr_pop |> select("id", "SAMPLINGDATE", "gfr_level"), by = c("id"), relationship = "many-to-many") |>
      filter(indexdate - years(5) <= SAMPLINGDATE & SAMPLINGDATE < indexdate) |>
      lazy_dt() |>
      summarise(
        d11 = any(str_detect(C_DIAG, "^D(N(0[01345]|1[79]|26)|I1(20|31|8[123]))") | gfr_level == "Stage 3"),
        .by = "id"
      ) |>
      as_tibble() |>
      right_join(basic_cohort, by = "id") |>
      mutate(across(starts_with("d"), \(x) replace_na(x, FALSE))) |>
      arrange(id) |>
      select("id", starts_with("d"), -c("dob")),
    by = "id"
  ) |>
  inner_join(
    basic_cohort |>
      inner_join(podiatry_cov |> select("PNR", "d_date", "KONTAKT"), by = c("pnr" = "PNR")) |>
      filter(indexdate - years(5) <= d_date & d_date < indexdate) |>
      lazy_dt() |>
      summarise(
        d21 = any(KONTAKT == 1),
        .by = "id"
      ) |>
      as_tibble() |>
      right_join(basic_cohort, by = "id") |>
      mutate(across(starts_with("d"), \(x) replace_na(x, FALSE))) |>
      arrange(id) |>
      select("id", starts_with("d"), -c("dob")),
    by = "id"
  ) |>
  inner_join(
    basic_cohort |>
      inner_join(all_diag |> select("PNR", "C_DIAG", "D_INDDTO"), by = c("pnr" = "PNR")) |>
      filter(indexdate - years(5) <= D_INDDTO & D_INDDTO < indexdate) |>
      inner_join(all_lmdb |> select("pnr", "eksd", "atc"), by = "pnr", relationship = "many-to-many") |> 
      # 1 year of lookback on drug use
      filter(indexdate - years(1) <= eksd & eksd < indexdate) |>
      lazy_dt() |>
      summarise(
        d22 = any(str_detect(C_DIAG, "^DF10") | str_detect(atc, "^N07BB0[13]")),
        .by = "id"
      ) |>
      as_tibble() |>
      right_join(basic_cohort, by = "id") |>
      mutate(across(starts_with("d"), \(x) replace_na(x, FALSE))) |>
      arrange(id) |>
      select("id", starts_with("d"), -c("dob")),
    by = "id"
  )
attr(covs_diag$d01, "label") <- "Heart failure"
attr(covs_diag$d02, "label") <- "Essential hypertension"
attr(covs_diag$d03, "label") <- "Secondary hypertension"
attr(covs_diag$d04, "label") <- "Ischemic heart disease"
attr(covs_diag$d05, "label") <- "Cerebrovascular disease"
attr(covs_diag$d06, "label") <- "Other CNS disorders"
attr(covs_diag$d07, "label") <- "Arrythmias"
attr(covs_diag$d08, "label") <- "Any malignancy"
attr(covs_diag$d09, "label") <- "Malignancy associated with hyponatremia"
attr(covs_diag$d10, "label") <- "Other malignancy"
attr(covs_diag$d11, "label") <- "Renal disorders"
attr(covs_diag$d12, "label") <- "Liver disease and peritonitis"
attr(covs_diag$d13, "label") <- "pancreatitis"
attr(covs_diag$d14, "label") <- "COPD"
attr(covs_diag$d15, "label") <- "Diabetes"
attr(covs_diag$d16, "label") <- "Dehydration"
attr(covs_diag$d17, "label") <- "Frail general health"
attr(covs_diag$d18, "label") <- "Physical impairment"
attr(covs_diag$d19, "label") <- "Mental impairment"
attr(covs_diag$d20, "label") <- "Rehabilitation contacts"
attr(covs_diag$d21, "label") <- "Podiatric contacts"
attr(covs_diag$d22, "label") <- "Alcohol abuse"
attr(covs_diag$d23, "label") <- "Drug abuse"
attr(covs_diag$d24, "label") <- "HIV"
attr(covs_diag$d25, "label") <- "Anorexia and primary polyclipsia"

saveRDS(covs_diag, paste0(cohort_path, "covs_diag.rds"))
rm(gfr_pop, podiatry_cov)

## table with covariates related to drug use
covs_drugs <- basic_cohort |>
  inner_join(all_lmdb |> select("pnr", "eksd", "atc"), by = "pnr") |>
  # 1 year of lookback on drug use
  filter(indexdate - years(1) <= eksd & eksd < indexdate) |> 
  lazy_dt() |>
  summarise(
    m01 = any(str_detect(atc, "^A02(A|B[AC])")),
    m02 = any(str_detect(atc, "^A0[67]")),
    m03 = any(str_detect(atc, "^A10B")),
    m04 = any(str_detect(atc, "^A10A")),
    m05 = any(str_detect(atc, "^B01A[AEF]")),
    m06 = any(str_detect(atc, "^B01AC06")),
    m07 = any(str_detect(atc, "^B01AC(04|22|24)")),
    m08 = any(str_detect(atc, "^C01DA")),
    m09 = any(str_detect(atc, "^C07")),
    m10 = any(str_detect(atc, "^C10")),
    m11 = any(str_detect(atc, "^H01BA02")),
    m12 = any(str_detect(atc, "^H03AA01")),
    m13 = any(str_detect(atc, "^M01A(?!X)")),
    m14 = any(str_detect(atc, "^N03A")),
    m15 = any(str_detect(atc, "^N02A")),
    m16 = any(str_detect(atc, "^N06A")),
    m17 = any(str_detect(atc, "^N05A")),
    m18 = any(str_detect(atc, "^R03")),
    n_drugs = length(unique(atc)),
    .by = "id"
  ) |>
  as_tibble() |>
  right_join(basic_cohort, by = "id") |>
  mutate(across(starts_with("m"), \(x) replace_na(x, FALSE))) |>
  mutate(n_drugs = replace_na(n_drugs, 0)) |>
  arrange(id) |>
  select("id", starts_with("m"), "n_drugs") |>
  mutate(
    n_drug_grp = factor(
      case_when(
        n_drugs < 1 ~ 0L,
        1 <= n_drugs & n_drugs < 4 ~ 1L,
        4 <= n_drugs & n_drugs < 7  ~ 2L,
        7 <= n_drugs & n_drugs < 10 ~ 3L,
        10 <= n_drugs ~ 4L,
      ),
      levels = 0:4
    )
  )
attr(covs_drugs$m01, "label") <- "PPI and anticides"
attr(covs_drugs$m02, "label") <- "Obstipation/diarrhea"
attr(covs_drugs$m03, "label") <- "Antidiabetics (not insulin)"
attr(covs_drugs$m04, "label") <- "Insulin"
attr(covs_drugs$m05, "label") <- "Oral anticoagulants"
attr(covs_drugs$m06, "label") <- "Aspirin"
attr(covs_drugs$m07, "label") <- "ADPi"
attr(covs_drugs$m08, "label") <- "Nitrates"
attr(covs_drugs$m09, "label") <- "Beta-blockers"
attr(covs_drugs$m10, "label") <- "Lipid lowering drugs"
attr(covs_drugs$m11, "label") <- "Desmopressin"
attr(covs_drugs$m12, "label") <- "Eltroxin"
attr(covs_drugs$m13, "label") <- "NSAIDs"
attr(covs_drugs$m14, "label") <- "Antiepileptics"
attr(covs_drugs$m15, "label") <- "Opioids"
attr(covs_drugs$m16, "label") <- "Antidepressives"
attr(covs_drugs$m17, "label") <- "Antipsychotics"
attr(covs_drugs$m18, "label") <- "Agents for COPD"
attr(covs_drugs$n_drugs, "label") <- "No. drugs in past yr"

saveRDS(covs_drugs, paste0(cohort_path, "covs_drugs.rds"))

rm(all_lmdb)

## Health care
out_diag <- all_diag |>
  filter(C_PATTYPE == "2") |>
  mutate(`_year` = year(D_INDDTO)) |>
  filter(`_year` < 2014 | (`_year` >= 2014 & C_INDM == 2))

# table with primary contacts
pri_contacts <- sssy |>
  filter(2013 <= SSSY_year & SSSY_year <= 2020) |>
  filter(KONTAKT == 1 & YDTYP %in% c("01", "02", "03", "04", "05")) |> # contact with primary care doctor
  filter(str_detect(SPEC2, "7[0-9]|80")) |> # speciality in general medicine
  arrange(PNR, HONUGE)

rm(sssy)

# Inpatient (days of inpatient hospitalization in year before indexdate (dd-1/mm/yyyy-1 - dd/mm/yyyy))
# lpr2:
cov_inp_days_lpr2 <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(t_adm, by = c("pnr" = "PNR")) |>
  filter(
    (indexdate - years(1) + 1 <= D_INDDTO & D_INDDTO < indexdate) |
      (indexdate - years(1) + 1 <= D_UDDTO & D_UDDTO < indexdate) |
      (D_INDDTO < indexdate - years(1) + 1 & indexdate < D_UDDTO)
  ) |>
  filter(C_PATTYPE == "0") |>
  arrange(id, D_INDDTO) |>
  mutate(# Only count each day once
    D_INDDTO = pmax(D_INDDTO, lag(D_UDDTO + 1), na.rm = TRUE), 
    .by = "id"
  ) |>
  filter(D_INDDTO <= D_UDDTO) |>
  summarise(# count both in date and out date as day of inpatient hospitalization
    n_hosp = n(),
    n_hosp_days = sum(
      1L + 
        as.integer(difftime(pmin(D_UDDTO, indexdate), pmax(D_INDDTO, indexdate - years(1) + 1), unit = "days"))
    ),
    .by = "id"
  ) |>
  arrange(id)

# lpr3:
# Define inpatient proxy as hospital contacts lasting more than 12 hours in a given date or stay overnight. 
# For contacts lasting multiple days, both in and out dates are counted as inpatient days.
# Define outpatient proxy as hospital contacts lasting less than 12 hours in a given date (at least one contact)
# contacts lasting multiple days are not counted towards outpatient contacts.
cov_hosp_days_lpr3 <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(
    t_kontakt_lpr3 |> 
      select("PNR", "enhedstype_ans", "dato_start", "tidspunkt_start", "dato_slut", "tidspunkt_slut", "kontakttype"), 
    by = c("pnr" = "PNR")
  ) |>
  filter(enhedstype_ans %in% c("klinisk enhed", "hospital")) |> # filter on hospital contacts
  filter(kontakttype == "ALCA00") |>
  filter(
    (indexdate - years(1) + 1 <= dato_start & dato_start < indexdate) |
      (indexdate - years(1) + 1 <= dato_slut & dato_slut < indexdate) |
      (dato_start < indexdate - years(1) + 1 & indexdate < dato_slut)
  ) |>
  filter(!is.na(dato_start) & !is.na(dato_slut) & !is.na(tidspunkt_start) & !is.na(tidspunkt_slut)) |>
  # limit contacts to within 1 year of index date:
  mutate(
    dato_start = pmax(dato_start, indexdate - years(1) + 1),
    tidspunkt_start = ifelse(
      indexdate - years(1) + 1 <= dato_start, 
      tidspunkt_start, 
      0
    ) |>
      structure(class = c("hms", "difftime"), units = "secs"),
    dato_slut = pmin(dato_slut, indexdate),
    tidspunkt_slut = ifelse(
      dato_slut < indexdate, 
      tidspunkt_slut, 
      0
    ) |>
      structure(class = c("hms", "difftime"), units = "secs")
  ) |>
  reframe(
    date_start_time = (\(datetime_start, datetime_end) {
      tbl <- tibble(datetime_start = datetime_start, datetime_end = datetime_end) |>
        arrange(datetime_start)
      n <- 1
      out <- dmy_hms("01-01-2024 00:00:00")
      for(i in seq_len(nrow(tbl))) {
        # Check for overlapping contacts and remove any overlaps. If a contact fully overlaps a previous contact
        # move to next contact.
        if (i > 1) {
          if (tbl$datetime_start[i] < tbl$datetime_end[i-1]) {
            if (tbl$datetime_end[i-1] < tbl$datetime_end[i]) {
              tbl$datetime_start[i] <- tbl$datetime_end[i-1]
            } else {
              next
            }
          }
        } 
        # If contact spans multiple days, create a separate row for each date, starting at midnight.
        # Note: this may create overlap with a previous contact. however, since the date is an inpatient day,
        # added contact time will not matter.
        if (date(tbl$datetime_start[i]) < date(tbl$datetime_end)[i]) {
          for (j in seq_len(as.numeric(date(tbl$datetime_end[i]) - date(tbl$datetime_start[i])) + 1)) {
            out[n] <- as_datetime(date(tbl$datetime_start[i]) + j - 1)
            n <- n + 1
          }
          tbl$datetime_end[i] <- max(tbl$datetime_end[i], as_datetime(date(tbl$datetime_end[i])) + hours(12))
        } else { 
          # If contact spans single date, start is latest of contact start time or previous contacts end time
          out[n] <- tbl$datetime_start[i]
          n <- n + 1
        }
      }
      return(out)
    })(datetime_start = ymd_hms(paste0(dato_start, tidspunkt_start)), datetime_end = ymd_hms(paste0(dato_slut, tidspunkt_slut))),
    date_end_time = (\(datetime_start, datetime_end) {
      tbl <- tibble(datetime_start = datetime_start, datetime_end = datetime_end) |>
        arrange(datetime_start)
      n <- 1
      out <- dmy_hms("01-01-2024 00:00:00")
      for(i in seq_len(nrow(tbl))) {
        # Check for overlapping contacts and remove any overlaps. If a contact fully overlaps a previous contact
        # move to next contact.
        if (i > 1) {
          if (tbl$datetime_start[i] < tbl$datetime_end[i-1]) {
            if (tbl$datetime_end[i-1] < tbl$datetime_end[i]) {
              tbl$datetime_start[i] <- tbl$datetime_end[i-1]
            } else {
              next
            }
          }
        } 
        # If contact spans multiple days, create a separate row for each date, ending earliest at noon for the out day.
        # Note: this may create overlap with a previous contact. however, since the date is an inpatient day,
        # added contact time will not matter.
        if (date(tbl$datetime_start[i]) < date(tbl$datetime_end)[i]) {
          for (j in seq_len(as.numeric(date(tbl$datetime_end[i]) - date(tbl$datetime_start[i])) + 1)) {
            if (j == (as.numeric(date(tbl$datetime_end[i]) - date(tbl$datetime_start[i])) + 1)) {
              out[n] <- max(
                tbl$datetime_end[i], 
                as_datetime(date(tbl$datetime_end[i])) + hours(12)
              )
            } else {
              out[n] <- as_datetime(date(tbl$datetime_start[i]) + j)
            }
            n <- n + 1
          }
          tbl$datetime_end[i] <- max(tbl$datetime_end[i], as_datetime(date(tbl$datetime_end[i])) + hours(12))
        } else { 
          # If contact spans single date, end is contact end time
          # all contacts start before index date, so contacts spanning only one date cannot end on or after index date
          out[n] <- tbl$datetime_end[i]
          n <- n + 1
        }
      }
      return(out)
    })(
      datetime_start = ymd_hms(paste0(dato_start, tidspunkt_start)), 
      datetime_end = ymd_hms(paste0(dato_slut, tidspunkt_slut))
    ),
    .by = "id"
  ) |>
  mutate(
    date_start = date(date_start_time)
  ) |>
  summarise(
    hosp_time = as.numeric(sum(difftime(date_end_time, date_start_time, units = "hours"))),
    .by = c("id", "date_start")
  ) |>
  summarise(
    n_hosp_days = sum(hosp_time >= 12),
    n_c02 = sum(hosp_time < 12),
    n_hosp = NA,
    .by = "id"
  )

# combine:
cov_inp_days <- bind_rows(
  cov_inp_days_lpr2,
  cov_hosp_days_lpr3 |> select("id", "n_hosp", "n_hosp_days")
) |>
  summarise(
    n_hosp = sum(n_hosp, na.rm = TRUE),
    n_hosp_days = sum(n_hosp_days),
    .by = "id"
  ) |>
  right_join(basic_cohort |> select("id"), by = c("id")) |>
  mutate(
    n_hosp = replace_na(n_hosp, 0),
    n_hosp_days = replace_na(n_hosp_days, 0),
  ) |>
  arrange(id)


# outpatient (number of days with outpatient contacts)
cov_out <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(out_diag, by = c("pnr" = "PNR")) |>
  filter(indexdate - years(1) + 1 <= D_INDDTO & D_INDDTO < indexdate) |>
  summarise(
    n_c02 = length(unique(D_INDDTO)),
    .by = "id"
  ) |>
  bind_rows(
    cov_hosp_days_lpr3 |> select("id", "n_c02")
  ) |>
  summarise(
    n_c02 = sum(n_c02),
    .by = "id"
  ) |>
  right_join(basic_cohort |> select("id"), by = c("id")) |>
  mutate(n_c02 = replace_na(n_c02, 0)) |>
  arrange(id)
attr(cov_out$n_c02, "label") <- "No of outpatient contacts"

# Primary care (number of weeks with primary care contacts)
cov_pri <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(pri_contacts, by = c("pnr" = "PNR")) |>
  mutate(
    week = as.numeric(str_sub(HONUGE, 3, 4)),
    d_date = as_date(
      paste0(
        SSSY_year, 
        case_when(
          1 + (week - 1) * 7 < 10 ~ "00", 
          1 + (week - 1) * 7 < 100 ~ "0",
          TRUE ~ ""
        ), 
        1 + (week - 1) * 7), 
      format = "%Y%j"
    ) # format date from year day format: yyyydd
  ) |>
  filter(indexdate - years(1) + 1 <= d_date & d_date < indexdate) |>
  summarise(
    n_c03 = length(unique(d_date)),
    .by = "id"
  ) |>
  right_join(basic_cohort |> select("id"), by = c("id")) |>
  mutate(n_c03 = replace_na(n_c03, 0)) |>
  arrange(id)
attr(cov_pri$n_c03, "label") <- "No of primary care contacts"

# combine health care tables
covs_utilization <- cov_inp_days |>
  inner_join(cov_out, by = "id") |>
  inner_join(cov_pri, by = "id") |>
  rename(
    "n_c09" = "n_hosp",
    "n_c01" = "n_hosp_days"
  ) 
attr(covs_utilization$n_c01, "label") <- "No of days of inpatient hospitalization"
attr(covs_utilization$n_c09, "label") <- "No of inpatient contacts"

saveRDS(covs_utilization, paste0(cohort_path, "covs_utilization.rds"))
rm(out_diag, pri_contacts, cov_inp_days_lpr2, cov_hosp_days_lpr3, cov_inp_days, cov_out, cov_pri)

## Demographics
# Education
udda <- readRDS(paste0(data_path, "UDDF202009.rds")) |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active")

hfaudd_niveau <- read.delim(paste0(cohort_path, "hfaudd_niveau.txt", sep = ";"))

uddannelse <- udda |>
  filter(HF_VFRA > dmy("01/01/1900")) |>
  left_join(hfaudd_niveau, by = "HFAUDD") |>
  right_join(basic_cohort, by = c("PNR" = "pnr")) |>
  filter(HF_VFRA < indexdate) |>
  select("pnr" = "PNR", "udd_level") |>
  distinct() |>
  arrange(pnr, udd_level) |>
  summarise(udd_level = last(udd_level), .by = "pnr") |>
  right_join(basic_cohort, by = "pnr") |>
  mutate(
    pnr = pnr,
    HAE = case_when(
      str_detect(udd_level, "^[12]0") ~ "01 Primary Edu",
      str_detect(udd_level, "^30") ~ "02 Secondary Edu",
      str_detect(udd_level, "^[56]0") ~ "03 Bachelor",
      str_detect(udd_level, "^[78]0") ~ "04 Master or >",
      TRUE ~ "05 Missing"
    ),
    .keep = "none"
  )

saveRDS(uddannelse, paste0(cohort_path, "uddannelse.rds"))
rm(udda, hfaudd_niveau)

# Income
family <- readRDS(paste0(data_path, "BEF1985_2020.rds")) |>
  select("PNR", "FAMILIE_ID", "year" = "bef_year") |> 
  distinct() |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active")

saveRDS(family, paste0(cohort_path, "family.rds"))

income <- readRDS(paste0(data_path, "IND1980_2020.rds")) |>
  distinct() |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active")

# get family disposable income by linking with family id and summing over disposable income of each family member
fam_income <- as_tibble(
  as.data.table(family)[
    as.data.table(select(income, "PNR", "DISPON_13", "year" = "Ind_year"))[
      , 
      list(DISPON_13 = sum(DISPON_13, na.rm = TRUE)), 
      by = list(PNR, year)
    ],
    on = list(PNR, year)
  ][
    ,
    list(famdisponibel_13 = sum(DISPON_13, na.rm = TRUE)),
    by = list(FAMILIE_ID, year)
  ]
) |>
  filter(!is.na(FAMILIE_ID))

# create table with one year and one person in rows with family id, birth date and that years disposable family income
fam_income_2 <- family |>
  filter(year > 2002) |>
  left_join(select(t_person, "PNR_active", "D_FODDATO"),  by = c("PNR" = "PNR_active")) |>
  left_join(fam_income, by = c("FAMILIE_ID", "year"))

rm(t_person, fam_income)

# table with number of adults in each family in each year, the family disposable income for each year and family income
# weighted by the number of family members (0.5 + 0.5 * n_>14 + 0.3 * n_<15)
fam_income_3 <- as_tibble(
  as.data.table(fam_income_2)[
    ,
    list(
      n_adult = sum((year - year(D_FODDATO) - 1) > 17),
      yr_eq_in = mean(famdisponibel_13) / 
        (0.5 + sum((year - year(D_FODDATO) - 1) > 14) * 0.5 + sum((year - year(D_FODDATO) - 1) < 15) * 0.3)
    ),
    by = list(FAMILIE_ID, year)
  ]
)

saveRDS(fam_income_3, paste0(cohort_path, "fam_income_3.rds"))
rm(fam_income_2)

# table with five year rolling sum of family income and number of adult-years earning that income within those years
# restrict to years after 2007
# using data.table dcast and rowSums
fam_income_4 <- completeDT(as.data.table(fam_income_3), cols = c("FAMILIE_ID", "year"))
fam_income_5 <- dcast(
  fam_income_4, 
  FAMILIE_ID ~ year,
  drop = FALSE,
  value.var = "yr_eq_in"
)
fam_income_5_voks <- dcast(
  fam_income_4, 
  FAMILIE_ID ~ year,
  drop = FALSE,
  value.var = "n_adult"
)

# map over 5 years at a time and take rowsums of income and number of adults
future::plan(sequential)
fam_income_6 <- furrr::future_map(
  .x = 2003:2020,
  .options = furrr_options(
    globals = c("fam_income_5"),
    packages = c("data.table")
  ),
  .f = \(year) {
    year_offsets <- seq(min(year - 2003, 4), 0, -1)
    data <- fam_income_5[, j, env = list(j = as.list(as.character(year - year_offsets)))]
    return(rowMeans(data, na.rm = TRUE))
  }
) |>
  structure(names = 2003:2020) |>
  as_tibble() |>
  mutate(
    FAMILIE_ID = fam_income_5[, FAMILIE_ID]
  ) |>
  as.data.table() |>
  melt(
    id.vars = "FAMILIE_ID",
    measure.vars = as.character(2003:2020),
    variable.name = "year",
    value.name = "income"
  ) |>
  as_tibble()
fam_income_6_voks <- furrr::future_map(
  .x = 2003:2020,
  .options = furrr_options(
    globals = c("fam_income_5_voks"),
    packages = c("data.table")
  ),
  .f = \(year) {
    year_offsets <- seq(min(year - 2003, 4), 0, -1)
    data <- fam_income_5_voks[, j, env = list(j = as.list(as.character(year - year_offsets)))]
    return(rowMeans(data, na.rm = TRUE))
  }
) |>
  structure(names = 2003:2020) |>
  as_tibble() |>
  mutate(
    FAMILIE_ID = fam_income_5_voks[, FAMILIE_ID]
  ) |>
  as.data.table() |>
  melt(
    id.vars = "FAMILIE_ID",
    measure.vars = as.character(2003:2020),
    variable.name = "year",
    value.name = "n_voks"
  ) |>
  as_tibble()

# add average income per adult per year in last 5 years in household
fam_income_7 <- fam_income_6 |>
  inner_join(fam_income_6_voks, by = c("FAMILIE_ID", "year")) |>
  mutate(
    year = as.character(year),
    income_std = ifelse(is.nan(income), NA, income / n_voks) # normalize to income per adult-year
  ) 

# calculate quintiles of income for each year after 2007
incomes <- fam_income_7 |>
  filter(year > 2007) |>
  mutate(
    income_quintiles = factor(case_when(
      income_std < quantile(income_std, 0.2, na.rm = TRUE) ~ 1L,
      income_std < quantile(income_std, 0.4, na.rm = TRUE) ~ 2L,
      income_std < quantile(income_std, 0.6, na.rm = TRUE) ~ 3L,
      income_std < quantile(income_std, 0.8, na.rm = TRUE) ~ 4L,
      income_std < max(income_std, na.rm = TRUE) ~ 5L,
      TRUE ~ NA_integer_
    )),
    .by = "year"
  ) |>
  mutate(year = as.integer(year))

# on year of indexdate, determine which quintile of family income each individual belongs to.
equi_incomes <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  mutate(indexyear = year(indexdate)) |>
  left_join(family, by = c("pnr" = "PNR", "indexyear" = "year")) |>
  left_join(incomes, by = c("FAMILIE_ID", "indexyear" = "year"))

saveRDS(equi_incomes, paste0(cohort_path, "equi_incomes.rds"))
rm(
  family, income, fam_income_3, fam_income_4, fam_income_5, fam_income_5_voks, 
  fam_income_6, fam_income_6_voks, fam_income_7, incomes
)

# Address
cohort_address <- basic_cohort |> 
  select("id", "pnr", "indexdate") |>
  left_join(address, by = c("pnr" = "PNR")) |>
  filter(D_TILFLYT_DATO <= indexdate | is.na(D_TILFLYT_DATO)) |>
  slice_max(tibble(as.numeric(D_TILFLYT_DATO), -as.numeric(D_FRAFLYT_DATO)), n = 1, by = "id", with_ties = FALSE) |>
  right_join(basic_cohort |> select("id", "pnr", "indexdate"), by = c("id", "pnr", "indexdate"))

saveRDS(cohort_address, paste0(cohort_path, "cohort_address.rds"))
rm(address)

# combine demographic information
covs_demos <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  left_join(cohort_address |> select("id", "p01" = "region"), by = "id") |>
  left_join(uddannelse |> select("pnr", "p02" = "HAE"), by = "pnr") |>
  left_join(equi_incomes |> select("id", "p04" = "income_std",  "p03" = "income_quintiles"), by = "id") |>
  mutate(
    y01 = (year(indexdate) == 2014),
    y02 = (year(indexdate) == 2015),
    y03 = (year(indexdate) == 2016),
    y04 = (year(indexdate) == 2017),
    y05 = (year(indexdate) == 2018),
    y06 = (year(indexdate) == 2019),
    y07 = (year(indexdate) == 2020)
  )
attr(covs_demos$p01, "label") <- "Region of residence"
attr(covs_demos$p02, "label") <- "Education"
attr(covs_demos$p03, "label") <- "Household income quintile"
attr(covs_demos$p04, "label") <- "Household income numeric"
attr(covs_demos$y01, "label") <- "Calender year: 2014"
attr(covs_demos$y02, "label") <- "Calender year: 2015"
attr(covs_demos$y03, "label") <- "Calender year: 2016"
attr(covs_demos$y04, "label") <- "Calender year: 2017"
attr(covs_demos$y05, "label") <- "Calender year: 2018"
attr(covs_demos$y06, "label") <- "Calender year: 2019"
attr(covs_demos$y07, "label") <- "Calender year: 2020"

saveRDS(covs_demos, paste0(cohort_path, "covs_demos.rds"))
rm(cohort_address, uddannelse, equi_incomes)

## Lab covariates
# load lab-data and select distinct combinations of relevant columns
future::plan("multisession", workers = 10)
lab <- furrr::future_map(
  .x = 1:10,
  .options = furrr_options(
    globals = c("data_path"),
    packages = c("glue", "dplyr")
  ),
  .f = \(i) {
    data <- readRDS(glue("{data_path}LAB_dm_forsker_{i}.rds")) |>
      filter(
        ANALYSISCODE %in% 
          c(
            "NPU03429", "DNK35131", "DNK35302", "NPU04998", "NPU18016", "NPU03230", "NPU02319", "NPU27300",
            "NPU19651", "NPU27783", "NPU19655", "NPU19673", "NPU01132", "NPU01446", "NPU04144", "NPU02647",
            "NPU03096", "NPU19658", "NPU01459", "NPU03568", "NPU02593", "NPU19763", "NPU10267", "NPU01423",
            "NPU19748", "NPU01566", "NPU18412", "NPU01568", "NPU10171"
          )
      )
    return(data)
  }
) |>
  list_rbind() |>
  left_join(pnr_active, by = "PNR") |>
  select(-c("PNR")) |>
  rename("PNR" = "PNR_active") |>
  distinct(PNR, SAMPLINGDATE, SAMPLINGTIME, LABORATORIUM_IDCODE, ANALYSISCODE, VALUE) |>
  mutate(
    value_num = as.numeric(VALUE),
    sampling_year = year(SAMPLINGDATE)
  )
future::plan("sequential")

# join with basic cohort and character column with type of lab-test
tmp_lab_pop <- basic_cohort |>
  select("id", "pnr", "indexdate") |>
  inner_join(lab, by = c("pnr" = "PNR")) |>
  filter(indexdate - years(5) <= SAMPLINGDATE & SAMPLINGDATE < indexdate) |>
  mutate(
    labtest = case_when(
      str_detect(ANALYSISCODE, "^NPU03429") ~ "Sodium",
      str_detect(ANALYSISCODE, "^DNK35131") ~ "mdrd_eGFR",
      str_detect(ANALYSISCODE, "^DNK35302") ~ "ckdepi_eGFR",
      str_detect(ANALYSISCODE, "^NPU(04998|18016)") ~ "kreatinin",
      str_detect(ANALYSISCODE, "^NPU03230") ~ "Potassium",
      str_detect(ANALYSISCODE, "^NPU02319") ~ "Hemoglobin",
      str_detect(ANALYSISCODE, "^NPU27300") ~ "HBA1c",
      str_detect(ANALYSISCODE, "^NPU19651") ~ "ALAT",
      str_detect(ANALYSISCODE, "^NPU(27783|19655)") ~ "BASP",
      str_detect(ANALYSISCODE, "^NPU(19673|01132)") ~ "Albumin",
      str_detect(ANALYSISCODE, "^NPU(01446|04144)") ~ "Calcium",
      str_detect(ANALYSISCODE, "^NPU02647") ~ "Magnesium",
      str_detect(ANALYSISCODE, "^NPU03096") ~ "Phosphate",
      str_detect(ANALYSISCODE, "^NPU19658") ~ "LDH",
      str_detect(ANALYSISCODE, "^NPU01459") ~ "Carbamide",
      str_detect(ANALYSISCODE, "^NPU03568") ~ "Thrombocytes",
      str_detect(ANALYSISCODE, "^NPU02593") ~ "Leukocytes",
      str_detect(ANALYSISCODE, "^NPU19763") ~ "Ferritin",
      str_detect(ANALYSISCODE, "^NPU10267") ~ "Vitamin_D",
      str_detect(ANALYSISCODE, "^NPU(01423|19748)") ~ "CRP",
      str_detect(ANALYSISCODE, "^NPU(01566|18412)") ~ "Cholesterol",
      str_detect(ANALYSISCODE, "^NPU(01568|10171)") ~ "LDL",
      TRUE ~ NA_character_
    )
  ) |>
  # Remove or replace relational or missing values (Still a few missing)
  mutate(
    VALUE = str_replace(VALUE, ",", ".")
  ) |>
  filter(
    !(VALUE %in% c("<12", "Negativ", "?", "<2", "<3")) 
  ) |>
  mutate(
    VALUE = case_when(
      VALUE == "<1.8" ~ "1.8",
      VALUE == "<10" ~ "0",
      VALUE == "<4" ~ "0",
      VALUE == ">300" ~ "300",
      TRUE ~ VALUE
    ),
    value_num = case_when(
      is.na(value_num) ~ as.numeric(VALUE),
      TRUE ~ value_num
    )
  ) |>
  # take latest lab-test taken before index date 
  slice_max(SAMPLINGDATE, by = c("id", "labtest"))

# pivot wide on value and date of lab test
tmp_lab_pop_2_1 <- tmp_lab_pop |>
  select("id", "labtest", "value_num") |>
  pivot_wider(names_from = "labtest", values_from = "value_num", values_fn = mean)

tmp_lab_pop_2_2 <- tmp_lab_pop |>
  select("id", "labtest", "SAMPLINGDATE") |>
  pivot_wider(names_from = "labtest", values_from = "SAMPLINGDATE", names_glue = "date_{labtest}", values_fn = first)

# create table with kidney test values (eGFR)
`kidney_` <- basic_cohort |>
  select("id", "sex", "indexdate", "dob") |>
  left_join(tmp_lab_pop_2_1 |> select("id", "kreatinin", "mdrd_eGFR", "ckdepi_eGFR"), by = "id") |>
  left_join(tmp_lab_pop_2_2 |> select("id", "date_kreatinin", "date_mdrd_eGFR", "date_ckdepi_eGFR"), by = "id") |>
  mutate(
    female = (sex == "K"),
    age = as.numeric(difftime(indexdate, dob, unit = "days")) / 365.25
  )

kidney4 <- `kidney_` |>
  # add egfr
  mutate(
    egfr = case_when(
      is.na(kreatinin) & is.na(mdrd_eGFR) & is.na(ckdepi_eGFR) ~ NA_real_,
      date_kreatinin == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~
        female       * (kreatinin <= 62) * (144 * (kreatinin / (0.7 * 88.4))^(-0.329) * 0.993^age) +
        female       * (kreatinin > 62)  * (144 * (kreatinin / (0.7 * 88.4))^(-1.209) * 0.993^age) +
        (1 - female) * (kreatinin <= 80) * (141 * (kreatinin / (0.9 * 88.4))^(-0.411) * 0.993^age) +
        (1 - female) * (kreatinin > 80)  * (141 * (kreatinin / (0.9 * 88.4))^(-1.209) * 0.993^age),
      date_ckdepi_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ ckdepi_eGFR,
      date_mdrd_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~
        (\(creat_est) {
          (mdrd_eGFR == 90) * 90 +
            (mdrd_eGFR != 90) * female       * (creat_est <= 62) * (144 * (creat_est / (0.7 * 88.4))^(-0.329) * 0.993^age) +
            (mdrd_eGFR != 90) * female       * (creat_est > 62)  * (144 * (creat_est / (0.7 * 88.4))^(-1.209) * 0.993^age) +
            (mdrd_eGFR != 90) * (1 - female) * (creat_est <= 80) * (141 * (creat_est / (0.9 * 88.4))^(-0.411) * 0.993^age) +
            (mdrd_eGFR != 90) * (1 - female) * (creat_est > 80)  * (141 * (creat_est / (0.9 * 88.4))^(-1.209) * 0.993^age)
        })(creat_est = 88.4 * ((175 * age^(-0.203) * (1 + (-1 + 0.742) * female)) / mdrd_eGFR)^(1/1.154))
    ),
    egfr_90 = case_when(
      is.na(kreatinin) & is.na(mdrd_eGFR) & is.na(ckdepi_eGFR) ~ NA_real_,
      date_kreatinin == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ NA_real_,
      date_ckdepi_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ ckdepi_eGFR == 90,
      date_mdrd_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ mdrd_eGFR == 90
    )
  ) |>
  # add renal function (low egfr or normal)
  mutate(
    renal_function = case_when(
      is.na(egfr) ~ "Missing",
      egfr < 30 ~ "Low_lt30",
      TRUE ~ "Normal"
    )
  ) |>
  # add egfr type
  mutate(
    date_egfr = case_when(
      is.na(date_kreatinin) & is.na(date_mdrd_eGFR) & is.na(date_ckdepi_eGFR) ~ NA_Date_,
      TRUE ~ pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE)
    ),
    egfr_type = case_when(
      is.na(date_kreatinin) & is.na(date_mdrd_eGFR) & is.na(date_ckdepi_eGFR) ~ "Missing",
      date_kreatinin == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ "creat",
      date_ckdepi_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ "ckd_epi",
      date_mdrd_eGFR == pmax(date_kreatinin, date_mdrd_eGFR, date_ckdepi_eGFR, na.rm = TRUE) ~ "mdrd"
    )
  )

cohort_kidney <- basic_cohort |>
  left_join(kidney4 |> select("id", "egfr", "renal_function", "date_egfr", "egfr_type"), by = "id") |>
  mutate(egfr = round(egfr))

saveRDS(cohort_kidney, paste0(cohort_path, "cohort_kidney.rds"))
rm(tmp_lab_pop, tmp_lab_pop_2_1, tmp_lab_pop_2_2, `kidney_`, kidney4)

# combine lab test covariates
covs_lab <- cohort_kidney |>
  select("id", "egfr", "renal_function", "date_egfr", "egfr_type") |>
  left_join(
    tmp_lab_pop_4_1 |>
      select(
        "id",
        "l01" = "Sodium", 
        "l03" = "Potassium",
        "l04" = "Hemoglobin",
        "l05" = "ALAT",
        "l06" = "BASP",
        "l07" = "Albumin",
        "l08" = "LDH",
        "l09" = "Carbamide",
        "l10" = "Thrombocytes",
        "l11" = "Leukocytes",
        "l12" = "CRP",
        "l13" = "Cholesterol",
        "l14" = "LDL",
      ),
    by = "id"
  ) |>
  left_join(
    tmp_lab_pop_4_2 |>
      select(
        "id",
        "date_l01" = "date_Sodium", 
        "date_l03" = "date_Potassium",
        "date_l04" = "date_Hemoglobin",
        "date_l05" = "date_ALAT",
        "date_l06" = "date_BASP",
        "date_l07" = "date_Albumin",
        "date_l08" = "date_LDH",
        "date_l09" = "date_Carbamide",
        "date_l10" = "date_Thrombocytes",
        "date_l11" = "date_Leukocytes",
        "date_l12" = "date_CRP",
        "date_l13" = "date_Cholesterol",
        "date_l14" = "date_LDL",
      ), 
    by = "id"
  ) |>
  arrange(id)
attr(covs_lab$l01, "label") <- attr(covs_lab$date_l01, "label") <-"Sodium"
attr(covs_lab$l03, "label") <- attr(covs_lab$date_l03, "label") <-"Potassium"
attr(covs_lab$l04, "label") <- attr(covs_lab$date_l04, "label") <-"Hemoglobin"
attr(covs_lab$l05, "label") <- attr(covs_lab$date_l05, "label") <-"ALAT"
attr(covs_lab$l06, "label") <- attr(covs_lab$date_l06, "label") <-"Alkaline phosphatase"
attr(covs_lab$l07, "label") <- attr(covs_lab$date_l07, "label") <-"Albumin"
attr(covs_lab$l08, "label") <- attr(covs_lab$date_l08, "label") <-"LDH"
attr(covs_lab$l09, "label") <- attr(covs_lab$date_l09, "label") <-"Carbamide"
attr(covs_lab$l10, "label") <- attr(covs_lab$date_l10, "label") <-"Thrombocytes"
attr(covs_lab$l11, "label") <- attr(covs_lab$date_l11, "label") <-"Leukocytes"
attr(covs_lab$l12, "label") <- attr(covs_lab$date_l12, "label") <-"CRP"
attr(covs_lab$l13, "label") <- attr(covs_lab$date_l13, "label") <-"Cholesterol"
attr(covs_lab$l14, "label") <- attr(covs_lab$date_l14, "label") <-"LDL-C"

saveRDS(covs_lab, paste0(cohort_path, "covs_lab.rds"))
rm(cohort_kidney)

# combine all covariates in one table
all_covs <- select(basic_cohort, "id", "pnr", "indexdrug", "sex", "reg_midt_censor_date", "indexdate", "dob") |>
  mutate(age = as.numeric(difftime(indexdate, dob, unit = "days")) / 365.25) |>
  select(-c("indexdate", "dob")) |>
  left_join(select(covs_demos, "id", "reg" = "p01", "edu" = "p02", "inc" = "p03", "inc_num" = "p04", starts_with("y0")), by = "id") |>
  left_join(select(covs_diag, "id", starts_with("d")), by = "id") |>
  left_join(select(x_hyponatremia, "id", "d26"), by = "id") |>
  left_join(select(covs_drugs, "id", starts_with("m"), "n_drugs", "n_drug_grp"), by = "id") |>
  left_join(covs_lab, by = "id") |>
  left_join(select(covs_utilization, "id", "n_c01", "n_c02", "n_c03"), by = "id") |>
  left_join(select(add_on_hyp_start_stop, -c("indexdrug")), by = "id") |>
  left_join(
    select(
      tmp_basic_cohort_1, 
      "id", "atc", "drugcat", "drugcat2", "drugcat3", "drugcat4", "indo", "strnum", "strunit", "x_multiple",
      "x_region_midt", "x_no_address", "x01", "x02", "x03", "x04"
    ), 
    by = "id"
  )
attr(all_covs$age, "label") <- "Age"
attr(all_covs$reg, "label") <- "Region"
attr(all_covs$inc_num, "label") <- "Income numeric"
attr(all_covs$d26, "label") <- "Previous hyponatremia"

# missingness (mode impute education and income)
all_covs <- all_covs |>
  mutate(
    edu_missing = ifelse(str_detect(edu, "^05"), "missing", "not missing"),
    inc_missing = ifelse(is.na(inc), "missing", "not missing"),
    edu = ifelse(edu_missing == "missing", "02 Secondary Edu", edu),
    inc = factor(
      ifelse(inc_missing == "missing", 5, inc),
      levels = 1:5,
      labels = paste0("Q", 1:5)
    ) 
  )
attr(all_covs$edu, "label") <- "Education"
attr(all_covs$inc, "label") <- "Income quintile"
attr(all_covs$edu_missing, "label") <- "Missing education"
attr(all_covs$inc_missing, "label") <- "Missing Income"

# group health care utilization
all_covs <- all_covs |>
  mutate(
    h01 = factor(
      case_when(
        n_c01 == 0 ~ 0L,
        n_c01 >= 1 & n_c01 <= 7 ~ 1L,
        n_c01 >= 8 & n_c01 <= 14 ~ 2L,
        n_c01 >= 15 ~ 3L,
      ),
      levels = 0:3,
      labels = c("0", "1-7", "8-14", "\u226515")
    ),
    h02 = factor(
      case_when(
        n_c02 == 0 ~ 0L,
        n_c02 >= 1 & n_c02 <= 3 ~ 1L,
        n_c02 >= 4 & n_c02 <= 6 ~ 2L,
        n_c02 >= 7 ~ 3L,
      ),
      levels = 0:3,
      labels = c("0", "1-3", "4-6", "\u22657")
    ),
    h03 = factor(
      case_when(
        n_c03 == 0 ~ 0L,
        n_c03 >= 1 & n_c03 <= 3 ~ 1L,
        n_c03 >= 4 & n_c03 <= 6 ~ 2L,
        n_c03 >= 7 ~ 3L,
      ),
      levels = 0:3,
      labels = c("0", "1-3", "4-6", "\u22657")
    )
  )
attr(all_covs$h01, "label") <- "No of days of inpatient hospitalization"
attr(all_covs$h02, "label") <- "No of outpatient contacts"
attr(all_covs$h03, "label") <- "No of primary care contacts"

saveRDS(all_covs, paste0(cohort_path, "all_covs.rds"))
rm(
  x_hyponatremia, tmp_basic_cohort_1, basic_cohort, add_on_hyp_start_stop, covs_diag, covs_drugs,
  covs_utilization, covs_demos, covs_lab
)

## Outcomes ------------------------------------------------------------------------------------------------------- ----
# Select sodium measurements from lab data
sodium_x <- lab |>
  filter(ANALYSISCODE == "NPU03429") |>
  filter(!(VALUE %in% c("?", "Negativ", "Positiv"))) |>
  mutate(
    hypo_na = case_when(
      value_num >= 130 & value_num < 135  ~ "01_mild",
      value_num >= 125 & value_num < 130  ~ "02_mode",
      value_num < 125                     ~ "03_seve",
      value_num >= 135 & value_num <= 145 ~ "04_norm",
      value_num > 145 & value_num < 150   ~ "05_high",
      value_num >= 150                    ~ "06_extr",
      TRUE ~ NA_character_
    ),
    any_na = "007_any"
  )

# join sodium measurements in follow-up period with cohort
tmp_sodium_pop <- basic_cohort |>
  select("id", "pnr", "entry", "exit") |>
  inner_join(
    sodium_x |> 
      select("PNR", "SAMPLINGDATE", "ANALYSISCODE", "value_num", "hypo_na", "any_na"), 
    by = c("pnr" = "PNR")
  ) |>
  filter(entry < SAMPLINGDATE & SAMPLINGDATE <= exit) |>
  select("id", "pnr", "SAMPLINGDATE", "ANALYSISCODE", "value_num", "hypo_na", "any_na")

# inpatient admissions from LPR2
tmp_hosp_in <- t_adm |>
  filter(year(D_INDDTO) > 2013 & (C_PATTYPE == "0" | (C_PATTYPE == "2" & C_INDM == "1"))) |>
  mutate(types = ifelse(C_PATTYPE == "0", "in", "ed"))

# combine inpatient admissions and sodium measurements on same date
hosp_in_all_lpr2 <- tmp_sodium_pop |>
  inner_join(tmp_hosp_in, by = c("pnr" = "PNR"), relationship = "many-to-many") |>
  filter(SAMPLINGDATE == D_INDDTO & value_num < 130) |>
  select("id", "pnr", "SAMPLINGDATE", "value_num", "hypo_na", "types") 

# inpatient proxy from LPR3
hosp_in_all_lpr3 <- t_kontakt_lpr3 |> 
  select("PNR", "enhedstype_ans", "dato_start", "tidspunkt_start", "dato_slut", "tidspunkt_slut", "kontakttype") |>
  filter(enhedstype_ans %in% c("klinisk enhed", "hospital")) |>
  filter(kontakttype == "ALCA00") |>
  inner_join(
    tmp_sodium_pop |> select("id", "pnr", "value_num", "hypo_na", "SAMPLINGDATE"), 
    by = c("PNR" = "pnr"), 
    relationship = "many-to-many"
  ) |>
  filter(!is.na(dato_start) & !is.na(dato_slut) & !is.na(tidspunkt_start) & !is.na(tidspunkt_slut)) |>
  filter(SAMPLINGDATE == dato_start & value_num < 130) |>
  mutate(
    datetime_start = ymd_hms(paste0(dato_start, tidspunkt_start)), 
    datetime_end = ymd_hms(paste0(dato_slut, tidspunkt_slut))
  ) |>
  mutate(
    datetime_end = as_datetime(ifelse(dato_start < dato_slut, as_datetime(dato_start + 1), datetime_end))
  ) |>
  reframe(
    pnr = dplyr::first(PNR),
    date_start_time = (\(datetime_start, datetime_end) {
      tbl <- tibble(datetime_start = datetime_start, datetime_end = datetime_end) |>
        arrange(datetime_start)
      n <- 1
      out <- dmy_hms("01-01-2024 00:00:00")
      for(i in seq_len(nrow(tbl))) {
        # Check for overlapping contacts and remove any overlaps. If a contact fully overlaps a previous contact
        # move to next contact.
        if (i > 1) {
          if (tbl$datetime_start[i] < tbl$datetime_end[i-1]) {
            if (tbl$datetime_end[i-1] >= tbl$datetime_end[i]) {
              next
            }
            tbl$datetime_start[i] <- tbl$datetime_end[i-1]
          }
        } 
        
        out[n] <- tbl$datetime_start[i]
        n <- n + 1
      }
      return(out)
    })(datetime_start = datetime_start, datetime_end = datetime_end),
    date_end_time = (\(datetime_start, datetime_end) {
      tbl <- tibble(datetime_start = datetime_start, datetime_end = datetime_end) |>
        arrange(datetime_start)
      n <- 1
      out <- dmy_hms("01-01-2024 00:00:00")
      for(i in seq_len(nrow(tbl))) {
        # Check for overlapping contacts and remove any overlaps. If a contact fully overlaps a previous contact
        # move to next contact.
        if (i > 1) {
          if (tbl$datetime_start[i] < tbl$datetime_end[i-1]) {
            if (tbl$datetime_end[i-1] >= tbl$datetime_end[i]) {
              next
            }
          }
        } 
        
        out[n] <- tbl$datetime_end[i]
        n <- n + 1
      }
      return(out)
    })(datetime_start = datetime_start, datetime_end = datetime_end),
    value_num = dplyr::first(value_num),
    hypo_na = dplyr::first(hypo_na),
    .by = c("id", "SAMPLINGDATE")
  ) |>
  summarise(
    pnr = dplyr::first(pnr),
    value_num = dplyr::first(value_num),
    hypo_na = dplyr::first(hypo_na),
    hosp_time = as.numeric(sum(difftime(date_end_time, date_start_time, units = "hours"))),
    .by = c("id", "SAMPLINGDATE")
  ) |>
  filter(hosp_time >= 12) |>
  mutate(types = "in") |>
  select("id", "pnr", "SAMPLINGDATE", "value_num", "hypo_na", "types")

# combine lpr2 and lpr3
hosp_in_all <- bind_rows(hosp_in_all_lpr2, hosp_in_all_lpr3)

# table with outcomes from each sodium measurement in the follow-up period
fx_outcomes <- basic_cohort |>
  select("id", "pnr", "entry", "exit") |>
  left_join(tmp_sodium_pop, by = c("id", "pnr")) |>
  mutate(
    mild = str_detect(hypo_na, "^01"),
    mode = str_detect(hypo_na, "^02"),
    seve = str_detect(hypo_na, "^03"),
    hosp = FALSE,
    norm = str_detect(hypo_na, "^04"),
    high = str_detect(hypo_na, "^05"),
    extr = str_detect(hypo_na, "^06"),
    any =  str_detect(any_na, "^007")
  ) |>
  select(
    "id", "pnr", "SAMPLINGDATE", "value_num", "hypo_na", "mild", "mode", "seve", "hosp", "norm", "high", "extr", "any"
  ) |>
  filter(mild | mode | seve | hosp | norm | high | extr) |>
  bind_rows(
    basic_cohort |>
      select("id", "pnr", "entry", "exit") |>
      left_join(hosp_in_all, by = c("id", "pnr")) |>
      filter(entry < SAMPLINGDATE & SAMPLINGDATE < exit) |>
      mutate(
        mild = FALSE,
        mode = FALSE,
        seve = FALSE,
        hosp = str_detect(hypo_na, "^0[23]"),
        norm = FALSE,
        high = FALSE,
        extr = FALSE,
        any =  TRUE
      ) |>
      select(
        "id", "pnr", "SAMPLINGDATE", "value_num", "hypo_na", "mild", 
        "mode", "seve", "hosp", "norm", "high", "extr", "any"
      ) |>
      filter(mild | mode | seve | hosp | norm | high | extr)
  )

# Table with death in cohort in follow-up period
death <- basic_cohort |>
  select("id", "pnr", "entry", "exit") |>
  inner_join(t_person |> select("pnr" = "PNR_active", "eventdate" = "D_STATUS_HEN_START", "C_STATUS"), by = "pnr") |>
  filter(entry < eventdate & eventdate < exit & C_STATUS == "90") |>
  mutate(eventtype = "07_death") |>
  select(-c("C_STATUS"))

# Table with all earliest events (columns eventtype and eventdate)
all_outcomes <- bind_rows(
  fx_outcomes |>
    filter(mode | seve) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "01_hypo_na_primary",
      .by = "id"
    ),
  fx_outcomes |>
    filter(mild | mode | seve) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "02_hypo_na_secondary",
      .by = "id"
    ),
  fx_outcomes |>
    filter(seve) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "03_severe_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(mode) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "04_moderate_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(mild) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "05_mild_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(hosp) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "06_hosp_na",
      .by = "id"
    ),
  death |>
    summarise(
      pnr = first(pnr),
      eventdate = min(eventdate),
      eventtype = eventtype,
      .by = "id"
    ),
  fx_outcomes |>
    filter(norm) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "08_norm_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(high | extr) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "09_high_na", # NOTE: includes both high and extreme sodium levels
      .by = "id"
    ),
  fx_outcomes |>
    filter(extr) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "10_extr_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(any) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "11_any_na",
      .by = "id"
    ),
  fx_outcomes |>
    filter(norm | high | extr) |>
    summarise(
      pnr = first(pnr),
      eventdate = min(SAMPLINGDATE),
      eventtype = "12_\u2265135_na",
      .by = "id"
    )
)
saveRDS(all_outcomes, paste0(cohort_path, "all_outcomes.rds"))

# Table with first sodium measurement after exposure drug initiation
sodium_continuous <- tmp_sodium_pop |>
  select("id", "value_num", "SAMPLINGDATE") |>
  slice_min(SAMPLINGDATE, n = 1, by = "id", with_ties = FALSE) |>
  rename(
    "first_on_drug_sodium" = "value_num",
    "date_first_on_drug_sodium" = "SAMPLINGDATE"
  )
saveRDS(sodium_continuous, paste0(cohort_path, "sodium_continuous.rds"))

## Study cohort --------------------------------------------------------------------------------------------------- ----
tmp_study_cohort <- basic_cohort |>
  select("id", "indexdate", "indexdrug", "entry", "exit", "stat", "stdt") |>
  left_join(all_covs |> select(,"id", "hyp_stop", "add_on_date_any_hyp", "reg_midt_censor_date"), by = "id")

study_cohort <- tmp_study_cohort |>
  mutate(
    exit = pmin(exit, entry + years(5), hyp_stop, add_on_date_any_hyp, reg_midt_censor_date, na.rm = TRUE),
    dead_exit = ifelse(!is.na(stdt) & stat == "90" & exit == (stdt + days(1)), 1, 0),
    emigrated_exit = ifelse(!is.na(stdt) & stat == "80" & exit == (stdt + days(1)), 1, 0),
    hyp_stop_add_on_exit = ifelse((!is.na(add_on_date_any_hyp) & exit == add_on_date_any_hyp) | exit == hyp_stop, 1, 0),
    reg_midt_censor_date_exit = ifelse(!is.na(reg_midt_censor_date) & exit == reg_midt_censor_date, 1, 0)
  ) |>
  filter(entry < exit) |>
  select(-c("stat", "stdt", "hyp_stop", "add_on_date_any_hyp", "reg_midt_censor_date"))

outcome_input_01 <- all_outcomes |>
  inner_join(study_cohort, by = "id") |>
  filter(str_detect(eventtype, "^01") & entry < eventdate & eventdate < exit) |>
  select("id", "pnr", "eventdate", "eventtype")

outcome_input_03 <- all_outcomes |>
  inner_join(study_cohort, by = "id") |>
  filter(str_detect(eventtype, "^03") & entry < eventdate & eventdate < exit) |>
  select("id", "pnr", "eventdate", "eventtype")

# time to outcome with <130mM events or censoring
out_01 <- study_cohort |>
  left_join(outcome_input_01 |> select("id", "eventdate"), by = "id") |>
  mutate(
    event_01 = ifelse(!is.na(eventdate) & eventdate < exit, 1, 0),
    t_01 = as.integer(difftime(pmin(eventdate, exit, na.rm = TRUE), entry, units = "days"))
  ) |>
  select("id", "indexdrug", "indexdate", "t_01", "event_01")

# time to outcome with <125mM events or censoring
out_03 <- study_cohort |>
  left_join(outcome_input_03 |> select("id", "eventdate"), by = "id") |>
  mutate(
    event_03 = ifelse(!is.na(eventdate) & eventdate < exit, 1, 0),
    t_03 = as.integer(difftime(pmin(eventdate, exit, na.rm = TRUE), entry, units = "days"))
  ) |>
  select("id", "indexdrug", "indexdate", "t_03", "event_03")

# Table with first sodium measurement after exposure drug initiation in study cohort
sodium_continuous_out <- sodium_continuous |>
  inner_join(study_cohort, by = "id") |>
  filter(entry < date_first_on_drug_sodium & date_first_on_drug_sodium < exit) |>
  select("id", "first_on_drug_sodium", "date_first_on_drug_sodium")

# Current study requires <130mM events and <125mM events
cohort_stw <- study_cohort |>
  select("id", "indexdate", "dead_exit", "emigrated_exit", "hyp_stop_add_on_exit", "reg_midt_censor_date_exit") |>
  left_join(out_01 |> select("id", "event_130" = "event_01", "t_130" = "t_01"), by = "id") |>
  left_join(out_03 |> select("id", "event_125" = "event_03", "t_125" = "t_03"), by = "id") |>
  left_join(all_covs, by = "id") |>
  left_join(sodium_continuous_out, by = "id") |>
  mutate(
    sex = ifelse(sex == "K", "Female", "Male"),
    low_egfr = !is.na(egfr) & egfr < 30
  )
attr(cohort_stw$event_130, "label") <- "Hyponatremia (<130mM) event"
attr(cohort_stw$`t_130`, "label") <- "Follow-up for (<130mM) event"
attr(cohort_stw$event_125, "label") <- "Hyponatremia (<125mM) event"
attr(cohort_stw$`t_125`, "label") <- "Follow-up for (<125mM) event"
saveRDS(cohort_stw, paste0(cohort_path, "cohort_all.rds"))

# filter on exclusion criteria
study_cohort <- cohort_stw |>
  filter(!(x_multiple | x_region_midt | x_no_address | x01 | x02 | x03 | x04 | low_egfr)) |>
  select(-c("low_egfr"))
saveRDS(study_cohort, paste0(cohort_path, "study_cohort.rds"))

# clean up R session
rm(list = ls())
gc()
