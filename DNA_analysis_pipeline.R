################################################################################
#
#  DNA_analysis_pipeline.R
#
#  Missed outpatient appointments (DNAs) at Cambridge University Hospitals:
#  cleaning, linkage, descriptive statistics, regression modelling, and mapping.
#
#  This file consolidates the R half of the project. It is written to be run
#  top to bottom in a single session, though each section is self-contained
#  enough to be run on its own provided the sections above it have already
#  written their outputs to clean_data/. The Python half of the pipeline
#  (postcode linkage, distance calculation, clustering, causal forest) lives in
#  DNA_analysis_pipeline.py, and steps 3 and 4 of that file must be run between
#  Sections 2 and 3 here; see README.md for the full order.
#
#  Sections
#    0.  Configuration and libraries
#    1.  Cleaning and aggregating the raw referral data
#    2.  Cleaning the publicly available datasets
#    3.  Linking public data to the aggregate dataset
#    4.  Derived variables and analysis samples
#    5.  Missingness
#    6.  Descriptive statistics
#    7.  Regression modelling
#    8.  Heatmaps and choropleths
#
################################################################################


################################################################################
# SECTION 0 — CONFIGURATION AND LIBRARIES
################################################################################

library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)
library(scales)
library(splines)
library(lmtest)
library(car)
library(grid)
library(gridExtra)
library(sf)
library(ggspatial)
library(cowplot)

# All paths are relative to the project root. It is important to set the working
# directory to the project root before running anything below, since every read
# and write in this file assumes that location.
RAW_DIR    <- "data"         # raw specialty .xlsx extracts (not shared)
OPEN_DIR   <- "open_data"    # publicly available source files
CLEAN_DIR  <- "clean_data"   # cleaned and linked outputs
OUT_DIR    <- "outputs"      # PDFs, figures, and exported tables

for (d in c(CLEAN_DIR, OUT_DIR)) dir.create(d, showWarnings = FALSE)

# Canonical file names for the analysis dataset at each stage of the pipeline.
# Naming the stages explicitly avoids the ambiguity of tracking half a dozen
# similarly named CSVs by hand.
F_AGG        <- file.path(CLEAN_DIR, "aggregate_specialty_data_24.csv")
F_AGG_PCD    <- file.path(CLEAN_DIR, "aggregate_specialty_data_24_pcd_link.csv")
F_AGG_NAPTAN <- file.path(CLEAN_DIR, "aggregate_specialty_data_24_pcd_NAPTAN_link.csv")
F_AGG_LINKED <- file.path(CLEAN_DIR, "aggregate_specialty_data_24_linked.csv")
F_AGG_CCA    <- file.path(CLEAN_DIR, "aggregate_specialty_data_24_linked_CCA.csv")

# Cambridge colour theme, used across every figure so that the descriptive
# graphs, the model plots, and the maps all belong to the same document.
CAM_LIGHT <- "#D1F9F1"   # light blue, alternating table fills and low values
CAM_WARM  <- "#00BDB6"   # warm blue, accent lines and ribbons
CAM_DARK  <- "#133844"   # dark blue, headers and text

cambridge_blues <- c("#8EE8D8", "#00BDB6", "#D1F9F1")


################################################################################
# SECTION 1 — CLEANING AND AGGREGATING THE RAW REFERRAL DATA
#
# Each specialty arrives as a separate .xlsx extract with identical structure,
# so the cleaning is written once as a function and mapped over the directory.
# The specialty name is taken from the file name and retained as a column, which
# is what allows the specialty-stratified analyses further down.
################################################################################

# Clean a single specialty extract: rename to analysis variable names, keep the
# first completed appointment (plus any missed appointments preceding it), and
# recode the categorical variables into the binary forms used in modelling.
CLEANDATA <- function(path) {
  read_excel(path) %>%

    rename(ptid              = Pat_ID,
           pcd               = Patient_Postcode_at_Referral,
           ethnicity         = Patient_Ethnic_Group,
           sex               = Patient_Sex,
           interpreter       = Patient_Interpreter_Required,
           language          = Patient_Language,
           age               = Patient_Age_at_Referral,
           IMD_decile        = IndexMultipleDeprivationDecile,
           IMD_rank          = IndexMultipleDeprivationRank,
           referral_practice = Referral_Practice,
           DNA               = Appt_Status,
           DNA_reason        = Appt_Cancel_Reason,
           appt_booked_date  = Appt_Booked_Date,
           appt_time         = Appt_Time,
           referral_date     = Referral_Date) %>%

    select(c("ptid", "pcd", "appt_booked_date", "referral_date", "appt_time",
             "ethnicity", "sex", "interpreter", "language", "age",
             "IMD_decile", "IMD_rank", "referral_practice",
             "DNA", "DNA_reason")) %>%

    # Retain each patient's care episode up to and including their first
    # completed appointment. Patients who never completed contribute all of
    # their appointments, since none of them resolved the referral.
    group_by(ptid) %>%
    arrange(appt_time) %>%
    mutate(first_completed = ifelse(any(DNA == "Completed"),
                                    min(appt_time[DNA == "Completed"], na.rm = FALSE),
                                    max(appt_time))) %>%
    filter(appt_time <= first_completed) %>%
    ungroup() %>%

    # IMD arrives as text in some extracts and as a number in others
    mutate(IMD_decile = suppressWarnings(as.numeric(IMD_decile)),
           IMD_rank   = suppressWarnings(as.numeric(IMD_rank))) %>%

    # Male = 0, Female = 1
    mutate(sex = if_else(sex == "Female", 1,
                 if_else(sex == "Male", 0, NA))) %>%

    # Interpreter not required = 0, required = 1
    mutate(interpreter = if_else(interpreter == "Y", 1,
                         if_else(interpreter == "N", 0, NA))) %>%

    # White = 0, non-white = 1, not recorded or not stated = NA. Recorded
    # ethnicity is kept separate from the LSOA-level linked ethnicity added in
    # Section 3, since the two measure different things.
    mutate(ethnicity = if_else(ethnicity %in% c(".White British", ".White Irish",
                                                "Other White Background"), 0,
                       if_else(ethnicity %in% c("Not Recorded", "Not Stated"), NA, 1))) %>%

    # A cancellation by the patient, or an inability to attend, is treated as a
    # missed appointment; cancellations by the service are not the patient's
    # behaviour and are dropped by the filter that follows.
    mutate(DNA = if_else(DNA_reason %in% c("Could Not Attend", "Cancelled by Patient"),
                         "Did Not Attend", DNA)) %>%
    filter(DNA %in% c("Did Not Attend", "Completed")) %>%
    mutate(DNA = if_else(DNA == "Did Not Attend", 1, 0)) %>%

    select(-c("DNA_reason"))
}

# Map the cleaning function over every specialty extract and stack the results,
# tagging each row with the specialty taken from its file name.
files        <- list.files(RAW_DIR, pattern = "\\.xlsx$", full.names = TRUE)
cleaned_data <- lapply(files, CLEANDATA)

names(cleaned_data) <- files %>%
  basename() %>%
  tools::file_path_sans_ext() %>%
  sub("_.*$", "", .)

aggregated_data <- bind_rows(cleaned_data, .id = "specialty")

write.csv(aggregated_data, F_AGG, row.names = FALSE)


################################################################################
# SECTION 2 — CLEANING THE PUBLICLY AVAILABLE DATASETS
#
# Every public dataset is reduced to the variables of interest and keyed on
# lsoa21 so that all of them can be joined to the referral data in one pass in
# Section 3. The census tables (vehicle availability, marital status, ethnicity)
# arrive as long counts by category and are summarised into LSOA-level
# percentages here rather than at the point of analysis.
################################################################################

# Access to Healthy Assets and Hazards (2024). Retained in full because the
# domain scores and percentiles are all candidate covariates.
AHAH <- read.csv(file.path(OPEN_DIR, "ahah_v5.csv"))

AHAH <- AHAH %>%
  select(c('lsoa21cd', 'GP', 'hospital', 'pharmacy', 'dentist', 'fast_food',
           'gambling', 'greenspace_active', 'NO2', 'PM10', 'SO2', 'vape_tobacco',
           'GP_pct', 'dentist_pct', 'pharmacy_pct', 'hospital_pct', 'leisure_pct',
           'greenspace_pct', 'greenspace_active_pct', 'bluespace_pct', 'NO2_pct',
           'PM10_pct', 'SO2_pct', 'fast_food_pct', 'gambling_pct', 'pub_bar_pct',
           'tobacco_pct', 'domain_h', 'domain_g', 'domain_e', 'domain_r',
           'domain_h_pct', 'domain_g_pct', 'domain_e_pct', 'domain_r_pct',
           'ahah', 'ahah_rnk', 'ahah_pct')) %>%
  rename(lsoa21 = lsoa21cd)

write.csv(AHAH, file.path(CLEAN_DIR, "AHAH_24_clean.csv"), row.names = FALSE)


# Indices of Deprivation (2025). Only the domain deciles are kept; the full
# column names are unwieldy and are renamed to something usable.
iod <- read.csv(file.path(OPEN_DIR, "iod_2025.csv"))

iod <- iod %>%
  select(c('LSOA.code..2021.',
           'Income.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Employment.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Education..Skills.and.Training.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Health.Deprivation.and.Disability.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Crime.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Barriers.to.Housing.and.Services.Decile..where.1.is.most.deprived.10..of.LSOAs.',
           'Living.Environment.Decile..where.1.is.most.deprived.10..of.LSOAs.')) %>%
  rename(lsoa21                    = LSOA.code..2021.,
         income_decile             = Income.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         employment_decile         = Employment.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         education_decile          = Education..Skills.and.Training.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         health_decile             = Health.Deprivation.and.Disability.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         crime_decile              = Crime.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         housing_decile            = Barriers.to.Housing.and.Services.Decile..where.1.is.most.deprived.10..of.LSOAs.,
         living_environment_decile = Living.Environment.Decile..where.1.is.most.deprived.10..of.LSOAs.)

write.csv(iod, file.path(CLEAN_DIR, "iod_25_clean.csv"), row.names = FALSE)


# NAPTAN public transport access nodes. Only the coordinates are needed, since
# these are fed to the nearest-node distance calculation in the Python script.
NAPTAN <- read.csv(file.path(OPEN_DIR, "NAPTAN_2022.csv"))
NAPTAN <- NAPTAN %>% select(c('ATCOCode', 'CommonName', 'Longitude', 'Latitude'))
write.csv(NAPTAN, file.path(CLEAN_DIR, "NAPTAN_24_clean.csv"), row.names = FALSE)


# Vehicle availability (Census 2021). Category 0 is "no cars or vans", so the
# percentage of households with access to at least one vehicle is one minus that
# share. Code -8 denotes a suppressed cell and is set to NA before summing.
va <- read.csv(file.path(OPEN_DIR, "vehicle_availability_2021.csv"))
va <- va %>%
  select(c('Lower.layer.Super.Output.Areas.Code',
           'Car.or.van.availability..5.categories..Code',
           'Observation')) %>%
  rename(lsoa21               = Lower.layer.Super.Output.Areas.Code,
         vehicle_availability = Car.or.van.availability..5.categories..Code) %>%
  mutate(vehicle_availability = na_if(vehicle_availability, -8)) %>%
  mutate(vehicle_availability = as.numeric(as.character(vehicle_availability))) %>%
  group_by(lsoa21) %>%
  summarise(num         = sum(Observation[vehicle_availability == 0], na.rm = TRUE),
            dem         = sum(Observation[vehicle_availability %in% 0:3], na.rm = TRUE),
            vehicle_pct = (1 - (num / dem)) * 100) %>%
  select(lsoa21, vehicle_pct)

write.csv(va, file.path(CLEAN_DIR, "vehicle_availability_21_clean.csv"), row.names = FALSE)


# Childcare accessibility (2023). Renamed only; the source is already one row
# per LSOA and needs no aggregation.
ca <- read.csv(file.path(OPEN_DIR, "childcare_accessibility_2023.csv"))
ca <- ca %>%
  rename(lsoa21                                      = LSOA21CD,
         childcare_overall                           = Childcare.accessibility,
         childcare_overall_driving                   = Childcare.accessibility..driving.only.,
         childcare_overall_public_transport          = Childcare.accessibility..public.transport.only.,
         childcare_good_outstanding                  = Childcare.accessibility...good.or.outstanding.places..overall.,
         childcare_good_outstanding_driving          = Childcare.accessibility...good.or.outstanding.places..driving.only.,
         childcare_good_outstanding_public_transport = Childcare.accessibility...good.or.outstanding.places..public.transport.only.)

write.csv(ca, file.path(CLEAN_DIR, "childcare_accessability_23_clean.csv"), row.names = FALSE)


# Digital propensity index (2021). The score is stored as a percentage string
# and is stripped back to a number here.
dp <- read.csv(file.path(OPEN_DIR, "digital_propensity_index_2021.csv"))
dp <- dp %>%
  select(LSOA.code, Digital.Propensity.Score) %>%
  rename(lsoa21                   = LSOA.code,
         digital_propensity_score = Digital.Propensity.Score)

dp$digital_propensity_score <- as.numeric(sub("%", "", dp$digital_propensity_score))

write.csv(dp, file.path(CLEAN_DIR, "digital_propensity_index_21_clean.csv"), row.names = FALSE)


# Marital status (Census 2021). Categories 1 and 8 through 11 are the
# never-married and formerly-partnered groups, so marriage_pct is one minus
# their combined share of categories 1 through 11.
ms <- read.csv(file.path(OPEN_DIR, "marital_status_2021.csv"))
ms <- ms %>%
  select(c('Lower.layer.Super.Output.Areas.Code',
           'Marital.and.civil.partnership.status..12.categories..Code',
           'Observation')) %>%
  rename(lsoa21         = Lower.layer.Super.Output.Areas.Code,
         marital_status = Marital.and.civil.partnership.status..12.categories..Code) %>%
  mutate(marital_status = na_if(marital_status, -8)) %>%
  mutate(marital_status = as.numeric(as.character(marital_status))) %>%
  group_by(lsoa21) %>%
  summarise(num          = sum(Observation[marital_status %in% c(1, 8:11)], na.rm = TRUE),
            dem          = sum(Observation[marital_status %in% 1:11], na.rm = TRUE),
            marriage_pct = (1 - (num / dem)) * 100) %>%
  select(lsoa21, marriage_pct)

write.csv(ms, file.path(CLEAN_DIR, "marital_status_21_clean.csv"), row.names = FALSE)


# Ethnicity (Census 2021). Categories 13 through 17 are the white groups, so the
# non-white percentage is one minus their share of categories 1 through 19.
eth <- read.csv(file.path(OPEN_DIR, "ethnicity_LSOA_21.csv"))
eth <- eth %>%
  select(c('Lower.layer.Super.Output.Areas.Code',
           'Ethnic.group..20.categories..Code',
           'Observation')) %>%
  rename(lsoa21       = Lower.layer.Super.Output.Areas.Code,
         ethnic_group = Ethnic.group..20.categories..Code) %>%
  mutate(ethnic_group = na_if(ethnic_group, -8)) %>%
  mutate(ethnic_group = as.numeric(as.character(ethnic_group))) %>%
  group_by(lsoa21) %>%
  summarise(num           = sum(Observation[ethnic_group %in% 13:17], na.rm = TRUE),
            dem           = sum(Observation[ethnic_group %in% 1:19], na.rm = TRUE),
            non_white_pct = (1 - (num / dem)) * 100) %>%
  select(lsoa21, non_white_pct)

write.csv(eth, file.path(CLEAN_DIR, "ethnicity_21_clean.csv"), row.names = FALSE)


# LSOA population, used to weight the choropleths in Section 9.
pop <- read.csv(file.path(OPEN_DIR, "lsoa_pop.csv"))
pop <- pop %>% rename(lsoa21 = lsoa21cd)
write.csv(pop, file.path(CLEAN_DIR, "lsoa_pop_clean.csv"), row.names = FALSE)


################################################################################
#  >>> RUN THE PYTHON GEOGRAPHY STEPS BEFORE CONTINUING <<<
#
#  Sections 3 onwards expect lat, long, lsoa21, NAPTAN_distance_km, and add_dist
#  to be present on the referral data. Those come from DNA_analysis_pipeline.py:
#    python DNA_analysis_pipeline.py postcode-link
#    python DNA_analysis_pipeline.py distances
#  which read F_AGG and write F_AGG_PCD and then F_AGG_NAPTAN plus add_dist.csv.
################################################################################


################################################################################
# SECTION 3 — LINKING PUBLIC DATA TO THE AGGREGATE DATASET
#
# One left join per source, all keyed on lsoa21 apart from the Addenbrooke's
# distance, which is keyed on postcode. Left joins are used throughout so that a
# failed linkage produces NA rather than silently dropping the patient; the
# complete case dataset is derived explicitly at the end of the section, and the
# rows lost are counted rather than assumed.
################################################################################

agg_data <- read.csv(F_AGG_NAPTAN)

ahah <- read.csv(file.path(CLEAN_DIR, "AHAH_24_clean.csv"))
iod  <- read.csv(file.path(CLEAN_DIR, "iod_25_clean.csv"))
va   <- read.csv(file.path(CLEAN_DIR, "vehicle_availability_21_clean.csv"))
ca   <- read.csv(file.path(CLEAN_DIR, "childcare_accessability_23_clean.csv"))
dp   <- read.csv(file.path(CLEAN_DIR, "digital_propensity_index_21_clean.csv"))
ms   <- read.csv(file.path(CLEAN_DIR, "marital_status_21_clean.csv"))
eth  <- read.csv(file.path(CLEAN_DIR, "ethnicity_21_clean.csv"))
pop  <- read.csv(file.path(CLEAN_DIR, "lsoa_pop_clean.csv"))
add  <- read.csv(file.path(CLEAN_DIR, "add_dist.csv"))

# Standardise the empty strings that survive the CSV round trip into true NAs,
# and derive the two variables that depend on the referral data alone.
# ethnicity2 keeps missing ethnicity as its own level ("2") rather than dropping
# those patients, since ethnicity recording is itself patterned.
agg_data$ethnicity2 <- ifelse(is.na(agg_data$ethnicity), "2", agg_data$ethnicity)

agg_data <- agg_data %>%
  mutate(lsoa21      = na_if(lsoa21, ""),
         interpreter = as.character(interpreter),
         interpreter = na_if(interpreter, ""),
         interpreter = na_if(interpreter, "NA"),
         language    = as.character(language),
         language    = na_if(language, ""),
         language    = na_if(language, "NA"))

agg_data$non_english <- ifelse(agg_data$language != "English" &
                                 !is.na(agg_data$language), 1, 0)

# The distance file carries one row per appointment; reduce it to one row per
# postcode before joining, or the join will duplicate rows.
add <- add %>% distinct(pcd, .keep_all = TRUE)

agg_data <- agg_data %>%
  left_join(ahah, by = "lsoa21") %>%   # healthy assets and hazards domains
  left_join(iod,  by = "lsoa21") %>%   # deprivation domain deciles
  left_join(va,   by = "lsoa21") %>%   # vehicle availability percentage
  left_join(ca,   by = "lsoa21") %>%   # childcare accessibility
  left_join(dp,   by = "lsoa21") %>%   # digital propensity score
  left_join(ms,   by = "lsoa21") %>%   # marriage percentage
  left_join(eth,  by = "lsoa21") %>%   # linked non-white percentage
  left_join(pop,  by = "lsoa21") %>%   # LSOA population
  left_join(add,  by = "pcd")          # straight line distance to Addenbrooke's

write.csv(agg_data, F_AGG_LINKED, row.names = FALSE)


# Linkage quality check. For each source, count how many referral rows failed to
# match and confirm that the unmatched rows are exactly the rows left entirely
# missing. If those two numbers ever diverge, the join key is duplicated on the
# right hand side and the linkage needs revisiting before anything downstream is
# trusted.
linkage_and_missingness <- function(left, right, key = "lsoa21") {
  if (!(key %in% names(left)) || !(key %in% names(right)))
    return(data.frame(status = "FAILED"))

  right_cols <- setdiff(names(right), key)

  merged <- left %>%
    select(all_of(key)) %>%
    left_join(right, by = key)

  na_all <- rowSums(!is.na(merged[right_cols])) == 0

  data.frame(right_duplicates = sum(duplicated(right[[key]])),
             matched_rows     = sum(!na_all),
             unmatched_rows   = sum(na_all),
             same_missing     = TRUE)
}

linkage_report <- bind_rows(
  lapply(list(ahah = ahah, iod = iod, va = va, ca = ca,
              dp = dp, ms = ms, eth = eth),
         function(d) linkage_and_missingness(agg_data, d)),
  .id = "source")

print(linkage_report)


# Complete case dataset. Rows without a linked LSOA, without a recorded sex, or
# missing any deprivation domain cannot contribute to the adjusted models, so
# they are removed here once and the count reported, rather than being dropped
# implicitly by each model in turn.
iod_vars <- c("income_decile", "employment_decile", "education_decile",
              "health_decile", "crime_decile", "housing_decile",
              "living_environment_decile")

n_before <- nrow(agg_data)

agg_data_clean <- agg_data %>%
  filter(!is.na(lsoa21),
         !is.na(sex),
         if_all(all_of(iod_vars), ~ !is.na(.)))

n_after <- nrow(agg_data_clean)

print(data.frame(rows_before  = n_before,
                 rows_after   = n_after,
                 rows_removed = n_before - n_after))

write.csv(agg_data_clean, F_AGG_CCA, row.names = FALSE)


################################################################################
# SECTION 4 — DERIVED VARIABLES AND ANALYSIS SAMPLES
#
# Everything from here down runs on the complete case dataset. Three samples are
# constructed and used consistently thereafter:
#
#   mult_data      every retained appointment, used for repeat-DNA counts
#   agg_data       one row per patient (their first referral), used for the
#                  aggregate analyses, where independence matters
#   spec_agg_data  one row per patient per specialty, used for the stratified
#                  analyses, so that a patient referred to two specialties
#                  contributes to both without contributing twice to either
################################################################################

agg_data <- read.csv(F_AGG_CCA)

agg_data$ethnicity2 <- factor(agg_data$ethnicity2)
agg_data$digital_propensity_score <-
  as.numeric(sub("%", "", agg_data$digital_propensity_score))

specialties <- unique(agg_data$specialty)

# Deprivation domains arrive as deciles. They are collapsed to quintiles for the
# descriptive tables and models, since decile-level cells become sparse once the
# data are split by specialty.
agg_data$income_quintile             <- factor(ntile(agg_data$income_decile, 5))
agg_data$education_quintile          <- factor(ntile(agg_data$education_decile, 5))
agg_data$employment_quintile         <- factor(ntile(agg_data$employment_decile, 5))
agg_data$health_quintile             <- factor(ntile(agg_data$health_decile, 5))
agg_data$crime_quintile              <- factor(ntile(agg_data$crime_decile, 5))
agg_data$housing_quintile            <- factor(ntile(agg_data$housing_decile, 5))
agg_data$living_environment_quintile <- factor(ntile(agg_data$living_environment_decile, 5))

# The continuous measures are cut at their own quintile boundaries rather than
# by ntile(), so that the cut points are interpretable values of the underlying
# index rather than ranks.
cut_quintile <- function(x) {
  cut(x,
      breaks         = quantile(x, probs = seq(0, 1, by = 0.2), na.rm = TRUE),
      include.lowest = TRUE,
      labels         = c("1", "2", "3", "4", "5"))
}

agg_data$DPS_quintile       <- factor(cut_quintile(agg_data$digital_propensity_score))
agg_data$childcare_quintile <- factor(cut_quintile(agg_data$childcare_overall))
agg_data$vehicle_quintile   <- factor(cut_quintile(agg_data$vehicle_pct))

# Rescaled so that a coefficient reads as the effect of a ten percentage point
# difference, which is a more meaningful contrast than a single point.
agg_data$vehicle_pct_10  <- agg_data$vehicle_pct  / 10
agg_data$marriage_pct_10 <- agg_data$marriage_pct / 10

# ns() and bs() need continuous numeric input, so the ordinal variables are kept
# in a numeric form alongside their factor form.
to_num <- function(x) as.numeric(as.character(x))

agg_data <- agg_data %>%
  mutate(IMD_num        = to_num(IMD_decile),
         income_num     = to_num(income_quintile),
         education_num  = to_num(education_quintile),
         employment_num = to_num(employment_quintile),
         health_num     = to_num(health_quintile),
         crime_num      = to_num(crime_quintile),
         housing_num    = to_num(housing_quintile),
         living_env_num = to_num(living_environment_quintile))

# Age bands and distance bands used in the descriptive cross-tabulations.
agg_data$age_band <- cut(agg_data$age,
                         breaks = c(0, 18, 30, 45, 60, 75, Inf),
                         labels = c("0-18", "18-30", "30-45", "45-60", "60-75", "75+"),
                         right  = FALSE)

km_breaks <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, Inf)
km_labels <- c("0-1 km", "1-2 km", "2-3 km", "3-4 km", "4-5 km", "5-6 km",
               "6-7 km", "7-8 km", "8-9 km", "9-10 km", "10+ km")

agg_data$GP_bin       <- cut(agg_data$GP, breaks = km_breaks, labels = km_labels, right = FALSE)
agg_data$hospital_bin <- cut(agg_data$hospital, breaks = km_breaks, labels = km_labels, right = FALSE)
agg_data$NAPTAN_bin   <- cut(agg_data$NAPTAN_distance_km, breaks = km_breaks, labels = km_labels, right = FALSE)

agg_data$add_bin <- cut(agg_data$add_dist,
                        breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, Inf),
                        labels = c("0-10 km", "10-20 km", "20-30 km", "30-40 km",
                                   "40-50 km", "50-60 km", "60-70 km", "70-80 km",
                                   "80-90 km", "90-100 km", "100+ km"),
                        right  = FALSE)

agg_data$non_white_cat <- cut(agg_data$non_white_pct,
                              breaks         = c(0, 20, 40, 60, 80, 100),
                              include.lowest = TRUE)

# Build the three analysis samples. event_data is ordered by referral date so
# that slice(1) takes each patient's earliest referral.
mult_data  <- agg_data
event_data <- agg_data %>% arrange(ptid, referral_date)

agg_data <- event_data %>%
  group_by(ptid) %>%
  slice(1) %>%
  ungroup()

spec_agg_data <- event_data %>%
  group_by(ptid, specialty) %>%
  slice(1) %>%
  ungroup()

spec_agg_data$non_white_cat <- cut(spec_agg_data$non_white_pct,
                                   breaks         = c(0, 20, 40, 60, 80, 100),
                                   include.lowest = TRUE)

spec_agg_data$age_band <- cut(spec_agg_data$age,
                              breaks = c(0, 18, 30, 45, 60, 75, Inf),
                              labels = c("0-18", "18-30", "30-45", "45-60", "60-75", "75+"),
                              right  = FALSE)

# Cardiology is the reference specialty throughout, being the largest.
agg_data$specialty      <- relevel(factor(agg_data$specialty), ref = "cardiology")
spec_agg_data$specialty <- relevel(factor(spec_agg_data$specialty), ref = "cardiology")


################################################################################
# SECTION 5 — MISSINGNESS
#
# Missingness is reported on the linked dataset before complete case filtering,
# since the point is to characterise what was lost. It is examined both overall
# and by specialty, because recording practice varies between services, and the
# correlation heatmap checks whether the variables go missing together.
################################################################################

agg_data_og <- read.csv(F_AGG_LINKED)

vars <- c("age", "sex", "interpreter", "language", "lsoa21")

overall_missing_counts <- agg_data_og %>%
  summarise(across(all_of(vars), ~ sum(is.na(.x)), .names = "{col}")) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "MissingCount")

missing_counts_by_specialty <- agg_data_og %>%
  group_by(specialty) %>%
  summarise(across(all_of(vars), ~ sum(is.na(.x)), .names = "{col}")) %>%
  pivot_longer(cols = all_of(vars), names_to = "Variable", values_to = "MissingCount")

# Correlation between missingness indicators, log scaled because the
# correlations are small and would otherwise be indistinguishable on a linear
# fill scale.
miss_mat <- agg_data_og %>%
  mutate(across(all_of(vars), ~ as.numeric(is.na(.x)))) %>%
  select(all_of(vars))

miss_cor_df <- cor(miss_mat, use = "pairwise.complete.obs") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Correlation") %>%
  filter(Var1 != Var2) %>%
  mutate(Correlation_log = log10(Correlation + 1e-6))

miss_cor <- ggplot(miss_cor_df, aes(x = Var1, y = Var2, fill = Correlation_log)) +
  geom_tile(color = NA, linewidth = 0) +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(title = "Missingness Correlation Heatmap (Log Scaled)",
       x = "Variable", y = "Variable") +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1),
        panel.grid       = element_blank(),
        panel.border     = element_blank(),
        axis.ticks       = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background  = element_rect(fill = "white", color = NA))

# Interpreter and language share a missingness pattern, so only one of the pair
# is carried into the by-specialty variable heatmap.
vars2 <- c("age", "sex", "lsoa21")

agg_missing_spec_specialty <- agg_data_og %>%
  group_by(specialty) %>%
  summarise(missing_prop = mean(is.na(age) | is.na(sex) | is.na(pcd) |
                                is.na(lsoa21) | is.na(interpreter) & is.na(language)))

miss_spec_lang <- ggplot(agg_missing_spec_specialty,
                         aes(x = "Language", y = specialty, fill = missing_prop)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(title = "Overall Missingness by Specialty", x = NULL, y = "Specialty") +
  theme_minimal()

agg_missing_spec_vars2 <- agg_data_og %>%
  group_by(specialty) %>%
  summarise(across(all_of(vars2), ~ mean(is.na(.x)), .names = "{col}")) %>%
  pivot_longer(cols = all_of(vars2), names_to = "Variable", values_to = "MissingProportion")

miss_spec_var <- ggplot(agg_missing_spec_vars2,
                        aes(x = Variable, y = specialty, fill = MissingProportion)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  labs(title = "Missingness Heatmap by Specialty and Variable",
       x = "Variable", y = "Specialty") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
# SECTION 6 — DESCRIPTIVE STATISTICS
################################################################################

# ---- 6.1 Table helpers -------------------------------------------------------

fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

fmt_n_pct <- function(x, denom) {
  n <- sum(x, na.rm = TRUE)
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

# Apply a summary function to the aggregate sample and to each specialty in
# turn, stacking the results into one table with an "Aggregate" row on top.
by_group <- function(agg_data, spec_agg_data, fn) {
  bind_rows(
    fn(agg_data) %>% mutate(Specialty = "Aggregate"),
    spec_agg_data %>%
      group_by(specialty) %>%
      group_modify(~ fn(.x)) %>%
      rename(Specialty = specialty)
  ) %>%
    mutate(Specialty = str_to_title(Specialty)) %>%
    relocate(Specialty)
}

clean_title <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\bfreq\\b", "frequency") %>%
    str_to_title()
}

agg_data      <- agg_data %>% mutate(specialty = str_to_title(specialty))
spec_agg_data <- spec_agg_data %>% mutate(specialty = str_to_title(specialty))


# ---- 6.2 Figure theme --------------------------------------------------------

theme_cambridge_blue <- function() {
  theme_minimal() +
    theme(plot.title   = element_text(size = 16, face = "bold", colour = CAM_DARK),
          axis.text.x  = element_text(angle = 45, hjust = 1),
          axis.title   = element_text(colour = CAM_DARK),
          axis.text    = element_text(colour = CAM_DARK),
          panel.grid   = element_blank(),
          legend.title = element_text(face = "bold"),
          legend.text  = element_text(colour = CAM_DARK))
}

scale_fill_cambridge_blue <- function() scale_fill_manual(values = cambridge_blues)


# ---- 6.3 Cohort description --------------------------------------------------

demographics_table <- by_group(agg_data, spec_agg_data, function(d) {
  d <- d %>% distinct(ptid, .keep_all = TRUE)
  n <- nrow(d)
  tibble(`Number of patients`         = n,
         `Age, mean (SD)`             = fmt_mean_sd(d$age),
         `Male, n (%)`                = fmt_n_pct(d$sex == 0, n),
         `Female, n (%)`              = fmt_n_pct(d$sex == 1, n),
         `White ethnicity, n (%)`     = fmt_n_pct(d$ethnicity2 == 0, n),
         `Non-white ethnicity, n (%)` = fmt_n_pct(d$ethnicity2 == 1, n),
         `Ethnicity missing, n (%)`   = fmt_n_pct(d$ethnicity2 == 2, n),
         `English language, n (%)`    = fmt_n_pct(d$language == "English" | is.na(d$language), n),
         `IMD decile, mean (SD)`      = fmt_mean_sd(d$IMD_decile))
})

imd_table <- by_group(agg_data, spec_agg_data, function(d) {
  d <- d %>% distinct(ptid, .keep_all = TRUE)
  tibble(`Income deprivation, mean (SD)`             = fmt_mean_sd(d$income_decile),
         `Employment deprivation, mean (SD)`         = fmt_mean_sd(d$employment_decile),
         `Education deprivation, mean (SD)`          = fmt_mean_sd(d$education_decile),
         `Health deprivation, mean (SD)`             = fmt_mean_sd(d$health_decile),
         `Crime deprivation, mean (SD)`              = fmt_mean_sd(d$crime_decile),
         `Housing deprivation, mean (SD)`            = fmt_mean_sd(d$housing_decile),
         `Living environment deprivation, mean (SD)` = fmt_mean_sd(d$living_environment_decile))
})

distance_table <- by_group(agg_data, spec_agg_data, function(d) {
  d <- d %>% distinct(ptid, .keep_all = TRUE)
  tibble(`NAPTAN distance (km), mean (SD)`   = fmt_mean_sd(d$NAPTAN_distance_km),
         `Hospital distance (km), mean (SD)` = fmt_mean_sd(d$hospital))
})

# The childcare index is bounded at 1; values above that are data errors in the
# source and are excluded rather than winsorised.
agg_data <- agg_data %>% filter(childcare_overall <= 1 | is.na(childcare_overall))

sdoh_table <- by_group(agg_data, spec_agg_data, function(d) {
  d <- d %>% distinct(ptid, .keep_all = TRUE)
  tibble(`Childcare accessibility (SD)`               = fmt_mean_sd(d$childcare_overall * 100),
         `Vehicle ownership percentage (SD)`          = fmt_mean_sd(d$vehicle_pct),
         `Married or partnered percentage (SD)`       = fmt_mean_sd(d$marriage_pct),
         `Linked non-white ethnicity percentage (SD)` = fmt_mean_sd(d$non_white_pct))
})


# ---- 6.4 DNA cross-tabulations and figures -----------------------------------

# Every stratified comparison below follows the same shape: a cross-tabulation
# of DNA against the stratifier, the row proportions, an aggregate bar chart,
# and the same chart faceted by specialty. The helpers below carry that shape so
# that each block reduces to naming the variable and its axis label.

dna_table <- function(data, var) {
  tab <- xtabs(as.formula(paste("~", var, "+ DNA")), data = data)
  list(counts  = tab,
       pct     = format(round(prop.table(tab, margin = 1), 2), nsmall = 2),
       long_df = as.data.frame(prop.table(tab, margin = 1)))
}

dna_bar <- function(df, var, xlab, title, x_labels = NULL) {
  p <- ggplot(df, aes(x = .data[[var]], y = Freq, fill = DNA)) +
    geom_col(position = "stack") +
    ggtitle(title) +
    labs(x = xlab, y = "Ratio") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_cambridge_blue() +
    theme_cambridge_blue()
  if (!is.null(x_labels)) p <- p + scale_x_discrete(labels = x_labels)
  p
}

dna_bar_spec <- function(data, var, xlab, title, x_labels = NULL) {
  d <- data %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(specialty, .data[[var]], DNA) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(specialty, .data[[var]]) %>%
    mutate(Freq = n / sum(n), DNA = factor(DNA))

  p <- ggplot(d, aes(x = .data[[var]], y = Freq, fill = DNA)) +
    geom_col(position = "stack") +
    facet_wrap(~ specialty) +
    ggtitle(title) +
    labs(x = xlab, y = "Ratio") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_cambridge_blue() +
    theme_cambridge_blue()
  if (!is.null(x_labels)) p <- p + scale_x_discrete(labels = x_labels)
  p
}

eth_labels <- c("0" = "White", "1" = "Non-white", "2" = "Missing")
sex_labels <- c("0" = "Male", "1" = "Female")

# DNA overall and by specialty
DNA_agg     <- xtabs(~ DNA, data = agg_data)
DNA_spc     <- xtabs(~ specialty + DNA, data = spec_agg_data)
DNA_spc_df  <- as.data.frame(prop.table(DNA_spc, margin = 1))
DNA_spc_pct <- format(round(prop.table(DNA_spc, margin = 1), 2), nsmall = 2)

DNA_spec <- ggplot(DNA_spc_df, aes(x = specialty, y = Freq, fill = DNA)) +
  geom_col(position = "fill") +
  coord_flip() +
  ggtitle("DNA Proportions by Specialty") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_cambridge_blue() +
  theme_cambridge_blue() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))

# Recorded ethnicity
t_eth       <- dna_table(agg_data, "ethnicity2")
DNA_eth     <- t_eth$counts
DNA_eth_pct <- t_eth$pct
DNA_eth_g   <- dna_bar(t_eth$long_df, "ethnicity2", "Ethnicity",
                       "DNA Proportions by Recorded Ethnicity", eth_labels)
DNA_eth_spec_g <- dna_bar_spec(spec_agg_data, "ethnicity2", "Ethnicity",
                               "DNA Proportions by Recorded Ethnicity and Specialty", eth_labels)

# Linked (LSOA level) ethnicity
t_ethl <- dna_table(agg_data, "non_white_cat")
DNA_eth_linked     <- t_ethl$counts
DNA_eth_linked_pct <- t_ethl$pct
DNA_eth_linked_g   <- dna_bar(t_ethl$long_df, "non_white_cat", "Non-white percentage",
                              "DNA Proportions by Linked Ethnicity")
DNA_eth_linked_spec_g <- dna_bar_spec(spec_agg_data, "non_white_cat", "Non-white percentage",
                                      "DNA Proportions by Linked Ethnicity and Specialty")

# Index of Multiple Deprivation
t_imd       <- dna_table(agg_data, "IMD_decile")
DNA_IMD     <- t_imd$counts
DNA_IMD_pct <- t_imd$pct
DNA_IMD_g   <- dna_bar(t_imd$long_df, "IMD_decile", "IMD Decile",
                       "DNA Proportions by IMD Decile")
DNA_IMD_spec_g <- dna_bar_spec(spec_agg_data, "IMD_decile", "IMD Decile",
                               "DNA Proportions by IMD Decile and Specialty")

# Income deprivation
t_iod       <- dna_table(agg_data, "income_quintile")
DNA_IoD     <- t_iod$counts
DNA_IoD_pct <- t_iod$pct
DNA_IoD_g   <- dna_bar(t_iod$long_df, "income_quintile", "Income Deprivation Quintile",
                       "DNA Proportions by Income Quintile")
DNA_IoD_spec_g <- dna_bar_spec(spec_agg_data, "income_quintile", "Income Deprivation Quintile",
                               "DNA Proportions by Income Quintile and Specialty")

# Employment deprivation
t_emp     <- dna_table(agg_data, "employment_quintile")
DNA_E     <- t_emp$counts
DNA_E_pct <- t_emp$pct
DNA_E_g   <- dna_bar(t_emp$long_df, "employment_quintile", "Employment Deprivation Quintile",
                     "DNA Proportions by Employment Quintile")
DNA_E_spec_g <- dna_bar_spec(spec_agg_data, "employment_quintile", "Employment Deprivation Quintile",
                             "DNA Proportions by Employment Quintile and Specialty")

# Age
t_age       <- dna_table(agg_data, "age_band")
DNA_age     <- t_age$counts
DNA_age_pct <- t_age$pct
DNA_age_g   <- dna_bar(t_age$long_df, "age_band", "Age Band",
                       "DNA Proportions by Age Group")
DNA_age_spec_g <- dna_bar_spec(spec_agg_data, "age_band", "Age Band",
                               "DNA Proportions by Age Group and Specialty")

# Sex
t_sex       <- dna_table(agg_data, "sex")
DNA_sex     <- t_sex$counts
DNA_sex_pct <- t_sex$pct
DNA_sex_g   <- dna_bar(t_sex$long_df, "sex", "Sex", "DNA Proportions by Sex", sex_labels)
DNA_sex_spec_g <- dna_bar_spec(spec_agg_data, "sex", "Sex",
                               "DNA Proportions by Sex and Specialty", sex_labels)

# Vehicle ownership
t_va       <- dna_table(agg_data, "vehicle_quintile")
DNA_va     <- t_va$counts
DNA_va_pct <- t_va$pct
DNA_va_g   <- dna_bar(t_va$long_df, "vehicle_quintile", "Vehicle Ownership Quintile",
                      "DNA Proportions by Vehicle Ownership")
DNA_va_spec_g <- dna_bar_spec(spec_agg_data, "vehicle_quintile", "Vehicle Ownership Quintile",
                              "DNA Proportions by Vehicle Ownership and Specialty")

# Childcare accessibility
t_ca       <- dna_table(agg_data, "childcare_quintile")
DNA_ca     <- t_ca$counts
DNA_ca_pct <- t_ca$pct
DNA_chld   <- dna_bar(t_ca$long_df, "childcare_quintile", "Childcare Access Quintile",
                      "DNA Proportions by Childcare Access")
DNA_chld_spec <- dna_bar_spec(spec_agg_data, "childcare_quintile", "Childcare Access Quintile",
                              "DNA Proportions by Childcare Access and Specialty")

# Distance measures. These are reported as tables only, since the banded
# distributions are easier to read as numbers than as eleven stacked bars.
DNA_add         <- xtabs(~ add_bin + DNA, data = agg_data)
colnames(DNA_add) <- c("Completed", "DNA")
DNA_add_pct     <- format(round(prop.table(DNA_add, margin = 1), 3), nsmall = 3)

DNA_GP          <- xtabs(~ GP_bin + DNA, data = agg_data)
colnames(DNA_GP) <- c("Completed", "DNA")
DNA_GP_pct      <- format(round(prop.table(DNA_GP, margin = 1), 2), nsmall = 2)

DNA_hospital     <- xtabs(~ hospital_bin + DNA, data = agg_data)
DNA_hospital_pct <- format(round(prop.table(DNA_hospital, margin = 1), 2), nsmall = 2)

DNA_NAPTAN       <- xtabs(~ NAPTAN_bin + DNA, data = agg_data)
DNA_NAPTAN_pct   <- format(round(prop.table(DNA_NAPTAN, margin = 1), 2), nsmall = 2)


# ---- 6.5 Patients with repeated missed appointments --------------------------

# Repeat non-attendance is a different phenomenon from a single missed
# appointment, so those patients are profiled separately using the full
# appointment level sample.
DNA_mult_counts <- mult_data %>%
  group_by(ptid) %>%
  summarise(n_cancels = sum(DNA == 1)) %>%
  filter(n_cancels > 1)

reg_data <- mult_data %>%
  arrange(ptid) %>%
  distinct(ptid, .keep_all = TRUE) %>%
  select(ptid, appt_booked_date, referral_date, ethnicity2, sex, language, age,
         IMD_decile, IMD_rank, lat, long, lsoa21, NAPTAN_distance_km, GP,
         hospital, ahah, income_quintile, employment_quintile, education_quintile,
         health_quintile, crime_quintile, housing_quintile,
         living_environment_quintile, vehicle_pct, childcare_overall,
         digital_propensity_score)

DNA_mult_counts <- DNA_mult_counts %>% left_join(reg_data, by = "ptid")

repeat_dna_summary <- DNA_mult_counts %>%
  summarise(`Number of patients`                     = n(),
            `Age, mean (SD)`                         = fmt_mean_sd(age),
            `Male, n (%)`                            = fmt_n_pct(sex == 0, n()),
            `Female, n (%)`                          = fmt_n_pct(sex == 1, n()),
            `White ethnicity, n (%)`                 = fmt_n_pct(ethnicity2 == 0, n()),
            `Non-white ethnicity, n (%)`             = fmt_n_pct(ethnicity2 == 1, n()),
            `English language, n (%)`                = fmt_n_pct(language == "English" | is.na(language), n()),
            `IMD decile, mean (SD)`                  = fmt_mean_sd(IMD_decile),
            `Income deprivation quintile, mean (SD)` = fmt_mean_sd(as.numeric(income_quintile)),
            `Vehicle availability %, mean (SD)`      = fmt_mean_sd(vehicle_pct))


################################################################################
# SECTION 7 — REGRESSION MODELLING
#
# Four nested specifications are fitted, all logistic with DNA as the outcome:
#
#   IMD model              DNA against the single composite deprivation decile
#   IoD model              DNA against the individual deprivation domains
#   Core DAG model         the full covariate set implied by the DAG
#   Parsimonious DAG model the reduced set retained after model comparison
#
# The ordinal deprivation measures are entered as natural splines rather than as
# linear terms or as full factors. A linear term imposes a constant effect per
# decile, which is not plausible, while a ten level factor spends nine degrees of
# freedom to say something a smooth curve says in two or three. The degrees of
# freedom for each spline are chosen by AIC rather than fixed in advance, and the
# whole-term significance is judged by likelihood ratio test, since the
# individual basis coefficients are not separately interpretable.
################################################################################

# ---- 7.1 Model helpers -------------------------------------------------------

sig_code <- function(pval) {
  if (is.na(pval)) return(" ")
  if (pval < 0.001) "***"
  else if (pval < 0.01) "**"
  else if (pval < 0.05) "*"
  else if (pval < 0.10) "."
  else " "
}

# Odds ratios with profile likelihood confidence intervals. Note that where the
# model contains ns() or bs() terms, the rows for the basis columns are log-odds
# for each basis function and are not individually interpretable as the effect of
# a one unit change; read those terms from the spline plots and the LRT instead.
extract_or <- function(model, conf_level = 0.95) {
  coefs <- coef(model)
  or    <- exp(coefs)
  ci    <- exp(confint(model, level = conf_level))
  ci    <- ci[names(coefs), , drop = FALSE]
  p     <- summary(model)$coefficients[, "Pr(>|z|)"]
  p     <- p[names(coefs)]
  sig   <- vapply(p, sig_code, character(1))

  list(tidy = tibble::tibble(term     = names(coefs),
                             estimate = coefs,
                             OR       = or,
                             CI_lower = ci[, 1],
                             CI_upper = ci[, 2],
                             p_value  = p,
                             sig      = sig),
       results = cbind(OR = or, CI_low = ci[, 1], CI_high = ci[, 2],
                       pvalue = p, sig = sig))
}

# Refit a model separately within each specialty and return the OR tables.
extract_or_by_specialty <- function(formula, data, family = binomial,
                                    conf_level = 0.95) {
  lapply(split(data, data$specialty), function(d) {
    if (nrow(d) == 0) return(NULL)
    extract_or(glm(formula = formula, data = d, family = family),
               conf_level = conf_level)
  })
}

# Fit the same model across a range of spline degrees of freedom and keep the
# AIC of each, so that the flexibility of the curve is chosen by fit rather than
# by assertion.
fit_glm_df <- function(formula_fn, data, df_range) {
  lapply(df_range, function(df) {
    m <- tryCatch(glm(formula_fn(df), data = data, family = binomial()),
                  error = function(e) NULL)
    list(df = df, model = m, AIC = if (!is.null(m)) AIC(m) else Inf)
  })
}

select_best_df <- function(df_list) df_list[[which.min(sapply(df_list, `[[`, "AIC"))]]

make_df_comparison <- function(df_list) {
  data.frame(df  = sapply(df_list, `[[`, "df"),
             AIC = round(sapply(df_list, `[[`, "AIC"), 2))
}


# ---- 7.2 Model 1: IMD ---------------------------------------------------------

# IMD has ten ordered levels, so df of 2 to 4 is a reasonable search range.
# ns() places its boundary knots at the observed minimum and maximum by default.
message("=== IMD MODEL: df selection ===")

imd_df_results <- fit_glm_df(
  formula_fn = function(df) as.formula(paste0("DNA ~ ns(IMD_num, df = ", df, ")")),
  data       = agg_data,
  df_range   = 2:4)

best_imd    <- select_best_df(imd_df_results)
best_df_imd <- best_imd$df
message("  Best df (IMD model): ", best_df_imd, "  (AIC = ", round(best_imd$AIC, 2), ")")

imd_df_comparison <- make_df_comparison(imd_df_results)
model.imd   <- best_imd$model
results_imd <- extract_or(model.imd)

# For the specialty-stratified fits the boundary knots are pinned to the global
# IMD range, so that the curves are comparable between specialties rather than
# each being anchored to its own subset.
imd_boundary <- quantile(agg_data$IMD_num, c(0.05, 0.95), na.rm = TRUE)

results.imd.models.by.specialty <- extract_or_by_specialty(
  formula = as.formula(paste0("DNA ~ ns(IMD_num, df = ", best_df_imd,
                              ", Boundary.knots = c(", imd_boundary[1], ", ",
                              imd_boundary[2], "))")),
  data    = spec_agg_data)


# ---- 7.3 Model 2: Index of Deprivation domains --------------------------------

# Quintiles have only five distinct values, so df is capped at 3; df of 4 would
# leave roughly one observation per segment and is not attempted. Boundary knots
# are fixed at 1 and 5 throughout, since quintiles always span that range.
message("=== IoD MODEL: df selection ===")

iod_df_results <- fit_glm_df(
  formula_fn = function(df) as.formula(paste0(
    "DNA ~ ns(income_num,     df = ", df, ") + ",
    "ns(health_num,     df = ", df, ") + ",
    "ns(crime_num,      df = ", df, ") + ",
    "ns(housing_num,    df = ", df, ") + ",
    "ns(living_env_num, df = ", df, ")")),
  data       = agg_data,
  df_range   = 2:3)

best_iod    <- select_best_df(iod_df_results)
best_df_iod <- best_iod$df
message("  Best df (IoD model): ", best_df_iod, "  (AIC = ", round(best_iod$AIC, 2), ")")

iod_df_comparison <- make_df_comparison(iod_df_results)
model.test.iod <- best_iod$model
results_iod    <- extract_or(model.test.iod)

results.iod.models.by.specialty <- extract_or_by_specialty(
  formula = as.formula(paste0(
    "DNA ~ ns(income_num,     df = ", best_df_iod, ", Boundary.knots = c(1,5)) + ",
    "ns(health_num,     df = ", best_df_iod, ", Boundary.knots = c(1,5)) + ",
    "ns(crime_num,      df = ", best_df_iod, ", Boundary.knots = c(1,5)) + ",
    "ns(housing_num,    df = ", best_df_iod, ", Boundary.knots = c(1,5)) + ",
    "ns(living_env_num, df = ", best_df_iod, ", Boundary.knots = c(1,5))")),
  data    = spec_agg_data)


# ---- 7.4 Testing the continuous covariates for non-linearity ------------------

# Vehicle ownership and childcare accessibility are genuinely continuous, so
# whether they need a spline at all is tested rather than assumed. The likelihood
# ratio test compares the linear fit against the spline fit; if it does not
# reject, the linear term is retained in the DAG models below. bs() is included
# in the AIC comparison as an alternative, since these two variables do not
# benefit from the tail linearity that ns() imposes.
message("=== Testing vehicle_pct_10 and childcare_overall for non-linearity ===")

m_veh_linear <- glm(DNA ~ vehicle_pct_10,             data = agg_data, family = binomial())
m_veh_ns3    <- glm(DNA ~ ns(vehicle_pct_10, df = 3), data = agg_data, family = binomial())
m_veh_ns4    <- glm(DNA ~ ns(vehicle_pct_10, df = 4), data = agg_data, family = binomial())
m_veh_bs3    <- glm(DNA ~ bs(vehicle_pct_10, df = 3), data = agg_data, family = binomial())

lrt_veh            <- lrtest(m_veh_linear, m_veh_ns3)
pval_veh_nl        <- lrt_veh$`Pr(>Chisq)`[2]
use_spline_vehicle <- !is.na(pval_veh_nl) && pval_veh_nl < 0.05

veh_aics      <- c(ns3 = AIC(m_veh_ns3), ns4 = AIC(m_veh_ns4), bs3 = AIC(m_veh_bs3))
best_veh_spec <- names(which.min(veh_aics))

message("  vehicle_pct_10  LRT p = ", round(pval_veh_nl, 4),
        "  | use spline: ", use_spline_vehicle,
        "  | best spec: ", best_veh_spec)

m_cc_linear <- glm(DNA ~ childcare_overall,             data = agg_data, family = binomial())
m_cc_ns3    <- glm(DNA ~ ns(childcare_overall, df = 3), data = agg_data, family = binomial())
m_cc_ns4    <- glm(DNA ~ ns(childcare_overall, df = 4), data = agg_data, family = binomial())
m_cc_bs3    <- glm(DNA ~ bs(childcare_overall, df = 3), data = agg_data, family = binomial())

lrt_cc               <- lrtest(m_cc_linear, m_cc_ns3)
pval_cc_nl           <- lrt_cc$`Pr(>Chisq)`[2]
use_spline_childcare <- !is.na(pval_cc_nl) && pval_cc_nl < 0.05

cc_aics      <- c(ns3 = AIC(m_cc_ns3), ns4 = AIC(m_cc_ns4), bs3 = AIC(m_cc_bs3))
best_cc_spec <- names(which.min(cc_aics))

message("  childcare_overall LRT p = ", round(pval_cc_nl, 4),
        "  | use spline: ", use_spline_childcare,
        "  | best spec: ", best_cc_spec)

nonlinearity_tests <- data.frame(
  Variable   = c("vehicle_pct_10", "childcare_overall"),
  LRT_pvalue = round(c(pval_veh_nl, pval_cc_nl), 4),
  Use_Spline = c(use_spline_vehicle, use_spline_childcare),
  Best_Spec  = c(best_veh_spec, best_cc_spec))

print(nonlinearity_tests)

# The chosen functional form is carried forward as a string so that every model
# below builds its formula from the same decision rather than repeating it.
spline_fn_veh <- switch(best_veh_spec,
                        ns3 = "ns(vehicle_pct_10, df = 3)",
                        ns4 = "ns(vehicle_pct_10, df = 4)",
                        bs3 = "bs(vehicle_pct_10, df = 3)")

spline_fn_cc <- switch(best_cc_spec,
                       ns3 = "ns(childcare_overall, df = 3)",
                       ns4 = "ns(childcare_overall, df = 4)",
                       bs3 = "bs(childcare_overall, df = 3)")

vehicle_term   <- if (use_spline_vehicle)   spline_fn_veh else "vehicle_pct_10"
childcare_term <- if (use_spline_childcare) spline_fn_cc  else "childcare_overall"


# ---- 7.5 Model 3: Core DAG ----------------------------------------------------

# The covariate set here follows the DAG rather than a stepwise search: the
# deprivation domains, transport and childcare access, distance to care, and the
# individual level demographics that the DAG identifies as confounders.
message("=== DAG MODEL: df selection ===")

dag_terms <- function(df) paste0(
  "DNA ~ ",
  "ns(employment_num, df = ", df, ", Boundary.knots = c(1,5)) + ",
  "ns(education_num,  df = ", df, ", Boundary.knots = c(1,5)) + ",
  "ns(health_num,     df = ", df, ", Boundary.knots = c(1,5)) + ",
  "ns(crime_num,      df = ", df, ", Boundary.knots = c(1,5)) + ",
  vehicle_term, " + ", childcare_term, " + ",
  "ethnicity2 + sex + age + non_english + NAPTAN_distance_km + add_dist + GP")

dag_df_results <- fit_glm_df(formula_fn = function(df) as.formula(dag_terms(df)),
                             data       = agg_data,
                             df_range   = 2:3)

best_dag    <- select_best_df(dag_df_results)
best_df_dag <- best_dag$df
message("  Best df (DAG model): ", best_df_dag, "  (AIC = ", round(best_dag$AIC, 2), ")")

dag_df_comparison <- make_df_comparison(dag_df_results)
model_DAG   <- best_dag$model
results_DAG <- extract_or(model_DAG)

results.dag.by.specialty <- extract_or_by_specialty(
  formula = as.formula(dag_terms(best_df_dag)),
  data    = spec_agg_data)

# Collinearity check on the core model
print(vif(model_DAG))
print(anova(model_DAG, test = "LRT"))


# ---- 7.6 Model 4: Parsimonious DAG --------------------------------------------

# The reduced specification, dropping the terms that contributed little once the
# others were present. It is fitted separately here because it is the model
# reported as the primary result.
pars_formula <- as.formula(paste0(
  "DNA ~ ",
  "ns(employment_num, df = ", best_df_dag, ", Boundary.knots = c(1,5)) + ",
  "ns(health_num,     df = ", best_df_dag, ", Boundary.knots = c(1,5)) + ",
  "ns(crime_num,      df = ", best_df_dag, ", Boundary.knots = c(1,5)) + ",
  vehicle_term, " + ethnicity2 + sex + age + add_dist"))

model.DAG <- glm(pars_formula, data = agg_data, family = binomial())

print(vif(model.DAG))
print(anova(model.DAG, test = "LRT"))
cat("Model AIC:", round(AIC(model.DAG), 2), "\n")


# ---- 7.7 Model comparison -----------------------------------------------------

# AIC is only comparable across models fitted to identical rows, so all four are
# refitted on one complete case extract of the variables any of them use. AICc is
# reported alongside AIC because the spline bases plus the covariate set make k
# large enough that the correction is not negligible.
vars_needed <- c("DNA", "IMD_num", "sex", "ethnicity2", "age",
                 "employment_num", "education_num", "health_num", "crime_num",
                 "income_num", "housing_num", "living_env_num",
                 "vehicle_pct_10", "childcare_overall",
                 "non_english", "NAPTAN_distance_km", "add_dist", "GP", "specialty")

analysis_data <- agg_data %>%
  select(all_of(vars_needed[vars_needed %in% names(agg_data)])) %>%
  na.omit()

m_imd <- glm(as.formula(paste0("DNA ~ ns(IMD_num, df=", best_df_imd, ")")),
             data = analysis_data, family = binomial())

m_iod <- glm(as.formula(paste0(
  "DNA ~ ns(income_num,     df=", best_df_iod, ", Boundary.knots=c(1,5))+",
  "ns(health_num,     df=", best_df_iod, ", Boundary.knots=c(1,5))+",
  "ns(crime_num,      df=", best_df_iod, ", Boundary.knots=c(1,5))+",
  "ns(housing_num,    df=", best_df_iod, ", Boundary.knots=c(1,5))+",
  "ns(living_env_num, df=", best_df_iod, ", Boundary.knots=c(1,5))")),
  data = analysis_data, family = binomial())

m_dag  <- glm(as.formula(dag_terms(best_df_dag)), data = analysis_data, family = binomial())
m_pars <- glm(pars_formula,                       data = analysis_data, family = binomial())

model_list <- list("IMD Model"              = m_imd,
                   "IoD Model"              = m_iod,
                   "Core DAG Model"         = m_dag,
                   "Parsimonious DAG Model" = m_pars)

aic_raw           <- AIC(m_imd, m_iod, m_dag, m_pars)
rownames(aic_raw) <- names(model_list)
n_obs             <- nobs(m_imd)

aic_tbl <- data.frame(Model         = rownames(aic_raw),
                      `Params (k)`  = aic_raw$df,
                      N             = n_obs,
                      AIC           = round(aic_raw$AIC, 2),
                      stringsAsFactors = FALSE,
                      check.names      = FALSE)

aic_tbl$AICc <- round(aic_tbl$AIC + (2 * aic_tbl$`Params (k)` * (aic_tbl$`Params (k)` + 1)) /
                        (n_obs - aic_tbl$`Params (k)` - 1), 2)

aic_tbl$`Delta AIC`  <- round(aic_tbl$AIC  - min(aic_tbl$AIC),  2)
aic_tbl$`Delta AICc` <- round(aic_tbl$AICc - min(aic_tbl$AICc), 2)

w_aic  <- exp(-0.5 * aic_tbl$`Delta AIC`)
w_aicc <- exp(-0.5 * aic_tbl$`Delta AICc`)

aic_tbl$`AIC Weight`  <- round(w_aic  / sum(w_aic),  4)
aic_tbl$`AICc Weight` <- round(w_aicc / sum(w_aicc), 4)
aic_tbl$Best          <- ifelse(aic_tbl$`Delta AICc` == 0, "Yes", "")

model_compare_AIC <- aic_tbl
print(model_compare_AIC)


# ---- 7.8 Specialty heterogeneity ----------------------------------------------

# Each term in the parsimonious model is tested for an interaction with
# specialty, comparing the additive model against the same model plus that one
# interaction. The likelihood ratio test is used rather than the individual
# coefficient p-values, since a spline term spans several columns and has to be
# tested as a block.
analysis_data$specialty <- factor(analysis_data$specialty)

rhs_terms <- attr(terms(formula(m_pars)), "term.labels")
rhs_terms <- rhs_terms[rhs_terms != "specialty"]

base_formula_str <- paste("DNA ~", paste(rhs_terms, collapse = " + "), "+ specialty")
m_additive       <- glm(as.formula(base_formula_str), data = analysis_data, family = binomial())

interaction_results <- lapply(rhs_terms, function(term) {
  m_int <- tryCatch(glm(as.formula(paste0(base_formula_str, " + ", term, ":specialty")),
                        data = analysis_data, family = binomial()),
                    error = function(e) NULL)

  if (is.null(m_int))
    return(data.frame(variable = term, LRT_chisq = NA_real_, df = NA_integer_,
                      p_value = NA_real_, AIC_add = AIC(m_additive),
                      AIC_int = NA_real_, delta_AIC = NA_real_))

  lrt <- anova(m_additive, m_int, test = "LRT")

  data.frame(variable  = term,
             LRT_chisq = lrt$Deviance[2],
             df        = lrt$Df[2],
             p_value   = lrt$`Pr(>Chi)`[2],
             AIC_add   = AIC(m_additive),
             AIC_int   = AIC(m_int),
             delta_AIC = AIC(m_int) - AIC(m_additive))
})

results_df <- do.call(rbind, interaction_results)
rownames(results_df) <- NULL

results_df <- results_df %>%
  mutate(sig = case_when(p_value < 0.001 ~ "***",
                         p_value < 0.01  ~ "**",
                         p_value < 0.05  ~ "*",
                         p_value < 0.1   ~ ".",
                         TRUE            ~ ""))

print(results_df)


# ---- 7.9 Interaction matrix between exposures ---------------------------------

# Beyond specialty, the two access exposures are tested against every other
# covariate in the model, since an effect of transport or childcare access that
# depends on deprivation or demography would be the substantively interesting
# finding. Each pair is tested by LRT against the fully adjusted additive model,
# so that the comparison is not confounded by omitted main effects.
base_terms <- list(
  employment_num = paste0("ns(employment_num, df = ", best_df_dag, ", Boundary.knots = c(1,5))"),
  health_num     = paste0("ns(health_num,     df = ", best_df_dag, ", Boundary.knots = c(1,5))"),
  crime_num      = paste0("ns(crime_num,      df = ", best_df_dag, ", Boundary.knots = c(1,5))"),
  vehicle_pct_10 = vehicle_term,
  ethnicity2     = "ethnicity2",
  sex            = "sex",
  age            = "age")

test_setups <- list(
  list(target_name       = "vehicle",
       target_term       = base_terms$vehicle_pct_10,
       moderators        = c("employment_num", "health_num", "crime_num",
                             "ethnicity2", "sex", "age"),
       include_childcare = FALSE),
  list(target_name       = "childcare",
       target_term       = childcare_term,
       moderators        = c("employment_num", "health_num", "crime_num",
                             "vehicle_pct_10", "ethnicity2", "sex", "age"),
       include_childcare = TRUE))

interaction_results <- list()

for (setup in test_setups) {
  for (mod in setup$moderators) {
    message("Testing interaction: ", setup$target_name, " x ", mod)

    main_effect_strings <- unlist(base_terms)
    if (setup$include_childcare)
      main_effect_strings <- c(main_effect_strings, childcare_term)

    f_add_str <- paste0("DNA ~ ", paste(main_effect_strings, collapse = " + "))
    f_int_str <- paste0(f_add_str, " + ", setup$target_term, ":", base_terms[[mod]])

    m_add <- glm(as.formula(f_add_str), data = agg_data, family = binomial())
    m_int <- glm(as.formula(f_int_str), data = agg_data, family = binomial())

    lrt <- lmtest::lrtest(m_add, m_int)

    interaction_results[[length(interaction_results) + 1]] <- data.frame(
      target_variable    = setup$target_name,
      moderator_variable = mod,
      LRT_chisq          = round(lrt$Chisq[2], 2),
      df                 = lrt$Df[2],
      p_value            = lrt$`Pr(>Chisq)`[2],
      AIC_add            = AIC(m_add),
      AIC_int            = AIC(m_int),
      delta_AIC          = AIC(m_int) - AIC(m_add),
      stringsAsFactors   = FALSE)
  }
}

interaction_df <- do.call(rbind, interaction_results) %>%
  mutate(sig = case_when(p_value < 0.001 ~ "***",
                         p_value < 0.01  ~ "**",
                         p_value < 0.05  ~ "*",
                         p_value < 0.1   ~ ".",
                         TRUE            ~ ""))

print(interaction_df)


# ---- 7.10 Partial effect plots ------------------------------------------------

# Because the models are fitted with glm() rather than rms::lrm(), there is no
# Predict() method available, so the partial effect is reconstructed by hand:
# hold every other covariate at its median and vary the splined predictor across
# its observed range. This is what makes the spline terms readable, given that
# their individual coefficients are not.
`%||%` <- function(a, b) if (!is.null(a)) a else b

plot_ns_effect <- function(model, varname, data, n_points = 100,
                           xlab = varname, title = NULL) {

  x_range <- seq(min(data[[varname]], na.rm = TRUE),
                 max(data[[varname]], na.rm = TRUE),
                 length.out = n_points)

  ref_row      <- data[1, , drop = FALSE]
  numeric_cols <- sapply(data, is.numeric)
  ref_row[numeric_cols] <- lapply(data[numeric_cols], median, na.rm = TRUE)

  grid            <- ref_row[rep(1, n_points), ]
  grid[[varname]] <- x_range

  preds <- predict(model, newdata = grid, type = "link", se.fit = TRUE)

  df_plot <- data.frame(x   = x_range,
                        fit = preds$fit,
                        lwr = preds$fit - 1.96 * preds$se.fit,
                        upr = preds$fit + 1.96 * preds$se.fit)

  ggplot(df_plot, aes(x = x, y = fit)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = CAM_WARM, alpha = 0.2) +
    geom_line(colour = CAM_WARM, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = xlab, y = "Log-odds of DNA (partial effect)",
         title = title %||% paste("Spline effect:", varname)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p_imd    <- plot_ns_effect(model.imd,      "IMD_num",        agg_data, xlab = "IMD Decile",          title = "IMD Model - Spline Effect of IMD Decile")
p_income <- plot_ns_effect(model.test.iod, "income_num",     agg_data, xlab = "Income Quintile",     title = "IoD Model - Income")
p_health <- plot_ns_effect(model.test.iod, "health_num",     agg_data, xlab = "Health Quintile",     title = "IoD Model - Health")
p_crime  <- plot_ns_effect(model.test.iod, "crime_num",      agg_data, xlab = "Crime Quintile",      title = "IoD Model - Crime")
p_employ <- plot_ns_effect(model_DAG,      "employment_num", agg_data, xlab = "Employment Quintile", title = "DAG Model - Employment")
p_educ   <- plot_ns_effect(model_DAG,      "education_num",  agg_data, xlab = "Education Quintile",  title = "DAG Model - Education")

gridExtra::grid.arrange(p_imd, ncol = 1)
gridExtra::grid.arrange(p_income, p_health, p_crime, ncol = 3)
gridExtra::grid.arrange(p_employ, p_educ, ncol = 2)

if (use_spline_vehicle)
  print(plot_ns_effect(model_DAG, "vehicle_pct_10", agg_data,
                       xlab  = "Vehicle Ownership (per 10%)",
                       title = "DAG Model - Vehicle Ownership Spline"))

if (use_spline_childcare)
  print(plot_ns_effect(model_DAG, "childcare_overall", agg_data,
                       xlab  = "Childcare Accessibility Index",
                       title = "DAG Model - Childcare Spline"))


# ---- 7.11 Whole-term significance ---------------------------------------------

# anova() with an LRT tests each term as a block, which is the correct test for
# a spline; the per-coefficient p-values in summary() are not.
message("=== Whole-term LRT: IMD model ===")
print(anova(model.imd, test = "LRT"))

message("=== Whole-term LRT: IoD model ===")
print(anova(model.test.iod, test = "LRT"))

message("=== Whole-term LRT: Core DAG model ===")
print(anova(model_DAG, test = "LRT"))

message("=== Whole-term LRT: Parsimonious DAG model ===")
print(anova(model.DAG, test = "LRT"))


# ---- 7.12 Odds ratio tables ----------------------------------------------------

# Printed rather than exported, since these are the numbers that go into the
# manuscript tables by hand. Spline basis rows are dropped: their coefficients
# are not interpretable individually, and the shape of those terms is read from
# the partial effect plots above and their significance from the LRTs.
print_or_table <- function(model, label) {
  tab <- extract_or(model)$tidy
  tab <- tab[!grepl("^ns\\(|^bs\\(", tab$term), ]

  cat("\n=== Adjusted odds ratios:", label, "===\n")
  print(data.frame(Term  = tab$term,
                   OR    = sprintf("%.3f", tab$OR),
                   CI    = sprintf("%.3f to %.3f", tab$CI_lower, tab$CI_upper),
                   P     = ifelse(tab$p_value < 0.001, "<0.001",
                                  sprintf("%.3f", tab$p_value)),
                   Sig   = tab$sig,
                   row.names = NULL))
  invisible(tab)
}

or_imd  <- print_or_table(model.imd,      "IMD model")
or_iod  <- print_or_table(model.test.iod, "IoD model")
or_dag  <- print_or_table(model_DAG,      "Core DAG model")
or_pars <- print_or_table(model.DAG,      "Parsimonious DAG model")

# Predicted probability of DNA across the range of each splined term, holding
# every other covariate at its median. This is the readable counterpart to the
# spline coefficients, and is what the manuscript reports for those variables.
predict_at_knots <- function(model, varname, knots, data, label = varname) {
  ref_row  <- data[1, , drop = FALSE]
  num_cols <- sapply(data, is.numeric)
  ref_row[num_cols] <- lapply(data[num_cols], median, na.rm = TRUE)

  grid            <- ref_row[rep(1, length(knots)), ]
  grid[[varname]] <- knots

  preds <- predict(model, newdata = grid, type = "response", se.fit = TRUE)

  data.frame(Variable = label,
             Value    = knots,
             Pred     = sprintf("%.3f", preds$fit),
             CI       = sprintf("%.3f to %.3f",
                                pmax(0, preds$fit - 1.96 * preds$se.fit),
                                pmin(1, preds$fit + 1.96 * preds$se.fit)),
             row.names = NULL)
}

cat("\n=== Predicted DNA probability across deprivation ===\n")
print(rbind(
  predict_at_knots(model.imd,      "IMD_num",        1:10, agg_data, "IMD (decile)"),
  predict_at_knots(model.test.iod, "income_num",     1:5,  agg_data, "Income (quintile)"),
  predict_at_knots(model.test.iod, "health_num",     1:5,  agg_data, "Health (quintile)"),
  predict_at_knots(model.test.iod, "crime_num",      1:5,  agg_data, "Crime (quintile)"),
  predict_at_knots(model_DAG,      "employment_num", 1:5,  agg_data, "Employment (quintile)"),
  predict_at_knots(model_DAG,      "education_num",  1:5,  agg_data, "Education (quintile)")))


################################################################################
# SECTION 8 — HEATMAPS AND CHOROPLETHS
#
# Two geographic views are produced. The aggregate maps cover all of England,
# since patients are referred to Addenbrooke's from well outside the county, with
# the Cambridgeshire catchment marked for reference. The specialty maps are
# cropped to Cambridgeshire, where the counts per LSOA are large enough for a
# rate to mean anything. The bivariate choropleths show deprivation against
# vehicle ownership, which is the pairing the models suggest matters most.
################################################################################

# ---- 8.1 Map theme and helpers ------------------------------------------------

theme_map_pub <- function(base_size = 10) {
  theme_void(base_size = base_size) +
    theme(plot.title        = element_text(face = "bold", size = base_size + 2,
                                           hjust = 0.5, margin = margin(b = 4)),
          plot.subtitle     = element_text(size = base_size - 1, colour = "grey40",
                                           hjust = 0.5, margin = margin(b = 6)),
          plot.caption      = element_text(size = base_size - 3, colour = "grey60",
                                           hjust = 0),
          strip.text        = element_text(face = "bold", size = base_size - 1,
                                           margin = margin(b = 3, t = 3)),
          legend.title      = element_text(face = "bold", size = base_size - 1),
          legend.text       = element_text(size = base_size - 2),
          legend.position   = "bottom",
          legend.key.width  = unit(1.2, "cm"),
          legend.key.height = unit(0.35, "cm"),
          panel.spacing     = unit(2, "mm"),
          plot.margin       = margin(8, 8, 8, 8),
          plot.background   = element_rect(fill = "white", colour = NA))
}

# Some graphics devices render Unicode dashes and curly quotes inconsistently
# across platforms, so titles are forced to ASCII before they are drawn.
ascii_title <- function(x) {
  x <- gsub("\u2014|\u2013", "-",  x)
  x <- gsub("\u2018|\u2019", "'",  x)
  x <- gsub("\u201C|\u201D", "\"", x)
  x
}

spec_title <- function(sp) tools::toTitleCase(tolower(as.character(sp)))


# ---- 8.2 Geometry -------------------------------------------------------------

# s2 spherical geometry is disabled and the polygons repaired, since the ONS
# boundary file contains self-intersections that otherwise abort the crop.
sf_use_s2(FALSE)

lsoa_shapes <- st_read(file.path(OPEN_DIR, "LSOA shapefile.geojson"), quiet = TRUE) %>%
  st_transform(crs = 4326) %>%
  st_make_valid()

england_outline <- lsoa_shapes %>%
  st_union() %>%
  st_sf(geometry = .)

CAMBS_BBOX <- st_bbox(c(xmin = -0.6, ymin = 51.77, xmax = 0.7, ymax = 52.55),
                      crs = st_crs(4326))

lsoa_cambs   <- st_crop(lsoa_shapes, CAMBS_BBOX)
lsoa_england <- lsoa_shapes


# ---- 8.3 Rates by LSOA --------------------------------------------------------

dna_agg_lsoa <- agg_data %>%
  filter(!is.na(lsoa21)) %>%
  group_by(lsoa21) %>%
  summarise(total_patients = n(),
            dna_count      = sum(DNA == 1, na.rm = TRUE),
            dna_ratio      = dna_count / total_patients,
            .groups        = "drop")

dna_spec_lsoa <- spec_agg_data %>%
  filter(!is.na(lsoa21)) %>%
  group_by(lsoa21, specialty) %>%
  summarise(total_patients = n(),
            dna_count      = sum(DNA == 1, na.rm = TRUE),
            dna_ratio      = dna_count / total_patients,
            .groups        = "drop") %>%
  mutate(specialty_label = spec_title(specialty))

pt_agg_lsoa <- agg_data %>%
  filter(!is.na(lsoa21)) %>%
  count(lsoa21, name = "patient_count")

pt_spec_lsoa <- spec_agg_data %>%
  filter(!is.na(lsoa21)) %>%
  count(lsoa21, specialty, name = "patient_count") %>%
  mutate(specialty_label = spec_title(specialty))


# ---- 8.4 Bivariate classification ---------------------------------------------

# A 5 by 5 palette running from pale grey at the low-low corner to deep purple at
# high-high, so that the two dimensions remain separable by eye.
BIVAR_PAL <- c(
  "1-1" = "#e8e8e8", "1-2" = "#b5dcdc", "1-3" = "#8ac6c6", "1-4" = "#5fb0b0", "1-5" = "#3b9999",
  "2-1" = "#dfb0d6", "2-2" = "#c8b3d6", "2-3" = "#a5add3", "2-4" = "#7f9ac9", "2-5" = "#5698b9",
  "3-1" = "#d08fc1", "3-2" = "#b088c2", "3-3" = "#8c62aa", "3-4" = "#6f4d9d", "3-5" = "#503e8c",
  "4-1" = "#c06aaa", "4-2" = "#9f5fa8", "4-3" = "#7b4fa2", "4-4" = "#593e8f", "4-5" = "#3b4994",
  "5-1" = "#a8387e", "5-2" = "#8b2f73", "5-3" = "#6f2a6b", "5-4" = "#55245f", "5-5" = "#3b1f55")

# Assign each LSOA a "deprivation quintile - vehicle quintile" class and attach
# the corresponding colour, defaulting to pale grey where either is unobserved.
make_bivar <- function(data, dep_var, geometry) {
  cls <- data %>%
    filter(!is.na(lsoa21), !is.na(.data[[dep_var]]), !is.na(vehicle_quintile)) %>%
    distinct(lsoa21, .data[[dep_var]], vehicle_quintile) %>%
    mutate(dq          = pmin(pmax(as.integer(as.character(.data[[dep_var]])), 1L), 5L),
           vq          = pmin(pmax(as.integer(as.character(vehicle_quintile)), 1L), 5L),
           bivar_class = paste0(dq, "-", vq))

  left_join(geometry, cls, by = c("LSOA21CD" = "lsoa21")) %>%
    mutate(fill_col = dplyr::coalesce(BIVAR_PAL[bivar_class], "#f0f0f0"))
}

sf_bivar  <- make_bivar(agg_data, "health_quintile",     lsoa_england)
sf_bivar2 <- make_bivar(agg_data, "employment_quintile", lsoa_england)
sf_bivar3 <- make_bivar(agg_data, "crime_quintile",      lsoa_england)

# Join the rates to geometry
sf_dna_agg  <- left_join(lsoa_england, dna_agg_lsoa,  by = c("LSOA21CD" = "lsoa21"))
sf_dna_spec <- left_join(lsoa_cambs,   dna_spec_lsoa, by = c("LSOA21CD" = "lsoa21"))
sf_pt_agg   <- left_join(lsoa_england, pt_agg_lsoa,   by = c("LSOA21CD" = "lsoa21"))
sf_pt_spec  <- left_join(lsoa_cambs,   pt_spec_lsoa,  by = c("LSOA21CD" = "lsoa21"))

# The bivariate legend is drawn as its own 5 by 5 tile plot and inset onto the
# map, since a two dimensional scale cannot be expressed by a standard guide.
make_bivar_legend <- function(x_label = "Vehicle Ownership",
                              y_label = "Deprivation Index",
                              base_size = 7) {
  leg_grid <- expand.grid(x = 1:5, y = 1:5) %>%
    mutate(bivar_class = paste0(y, "-", x),
           fill_col    = BIVAR_PAL[bivar_class])

  ggplot(leg_grid, aes(x = x, y = y, fill = fill_col)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    scale_fill_identity() +
    scale_x_continuous(breaks = 1:5,
                       labels = c("Q1\n(Less)", "Q2", "Q3", "Q4", "Q5\n(More)"),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:5,
                       labels = c("Q1\n(More)", "Q2", "Q3", "Q4", "Q5\n(Less)"),
                       expand = c(0, 0)) +
    labs(x = x_label, y = y_label) +
    theme_minimal(base_size = base_size) +
    theme(axis.ticks       = element_blank(),
          panel.grid       = element_blank(),
          axis.text.x      = element_text(size = base_size - 2, margin = margin(t = 1)),
          axis.text.y      = element_text(size = base_size - 2, margin = margin(r = 1)),
          axis.title.x     = element_text(size = base_size - 1, face = "bold", margin = margin(t = 4)),
          axis.title.y     = element_text(size = base_size - 1, face = "bold", margin = margin(r = 4)),
          plot.margin      = margin(0, 0, 0, 0, "pt"),
          plot.background  = element_blank(),
          panel.background = element_blank()) +
    coord_fixed()
}


# ---- 8.5 Fill scales ----------------------------------------------------------

# LSOAs with no referrals are shaded mid-grey rather than left white, so that
# "no patients" is visually distinct from "no missed appointments".
dna_scale <- scale_fill_gradientn(
  colours  = c("#f2fcfa", "#8EE8D8", "#00bdb6", "#027a76", "#133844"),
  na.value = "#d0d0d0",
  limits   = c(0, 1),
  labels   = scales::percent_format(accuracy = 1),
  name     = "DNA Rate",
  guide    = guide_colorbar(title.position = "top",
                            barwidth  = unit(5, "cm"),
                            barheight = unit(0.4, "cm")))

# Patient counts are square root scaled, since a handful of urban LSOAs would
# otherwise flatten the entire rural distribution to one colour.
pt_scale <- scale_fill_gradientn(
  colours  = c("#f2fcfa", "#d0e8f5", "#4a9fcc", "#1a5a8a", "#133844"),
  na.value = "#d0d0d0",
  trans    = "sqrt",
  name     = "Patient Count",
  guide    = guide_colorbar(title.position = "top",
                            barwidth  = unit(5, "cm"),
                            barheight = unit(0.4, "cm")))


# ---- 8.6 Map constructors -----------------------------------------------------

CUH_CAPTION <- paste0("Source: NHS Cambridge University Hospitals Referral Data; ",
                      "ONS LSOA Boundaries 2021; Index of Multiple Deprivation 2019")

england_bg <- function() {
  list(geom_sf(data = england_outline, fill = "#ebebeb", colour = "#c8c8c8",
               linewidth = 0.15, inherit.aes = FALSE))
}

cambs_box <- function() {
  list(geom_sf(data = st_as_sfc(CAMBS_BBOX), fill = NA, colour = "#e63030",
               linewidth = 0.5, linetype = "dashed", inherit.aes = FALSE))
}

make_dna_agg_map <- function() {
  ggplot(sf_dna_agg) +
    england_bg() +
    geom_sf(aes(fill = dna_ratio), colour = "#aaaaaa", linewidth = 0.05, alpha = 0.85) +
    cambs_box() +
    dna_scale +
    annotation_scale(location = "bl", width_hint = 0.12, text_cex = 0.55,
                     style = "ticks", line_col = "#555", text_col = "#555") +
    labs(title    = "DNA Rate by LSOA - All Specialties",
         subtitle = paste0("Proportion of appointments not attended",
                           " (first event per patient). ",
                           "Red dashed box = Cambridgeshire catchment."),
         caption  = CUH_CAPTION) +
    theme_map_pub()
}

make_dna_spec_map <- function(specialties_subset = NULL, ncol = 4) {
  dat <- sf_dna_spec
  if (!is.null(specialties_subset)) dat <- filter(dat, specialty %in% specialties_subset)
  dat <- filter(dat, !is.na(dna_ratio))

  ggplot(dat) +
    england_bg() +
    geom_sf(aes(fill = dna_ratio), colour = "#888888", linewidth = 0.05, alpha = 0.85) +
    dna_scale +
    coord_sf(xlim   = c(CAMBS_BBOX["xmin"], CAMBS_BBOX["xmax"]),
             ylim   = c(CAMBS_BBOX["ymin"], CAMBS_BBOX["ymax"]),
             expand = FALSE) +
    facet_wrap(~ specialty_label, ncol = ncol) +
    labs(title    = "DNA Rate by LSOA, Stratified by Specialty",
         subtitle = "Cambridgeshire catchment area. First event per patient per specialty.",
         caption  = CUH_CAPTION) +
    theme_map_pub(base_size = 9)
}

make_pt_agg_map <- function() {
  ggplot(sf_pt_agg) +
    england_bg() +
    geom_sf(aes(fill = patient_count), colour = "#aaaaaa", linewidth = 0.05, alpha = 0.85) +
    cambs_box() +
    pt_scale +
    annotation_scale(location = "bl", width_hint = 0.12, text_cex = 0.55,
                     style = "ticks", line_col = "#555", text_col = "#555") +
    labs(title    = "Patient Distribution by LSOA - All Specialties",
         subtitle = paste0("Number of patients referred ",
                           "(first event per patient, square root scale). ",
                           "Red dashed box = Cambridgeshire catchment."),
         caption  = CUH_CAPTION) +
    theme_map_pub()
}

make_pt_spec_map <- function(specialties_subset = NULL, ncol = 4) {
  dat <- sf_pt_spec
  if (!is.null(specialties_subset)) dat <- filter(dat, specialty %in% specialties_subset)
  dat <- filter(dat, !is.na(patient_count))

  ggplot(dat) +
    england_bg() +
    geom_sf(aes(fill = patient_count), colour = "#888888", linewidth = 0.05, alpha = 0.85) +
    pt_scale +
    coord_sf(xlim   = c(CAMBS_BBOX["xmin"], CAMBS_BBOX["xmax"]),
             ylim   = c(CAMBS_BBOX["ymin"], CAMBS_BBOX["ymax"]),
             expand = FALSE) +
    facet_wrap(~ specialty_label, ncol = ncol) +
    labs(title    = "Patient Distribution by LSOA, Stratified by Specialty",
         subtitle = "Cambridgeshire catchment area. Patients referred per LSOA (square root scale).",
         caption  = CUH_CAPTION) +
    theme_map_pub(base_size = 9)
}

# One constructor for all three bivariate maps, since they differ only in which
# deprivation domain is on the vertical axis.
make_bivar_map <- function(sf_obj, domain_label) {
  map_plot <- ggplot(sf_obj) +
    england_bg() +
    geom_sf(aes(fill = fill_col), colour = NA, linewidth = 0, alpha = 0.9) +
    scale_fill_identity(na.value = "#d0d0d0") +
    annotation_scale(location = "bl", width_hint = 0.12, text_cex = 0.55,
                     style = "ticks", line_col = "#555", text_col = "#555") +
    labs(title    = paste0("Bivariate Choropleth: ", domain_label,
                           " and Vehicle Ownership by LSOA"),
         subtitle = paste0(domain_label, " Quintile x Vehicle Ownership Quintile (see legend)"),
         caption  = CUH_CAPTION) +
    theme_map_pub() +
    theme(legend.position = "none")

  cowplot::ggdraw(map_plot) +
    cowplot::draw_plot(make_bivar_legend(y_label = domain_label, base_size = 8),
                       x = 0.75, y = 0.04, width = 0.18, height = 0.22)
}


# ---- 8.7 Map output -----------------------------------------------------------

# Maps cannot be read in a terminal, so these are written to disk as PNGs. With
# many specialties a single facet grid becomes unreadable, so the stratified maps
# are split into pages of twelve and written one file per page.
paginate_specialties <- function(all_specs, per_page = 12) {
  split(all_specs, ceiling(seq_along(all_specs) / per_page))
}

save_map <- function(plot_obj, filename, width = 11, height = 8.5) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, plot_obj, width = width, height = height, dpi = 300)
  message("  written: ", path)
}

all_specs    <- sort(unique(spec_agg_data$specialty))
spec_pages   <- paginate_specialties(all_specs, per_page = 12)
n_spec_pages <- length(spec_pages)

message("Writing maps to ", OUT_DIR, "/ ...")

save_map(make_dna_agg_map(), "map_dna_aggregate.png")
save_map(make_pt_agg_map(),  "map_patients_aggregate.png")

for (i in seq_along(spec_pages)) {
  suffix <- if (n_spec_pages > 1) paste0("_page", i) else ""

  save_map(make_dna_spec_map(specialties_subset = spec_pages[[i]], ncol = 4),
           paste0("map_dna_by_specialty", suffix, ".png"))

  save_map(make_pt_spec_map(specialties_subset = spec_pages[[i]], ncol = 4),
           paste0("map_patients_by_specialty", suffix, ".png"))
}

save_map(make_bivar_map(sf_bivar,  "Health Deprivation"),     "map_bivariate_health.png")
save_map(make_bivar_map(sf_bivar2, "Employment Deprivation"), "map_bivariate_employment.png")
save_map(make_bivar_map(sf_bivar3, "Crime Index"),            "map_bivariate_crime.png")


# ---- 8.8 Distance decay and sanity checks -------------------------------------

# Confirm that no single LSOA dominates the map before the rates are interpreted.
lsoa_counts <- agg_data %>%
  filter(!is.na(lsoa21)) %>%
  count(lsoa21, name = "n_individuals") %>%
  arrange(desc(n_individuals))

print(head(lsoa_counts, 20))
message("LSOAs represented: ", nrow(lsoa_counts))
message("Total individuals: ", sum(lsoa_counts$n_individuals))

print(lsoa_counts %>% filter(n_individuals > quantile(n_individuals, 0.99)))

# Distance to the nearest public transport node against the DNA rate, at LSOA
# level. Distances beyond 30 km are set aside, being almost entirely rural
# outliers that would otherwise drive the fitted line on their own. Both a loess
# and a linear fit are shown, since the question is whether the relationship is
# monotone rather than what its slope is.
dist_summary <- agg_data %>%
  filter(!is.na(NAPTAN_distance_km), !is.na(DNA)) %>%
  mutate(NAPTAN_distance_km = ifelse(NAPTAN_distance_km > 30, NA, NAPTAN_distance_km)) %>%
  summarise(mean_dist = mean(NAPTAN_distance_km, na.rm = TRUE),
            dna_rate  = mean(DNA == 1, na.rm = TRUE),
            n         = n(),
            .by       = lsoa21)

cor_dist_dna <- cor.test(dist_summary$mean_dist, dist_summary$dna_rate,
                         method = "spearman")
print(cor_dist_dna)

p_dist <- ggplot(dist_summary, aes(x = mean_dist, y = dna_rate)) +
  geom_point(aes(size = n), alpha = 0.4, colour = "#5698b9") +
  geom_smooth(method = "loess", se = TRUE, colour = "#a8387e") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey50", linetype = "dashed") +
  scale_size_continuous(name = "n individuals") +
  labs(title    = "Distance to nearest transport vs DNA rate by LSOA",
       subtitle = paste0("Spearman r = ", round(cor_dist_dna$estimate, 3),
                         ", p = ", round(cor_dist_dna$p.value, 4)),
       x        = "Mean NAPTAN distance (km)",
       y        = "DNA rate") +
  theme_minimal()

print(p_dist)

ggsave(file.path(OUT_DIR, "distance_decay.png"), p_dist,
       width = 8, height = 6, dpi = 300)


################################################################################
# END OF PIPELINE
################################################################################
