# Import Library
library(tidyverse)
library(data.table)


################################################################################# I ### Load data
load(paste0(here::here(),
            "/data/raw_data",
            "/CHinORIEN_AllTumorTypes_Files_20260730.RData"))

path_raw <- fs::path("", "Volumes", "Gillis_Research",
                     "Lab_Data", "CHinORIEN")

lung_patients <- 
  read.delim(paste0(path_raw,
                    "/ProcessedData/SampleList/AllLung_Samplelist_wRNA.txt"))

# ClinicalMolLinkage file in Normalized folder doesn't have WESid
# Use the file from first folder
ClinicalMolLinkage <- read_csv(paste0(
  path_raw,
  "/RawData",
  "/23PRJ127MCC_ClinicalData_20230824",
  "/23PRJ127MCC_20230620_ClinicalMolLinkage_V4.csv"))

parent_dir_path <- dirname(path_raw)
drug_class <- 
  read.csv(paste0(#here::here(), 
    # "/data/processed_data",
    parent_dir_path,
    "/SharedResources/BoltonDrugCategories",
    "/CHevolution_Updated_BoltonChemoDosing_20260713.csv"))

removed_drug <- 
  read.csv(paste0(#here::here(), 
    # "/data/processed_data",
    parent_dir_path,
    "/SharedResources/NotCancerDrugs",
    "/CHinOvary_NotCancerDrug_ToRemovedFromData_20260617.csv"))


################################################################################# II ### Data cleaning
ClinicalMolLinkage <- ClinicalMolLinkage %>% 
  select(ORIENAvatarKey, sample_tumor_Disease.Type = Disease.Type, 
         SpecimenSiteOfCollection, WES, RNASeq, Age.At.Specimen.Collection)

lung_patients1 <- lung_patients %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(-c(RNASeq, sample_tumor_Disease.Type), 
                     age_at_germline_collection = Age.At.Specimen.Collection, 
                     germline_SpecimenSiteOfCollection = SpecimenSiteOfCollection), 
            by = c("ORIENAvatarKey", "Germline" = "WES")) %>% 
  # rename() %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(-RNASeq, 
                     age_at_tunor_collection = Age.At.Specimen.Collection,
                     tumor_SpecimenSiteOfCollection = SpecimenSiteOfCollection), 
            by = c("ORIENAvatarKey", "Tumor_WES" = "WES")) %>% 
  rename(Germline_WES = Germline) %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(ORIENAvatarKey, RNASeq, age_at_rnaseq_collection = Age.At.Specimen.Collection) %>% 
              distinct(RNASeq, .keep_all = TRUE), 
            by = c("ORIENAvatarKey", "RNASeq")) %>% 
  mutate(age_at_germline_collection_is_more_than_90 = case_when(
    age_at_germline_collection == "Age 90 or older"                  ~ "Yes"
  )) |> 
  mutate(across(c("age_at_germline_collection",
                  "age_at_tunor_collection",
                  "age_at_rnaseq_collection"), ~ case_when(
                    . == "Age 90 or older"                           ~ 90,
                    . == "Unknown/Not Applicable"                    ~ NA_real_,
                    TRUE                                             ~ as.numeric(.)
                  )))

# Cytogenetic ----
CytogeneticAbnormalities <- CytogeneticAbnormalities %>%
  filter(CytogenAbnormResult == "Positive")

CytogeneticAbnormalities <-
  CytogeneticAbnormalities %>%
  select(AvatarKey, CytogenAbnormName, CytogenAbnormInd) %>%
  distinct() %>%
  pivot_wider(id_cols = c(AvatarKey),
              names_from = CytogenAbnormName,
              values_from = CytogenAbnormInd)

# Demographics ----
demographics <- PatientMaster %>%
  mutate(race_all = case_when(
    str_detect(Race, "American Indian")                              ~ "American Indian or Alaska Native",
    Race == "White"                                                  ~ "White",
    str_detect(Race, "Black")                                        ~ "Black",
    str_detect(Race, "Asian") |
      str_detect(Race, "Cambodian") |
      Race == "Chinese" |
      Race == "Filipino" |
      Race == "Japanese" |
      Race == "Korean" |
      Race == "Pakistani" |
      Race == "Thai" |
      Race == "Vietnamese" |
      Race == "Laotian"                                              ~ "Asian",
    str_detect(Race, "Hawaiian") |
      Race == "Micronesian, NOS" |
      Race == "Pacific Islander, NOS" |
      Race == "Polynesian, NOS" |
      Race == "Samoan" |
      Race == "Tongan"                                               ~ "Native Hawaiian or Other Pacific Islander",
    str_detect(Race, "Unknown")                                      ~ NA_character_,
    Race == "Some other race"                                        ~ "Other",
    TRUE                                                             ~ Race
  )) %>%
  mutate(race = case_when(
    race_all == "White"                                              ~ "White",
    race_all == "Black"                                              ~ "Black",
    race_all == "Asian"                                              ~ "Asian",
    !is.na(race_all)                                                 ~ "Other"
  )) |> 
  mutate(Ethnicity = case_when(
    Ethnicity == "Spanish surname only"                              ~ "Non-Hispanic",
    str_detect(Ethnicity, "Non-Hispanic")                            ~ "Non-Hispanic",
    str_detect(Ethnicity, "Unknown")                                 ~ NA_character_,
    TRUE                                                             ~ "Hispanic" # all other are Hispanic
  )) %>%
  mutate(race_eth = case_when(
    race == "White" &
      Ethnicity == "Non-Hispanic"                                    ~ "White, Non-Hispanic",
    race == "Black" &
      Ethnicity == "Non-Hispanic"                                    ~ "Black, Non-Hispanic",
    race == "Others" &
      Ethnicity == "Non-Hispanic"                                    ~ "Others, Non-Hispanic",
    Ethnicity == "Hispanic"                                          ~ "Hispanic, any race"
  ))# %>% 
# mutate(Race = factor(Race, levels = c("White", "Black", "Asian",
#                                       "American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander",
#                                       "Unknown")))

# Patient History----
PatientHistory <- PatientHistory %>%
  mutate(SmokingStatus = case_when(
    SmokingStatus == "Never"                                         ~ "Never",
    SmokingStatus == "Current" |
      str_detect(SmokingStatus, "Ever")                              ~ "Ever"
  ), SmokingStatus = factor(SmokingStatus, levels = c("Never", "Ever"))) %>%
  mutate(AlcoholUse = case_when(
    AlcoholUse == "Never"                                            ~ "Never",
    AlcoholUse == "Current" |
      str_detect(AlcoholUse, "Ever")                                 ~ "Ever"
  ), AlcoholUse = factor(AlcoholUse, levels = c("Never", "Ever")))

# Diagnosis ----
Diagnosis_save <- Diagnosis
Diagnosis <- Diagnosis_save

Diagnosis <- Diagnosis %>%
  filter(str_detect(AvatarKey, paste0(lung_patients1$ORIENAvatarKey, collapse = "|"))) %>%
  group_by(AvatarKey) %>% 
  mutate(number_of_dx = n(), .before = AgeAtDiagnosis) %>% 
  ungroup() |> 
  
  # Few missing age at dx
  # It is ok to use the age at first contact to fill missing age at dx
  mutate(not_real_dxage = case_when(
    AgeAtDiagnosis == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtDiagnosis", 
                  "AgeAtFirstContact"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>%
  mutate(dx_age_is_first_contact_age = case_when(
    !is.na(AgeAtDiagnosis)                ~ "Age is from AgeAtDiagnosis",
    is.na(AgeAtDiagnosis) &
      !is.na(AgeAtFirstContact)           ~ "Age is from AgeAtFirstContact"
  ), .after = AgeAtDiagnosis) %>% 
  mutate(AgeAtDiagnosis = coalesce(AgeAtDiagnosis, AgeAtFirstContact)) %>% 
  arrange(AvatarKey, AgeAtDiagnosis) %>% 
  # Create var specific to ovarian dx
  group_by(AvatarKey) %>% 
  mutate(diagnosis_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  ungroup() %>% 
  mutate(is_lung_dx = case_when(
    str_detect(PrimaryDiagnosisSite, "lung")    ~ "Yes",
    str_detect(PrimaryDiagnosisSite, "Lung")    ~ "Yes",
    str_detect(PrimaryDiagnosisSite, "Main bronchus")    ~ "Yes",
    PrimaryDiagnosisSite == "Head of pancreas" &
      AgeAtDiagnosis == 59.123                  ~ "Yes"
  )) |> 
  # group_by(AvatarKey) |> 
  # fill(is_lung_dx, .direction = "updown") |>
  # ungroup()
  group_by(AvatarKey, is_lung_dx) |> 
  mutate(number_of_lung_dx = case_when(
    is_lung_dx == "Yes"          ~ n()
  ), .before = AgeAtDiagnosis) %>% 
  ungroup() |> 
  filter(is_lung_dx == "Yes") |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  select(AvatarKey, number_of_dx, diagnosis_sequence,
         AgeAtFirstContact,
         AgeAtDiagnosis, not_real_dxage, AgeAtDiagnosisFlag, dx_age_is_first_contact_age,
         YearOfDiagnosis,
         PrimaryDiagnosisSiteCode : Histology,
         ClinGroupStage, PathGroupStage,
         GradeClinical, GradePathological,
         CurrentlySeenForPrimaryOrRecurr,
         PerformStatusAtDiagnosis, 
         OtherStagingSystem, OtherStagingValue, everything(),
         -c(is_lung_dx, AgeAtFirstContactFlag))

Diagnosis <- Diagnosis %>%
  mutate(ECOG = str_match(PerformStatusAtDiagnosis, "ECOG ([:digit:])")[,2], 
         .after = PerformStatusAtDiagnosis) %>%
  mutate(Karnofsky = str_match(PerformStatusAtDiagnosis, "Karnofsky ([:digit:].*)%")[,2], 
         .after = PerformStatusAtDiagnosis) %>%
  # Stage
  mutate(ClinGroupStage = case_when(
    str_detect(ClinGroupStage, "IV")                        ~ "4",
    str_detect(ClinGroupStage, "III")                       ~ "3",
    str_detect(ClinGroupStage, "II")                        ~ "2",
    str_detect(ClinGroupStage, "I")                         ~ "1",
    str_detect(ClinGroupStage, "0")                         ~ "0"
  )) %>% 
  mutate(PathGroupStage = case_when(
    str_detect(PathGroupStage, "IV")                        ~ "4",
    str_detect(PathGroupStage, "III")                       ~ "3",
    str_detect(PathGroupStage, "II")                        ~ "2",
    str_detect(PathGroupStage, "I")                         ~ "1",
    str_detect(PathGroupStage, "0")                         ~ "0"
  )) %>% 
  group_by(AvatarKey) %>% 
  mutate(stage = max(ClinGroupStage, PathGroupStage, na.rm = TRUE), .after = PathGroupStage) %>%
  ungroup()

VitalStatus <- VitalStatus %>%
  filter(str_detect(AvatarKey, paste0(lung_patients1$ORIENAvatarKey, collapse = "|"))) %>%
  select(AvatarKey, VitalStatus, AgeAtLastContact, AgeAtDeath, CauseOfDeath) %>%
  mutate(across(c("AgeAtLastContact", "AgeAtDeath"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

# Medications ----
Medications1 <- Medications %>% 
  filter(str_detect(AvatarKey, paste0(lung_patients1$ORIENAvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtMedStart",
                  "AgeAtMedStop"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  # Fix drug names I already know about
  mutate(Medication = str_remove(
    Medication, 
    " Hydrochloride|Liposomal | Sulfate| Phosphate| Citrate| Acetate| Camsylate| Disodium| Mesylate| Ditosylate| Tosylate| Tartrate| \\(Leucovorin\\)"), 
    Medication = case_when(
      str_detect(Medication, "Paclitaxel")        ~ "Paclitaxel",
      str_detect(Medication, "^Interferon")       ~ "Interferon",
      Medication == "Bevacizumab-maly"            ~ "Bevacizumab",
      Medication == "Bevacizumab-adcd"            ~ "Bevacizumab",
      TRUE                                        ~ Medication
    )) %>% 
  mutate(Medication = str_to_lower(Medication)) %>% 
  # Remove non cancer drugs from data so their date start will not count
  left_join(., removed_drug, by = c("Medication" = "drug_name")) %>% 
  filter(is.na(drug_type)) %>% select(-drug_type) %>% 
  # Add Bolton chemo cat
  left_join(., drug_class, by = c("Medication" = "drug_name")) %>% 
  # and code an overall chemotherapy variable
  mutate(is_chemotherapy = case_when(
    !is.na(narrow_drug_class_cytotoxic_only)        ~ "Yes",
    Medication == "chemo, nos"                      ~ "Yes"
  )) |> 
  select(AvatarKey, MedicationInd, is_chemotherapy,
         Medication, MedLineRegimen,
         AgeAtMedStart, AgeAtMedStop, 
         narrow_drug_class_cytotoxic_only, general_drug_class,
         MedContinuing, ChangeOfTreatment,
         everything())

# Separate patient who didn't receive drugs
Medications_never <- Medications1 %>% 
  filter(MedicationInd == "No") %>% 
  mutate(has_first_line = "Never received any drug", .after = MedicationInd) %>%
  mutate(is_chemotherapy =  "Never received any drug", .after = MedicationInd)

Medications_yes <- Medications1 %>% 
  # Filter patients who received drugs
  filter(MedicationInd == "Yes") %>% 
  # Fix age for a couple of patients
  mutate(AgeAtMedStart = coalesce(AgeAtMedStart, AgeAtMedStop)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_med_age = min(AgeAtMedStart, na.rm = TRUE),
         ever_first_med_age = na_if(ever_first_med_age, Inf)) %>% 
  ungroup() %>% 
  # recode line
  # There are regimen line are unknown 
  # but some have the same age as other row for which the line is known, use to fill it up
  mutate(MedLineRegimen = case_when( # m
    str_detect(MedLineRegimen, "Unknown")            ~ NA_character_,
    TRUE                                             ~ MedLineRegimen
  )) %>%
  group_by(AvatarKey, AgeAtMedStart) %>% 
  fill(MedLineRegimen, .direction = "updown") %>% 
  ungroup() %>% 
  # Do the reverse
  group_by(AvatarKey, MedLineRegimen) %>% 
  fill(AgeAtMedStart, .direction = "updown") %>% 
  ungroup() %>% 
  # Organize vars
  distinct() %>% 
  # Recode line with numbers to be able to sort and make dense rank
  mutate(regimen_line = case_when(
    is.na(MedLineRegimen)                            ~ 999,
    str_detect(MedLineRegimen, "First") &
      str_detect(MedLineRegimen, "Neoadjuvant")      ~ -1,
    str_detect(MedLineRegimen, "First") &
      str_detect(MedLineRegimen, "Adjuvant")         ~ 1,
    str_detect(MedLineRegimen, "Second")             ~ 2,
    str_detect(MedLineRegimen, "Third")              ~ 3,
    str_detect(MedLineRegimen, "Fourth")             ~ 4,
    str_detect(MedLineRegimen, "Fifth")              ~ 5,
    str_detect(MedLineRegimen, "Sixth")              ~ 6,
    str_detect(MedLineRegimen, "Seventh")            ~ 7,
    str_detect(MedLineRegimen, "Eighth")             ~ 8,
    str_detect(MedLineRegimen, "Ninth")              ~ 9,
    str_detect(MedLineRegimen, "Tenth")              ~ 10,
    str_detect(MedLineRegimen, "Eleventh")           ~ 11,
    str_detect(MedLineRegimen, "Twelfth")            ~ 12,
    MedLineRegimen == "Consolidation"                ~ 90,
    MedLineRegimen == "Maintenance"                  ~ 91,
    MedLineRegimen == "Palliative"                   ~ 92,
    TRUE                                             ~ 1000 # none?
  ), .after = MedLineRegimen) %>% 
  # Manually code line if thee are no line info at all
  arrange(AvatarKey, AgeAtMedStart) %>% 
  # mutate(no_line_info = case_when(
  #   regimen_line == 999                              ~ 999
  # )) %>% 
  # group_by(AvatarKey) %>% 
  # fill(no_line_info, .direction = "updown") %>% 
  # ungroup() %>% 
  group_by(AvatarKey) %>% 
  mutate(regimen_line1 = dense_rank(interaction(regimen_line, AgeAtMedStart)), .after = regimen_line) %>% 
  ungroup() |> 
  # # Removed this. I cheched the patient that have a line + 999. They are not abstracted well.
  # # It is better to use the coded regimen line.
  # mutate(has_some_sort_line_info = case_when(
  #   no_line_info == 999 &
  #     !is.na(MedLineRegimen) &
  #     ((regimen_line == 1 |
  #        regimen_line == -1) &
  #     regimen_line2 == 1)                            ~ "take regimen_line2"
  # )) %>% 
  # # select(AvatarKey, MedLineRegimen, regimen_line, no_line_info, has_some_sort_line_info, AgeAtMedStart) %>% 
  # group_by(AvatarKey) %>% 
  # fill(has_some_sort_line_info, .direction = "updown") %>% 
  # ungroup() %>% 
  # mutate(has_some_sort_line_info = case_when(
  #   has_some_sort_line_info == "take regimen_line2"  ~ "take regimen_line2",
  #   no_line_info == 999 &
  #     !is.na(MedLineRegimen)                         ~ "has_some_sort_line_info"
  # )) %>% 
  # group_by(AvatarKey) %>% 
  # fill(has_some_sort_line_info, .direction = "updown") %>% 
  # ungroup() %>% 
  # select(-c(no_line_info, has_some_sort_line_info, regimen_line2))
  arrange(AvatarKey, regimen_line, AgeAtMedStart, Medication)

Medications_final <- bind_rows(Medications_yes, Medications_never) %>% 
  distinct(AvatarKey, AgeAtMedStart, Medication, AgeAtMedStop, .keep_all = TRUE) %>% 
  mutate(ever_received_chemotherapy = case_when(
    is_chemotherapy == "Yes"               ~ "Ever"
  )) %>% 
  group_by(AvatarKey, ever_received_chemotherapy) %>% 
  mutate(year_first_chemo = case_when(
    is_chemotherapy == "Yes"               ~ first(YearOfMedStart)
  ), .after = YearOfMedStart) %>% 
  
  group_by(AvatarKey) %>% 
  fill(ever_received_chemotherapy, year_first_chemo, .direction = "updown") %>% 
  mutate(year_first_medication = first(YearOfMedStart), .after = YearOfMedStart) %>% 
  group_by(AvatarKey, is_chemotherapy) %>% 
  mutate(first_chemo_age = case_when(
    is_chemotherapy == "Yes"            ~ first(AgeAtMedStart)
  ), .after = is_chemotherapy) %>% 
  ungroup() %>% 
  group_by(AvatarKey) %>% 
  fill(first_chemo_age, .direction = "updown") %>% 
  ungroup()

# Create regimen - keep as "long" data
Medications_regimen <- Medications_final %>%
  # Same age start and stop
  group_by(AvatarKey, regimen_line, 
           AgeAtMedStart, MedicationInd, AgeAtMedStop, 
           ever_received_chemotherapy, first_chemo_age, year_first_chemo,
           ever_first_med_age, year_first_medication
  ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class), str_c, collapse = "; ") %>%
  # same age start
  group_by(AvatarKey, regimen_line, 
           AgeAtMedStart, MedicationInd, 
           ever_received_chemotherapy, first_chemo_age, year_first_chemo,
           ever_first_med_age, year_first_medication
  ) %>%
  summarise_at(vars(Medication, 
                    narrow_drug_class_cytotoxic_only, 
                    general_drug_class,
                    AgeAtMedStop), str_c, collapse = "; ") %>%
  mutate(general_drug_class = paste(unique(unlist(strsplit(general_drug_class, "; "))), 
                                    collapse = "; ")) %>% 
  mutate(narrow_drug_class_cytotoxic_only = paste(unique(unlist(strsplit(narrow_drug_class_cytotoxic_only, "; "))), 
                                                  collapse = "; ")) %>% 
  # # remove regimen line - not right if NAs
  # group_by(AvatarKey, AgeAtMedStart, MedicationInd, 
  #          ever_first_med_age, year_first_medication
  # ) %>%
  # summarise_at(vars(Medication, 
  #                   AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup() %>% 
  rename(regimen_drugname = Medication) %>% 
  arrange(AvatarKey, regimen_line, AgeAtMedStart)

# write_rds(Medications_regimen,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Regimen long format_",
#                  today(), ".rds"))
# write_rds(Medications_regimen,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Regimen long format_",
#                  today(), ".csv"))

Medications_wide <- dcast(setDT(Medications_regimen),
                          AvatarKey + MedicationInd + 
                            ever_received_chemotherapy + first_chemo_age + year_first_chemo +
                            ever_first_med_age + year_first_medication
                          ~ rowid(AvatarKey),
                          value.var = c("AgeAtMedStart", "regimen_drugname", "AgeAtMedStop",
                                        "narrow_drug_class_cytotoxic_only", "general_drug_class")) %>% 
  rename(drug_ever = MedicationInd) %>%
  mutate(has_medication_data = "Yes")

# write_rds(Medications_wide,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Regimen wide format_",
#                  today(), ".rds"))
# write_rds(Medications_wide,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Regimen wide format_",
#                  today(), ".csv"))

# radiation----
Radiation <- Radiation %>%
  filter(str_detect(AvatarKey, paste0(lung_patients1$ORIENAvatarKey, collapse = "|"))) %>% 
  mutate(across(c("AgeAtRadiationStart",
                  "AgeAtRadiationStop"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>% 
  mutate(RadDose = as.numeric(RadDose)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_rad_age = min(AgeAtRadiationStart, na.rm = TRUE),
         ever_first_rad_age = na_if(ever_first_rad_age, Inf)) %>% 
  mutate(year_first_radiation = first(YearOfRadiationStart), .after = YearOfRadiationStart) %>% 
  ungroup() %>% 
  distinct(AvatarKey, AgeAtRadiationStart, RadDose, AgeAtRadiationStop, .keep_all = TRUE)

# write_rds(Radiation,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Radiation long format_",
#                  today(), ".rds"))
# write_rds(Radiation,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Radiation long format_",
#                  today(), ".csv"))

Radiation_wide <- dcast(setDT(Radiation),
                        AvatarKey + RadiationTherapyInd + ever_first_rad_age + year_first_radiation
                        ~ rowid(AvatarKey),
                        value.var = c("AgeAtRadiationStart", "RadDose", "RadFractions", "AgeAtRadiationStop")) %>% 
  select(AvatarKey, radiation_ever = RadiationTherapyInd, ever_first_rad_age, 
         starts_with("AgeAtRadiationStart_"),
         starts_with("RadDose_"), starts_with("AgeAtMedStop_"), everything()) %>%
  mutate(has_radiation_data = "Yes")

# write_rds(Radiation_wide,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Radiation wide format_",
#                  today(), ".rds"))
# 
# write_rds(Radiation_wide,
#           paste0(here::here(),
#                  "/data/processed_data",
#                  "/Radiation wide format_",
#                  today(), ".rds"))

# PFS ----
Outcomes_save <- Outcomes
Outcomes <- Outcomes_save

Outcomes <- Outcomes %>%
  filter(str_detect(AvatarKey, paste0(lung_patients1$ORIENAvatarKey, collapse = "|"))) %>% 
  mutate(across(everything(), ~ 
                  str_replace_all(., 
                                  "Unknown/Not Applicable|Unknown/Not Performed|Age Unknown/Not Recorded",
                                  NA_character_))) |> 
  distinct() |> 
  purrr::keep(~!all(is.na(.))) |> 
  # mutate(is_ovary_outcome = case_when(
  #   OutcomesPrimaryDiagnosisSite %in% c(
  #     "Ovary", "Fallopian tube", 
  #     "Peritoneum, NOS", 
  #     "Specified parts of peritoneum",
  #     "Overlapping lesion of female genital organs")      ~ "Yes"
  # )) %>% 
  mutate(not_real_pfsage = case_when(
    AgeAtProgRecur == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  mutate(across(c("AgeAtProgRecur",
                  "AgeAtCurrentDiseaseStatus"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  ))) %>%
  select(AvatarKey, ProgRecurInd, 
         AgeAtProgRecur, 
         YearOfProgRecur, AgeAtCurrentDiseaseStatus, AgeAtPerformStatusMostRecent,
         OutcomesPrimaryDiagnosisSite, OutcomesPrimaryDiagnosisSiteCode, 
         everything())

Outcomes1 <- Outcomes |> 
  mutate(is_lung_pfs_info = case_when(
    str_detect(OutcomesPrimaryDiagnosisSite, "lung")    ~ "Yes",
    str_detect(OutcomesPrimaryDiagnosisSite, "Lung")    ~ "Yes",
    str_detect(OutcomesPrimaryDiagnosisSite, "Main bronchus")    ~ "Yes",
    OutcomesPrimaryDiagnosisSite == "Head of pancreas" &
      AgeAtProgRecur == 62.838                          ~ "Yes"
  )) %>%
  inner_join(., Diagnosis %>%
               select(AvatarKey, number_of_dx),
             by = "AvatarKey") |> 
  # 4 patients have outcomes info for missing site but they all have only 1 dx 
  # keep info into data and age meet dx age
  # group_by(AvatarKey) |>
  # fill(is_lung_pfs_info, .direction = "updown") |>
  # ungroup() |> 
  # Patients with oucomes data for another site but none for lung
  # Those outcomes are truely for the other site - exclude
  filter(is_lung_pfs_info == "Yes" | 
           (is.na(is_lung_pfs_info) & is.na(OutcomesPrimaryDiagnosisSite))|
           (is.na(is_lung_pfs_info) & number_of_dx == 1)) |> 
  select(-number_of_dx) |> 
  mutate(age_at_disease_check = coalesce(AgeAtProgRecur, AgeAtCurrentDiseaseStatus), 
         .after = AgeAtCurrentDiseaseStatus) %>% 
  arrange(AvatarKey, age_at_disease_check)

Outcomes2 <- Outcomes1 %>%
  filter(!is.na(ProgRecurInd)) |> 
  mutate(AgeAtProgRecur = coalesce(AgeAtProgRecur, AgeAtCurrentDiseaseStatus)) |> 
  arrange(AvatarKey, AgeAtProgRecur) |> 
  distinct(AvatarKey, .keep_all = TRUE)

pfs_free <- Outcomes1 %>%
  filter(is.na(ProgRecurInd)) |> 
  filter(!str_detect(AvatarKey, paste0(Outcomes2$AvatarKey, collapse = "|"))) |> 
  arrange(AvatarKey, age_at_disease_check) |> 
  distinct(AvatarKey, .keep_all = TRUE) |> 
  mutate(ProgRecurInd = "No")

Outcomes <- bind_rows(Outcomes2, pfs_free) |> 
  mutate(has_outcomes_data = "Yes") |> 
  distinct(AvatarKey, .keep_all = TRUE)
  

################################################################################# III ### Merging
data <- lung_patients1 %>% 
  rename(AvatarKey = ORIENAvatarKey) %>%
  left_join(., demographics, by = "AvatarKey") %>% 
  left_join(., Diagnosis, by = "AvatarKey") %>% 
  left_join(., VitalStatus, by = "AvatarKey") %>% 
  left_join(., Outcomes, by = "AvatarKey") %>% 
  left_join(., Medications_wide, by = "AvatarKey") %>% 
  left_join(., Radiation_wide, by = "AvatarKey") %>% 
  left_join(., PatientHistory, by = "AvatarKey") %>% 
  mutate(drug_radiation_ever = case_when(
    drug_ever == "Yes" &
      radiation_ever == "Yes"                              ~ "Drug+Radiation",
    drug_ever == "Yes"                                     ~ "Drug only",
    radiation_ever == "Yes"                                ~ "Radiation only",
    drug_ever == "No" &
      radiation_ever == "No"                               ~ "No drug or radiation",
    drug_ever == "No"                                      ~ "Not received drug",
    radiation_ever == "No"                                 ~ "Not received rad"
  )) %>% 
  
  mutate(blood_before_drug = case_when(
    age_at_germline_collection <= AgeAtMedStart_1          ~ "blood before drug",
    age_at_germline_collection > AgeAtMedStart_1           ~ "blood after drug",
    drug_ever == "No"                                      ~ "No drug received",
    is.na(age_at_germline_collection)                      ~ "no germline age",
    is.na(AgeAtMedStart_1)                                 ~ "missing drug age"
  ), .after = age_at_germline_collection) %>% 
  mutate(years_blood_before_drug = case_when(
    blood_before_drug == "blood before drug"               ~ age_at_germline_collection - AgeAtMedStart_1
  )) %>% 
  mutate(years_blood_after_drug = case_when(
    blood_before_drug == "blood after drug"                ~ AgeAtMedStart_1 - age_at_germline_collection
  )) %>% 
  mutate(blood_before_rad = case_when(
    age_at_germline_collection <= AgeAtRadiationStart_1    ~ "blood before rad",
    age_at_germline_collection > AgeAtRadiationStart_1     ~ "blood after rad",
    radiation_ever == "No"                                 ~ "No rad received",
    is.na(age_at_germline_collection)                      ~ "no germline age",
    is.na(AgeAtRadiationStart_1)                           ~ "missing rad age",
  ), .after = blood_before_drug) %>% 
  mutate(years_blood_before_rad = case_when(
    blood_before_rad == "blood before rad"                 ~ age_at_germline_collection - AgeAtRadiationStart_1
  )) %>% 
  mutate(years_blood_after_rad = case_when(
    blood_before_rad == "blood after rad"                  ~ AgeAtRadiationStart_1 - age_at_germline_collection
  )) %>% 
  mutate(blood_before_drugorrad = case_when(
    age_at_germline_collection <= AgeAtMedStart_1 &
      age_at_germline_collection <= AgeAtRadiationStart_1  ~ "blood before drug/rad",
    age_at_germline_collection > AgeAtMedStart_1 &
      age_at_germline_collection > AgeAtRadiationStart_1   ~ "blood after drug/rad",
    drug_ever == "No" &
      radiation_ever == "No"                               ~ "No drug/rad received",
    
    radiation_ever == "No" &
      age_at_germline_collection <= AgeAtMedStart_1        ~ "blood before drug/rad",
    radiation_ever == "No" &
      age_at_germline_collection > AgeAtMedStart_1         ~ "blood after drug/rad",
    drug_ever == "No" &
      age_at_germline_collection <= AgeAtRadiationStart_1  ~ "blood before drug/rad",
    drug_ever == "No" &
      age_at_germline_collection > AgeAtRadiationStart_1   ~ "blood after drug/rad",
    
    drug_ever == "Yes" &
      radiation_ever == "No" &
      is.na(AgeAtMedStart_1)                               ~ "missing drug age",
    drug_ever == "Yes" &
      radiation_ever == "No"                               ~ "Drug only",
    
    radiation_ever == "Yes" &
      drug_ever == "No" &
      is.na(AgeAtRadiationStart_1)                         ~ "missing rad age",
    radiation_ever == "Yes" &
      drug_ever == "No"                                    ~ "Radiation only",
    
    is.na(drug_ever)                                       ~ "don't know if ever received drug",
    is.na(radiation_ever)                                  ~ "don't know if ever received rad",
    
    drug_ever == "Yes" &
      radiation_ever == "Yes" &
      is.na(AgeAtMedStart_1) &
      is.na(AgeAtRadiationStart_1)                         ~ "missing drug and rad age",
    
    is.na(AgeAtMedStart_1)                                 ~ "missing drug age",
    is.na(AgeAtRadiationStart_1)                           ~ "missing rad age",
    
    age_at_germline_collection > AgeAtMedStart_1 &
      age_at_germline_collection <= AgeAtRadiationStart_1  ~ "drug/blood/rad",
    age_at_germline_collection <= AgeAtMedStart_1 &
      age_at_germline_collection > AgeAtRadiationStart_1  ~ "rad/blood/drug"
    
    
  ), .after = age_at_germline_collection) %>% 
  
  mutate(years_blood_before_drugrad = case_when(
    blood_before_drugorrad == "blood before drug/rad" &
      AgeAtMedStart_1 <= AgeAtRadiationStart_1       ~ age_at_germline_collection - AgeAtMedStart_1,
    blood_before_drugorrad == "blood before drug/rad" &
      AgeAtMedStart_1 > AgeAtRadiationStart_1        ~ age_at_germline_collection - AgeAtRadiationStart_1
  )) %>% 
  mutate(years_blood_after_drugrad = case_when(
    blood_before_drugorrad == "blood after drug/rad" &
      AgeAtMedStart_1 <= AgeAtRadiationStart_1       ~ age_at_germline_collection - AgeAtMedStart_1,
    blood_before_drugorrad == "blood after drug/rad" &
      AgeAtMedStart_1 > AgeAtRadiationStart_1        ~ age_at_germline_collection - AgeAtRadiationStart_1
  )) %>% 
  
  mutate(treatment_sequence = case_when(
    AgeAtMedStart_1 <= AgeAtRadiationStart_1               ~ "drug/rad",
    AgeAtMedStart_1 > AgeAtRadiationStart_1                ~ "rad/drug",
    drug_ever == "No" &
      radiation_ever == "No"                               ~ "No drug or radiation",
    drug_ever == "Yes" &
      radiation_ever == "No"                               ~ "Drug only",
    radiation_ever == "Yes" &
      drug_ever == "No"                                    ~ "Radiation only",
    is.na(AgeAtMedStart_1) &
      is.na(AgeAtRadiationStart_1)                         ~ "missing drug and rad age",
    is.na(AgeAtMedStart_1)                                 ~ "missing drug age",
    is.na(AgeAtRadiationStart_1)                           ~ "missing rad age"
  ), .after = blood_before_rad) %>% 
  select(AvatarKey : treatment_sequence, contains("ever_first"), everything())


# PFS and OS
data <- data |> 
  # Age at first treatment
  mutate(first_treatment = case_when(
      (AgeAtRadiationStart_1 < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(AgeAtRadiationStart_1)))              ~ "Radiation",
    
    (AgeAtMedStart_1 < AgeAtRadiationStart_1 | 
       (is.na(AgeAtRadiationStart_1) &
          !is.na(AgeAtMedStart_1)))                      ~ "Drugs"
  )) %>% 
  mutate(age_at_first_treatment = case_when(
    first_treatment == "Radiation"                  ~ AgeAtRadiationStart_1,
    first_treatment == "Drugs"                      ~ AgeAtMedStart_1
  )) %>% 
  # OS
  mutate(os_event = case_when(
    VitalStatus == "Alive"                          ~ 0,
    VitalStatus == "Dead"                           ~ 1
  )) %>% 
  mutate(os_age = coalesce(AgeAtDeath, AgeAtLastContact)) %>% 
  mutate(os_time_from_dx_years = os_age - AgeAtDiagnosis) %>% 
  mutate(os_time_from_treatment_years = os_age - age_at_first_treatment) %>% 
  # PFS
  mutate(pfs_event = case_when(
    ProgRecurInd == "Progression"                   ~ 1,
    ProgRecurInd == "Recurrence"                    ~ 1,
    os_event == 1                                   ~ 1,
    ProgRecurInd == "No"                            ~ 0
  )) %>% 
  mutate(pfs_age = case_when(
    ProgRecurInd == "Progression"                   ~ AgeAtProgRecur,
    ProgRecurInd == "Recurrence"                    ~ AgeAtProgRecur,
    os_event == 1                                   ~ os_age,
    ProgRecurInd == "No"                            ~ os_age
  )) %>% 
  mutate(pfs_time_from_dx_years = pfs_age - AgeAtDiagnosis) %>% 
  mutate(pfs_time_from_treatment_years = pfs_age - age_at_first_treatment)


write_rds(data,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Lung_CHinORIEN_",
                 str_remove_all(today(), "-"), ".rds"))

write_csv(data,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Lung_CHinORIEN_",
                 str_remove_all(today(), "-"), ".csv"))

write_rds(data,
          paste0(path_raw,
                 "/ProcessedData",
                 "/Lung_CHinORIEN_",
                 str_remove_all(today(), "-"), ".rds"))

write_csv(data,
          paste0(path_raw,
                 "/ProcessedData",
                 "/Lung_CHinORIEN_",
                 str_remove_all(today(), "-"), ".csv"))


# END ----            

