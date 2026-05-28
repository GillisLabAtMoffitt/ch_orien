# Import Library
library(tidyverse)
library(data.table)
# library(lubridate)


################################################################################# I ### Load data
# # Load each file in its own dataframe
# path_raw <- fs::path("", "Volumes", "Gillis_Research",
#                      "Lab_Data", "CHinORIEN")
# filenames <- list.files(path = paste0(path_raw,
#                                       "/RawData",
#                                       "/23PRJ127MCC_NormalizedFiles_20240715"),
#                         pattern = "23PRJ127MCC_20240715_.*.csv")
# names <- str_match(filenames, "23PRJ127MCC_20240715_(.*)_V4.csv")[,2]
# 
# for(i in names){
#     filepath <- file.path(paste0(path_raw,
#                                  "/RawData",
#                                  "/23PRJ127MCC_NormalizedFiles_20240715/",
#                                  paste("23PRJ127MCC_20240715_", i, "_V4.csv",sep="")))
#     assign(i, read.csv(filepath, na.strings = ""))
# }
# 
# rm(filenames, names, filepath, i)
# # Saved environment image in a .Rdata
# save.image(file = paste0(path_raw,
#                          "/RawData",
#                          "/23PRJ127MCC_NormalizedFiles_20240715/",
#                          "CHinORIEN_AllTumorTypes_Files.RData"))
# load(paste0(path_raw,
#             "/RawData",
#             "/23PRJ127MCC_NormalizedFiles_20240715/",
#             "CHinORIEN_AllTumorTypes_Files.RData"))
load(paste0(here::here(),
            "/data/raw_data",
            "/CHinORIEN_AllTumorTypes_Files.RData"))

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

################################################################################# II ### Treatment data cleaning
ClinicalMolLinkage <- ClinicalMolLinkage %>% 
  select(ORIENAvatarKey, WES, RNASeq, Age.At.Specimen.Collection)

lung_patients1 <- lung_patients %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(-RNASeq), 
            by = c("ORIENAvatarKey", "Germline" = "WES")) %>% 
  rename(age_at_germline_collection = Age.At.Specimen.Collection) %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(-RNASeq), 
            by = c("ORIENAvatarKey", "Tumor_WES" = "WES")) %>% 
  rename(age_at_tunor_collection = Age.At.Specimen.Collection) %>% 
  left_join(., ClinicalMolLinkage %>% 
              select(-WES) %>% 
              distinct(RNASeq, .keep_all = TRUE), 
            by = c("ORIENAvatarKey", "RNASeq")) %>% 
  rename(age_at_rnaseq_collection = Age.At.Specimen.Collection) %>% 
  mutate(across(c("age_at_germline_collection",
                  "age_at_tunor_collection",
                  "age_at_rnaseq_collection"), ~ case_when(
                    . == "Age 90 or older"                ~ 90,
                    . == "Unknown/Not Applicable"         ~ NA_real_,
                    TRUE                                  ~ as.numeric(.)
                  )))

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
    " Hydrochloride|Liposomal | Sulfate| Phosphate| Citrate| Acetate| Camsylate| Disodium| Mesylate| Ditosylate| Tosylate| Tartrate"), 
    Medication = case_when(
      str_detect(Medication, "Paclitaxel")        ~ "Paclitaxel",
      Medication == "Interferon, NOS"             ~ "Interferon",
      Medication == "Bevacizumab-maly"            ~ "Bevacizumab",
      TRUE                                        ~ Medication
    )) %>% 
  mutate(Medication = str_to_lower(Medication))

# Separate patient who didn't receive drugs
Medications_never <- Medications1 %>% 
  filter(MedicationInd == "No") #%>% 
  # mutate(has_first_line = "Never received any drug", .after = MedicationInd) %>% 
  # mutate(is_chemotherapy =  "Never received any drug", .after = MedicationInd)
  
Medications_yes <- Medications1 %>% 
  # Fix age for a couple of patients
  mutate(AgeAtMedStart = coalesce(AgeAtMedStart, AgeAtMedStop)) %>% 
  group_by(AvatarKey) %>% 
  mutate(ever_first_med_age = min(AgeAtMedStart, na.rm = TRUE),
         ever_first_med_age = na_if(ever_first_med_age, Inf)) %>% 
  ungroup() %>% 
  # Recode line with numbers to be able to sort and make dense rank
  filter(MedicationInd == "Yes") %>% 
  mutate(regimen_line = case_when(
    str_detect(MedLineRegimen, "Unknown")            ~ 999,
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
    MedLineRegimen == "Maintenance"                  ~ 90,
    MedLineRegimen == "Palliative"                   ~ 91,
    TRUE                                             ~ 1000 # none
  ), .after = MedLineRegimen) %>% 
  # Manually code line if thee are no line info at all
  arrange(AvatarKey, AgeAtMedStart) %>% 
  mutate(no_line_info = case_when(
    str_detect(MedLineRegimen, "Unknown")            ~ 999
  )) %>% 
  group_by(AvatarKey) %>% 
  fill(no_line_info, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(has_some_sort_line_info = case_when(
    no_line_info == 999 &
      !str_detect(MedLineRegimen, "Unknown")         ~ "has_some_sort_line_info"
  )) %>% 
  # select(AvatarKey, MedLineRegimen, regimen_line, no_line_info, has_some_sort_line_info, AgeAtMedStart) %>% 
  group_by(AvatarKey) %>% 
  fill(has_some_sort_line_info, .direction = "updown") %>% 
  ungroup() %>% 
  group_by(AvatarKey) %>% 
  mutate(regimen_line2 = dense_rank(interaction(regimen_line, AgeAtMedStart)), .after = regimen_line) %>% 
  ungroup() %>% 
  mutate(regimen_line = case_when(
    no_line_info == 999 &
      is.na(has_some_sort_line_info)                 ~ regimen_line2,
    TRUE                                             ~ regimen_line
  )) %>% 
  select(-c(no_line_info, has_some_sort_line_info, regimen_line2))

Medications_final <- bind_rows(Medications_yes, Medications_never) %>% 
  distinct(AvatarKey, AgeAtMedStart, Medication, AgeAtMedStop, .keep_all = TRUE) %>% 
  group_by(AvatarKey) %>% 
  mutate(year_first_medication = first(YearOfMedStart), .after = YearOfMedStart) %>% 
  ungroup()

# Create regimen - keep as "long" data
Medications_regimen <- Medications_final %>%
  # Same age start and stop
  group_by(AvatarKey, regimen_line, AgeAtMedStart, MedicationInd, AgeAtMedStop, 
           ever_first_med_age, year_first_medication
  ) %>%
  summarise_at(vars(Medication), str_c, collapse = "; ") %>%
  # same age start
  group_by(AvatarKey, regimen_line, AgeAtMedStart, MedicationInd, 
           ever_first_med_age, year_first_medication
  ) %>%
  summarise_at(vars(Medication, 
                    AgeAtMedStop), str_c, collapse = "; ") %>%
  # # remove regimen line - not right if NAs
  # group_by(AvatarKey, AgeAtMedStart, MedicationInd, 
  #          ever_first_med_age, year_first_medication
  # ) %>%
  # summarise_at(vars(Medication, 
  #                   AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup() %>% 
  rename(regimen_drugs = Medication) %>% 
  arrange(AvatarKey, AgeAtMedStart)

write_rds(Medications_regimen,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen long format_",
                 today(), ".rds"))
write_rds(Medications_regimen,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen long format_",
                 today(), ".csv"))

Medications_wide <- dcast(setDT(Medications_regimen),
                          AvatarKey + MedicationInd + ever_first_med_age + year_first_medication
                          ~ rowid(AvatarKey),
                          value.var = c("AgeAtMedStart", "regimen_drugs", "AgeAtMedStop")) %>% 
  rename(drug_ever = MedicationInd) %>%
  mutate(has_medication_data = "Yes") %>% 
  mutate(had_prior_medication = case_when(
    AgeAtMedStart_1 > ever_first_med_age       ~ "Yes"
  ))

write_rds(Medications_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen wide format_",
                 today(), ".rds"))
write_rds(Medications_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Regimen wide format_",
                 today(), ".csv"))

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

write_rds(Radiation,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Radiation long format_",
                 today(), ".rds"))
write_rds(Radiation,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Radiation long format_",
                 today(), ".csv"))

Radiation_wide <- dcast(setDT(Radiation),
                        AvatarKey + RadiationTherapyInd + ever_first_rad_age + year_first_radiation
                        ~ rowid(AvatarKey),
                        value.var = c("AgeAtRadiationStart", "RadDose", "AgeAtRadiationStop")) %>% 
  select(AvatarKey, radiation_ever = RadiationTherapyInd, ever_first_rad_age, 
         starts_with("AgeAtRadiationStart_"),
         starts_with("RadDose_"), starts_with("AgeAtMedStop_"), everything()) %>%
  mutate(has_radiation_data = "Yes")

write_rds(Radiation_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Radiation wide format_",
                 today(), ".rds"))

write_rds(Radiation_wide,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Radiation wide format_",
                 today(), ".rds"))

data <- lung_patients1 %>% 
  left_join(., Medications_wide, by = c("ORIENAvatarKey" = "AvatarKey")) %>% 
  left_join(., Radiation_wide, by = c("ORIENAvatarKey" = "AvatarKey")) %>% 
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
    drug_ever == "Yes" &
      radiation_ever == "No"                               ~ "Drug only",
    radiation_ever == "Yes" &
      drug_ever == "No"                                    ~ "Radiation only",
    
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
  select(ORIENAvatarKey : treatment_sequence, contains("ever_first"), everything())

write_rds(data,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Lung_data_",
                 today(), ".rds"))

write_csv(data,
          paste0(here::here(),
                 "/data/processed_data",
                 "/Lung_data_",
                 today(), ".csv"))

write_rds(data,
          paste0(path_raw,
                 "/ProcessedData",
                 "/Lung_CHinORIEN_",
                 today(), ".rds"))

write_csv(data,
          paste0(path_raw,
                 "/ProcessedData",
                 "/Lung_CHinORIEN_",
                 today(), ".csv"))
            
            
            
            
            
            
            
            
            
            
            
            
            