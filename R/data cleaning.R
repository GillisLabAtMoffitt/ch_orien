library(tidyverse)
library(data.table)

# I. Load data ----
load("~/Documents/GitHub/Gillis/ch_orien/Rmd/R01 Gillis-Teng 2023-2024/initial_files.RData")

path <- fs::path("", "Volumes", "Lab_Gillis","Data", "ORIEN")
mm_samples <- 
  readxl:::read_xlsx(paste0(path, "/ORIEN_conference/data_from_YiHan/ORIEN_MM_sample_list.xlsx"))
MSI_marker <-
  readxl::read_xlsx(paste0(here::here(), 
                           "/Rmd/R01 Gillis-Teng 2023-2024/MSI_23PRJ127MCC_20230620.xlsx"))

# II. Clean data ----
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
  mutate(Race = case_when(
    str_detect(Race, "American Indian")        ~ "American Indian or Alaska Native",
    Race == "White"                            ~ "White",
    Race == "Black"                            ~ "Black",
    str_detect(Race, "Asian") |
      str_detect(Race, "Cambodian") |
      Race == "Chinese" |
      Race == "Filipino" |
      Race == "Japanese" |
      Race == "Korean" |
      Race == "Pakistani" |
      Race == "Thai" |
      Race == "Vietnamese" |
      Race == "Laotian"                        ~ "Asian",
    Race == "Hawaiian" |
      Race == "Micronesian, NOS" |
      Race == "Pacific Islander, NOS" |
      Race == "Polynesian, NOS" |
      Race == "Samoan" |
      Race == "Tongan ese"                     ~ "Native Hawaiian or Other Pacific Islander",
    str_detect(Race, "Unknown")                ~ "Unknown",
    Race == "Other"                            ~ "Unknown",
    TRUE                                       ~ Race
  )) %>%
  mutate(Ethnicity = case_when(
    Ethnicity == "Spanish surname only"        ~ "Non-Hispanic",
    str_detect(Ethnicity, "Non-Hispanic")      ~ "Non-Hispanic",
    str_detect(Ethnicity, "Unknown")           ~ "Unknown",
    TRUE                                       ~ "Hispanic"
  )) %>%
  mutate(race_eth = case_when(
    Race == "White" &
      Ethnicity == "Non-Hispanic"              ~ "White Non-Hispanic",
    Race == "Black" &
      Ethnicity == "Non-Hispanic"              ~ "Black",
    Race == "Others" &
      Ethnicity == "Non-Hispanic"              ~ "Others Non-Hispanic",
    Ethnicity == "Hispanic"                    ~ "Hispanic"
  ))# %>% 
  # mutate(Race = factor(Race, levels = c("White", "Black", "Asian",
  #                                       "American Indian or Alaska Native", "Native Hawaiian or Other Pacific Islander",
  #                                       "Unknown")))

# Patient History----
PatientHistory <- PatientHistory %>%
  mutate(SmokingStatus = case_when(
    SmokingStatus == "Never"                   ~ "Never",
    SmokingStatus == "Current" |
      str_detect(SmokingStatus, "Ever")        ~ "Ever"
  ), SmokingStatus = factor(SmokingStatus, levels = c("Never", "Ever")))

# Diagnosis ----
Diagnosis_save <- Diagnosis
Diagnosis <- Diagnosis_save

Diagnosis %>%
  select(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode) %>%
  distinct() %>%
  arrange(PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite)

Diagnosis1 <- Diagnosis %>%
  select(AvatarKey, AgeAtDiagnosis, YearOfDiagnosis,
         PrimaryDiagnosisSiteCode : Histology,
         ClinGroupStage,
         CurrentlySeenForPrimaryOrRecurr,
         PerformStatusAtDiagnosis, 
         OtherStagingSystem, OtherStagingValue) %>%
  # join and filter by patients with MM samples
  inner_join(mm_samples %>% 
               rename("AvatarKey" = "ORIENAvatarKey"), .,
             by = c("AvatarKey")) %>% 
  mutate(not_real_age = case_when(
    AgeAtDiagnosis == "Age 90 or older"   ~ "Age 90 or older"
  )) %>%
  # mutate(Age.At.Specimen.Collection = as.numeric(Age.At.Specimen.Collection)) %>% 
  mutate(across(c("AgeAtDiagnosis", 
                  "Age.At.Specimen.Collection"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(ECOG = str_match(PerformStatusAtDiagnosis, "ECOG ([:digit:])")[,2]) %>%
  mutate(Karnofsky = str_match(PerformStatusAtDiagnosis, "Karnofsky ([:digit:].*)%")[,2]) %>%
  mutate(iss_stage = case_when(
    OtherStagingSystem == "ISS"           ~ OtherStagingValue
  )) %>% 
  
  arrange(AvatarKey, AgeAtDiagnosis) %>%
  group_by(AvatarKey) %>%
  mutate(how_many_dx = n(), .after = AvatarKey) %>%
  mutate(cancer_sequence = row_number(AvatarKey), .after = AvatarKey) %>%
  mutate(had_MM_at_anytime = case_when(
    Histology == "Multiple myeloma"       ~ "Yes"
  )) %>% 
  group_by(AvatarKey, Histology) %>% 
  mutate(MM_sequence  = case_when(
    Histology == "Multiple myeloma"       ~ as.character(row_number(Histology)),
    TRUE                                  ~ paste0(Histology, " (",PrimaryDiagnosisSite, ")")
    ), MM_sequence = case_when(
      MM_sequence == "1"                  ~ "1st MM dx",
      MM_sequence == "2"                  ~ "2nd MM dx",
      MM_sequence == "3"                  ~ "3rd MM dx",
      TRUE                                ~ MM_sequence
    )) %>% 
  mutate(age_at_first_MM_dx = case_when(
    MM_sequence == "1st MM dx"            ~ AgeAtDiagnosis
  )) %>% 
  mutate(age_at_last_MM_dx = case_when(
    MM_sequence == "3rd MM dx"            ~ AgeAtDiagnosis,
    MM_sequence == "2nd MM dx"            ~ AgeAtDiagnosis,
    MM_sequence == "1st MM dx"            ~ AgeAtDiagnosis
  )) %>% 
  mutate(age_at_MGUS_dx = case_when(
    Histology == "Monoclonal gammopathy"  ~ AgeAtDiagnosis
  )) %>% 
  group_by(AvatarKey) %>% 
  fill(had_MM_at_anytime, age_at_first_MM_dx,
       age_at_last_MM_dx, age_at_MGUS_dx, .direction = "updown") %>% 
  ungroup() %>% 
  # mutate(suid = AvatarKey) %>% 
  # mutate(sample_age = Age.At.Specimen.Collection) %>% 
  mutate(is_patient_MM = case_when(
    age_at_last_MM_dx >= age_at_MGUS_dx   ~ "Yes",
    age_at_last_MM_dx < age_at_MGUS_dx    ~ "MGUS patient",
    had_MM_at_anytime == "Yes"            ~ "Yes",
    is.na(had_MM_at_anytime)              ~ "No"
  )) %>% 
  mutate(sample_is_for_MM_dx = case_when(
    is_patient_MM == "MGUS patient"       ~ "No",
    is_patient_MM == "No"                 ~ "No",
    !is.na(age_at_MGUS_dx) &
      Age.At.Specimen.Collection > 
      age_at_MGUS_dx                      ~ "Yes",
    !is.na(age_at_MGUS_dx) &
      age_at_last_MM_dx == 
      age_at_MGUS_dx                      ~ "Yes",
    had_MM_at_anytime == "Yes"            ~ "Yes"
  ))
  
  # ungroup()
  # 
  # mutate(time_diff_dx_blood_days = case_when(
  #   Histology == "Multiple myeloma"       ~ round(
  #     (Age.At.Specimen.Collection - AgeAtDiagnosis) * 365.25,
  #     0)
  # ), .after = AvatarKey)
#   select("AvatarKey", "how_many_dx", "MM_sequence", cancer_sequence,
#          "time_diff_dx_blood_days",
#          AgeAtDiagnosis, Age.At.Specimen.Collection, Histology) %>% 
#   group_by(AvatarKey) %>% 
#   mutate(min_time = min(time_diff_dx_blood_days), .after = time_diff_dx_blood_days) %>% 
#   ungroup() %>% 
#   mutate(which_dx_is_blood_from = case_when(
#     how_many_dx == 1 &
#       time_diff_dx_blood_days == min_time            ~ MM_sequence,
#     
#     how_many_dx > 1 &
#       time_diff_dx_blood_days == min_time            ~ MM_sequence
#   ), .after = min_time)
#   
# a <- Diagnosis1 %>% filter(how_many_dx > 1)
# b <- Diagnosis1 %>% filter(how_many_dx == 1)


  # mutate(blood_is_in_range = case_when(
  #   Histology == "Multiple myeloma" &
  #     time_diff_dx_blood_days <= -7                     ~ "Yes"
  # ), .after = time_diff_dx_blood_days) 
  # # mutate(MM_first_dx_age = case_when(
  # #   Histology == "Multiple myeloma" &
  # #     MM_sequence == 1                  ~ AgeAtDiagnosis
  # # )) %>% 
  # # mutate(is_blood_at_MMdx = case_when(
  # #   Age.At.Specimen.Collection * 365.25 >= 
  # #     (MM_first_dx_age * 365.25 - 7)            ~ "Yes"
  # # )) %>% 
  # mutate(is_blood_at_MMdx = case_when(
  #   blood_is_after == "Yes" &
  #     MM_sequence == 1 &
  #     time_diff_dx_blood_days == min_time            ~ "Yes"
  # )) %>% 
  # mutate(sample_closest_to = case_when(
  #   time_diff_dx_blood_days == min_time              ~ MM_sequence
  # )) %>% 
  # fill(is_blood_at_MMdx, which_dx_is_blood_from,
  #      sample_closest_to,
  #      had_MM_at_anytime, MM_first_dx_age,
  #      iss_stage, .direction = "updown") %>%
  # ungroup() %>% 
  
  # select("AvatarKey","how_many_dx","MM_sequence",cancer_sequence,"MM_first_dx_age",
  #        "time_diff_dx_blood_days","blood_is_after","min_time","is_blood_at_MMdx",
  #        AgeAtDiagnosis, Age.At.Specimen.Collection, Histology)
  # select(AvatarKey, MM_sequence, everything())
  


subsequent_diagnosis <- Diagnosis %>%
  arrange(AvatarKey, AgeAtDiagnosis) %>% 
  group_by(AvatarKey) %>%
  summarise_at(vars(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode,
                    Histology,
                    AgeAtDiagnosis), str_c, collapse = "; ") %>%
  ungroup() %>%
  separate_wider_delim(cols = PrimaryDiagnosisSite, delim = "; ",
                       names = c("PrimaryDiagnosisSite", "subsequent_cancer_PrimaryDiagnosisSite"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  separate_wider_delim(cols = PrimaryDiagnosisSiteCode, delim = "; ",
                       names = c("PrimaryDiagnosisSiteCode", "subsequent_cancer_PrimaryDiagnosisSiteCode"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  separate_wider_delim(cols = Histology, delim = "; ",
                       names = c("Histology", "subsequent_cancer_Histology"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  separate_wider_delim(cols = AgeAtDiagnosis, delim = "; ",
                       names = c("AgeAtDiagnosis", "subsequent_cancer_AgeAtDiagnosis"),
                       too_few = "align_start", too_many = "merge",
                       cols_remove = TRUE) %>%
  select(-c(PrimaryDiagnosisSite, PrimaryDiagnosisSiteCode,
            Histology,
            AgeAtDiagnosis))

Diagnosis <- Diagnosis1 %>%
  arrange(AvatarKey, AgeAtDiagnosis) %>% 
  distinct(AvatarKey, .keep_all = TRUE) %>%
  left_join(., subsequent_diagnosis,
            by = c("AvatarKey")) #%>% 
  # mutate(is_MM_first_dx = case_when(
  #   MM_dx_age == AgeAtDiagnosis           ~ "Yes",
  #   TRUE                                  ~ "No"
  # ))

# VitalStatus----
VitalStatus <- VitalStatus %>%
  select(AvatarKey, VitalStatus, AgeAtLastContact, AgeAtDeath) %>%
  mutate(across(c("AgeAtLastContact", "AgeAtDeath"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

# Outcomes----
Outcomes1 <- Outcomes %>%
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "OutcomesPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "OutcomesPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtProgRecur) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Outcomes2 <- Outcomes %>%
  # For the patients with no diagnosis site
  filter(OutcomesPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtProgRecur) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Outcomes <- Outcomes1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Outcomes2) %>%
  select(AvatarKey, OutcomesPrimaryDiagnosisSiteCode, OutcomesPrimaryDiagnosisSite,
         ProgRecurInd, AgeAtProgRecur, RelapseStatus) %>%
  mutate(across(c("AgeAtProgRecur"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_outcomes_data = "Yes")

# MetastaticDisease----
MetastaticDisease1 <- MetastaticDisease %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "MetastaticDiseaseSiteCode" = "PrimaryDiagnosisSiteCode",
                    "MetastaticDiseaseSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE)

MetastaticDisease2 <- MetastaticDisease %>%
  # For the patients with no diagnosis site
  filter(MetastaticDiseaseSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtMetastaticSite) %>%
  distinct(AvatarKey, .keep_all = TRUE)

MetastaticDisease <- MetastaticDisease1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., MetastaticDisease2) %>%
  select(AvatarKey, MetastaticDiseaseSiteCode, MetastaticDiseaseSite,
         had_metastasis = MetastaticDiseaseInd, AgeAtMetastaticSite, MetsDzPrimaryDiagnosisSite) %>%
  mutate(across(c("AgeAtMetastaticSite"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_metastasis_data = "Yes")

# Medications----
Medications <- Medications %>%
  filter(MedicationInd != "(Migrated) Cannot determine from available documentation")

Medications1 <- Medications %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "MedPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "MedPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtMedStart)

Medications2 <- Medications %>%
  # For the patients with no diagnosis site
  filter(MedPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtMedStart)

Medications <- Medications1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Medications2) %>%
  group_by(AvatarKey) %>%
  mutate(n = dense_rank(MedPrimaryDiagnosisSiteCode), .after = MedPrimaryDiagnosisSiteCode) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(AvatarKey, MedPrimaryDiagnosisSiteCode, MedPrimaryDiagnosisSite,
         MedicationInd, AgeAtMedStart, Medication, AgeAtMedStop) %>%
  mutate(across(c("AgeAtMedStart"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

Medications1 <- Medications %>%
  arrange(AvatarKey, AgeAtMedStart, Medication, AgeAtMedStop) %>%
  mutate(received_imids = case_when(
    str_detect(Medication, "domide")       ~ "Yes",
    Medication %in% c("Revlimid",
                      "Pomalyst")         ~ "Yes"
  )) %>% 
  group_by(AvatarKey) %>% 
  fill(received_imids, .direction = "updown") %>% 
  group_by(AvatarKey, AgeAtMedStart, MedicationInd, received_imids) %>%
  summarise_at(vars(Medication, AgeAtMedStop), str_c, collapse = "; ") %>%
  ungroup()

library(data.table)
Medications <- dcast(setDT(Medications1),
                     AvatarKey + received_imids + MedicationInd ~ rowid(AvatarKey),
                     value.var = c("AgeAtMedStart", "Medication", "AgeAtMedStop")) %>%
  select(AvatarKey, received_imids, drugs_ever = MedicationInd, AgeAtMedStart_1 : AgeAtMedStart_10,
         Medication_1 : Medication_10, AgeAtMedStop_1 : AgeAtMedStop_10) %>%
  mutate(has_medication_data = "Yes")

# StemCellTransplant----
StemCellTransplant1 <- StemCellTransplant %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "SCTPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "SCTPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtTransplant) %>%
  distinct(AvatarKey, .keep_all = TRUE)

StemCellTransplant2 <- StemCellTransplant %>%
  # For the patients with no diagnosis site
  filter(SCTPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtTransplant) %>%
  distinct(AvatarKey, .keep_all = TRUE)

StemCellTransplant <- StemCellTransplant1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., StemCellTransplant2) %>%
  select(AvatarKey, sct_ever = SCTInd, AgeAtTransplant,
         TransplantType, TransplantCellSource) %>%
  mutate(across(c("AgeAtTransplant"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_stemcell_data = "Yes")

# Radiation----
Radiation <- Radiation %>%
  filter(RadiationTherapyInd != "(Migrated) Cannot determine from available documentation")

Radiation1 <- Radiation %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "RadPrimaryDiagnosisSiteCode" = "PrimaryDiagnosisSiteCode",
                    "RadPrimaryDiagnosisSite" = "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtRadiationStart) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Radiation2 <- Radiation %>%
  # For the patients with no diagnosis site
  filter(RadPrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtRadiationStart) %>%
  distinct(AvatarKey, .keep_all = TRUE)

Radiation <- Radiation1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., Radiation2) %>%
  select(AvatarKey, RadPrimaryDiagnosisSiteCode, RadPrimaryDiagnosisSite,
         radiation_ever = RadiationTherapyInd, AgeAtRadiationStart) %>%
  mutate(across(c("AgeAtRadiationStart"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(has_radiation_data = "Yes")

# SurgeryBiopsy----
SurgeryBiopsy <- SurgeryBiopsy %>%
  filter(SurgeryBiopsyInd != "(Migrated) Cannot determine from available documentation")

SurgeryBiopsy1 <- SurgeryBiopsy %>%
  # For the patient with a diagnosis site
  inner_join(., Diagnosis %>%
               select(AvatarKey, PrimaryDiagnosisSiteCode, PrimaryDiagnosisSite),
             by = c("AvatarKey", "PrimaryDiagnosisSiteCode",
                    "PrimaryDiagnosisSite")) %>%
  arrange(AvatarKey, AgeAtSurgeryBiopsy) %>%
  distinct(AvatarKey, .keep_all = TRUE)

SurgeryBiopsy2 <- SurgeryBiopsy %>%
  # For the patients with no diagnosis site
  filter(PrimaryDiagnosisSite == "Unknown/Not Applicable") %>%
  arrange(AvatarKey, AgeAtSurgeryBiopsy) %>%
  distinct(AvatarKey, .keep_all = TRUE)

SurgeryBiopsy <- SurgeryBiopsy1 %>%
  # bind the data with a primary site known as first diagnosis
  # then the one without which show some No Metastasis
  bind_rows(., SurgeryBiopsy2) %>%
  select(AvatarKey, SurgPrimaryDiagnosisSiteCode = PrimaryDiagnosisSiteCode,
         SurgPrimaryDiagnosisSite = PrimaryDiagnosisSite,
         surgery_ever = SiteTherapeutic, AgeAtSurgeryBiopsy) %>%
  mutate(across(c("AgeAtSurgeryBiopsy"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  mutate(surgery_ever = case_when(
    surgery_ever == "Unknown/Not Applicable" &
      !is.na(AgeAtSurgeryBiopsy)          ~ "Yes",
    surgery_ever == "Unknown/Not Applicable" &
      is.na(AgeAtSurgeryBiopsy)           ~ "No",
    TRUE                                  ~ surgery_ever
  )) %>%
  mutate(has_surgery_data = "Yes")

# PhysicalAssessment----
PhysicalAssessment <- PhysicalAssessment %>%
  left_join(., Diagnosis %>%
              select(AvatarKey, AgeAtDiagnosis),
            by = c("AvatarKey")) %>%
  mutate(across(c("AgeAtPhysicalExam"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(int = abs(AgeAtPhysicalExam - AgeAtDiagnosis)) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  select(AvatarKey, AgeAtPhysicalExam, BodyWeight, BodyHeight, BMI) %>%
  mutate(across(c("BodyWeight", "BodyHeight", "BMI"), ~ as.numeric(.)
  )) %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5                  ~ "Underweight",
    BMI >= 18.5 &
      BMI < 25                  ~ "Healthy",
    BMI >= 25 &
      BMI < 30                  ~ "Overweight",
    BMI >= 30                   ~ "Obese"
  )) %>%
  mutate(bmi_cat = factor(bmi_cat, levels = c("Underweight", "Healthy", "Overweight", "Obese")))

# Labs for ISS----
serum_albumin <- Labs %>%
  filter(LabTest == "Albumin (Serum)") %>% 
  mutate(across(c("AgeAtLabResults"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  left_join(., Diagnosis %>%
              select(AvatarKey, AgeAtDiagnosis),
            by = c("AvatarKey")) %>%
  mutate(int = abs(AgeAtLabResults - AgeAtDiagnosis)) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>% 
  mutate(LabResults = as.numeric(LabResults)) %>% 
  select(AvatarKey,
         albumin_value = LabResults, 
         albumin_unit = LabUnits)

iss_b2 <- TumorMarker %>%
  mutate(across(c("AgeAtTumorMarkerTest"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  filter(TMarkerTest == "Beta-2-microglobulin (B2M) level") %>% 
  select(AvatarKey, AgeAtTumorMarkerTest, 
         TMarkerTest, TMarkerResultValue, TMarkerValueUOM)

iss_b2_ <- Labs %>%
  mutate(across(c("AgeAtLabResults"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  filter(LabTest == "Beta 2 Microglobulin") %>% 
  select(AvatarKey, AgeAtTumorMarkerTest = AgeAtLabResults, 
         TMarkerTest = LabTest, 
         TMarkerResultValue = LabResults, 
         TMarkerValueUOM = LabUnits)

b2_globulin <- bind_rows(iss_b2, iss_b2_) %>% 
  arrange(AvatarKey, AgeAtTumorMarkerTest) %>% 
  left_join(., Diagnosis %>%
              select(AvatarKey, AgeAtDiagnosis),
            by = c("AvatarKey")) %>%
  mutate(int = abs(AgeAtTumorMarkerTest - AgeAtDiagnosis)) %>%
  arrange(AvatarKey, int) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>% 
  select(AvatarKey,
         b2_value = TMarkerResultValue, 
         b2_unit = TMarkerValueUOM) %>% 
  mutate(b2_value = as.numeric(b2_value)) %>% 
  mutate(b2_value = case_when(
    b2_unit == "mg/L"                     ~ b2_value,
    b2_unit == "mg/dL"                    ~ b2_value * 10,
    b2_unit == "ng/mL"                    ~ b2_value / 1000,
    b2_unit == "ug/mL (mcg/mL)"           ~ b2_value,
    b2_unit == "Unknown"                  ~ NA_real_
  )) %>% 
  mutate(b2_unit = "mg/L")

iss_stage_calculated <- 
  full_join(serum_albumin, b2_globulin, by = "AvatarKey") %>% 
  mutate(iss_stage_calculated = case_when(
    b2_value < 3.5 &
      albumin_value >= 3.5                ~ "I",
    b2_value < 3.5 &
      albumin_value < 3.5                 ~ "II",
    b2_value > 3.5 &
      b2_value <= 5.5                     ~ "II",
    b2_value > 5.5                        ~ "III"
  )) %>% 
  filter(!is.na(iss_stage_calculated))

Diagnosis <- Diagnosis %>% 
  left_join(., iss_stage_calculated, 
            by = "AvatarKey") %>% 
  mutate(iss_stage = coalesce(iss_stage, iss_stage_calculated))


# TumorMarker----
TumorMarker <- TumorMarker %>%
  mutate(across(c("AgeAtTumorMarkerTest"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  mutate(across(where(is.character), ~ na_if(., "Unknown/Not Applicable"))) %>%
  rename(marker_test = TMarkerTest, marker_unit = TMarkerValueUOM,
         marker_value = TMarkerResultValue) %>%
  group_by(marker_test) %>%
  fill(TMarkerLowRange, TMarkerHighRange, .direction = "updown") %>%
  ungroup() %>%
  mutate(across(c("marker_value", "TMarkerLowRange",
                  "TMarkerHighRange"), ~ as.numeric(.)
  )) %>%
  mutate(marker_interpretation = case_when(
    marker_value >= TMarkerLowRange &
      marker_value <= TMarkerHighRange     ~ "WNL (Within Normal Limits)",
    marker_value < TMarkerLowRange         ~ "Low",
    marker_value > TMarkerHighRange        ~ "High"
  ), marker_interpretation = coalesce(TMarkerRangeIndicator, marker_interpretation)) %>%
  mutate(marker_interpretation = case_when(
    TMarkerResult == "Value"              ~ marker_interpretation,
    TRUE                                  ~ TMarkerResult
  )) %>%
  mutate(marker_test = case_when(
    str_detect(marker_test, "HER2")       ~ "HER2",
    TRUE                                  ~ marker_test
  )) %>% 
  group_by(AvatarKey, marker_test) %>%
  mutate(er_result = case_when(
    marker_test == "Estrogen receptor (ER)"        ~ first(TMarkerResult)
  )) %>%
  mutate(pr_result = case_when(
    marker_test == "Progesterone receptor (PR)"    ~ first(TMarkerResult)
  )) %>%
  mutate(clean_her2 = case_when(
    str_detect(marker_test, "HER2") & 
      (TMarkerResult == "Positive" | 
         TMarkerResult == "Negative")                ~ TMarkerResult,
    str_detect(marker_test, "HER2") & 
      TMarkerResult == "Amplified"                 ~ "Positive"
  ), .after = marker_test) %>% 
  mutate(her_result = case_when(
    str_detect(marker_test, "HER2")                ~ first(clean_her2)
  ), .after = marker_test) %>% 
  mutate(ER_PR_HER2_test = case_when(
    marker_test == "Estrogen receptor (ER)" |
      marker_test == "Progesterone receptor (PR)" |
      marker_test == "HER2"                        ~ "ER_PR_HER2_status"
  )) %>% 
  mutate(er_receptor = case_when(
    er_result == "Positive"         ~ "ER+",
    er_result == "Negative"         ~ "ER-"
  )) %>% 
  mutate(pr_receptor = case_when(
    pr_result == "Positive"         ~ "PR+",
    pr_result == "Negative"         ~ "PR-"
  )) %>% 
  mutate(her_receptor = case_when(
    her_result == "Positive"         ~ "HER2+",
    her_result == "Negative"         ~ "HER2-"
  )) %>% 
  group_by(AvatarKey, ER_PR_HER2_test) %>% 
  fill(er_receptor, pr_receptor, her_receptor, .direction = "downup") %>% 
  ungroup() %>% 
  unite(ER_PR_HER2_status, c(er_receptor, pr_receptor, her_receptor), sep = "/", remove = FALSE, na.rm = TRUE) %>% 
  # TMarkerPercentStainResultValue is not helpful
  mutate(marker_test = coalesce(ER_PR_HER2_test, marker_test)) %>% 
  mutate(marker_interpretation = coalesce(ER_PR_HER2_status, marker_interpretation)) %>% 
  select(AvatarKey, AgeAtTumorMarkerTest, marker_test,
         marker_interpretation, marker_value, marker_unit)

TumorMarker <- TumorMarker %>%
  arrange(AvatarKey, marker_test, AgeAtTumorMarkerTest) %>%
  distinct(AvatarKey, marker_test, .keep_all = TRUE) %>%
  pivot_wider(id_cols = AvatarKey,
              names_from = marker_test,
              values_from = marker_interpretation)


# MSI ----
MSI_marker <- MSI_marker %>% 
  mutate(MSI_high_score = case_when(
    `MSIsensor2 Score` >= 20       ~ "Yes",
    `MSIsensor2 Score` < 20        ~ "No"
  )) %>% 
  select(ORIENAvatarKey, `WES SLID`, MSI_high_score)

# FamilyHistory----
FamilyHistory <- FamilyHistory %>%
  arrange(AvatarKey, desc(CancerInFamilyInd)) %>%
  distinct(AvatarKey, .keep_all = TRUE) %>%
  select(AvatarKey, CancerInFamilyInd)

# last_date----
last_date <- bind_rows(
  Imaging %>% select(AvatarKey, age_at_lab = AgeAtImageScan),
  Labs %>% select(AvatarKey, age_at_lab = AgeAtLabResults)
) %>%
  mutate(across(c("age_at_lab"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  ))) %>%
  arrange(AvatarKey, desc(age_at_lab)) %>%
  distinct(AvatarKey, .keep_all = TRUE)

# TumorSequencing----
TumorSequencing <- TumorSequencing %>%
  filter(TumorSequencingInd == "Yes") %>%
  select(AvatarKey, AgeAtTumorSequencing) %>%
  mutate(across(c("AgeAtTumorSequencing"), ~ case_when(
    . == "Age 90 or older"                ~ 90,
    . == "Unknown/Not Applicable"         ~ NA_real_,
    TRUE                                  ~ as.numeric(.)
  )))

Tumor_wes <-
  ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Tumor") %>%
  left_join(., MSI_marker,
            by = c("ORIENAvatarKey", "WES" = "WES SLID")) %>% 
  select(ORIENAvatarKey, Tumor_WES = WES,
         tumor_collection_age = Age.At.Specimen.Collection,
         Disease.Type,
         tumor_SpecimenSiteOfCollection = SpecimenSiteOfCollection,
         MSI_high_score)
Germline_wes <-
  ClinicalMolLinkage %>%
  filter(Tumor.Germline == "Germline") %>%
  select(ORIENAvatarKey, Germline_WES = WES,
         germline_collection_age = Age.At.Specimen.Collection,
         germline_SpecimenSiteOfCollection = SpecimenSiteOfCollection)

ClinicalMolLinkage <- full_join(Germline_wes, Tumor_wes,
                                by = "ORIENAvatarKey")

# CH_status <- CH_status %>%
#   left_join(., ClinicalMolLinkage,
#             by = c("AvatarKey" = "ORIENAvatarKey", "Germline_WES", "Tumor_WES")) %>%
#   mutate(across(c("germline_collection_age",
#                   "tumor_collection_age"), ~ case_when(
#                     . == "Age 90 or older"                ~ 90,
#                     . == "Unknown/Not Applicable"         ~ NA_real_,
#                     TRUE                                  ~ as.numeric(.)
#                   )))




# III. Join data----
data <- 
  left_join(Diagnosis, demographics, by = "AvatarKey") %>%
  left_join(., PatientHistory %>%
              select(AvatarKey, SmokingStatus),
            by = "AvatarKey") %>%
  left_join(., VitalStatus, by = "AvatarKey") %>%
  left_join(., CytogeneticAbnormalities, by = "AvatarKey") %>%
  left_join(., TumorMarker, by = "AvatarKey") %>%
  left_join(., PhysicalAssessment, by = "AvatarKey") %>%
  left_join(., Medications, by = "AvatarKey") %>%
  left_join(., StemCellTransplant, by = "AvatarKey") %>%
  left_join(., Radiation, by = "AvatarKey") %>%
  left_join(., SurgeryBiopsy, by = "AvatarKey") %>%
  left_join(., Outcomes, by = "AvatarKey") %>%
  left_join(., MetastaticDisease, by = "AvatarKey") %>%
  left_join(., last_date, by = "AvatarKey") %>%
  left_join(., FamilyHistory, by = "AvatarKey")

# IV. Create new variables----
data <- data %>% 
  # Age at last contact
  mutate(is_agelastcontact_last_date = case_when(
    AgeAtLastContact >= age_at_lab                  ~ "Yes",
    AgeAtLastContact < age_at_lab                   ~ "No"
  )) %>% 
  mutate(AgeAtLastContact = case_when(
    is_agelastcontact_last_date == "Yes"            ~ AgeAtLastContact,
    is_agelastcontact_last_date == "No"             ~ age_at_lab
  )) %>% 
  # Age at first treatment
  mutate(first_treatment = case_when(
    (AgeAtRadiationStart < AgeAtSurgeryBiopsy |
       (is.na(AgeAtSurgeryBiopsy) &
          !is.na(AgeAtRadiationStart))) &
      
      (AgeAtRadiationStart < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(AgeAtRadiationStart))) &
      
      (AgeAtRadiationStart < AgeAtTransplant | 
         (is.na(AgeAtTransplant) &
            !is.na(AgeAtRadiationStart)))              ~ "Radiation",
    
    (AgeAtSurgeryBiopsy < AgeAtRadiationStart | 
       (is.na(AgeAtRadiationStart) &
          !is.na(AgeAtSurgeryBiopsy))) &
      
      (AgeAtSurgeryBiopsy < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(AgeAtSurgeryBiopsy))) &
      
      (AgeAtSurgeryBiopsy < AgeAtTransplant | 
         (is.na(AgeAtTransplant) &
            !is.na(AgeAtSurgeryBiopsy)))          ~ "Surgery",
    
    (AgeAtMedStart_1 < AgeAtRadiationStart | 
       (is.na(AgeAtRadiationStart) &
          !is.na(AgeAtMedStart_1))) &
      
      (AgeAtMedStart_1 < AgeAtSurgeryBiopsy |
         (is.na(AgeAtSurgeryBiopsy) &
            !is.na(AgeAtMedStart_1))) &
      
      (AgeAtMedStart_1 < AgeAtTransplant | 
         (is.na(AgeAtTransplant) &
            !is.na(AgeAtMedStart_1)))             ~ "Drugs",
    
    (AgeAtTransplant < AgeAtRadiationStart | 
       (is.na(AgeAtRadiationStart) &
          !is.na(AgeAtTransplant))) &
      
      (AgeAtTransplant < AgeAtSurgeryBiopsy |
         (is.na(AgeAtSurgeryBiopsy) &
            !is.na(AgeAtTransplant))) &
      
      (AgeAtTransplant < AgeAtMedStart_1 |
         (is.na(AgeAtMedStart_1) &
            !is.na(AgeAtTransplant)))             ~ "SCT"
    
  )) %>% 
  mutate(age_at_first_treatment = case_when(
    first_treatment == "Radiation"                  ~ AgeAtRadiationStart,
    first_treatment == "Surgery"                    ~ AgeAtSurgeryBiopsy,
    first_treatment == "Drugs"                      ~ AgeAtMedStart_1,
    first_treatment == "SCT"                        ~ AgeAtTransplant,
  )) %>% 
  mutate(received_imids = case_when(
    received_imids == "Yes"                         ~ "Yes",
    TRUE                                            ~ "No"
  )) %>% 
  # OS
  mutate(os_event = case_when(
    VitalStatus == "Alive"                          ~ 0,
    VitalStatus == "Lost to follow-up"              ~ 0,
    VitalStatus == "Dead"                           ~ 1
  )) %>% 
  mutate(os_age = coalesce(AgeAtDeath, AgeAtLastContact)) %>% 
  mutate(os_time_from_dx_years = os_age - AgeAtDiagnosis) %>% 
  mutate(os_time_from_treatment_years = os_age - age_at_first_treatment) %>% 
  # PFS
  mutate(pfs_event = case_when(
    ProgRecurInd == "No"                            ~ 0,
    ProgRecurInd == "Progression"                   ~ 1,
    ProgRecurInd == "Recurrence"                    ~ 1
  )) %>% 
  mutate(pfs_age = coalesce(AgeAtProgRecur, AgeAtLastContact)) %>% 
  mutate(pfs_time_from_dx_years = pfs_age - AgeAtDiagnosis) %>% 
  mutate(pfs_time_from_treatment_years = pfs_age - age_at_first_treatment) %>% 
  # Metastasis
  mutate(pfs_event = case_when(
    had_metastasis == "Yes"                         ~ 1,
    is.na(had_metastasis)                           ~ 0
  )) %>% 
  mutate(met_age = coalesce(AgeAtMetastaticSite, AgeAtLastContact)) %>% 
  mutate(met_time_from_dx_years = met_age - AgeAtDiagnosis) %>% 
  mutate(met_time_from_treatment_years = met_age - age_at_first_treatment)

# V. Save----
write_csv(data %>% select(AvatarKey, 
                          Sex, Race, Ethnicity, race_eth, 
                          drugs_ever, received_imids,
                          surgery_ever,
                          iss_stage, ClinGroupStage,
                          SmokingStatus,
                          os_event, os_time_from_dx_years,
                          os_time_from_treatment_years,
                          pfs_event, pfs_time_from_dx_years,
                          pfs_time_from_treatment_years,
                          is_patient_MM), 
          paste0("MM_data_", today() ,".csv"))
