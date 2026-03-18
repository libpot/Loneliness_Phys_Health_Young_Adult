############################################################################
### CHRONIC LONELINESS AND PRESCRIPTIONS FOR PHYSICAL HEALTH CONDITIONS ####
############################################################################



###########################
### DATA PRE-PROCESSING ###
###########################

# loading packages
library(tidyverse) # data manipulation
library(haven) # loading SPSS dataset

# loading data
NorPD_UiO_agg_orig <- read_sav("disaggregated_300316.sav")




#######################
### YOUNG IN NORWAY ###
#### SURVEY/PANEL ####
######################

### SELECTING VARIABLES ###

NorPD_UiO_agg <- NorPD_UiO_agg_orig %>% 
  distinct_at(vars(id), .keep_all = T) %>%   
  select(id, # identification number
         
  # BASELINE COVARIATES
         # socio-demographic factors
         Gende1_1,# gender (self-reported and register-based)
         Age1_1, Age2_1,  # age and birth date
         ParEd1_1:ParEd1_9,  # parental education (self-reported and register-based)
         Origi1p1, Origi1p2, # ethnic background
         Urban1_1, # urbanity
         # early family environment 
         PBI1_06:PBI1_10, # parental warmth/care
         Subst1p1, Subst1p2, Subst1p3, Subst3p4, Subst3p5, # parental alcohol use 
         Smok3_02, Smok3_03, # parental smoking
         # chronic health conditions
         Disa1_03, Disa1_04, # asthma/allergy
         Disa1_05, # physical disability
         Disa2_09, # diabetes m     
         
# LONELINESS
         starts_with("UCLA") & !ends_with("5") & !ends_with(c("n", "m")), # UCLA loneliness scale
         starts_with("UCLA") & ends_with("5") & !ends_with(c("n", "m")), # direct loneliness measure
         
# TIME-VARYING COVARIATES
         # socio-demographic factors
         Occup2_1, Occ2_1a, Occ2_1b, Occ2_1c, Occ2_1d,
         ParJo1_1, ParJo1_2, ParJo2_1, ParJo2_2, Occup3_2, # employment status
         # lifestyle factors
         starts_with("Subst") & ends_with(c("a1", "a2")), # alcohol use (frequency and quantity, respectively)
         starts_with("Subst") & ends_with("_1"), # binge drinking
         starts_with("Subst") & ends_with("_2"), starts_with("Subst") & ends_with("_4"), # illicit drugs
         starts_with("Smoke") & ends_with("_1"), # smoking
         starts_with("Weigh") & ends_with("_1"), # weight
         starts_with("Heigh") & ends_with("_1"), # height
         starts_with("Train") & ends_with(c("_1", "_2", "_3")), # physical activity
         starts_with("LeAc") & ends_with(c("03", "04","06")), # physical activity II
         # social isolation indicators
         House1_1, House2_1, House3_1, # living situation
         Girlf1_1, Partn3_1, Partn3y2, Partn3y3, Marit3_1, # relationship status
         SS2b08, SS2c08, SS2d08, SS2e04, 
         SS3b08, SS3c08, SS3d08, SS3e04, 
         SS4b08, SS4c08, SS4d08, SS4e04, 
         SS1b08, SS1c08, SS1d08,  
         SPPA1sa2, SPPA2sa2, SPPA3sa2, # social network size
         # depressive symptoms
         starts_with("SCL") & ends_with(c("_07", "_08", "_09", "_10", "_11", "_12")) & !ends_with("tot"), 
         # perceived social support (in hypothetical situations)
         SS1b01, SS2b01, SS3b01,  # educational choice
         SS1c01, SS2c01, SS3c01, # personal problem
         SS1d01, SS2d01, SS3d01, # done something illegal
         SS2e01, SS3e01) # feeling down



### CREATING VARIABLES ###

NorPD_UiO_agg <- NorPD_UiO_agg %>% 
  # reverse coding and re-coding some items
  mutate(across(c(starts_with("UCLA") & ends_with(c("1", "2")), PBI1_07, PBI1_08), # higher values indicate over-controlling and cold parenting
                ~ 5 - .x),
         Heigh2_1 = ifelse(Heigh2_1 < 100, Heigh2_1 + 100, Heigh2_1), # coding errors: some heights were missing 1 at the beginning (e.g., 68 instead of 168)
         Heigh3_1 = ifelse(Heigh3_1 < 130, (Heigh2_1+Heigh4_1)/2, Heigh3_1),
         across(c(SS1b01, SS1c01, SS1d01),
                ~ ifelse(.x %in% c(2,3), 0, .x)),
         Occup3_2 = case_when(Occup3_2 %in% c(1,2,3) ~ "employed",
                              Occup3_2 == 4 ~ "unemployed",
                              Occup3_2 %in% c(5,6,7) ~ "inactive"),
         Girlf2_1 = case_when(Partn3_1 == 1 | Partn3y2 > 94 ~ "single",
                              Partn3y2 <= 94 & Partn3y3 > 94 | Partn3y2 <= 94 & is.na(Partn3y3) ~ "in relationship",
                              TRUE ~ NA),
         Girlf2_supp = rowSums(select(., SS2b08, SS2c08, SS2d08, SS2e04), na.rm = F),
         Girlf1_supp = rowSums(select(., SS1b08, SS1c08, SS1d08), na.rm = F),
         Girlf3_supp = rowSums(select(., SS3b08, SS3c08, SS3d08, SS3e04), na.rm = F),
         across(c(SPPA1sa2, SPPA2sa2, SPPA3sa2),
                ~ case_when(.x %in% c(1,2) ~ "little",
                            .x %in% c(3,4) ~ "a lot"))) %>% 
# BASELINE COVARIATES
         # sociodemographic factors
  mutate(age = ifelse(is.na(Age1_1), Age2_1 - 2, Age1_1),
         age = ifelse(age > 25, Age2_1 - 2, age), # re-coding non-sensical values
         gender = Gende1_1,
         gender = case_when(gender == 1 ~ "female",
                            gender == 0 ~ "male",
                            TRUE ~ NA),
         ethnicity = case_when(Origi1p1 == 1 | Origi1p2 == 1 ~ "Norway-born",
                               Origi1p1 == 2 & Origi1p2 == 2 ~ "migrant"),
         urbanity = case_when(Urban1_1 %in% c(1,2) ~ "rural", 
                              Urban1_1 %in% c(3,4) ~ "urban", 
                              TRUE ~ NA),
         parental_education = case_when(ParEd1_7 %in% c(1,2,3) | ParEd1_6 %in% c(1,2,3) ~ "college/university",
                                        ParEd1_3 %in% c(1,2,3) | ParEd1_4 %in% c(1,2,3) | ParEd1_5 %in% c(1,2,3) | ParEd1_1 %in% c(1,2,3) | ParEd1_2 %in% c(1,2,3) ~ "lower",
                                        TRUE ~ NA),
         
         # early family environment 
         warm_parenting = rowMeans(select(., PBI1_06:PBI1_10), na.rm = T),
         parental_alcoholuse = case_when(Subst1p1 %in% c(4,5) | Subst1p2 %in% c(4,5) | Subst1p3 %in% c(4,5) | Subst3p4 %in% c(4,5) | Subst3p5 %in% c(4,5) ~ "heavy use",
                                         Subst1p1 == 1 & Subst1p2 == 1 & Subst1p3 == 1 ~ "abstinence",
                                         Subst1p1 %in% c(2,3) | Subst1p2 %in% c(2,3) | Subst1p3 %in% c(1,2,3)  ~ "occasional/low use",
                                         TRUE ~ NA),
         parental_smoking = case_when(Smok3_02 == 4 | Smok3_03 == 4 ~ "daily smoking",
                                      Smok3_02 %in% c(1,2,3) | Smok3_03 %in% c(1,2,3) ~ "non-smoking",
                                      TRUE ~ NA),
         # chronic conditions
         asthma_allergy = case_when(Disa1_03 == 1 | Disa1_04 == 1 ~ 1, 
                                    is.na(Disa1_03) & is.na(Disa1_04) ~ NA,
                                    TRUE ~ 0),
         diabetes = case_when(Disa2_09 == 1 ~ 1, 
                              is.na(Disa2_09) ~ NA,
                              TRUE ~ 0),
         phys_disability = case_when(Disa1_05 == 1 ~ 1,
                                     is.na(Disa1_05) ~ NA,
                                     TRUE ~ 0)) %>% 
  
  # TIME-VARYING COVARIATES
  # eomployment status
  mutate(work_status_1f = case_when(ParJo1_1 %in% c(1,2) ~ "employed",
                                    ParJo1_1 %in% c(3,4,5,6) ~ "unemployed",
                                    TRUE ~ NA),
         work_status_1m = case_when(ParJo1_2 %in% c(1,2) ~ "employed",
                                    ParJo1_2 %in% c(3,4,5,6) ~ "unemployed",
                                    TRUE ~ NA),
         employment_1 = case_when(work_status_1f == "employed"| work_status_1m == "employed" ~ "employed",
                                  is.na(work_status_1f) ~ work_status_1m,
                                  is.na(work_status_1m) ~ work_status_1f,
                                  work_status_1f == "unemployed" & work_status_1m == "unemployed" ~ "unemployed"),
         work_status_2f = case_when(ParJo2_1 %in% c(1,2) ~ "employed",
                                    ParJo2_1 %in% c(3,4,5,6) ~ "unemployed",
                                    TRUE ~ NA),
         work_status_2m = case_when(ParJo2_2 %in% c(1,2) ~ "employed",
                                    ParJo2_2 %in% c(3,4,5,6) ~ "unemployed",
                                    TRUE ~ NA),
         employment_2 = case_when(work_status_2f == "employed"| work_status_2m == "employed" ~ "employed",
                                  is.na(work_status_2f) ~ work_status_2m,
                                  is.na(work_status_2m) ~ work_status_2f,
                                  work_status_2f == "unemployed" & work_status_2m == "unemployed" ~ "unemployed"),
         employment_3 = Occup3_2,
         # social isolation indicators
         across(c(House1_1,House2_1), ~ case_when(.x %in% c(2,3,4,5,6,7,8) ~ "not with both parents",
                                                  .x == 1 ~ "with both parents")),
         House3_1 = case_when(House3_1 == 2 ~ "alone",
                              House3_1 %in% c(1,3,4,5) ~ "not alone"),
         living_situation_1 = House1_1,
         living_situation_2 = House2_1,
         living_situation_3 = House3_1,
         across(c(Girlf1_supp, Girlf2_supp, Girlf3_supp),
                ~ case_when(.x >= 1 ~ "in relationship",
                            .x == 0 ~ "single")),
         Marit3_1 = case_when(Marit3_1 == 1 ~ "single",
                              Marit3_1 %in% c(2,3) ~ "in relationship"),
         relationship_1 = case_when(Girlf1_1 == 1 ~ "in relationship",
                                    Girlf1_1 %in% c(2,3) ~ "single"),
         relationship_1 = ifelse(is.na(relationship_1), Girlf1_supp, relationship_1),
         relationship_2 = ifelse(is.na(Girlf2_supp), Girlf2_1, Girlf2_supp),
         relationship_3 = ifelse(is.na(Girlf3_supp), Marit3_1, Girlf3_supp),
         friends_1 = SPPA1sa2,
         friends_2 = SPPA2sa2,
         friends_3 = SPPA3sa2) %>% 
  # social support
  mutate(social_support_1 = rowSums(select(., SS1b01,SS1c01,SS1d01), na.rm = F),
         social_support_2 = rowSums(select(., SS2b01,SS2c01,SS2d01,SS2e01), na.rm = F),
         social_support_3 = rowSums(select(., SS3b01,SS3c01,SS3d01,SS3e01), na.rm = F),
         #across(c(social_support_1:social_support_3),
         #      ~ ifelse(.x == max(.x, na.rm = T), .x - 1, .x)), # too few cases with extremely low social support, merging the two lowest categories
         # depressive symptoms
         depression_1 = rowMeans(select(., SCL1_07:SCL1_12), na.rm = T),
         depression_2 = rowMeans(select(., SCL2_07:SCL2_12), na.rm = T),
         depression_3 = rowMeans(select(., SCL3_07:SCL3_12), na.rm = T)) %>% 
  
  # LONELINESS
  mutate(loneliness_1 = rowMeans(select(., UCLA1_1:UCLA1_4), na.rm = T),
         loneliness_2 = rowMeans(select(., UCLA2_1:UCLA2_4), na.rm = T),
         loneliness_3 = rowMeans(select(., UCLA3_1:UCLA3_4), na.rm = T),
         loneliness_4 = rowMeans(select(., UCLA4_1:UCLA4_4), na.rm = T),
         loneli_emo_1 = rowMeans(select(., UCLA1_3:UCLA1_4), na.rm = T),
         loneli_emo_2 = rowMeans(select(., UCLA2_3:UCLA2_4), na.rm = T),
         loneli_emo_3 = rowMeans(select(., UCLA3_3:UCLA3_4), na.rm = T),
         loneli_emo_4 = rowMeans(select(., UCLA4_3:UCLA4_4), na.rm = T),
         loneli_soc_1 = rowMeans(select(., UCLA1_1:UCLA1_2), na.rm = T),
         loneli_soc_2 = rowMeans(select(., UCLA2_1:UCLA2_2), na.rm = T),
         loneli_soc_3 = rowMeans(select(., UCLA3_1:UCLA3_2), na.rm = T),
         loneli_soc_4 = rowMeans(select(., UCLA4_1:UCLA4_2), na.rm = T),
         loneli_direct_1 = UCLA1_5,
         loneli_direct_2 = UCLA2_5,
         loneli_direct_3 = UCLA3_5,
         loneli_direct_4 = UCLA4_5) %>% 
  
  # HEALTH BEHAVIORS
  mutate(across(starts_with("Smoke") & ends_with("_1"),
                ~ case_when(.x %in% c(1,2,3) ~ "non-smoking",
                            .x %in% c(4,5) ~ "smoking",
                            TRUE ~ NA)),
         Train3_2 = ifelse(is.na(Train3_2), 0, Train3_2),
         Subst1a2 = ifelse(Subst1a1 == 0, 0, Subst1a2),
         Subst2a2 = ifelse(Subst2a1 == 0, 0, Subst2a2),
         Subst3a2 = ifelse(Subst3a1 == 0, 0, Subst3a2),
         illicit_drug_1 = case_when(Subst1_2 %in% c(2,3,4,5,6) | Subst1_4 %in% c(2,3,4,5,6) ~ 1,
                                    Subst1_2 == 1 &  Subst1_4 == 1 ~ 0,
                                    TRUE ~ NA),
         illicit_drug_2 = case_when(Subst2_2 %in% c(2,3,4,5,6) | Subst2_4 %in% c(2,3,4,5,6) ~ 1,
                                    Subst2_2 == 1 &  Subst2_4 == 1 ~ 0,
                                    TRUE ~ NA),
         illicit_drug_3 = case_when(Subst3_2 %in% c(2,3,4,5,6) | Subst3_4 %in% c(2,3,4,5,6) ~ 1,
                                    Subst3_2 == 1 &  Subst3_4 == 1 ~ 0,
                                    TRUE ~ NA),
         bmi_1 = Weigh1_1/((Heigh1_1/100)*(Heigh1_1/100)),
         bmi_2 = Weigh2_1/((Heigh2_1/100)*(Heigh2_1/100)),
         bmi_3 = Weigh3_1/((Heigh3_1/100)*(Heigh3_1/100))) %>% 
  mutate(phys_exer_1 = (LeAc1_03 + LeAc1_04 + LeAc1_06)*60,
         phys_exer_2 = (LeAc2_03 + LeAc2_04 + LeAc2_06)*60,
         phys_exer_3 = case_when(Train3_3 == 1 ~ 0,
                                 Train3_3 == 0 ~ Train3_1*60 + Train3_2,
                                 TRUE ~ NA),
         smoking_1 = Smoke1_1,
         smoking_2 = Smoke2_1,  
         smoking_3 = Smoke3_1,
         alcohol_use_1 = Subst1a1*Subst1a2,
         alcohol_use_2 = Subst2a1*Subst2a2,
         alcohol_use_3 = Subst3a1*Subst3a2) %>%
  
  # FIDELITY CHECK: replacing non-nonsensical values with NA
  mutate(across(starts_with("alcohol_use_"),
                ~ ifelse(.x > 300, NA, .x)),  # more than 300 drinks per 4 weeks (e.g., 30 times x 10 drinks) = implausible and coded as missing 
         across(starts_with("phys_exer"),
                ~ ifelse(.x > 2400, NA, .x))) %>%  # more than 2400 min (40h) of exercise per week = implausible and coded as missing 
  # coding variable type: binary/ordinal/numeric
  mutate(across(c(age, warm_parenting,
                  starts_with(c("social_support_", "depression_")),
                  starts_with(c("alcohol_use_", "bingedrink_", "phys_exer_", "bmi_")),
                  starts_with(c("loneliness", "loneli_direct", "UCLA"))),
                ~ as.numeric(.x)),
         across(c(gender, ethnicity, parental_education, urbanity,
                  parental_smoking, parental_alcoholuse,
                  asthma_allergy, diabetes, phys_disability,
                  starts_with(c("living_situation_", "relationship_", "friends_", "employment")),
                  starts_with(c("smoking", "illicit_drug"))),
                ~ as.factor(.x))) %>% 
  # standardizing (z-transforming) loneliness scores
  mutate(across(starts_with("loneliness_"),
                ~ (.x - mean(.x, na.rm = T))/sd(.x, na.rm = T))) %>% 
  
  # CREATING FINAL DATASET: selecting final variables
  select(id, 
         age, gender, ethnicity, parental_education, urbanity, # baseline socio-demographics
         warm_parenting, parental_alcoholuse, parental_smoking, # early family environment 
         asthma_allergy, diabetes, phys_disability, # chronic health conditions
         starts_with(c("social_support_", "depression_")), # social support and depression
         starts_with(c("living_situation_", "relationship_", "friends_", "employment")),  # socio-demographics
         starts_with(c("smoking_", "alcohol_use_", "illicit_drug", "phys_exer_", "bmi_")), # health factors
         starts_with("loneli")) # loneliness






######################################
#### REGISTER-BASED PRESCRIPTIONS ####
### FOR PHYSICAL HEALTH CONDITIONS ###
######################################

NorPD_tidy <- NorPD_UiO_agg_orig %>% 
  select(id, UtleveringsDato, ATCKode, OrdinasjonAntallDDD) %>% 
  mutate(across(c(UtleveringsDato, ATCKode), 
         ~ ifelse(.x == "", NA, .x))) %>% 
  transmute(id = id,
            date = as_date(UtleveringsDato),
            atc = ATCKode,
            ddd = as.numeric(gsub(",", ".", OrdinasjonAntallDDD)))

# 60 did not receive any prescription at all 
length(unique(filter(NorPD_tidy, if_any(date:atc, is.na))$id))

# 2,400 individuals (out of 2,062, 92.%) were prescribed at least one examined medication during the 10-year follow-up
length(unique(filter(NorPD_tidy, date >= '2006-09-1', str_detect(atc, all_atcs))$id))

NorPD_tidy %>% 
filter(str_detect(atc, all_atcs)) %>% 
       filter(if_any(date:ddd, is.na))
# 4 individuals with missing DDD, but no missing in date or atc codes

# counting overall number of prescriptions and prescriptions' DDD for analgesics and antibiotics
NorPD_acute_outcomes <- NorPD_tidy %>% 
  arrange(id, date) %>%               
  group_by(id) %>%    
  summarise(analgesics_nr = length(id[grepl("^M01A", atc) & date >= '2006-09-1']),
            analgesics_DDD = sum(ddd[grepl("^M01A", atc) & date >= '2006-09-1'], na.rm = T),
            antibiotics_nr = length(id[grepl("^J01", atc) & date >= '2006-09-1']),
            antibiotics_DDD = sum(ddd[grepl("^J01", atc) & date >= '2006-09-1'], na.rm = T)) %>% 
  ungroup() 

# creating a time-to-event dataset for chronic health conditions
NorPD_chronic_outcomes <- NorPD_tidy %>%
  filter(date >= '2006-09-1') %>% 
  mutate(med_out = as.factor(case_when(grepl("^C01|^C02|^C03|^C07|^C08|^C09|^C10", atc) ~ "C",
                                       grepl("^H03A", atc) ~ "H03A", 
                                       grepl("^A07A|^A07E", atc) ~ "A07", 
                                       grepl("^N05B", atc) ~ "N05B",
                                       grepl("^N06A", atc) ~ "N06A",
                                       grepl("^H02", atc) ~ "H02", 
                                       grepl("^L04A", atc) ~ "L04A",
                                       grepl("^N02CC", atc) ~ "N02CC"))) %>%
  filter(!is.na(med_out)) %>% 
  group_by(id, med_out) %>% 
  arrange(id, date, med_out) %>%    
  mutate(two_prescr_12m = as.numeric((date - lag(date)) <= 365), 
         ddd_100_12m = map_dbl(date, ~ as.numeric(sum(ddd[date >= .x - months(12) & date <= .x], na.rm = TRUE) >= 100)),
         event = as.numeric(two_prescr_12m | ddd_100_12m),
         time = ifelse(event == 1, as.numeric(date - as_date('2006-09-1'))/365.25, NA)) %>% 
  ungroup() %>% 
  mutate(event = ifelse(is.na(event), 0, event),
         time = ifelse(is.na(time), 9.33, time)) %>%  # 9.33 is the length of the follow-up period in years
  group_by(id, med_out) %>% 
  arrange(id, date, med_out) %>%    
  summarise(event = as.numeric(sum(event, na.rm = T) >= 1),
            time  = min(time, na.rm = T)) %>% 
  ungroup() %>% 
  filter(event == 1) %>% 
  pivot_wider(names_from = med_out, 
              values_from = c(event, time))

NorPD_antidep_washout <- NorPD_tidy %>%
  filter(grepl("^N06A", atc)) %>%
  arrange(id, date) %>%  
  group_by(id) %>% 
  mutate(pre_exp_med = ifelse(date <= '2004-12-31', TRUE, FALSE),
         no_preexp_med = all(pre_exp_med == FALSE)) %>%  # no pre-exposure antidepressants (1-year washout period) 
  filter(date >= '2006-09-1') %>% # outcome: post-exposure antidepressants
  mutate(two_prescr_12m = as.numeric((date - lag(date)) <= 365), 
         ddd_100_12m = map_dbl(date, ~ as.numeric(sum(ddd[date >= .x - months(12) & date <= .x], na.rm = TRUE) >= 100)),
         event = as.numeric(two_prescr_12m | ddd_100_12m),
         time = ifelse(event == 1, as.numeric(date - as_date('2006-09-1'))/365.25, NA)) %>% 
  ungroup() %>% 
  mutate(event = ifelse(is.na(event), 0, event),
         time = ifelse(is.na(time), 9.33, time)) %>%  # 9.33 is the length of the follow-up period in years
  group_by(id, no_preexp_med) %>% 
  arrange(id, date) %>%    
  summarise(event = as.numeric(sum(event, na.rm = T) >= 1),
            time  = min(time, na.rm = T)) %>% 
  ungroup() %>% 
  filter(event == 1)
write.csv(NorPD_antidep_washout, "antidep_washout")


### MERGING SURVEY AND REGISTER DATA ###

NorPD_UiN_analysis <- left_join(NorPD_UiO_agg, NorPD_chronic_outcomes, by = "id") %>% 
  mutate(across(starts_with("event_"), 
                ~ ifelse(is.na(.x), 0, .x)),
         across(starts_with("time_"), 
                ~ ifelse(is.na(.x), 9.33, .x))) %>%  #  9.33 is the length of the follow-up period in years
left_join(., NorPD_acute_outcomes, by = "id") 

NorPD_meds <- left_join(NorPD_UiO_agg[,"id"], NorPD_chronic_outcomes, by = "id") %>% 
  mutate(across(starts_with("event_"), 
                ~ ifelse(is.na(.x), 0, .x)),
         across(starts_with("time_"), 
                ~ ifelse(is.na(.x), 9.33, .x))) %>%  #  9.33 is the length of the follow-up period in years
  left_join(., NorPD_acute_outcomes, by = "id")   
  
write.csv(NorPD_UiN_analysis, "UiN_analysis")




###################################
### PRESCRIPTIONS: DESCRIPTIVES ###
###################################

atc_chronic <- c("H03A", "C01|^C02|^C03|^C07|^C08|^C09|^C10", "A07A|^A07E", "H02", "L04A", "N02CC")
all_atcs <- "^M01A|^J01|^H03A|^C01|^C02|^C03|^C07|^C08|^C09|^C10|^A07A|^A07E|^H02|^L04A|^N02CC"

# function for counting number of individuals with at least 2 prescriptions within 12 consecutive months or prescription with 100+ DDD covering three months
prescriptions_within_12_months <- function(diag) {
  sum <- NorPD_tidy %>%
    filter(date >= '2006-09-1',
           str_detect(atc, paste0("^", diag))) %>%
    arrange(id, date) %>%               
    group_by(id) %>% 
    mutate(ddd_100_12m = map_lgl(date, ~ sum(ddd[date >= .x - months(12) & date <= .x], na.rm = TRUE) >= 100)) %>%
    summarise(two_prescr = sum(diff(date) <= 365) > 0,
              sum_100ddd = any(ddd_100_12m),
              any = two_prescr == TRUE | sum_100ddd == TRUE) %>%  
    ungroup() %>% 
    summarise(two_prescr = sum(two_prescr),
              sum_100ddd = sum(sum_100ddd),
              sum_any = sum(any)) 
}

# function for counting median prescription DDD in individuals with at least 2 prescriptions within 12 consecutive months
DDD_mean <- function(diag) {
  sum <- NorPD_tidy %>%
    filter(date >= '2006-09-1', 
           str_detect(atc, paste0("^", diag))) %>%    
    arrange(id, date) %>%               
    group_by(id) %>%                       
    mutate(cond = sum(diff(date) <= 365) > 0) %>%  
    filter(cond == TRUE) %>% 
    ungroup() %>%            
    summarise(median_iqr = paste0(round(median(ddd, na.rm = T)), " (", round(quantile(ddd, 0.25, na.rm = T)), " to ", round(quantile(ddd, 0.75, na.rm = T)), ")"))
}

# function for counting median number of prescriptions in individuals with outcome
prescr_mean <- function(diag) {
  sum <- NorPD_tidy %>%
    filter(date >= '2006-09-1', 
           str_detect(atc, paste0("^", diag))) %>%    
    arrange(id, date) %>%               
    group_by(id) %>%                       
    mutate(ddd_100_12m = map_lgl(date, ~ sum(ddd[date >= .x - months(12) & date <= .x], na.rm = TRUE) >= 100),
           two_prescr = sum(diff(date) <= 365) > 0,
           sum_100ddd = any(ddd_100_12m),
           any = two_prescr == TRUE | sum_100ddd == TRUE) %>%  
    filter(any == TRUE) %>% 
    summarise(ind_prescr = sum(any)) %>% 
    ungroup() %>% 
    summarise(median_iqr = paste0(round(median(ind_prescr, na.rm = T)), " (", round(quantile(ind_prescr, 0.25, na.rm = T)), " to ", round(quantile(ind_prescr, 0.75, na.rm = T)), ")"))
}

# counting overall number of prescriptions and prescriptions' DDD
NorPD_acute_stats <- NorPD_acute_outcomes %>% 
  summarise(across(analgesics_prescr_nr:antibiotics_DDD,
                   ~ paste0(median(.x, na.rm = T), " (", min(.x, na.rm = T), " to ", round(max(.x, na.rm = T)), ")")))

# applying the function for each ATC code
chronic_cond <- data.frame(condition = atc_chronic,
                           sum = map_dfr(atc_chronic, prescriptions_within_12_months),
                           mean_DDD = map_dfr(atc_chronic, DDD_mean),
                           mean_prescr = map_dfr(atc_chronic, prescr_mean)) %>% 
  mutate(condition = case_when(condition == "C01|^C02|^C03|^C07|^C08|^C09|^C10" ~ "Cx",
                               condition == "A07A|^A07E" ~ "A07x",
                               TRUE ~ condition))
write.csv(chronic_cond, "descriptives_chronic_cond")

  
  

##############################
### SAMPLE CHARACTERISTICS ###
##############################

sample_characteristics_cat <- NorPD_UiO_agg %>% 
  select(gender, ethnicity, parental_education, urbanity, employment_1,
         parental_alcoholuse, parental_smoking,
         asthma_allergy, diabetes, phys_disability,
         living_situation_1, relationship_1, friends_1, 
         smoking_1, illicit_drug_1) %>% 
  mutate(across(everything(), as.factor)) %>% 
  pivot_longer(everything()) %>% 
  count(name, value) %>% 
  group_by(name) %>% 
  mutate(prop = formatC(round(n/sum(n)*100, 2), format = "f", digits = 2),
         stat = paste(n, paste0("(", prop, ")"))) %>% 
  filter(!is.na(value)) %>% 
  select(name, value, stat) 

sample_characteristics_conti <- NorPD_UiO_agg %>% 
  select(age,  
         warm_parenting,
         depression_1, social_support_1,
         alcohol_use_1, phys_exer_1, bmi_1,
         loneliness_1) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  summarise(stat = paste0(formatC(round(mean(value, na.rm = TRUE), 2), format = "f", digits = 2), " (", 
                          formatC(round(sd(value, na.rm = TRUE), 2), format = "f", digits = 2), ")"),
            median = formatC(round(median(value, na.rm = TRUE), 2), format = "f", digits = 2),
            min_max = paste0(formatC(round(min(value, na.rm = TRUE), 0), format = "f", digits = 0), " to ", 
                             formatC(round(max(value, na.rm = TRUE), 0), format = "f", digits = 0)))

sample_characteristics <- rbind(sample_characteristics_cat, sample_characteristics_conti) %>%
  mutate(name = factor(name, levels = c("gender", "age", "ethnicity", "urbanity", "parental_education", "employment_1",
                                        "warm_parenting", "parental_alcoholuse", "parental_smoking",
                                        "asthma_allergy", "diabetes", "phys_disability",
                                        "living_situation_1", "relationship_1", "friends_1", 
                                        "alcohol_use_1", "smoking_1", "illicit_drug_1", "phys_exer_1", "bmi_1",
                                        "depression_1", "social_support_1", 
                                        "loneliness_1"), ordered = TRUE)) %>% 
  arrange(name) %>% 
  select(-median) 
write.csv(sample_characteristics, "sample_characteristics")




################################
### MISSING DATA EXPLORATION ###
################################

missing_proportions_long <- NorPD_UiO_agg %>% 
  select(-id) %>% 
  summarise(across(everything(), ~ round(mean(is.na(.))*100))) %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "value")  

missing_proportions_long %>% arrange(desc(value))

missing_proportions_table <- missing_proportions_long %>% 
  mutate(variable = ifelse(grepl("\\d$", variable), variable, paste0(variable, "_1"))) %>%
  separate(variable, into = c("name", "wave"), sep = "_(?=\\d+$)") %>% 
  pivot_wider(names_from = "wave",
              values_from = "value") 
write.csv(missing_proportions_table, "missing_props")


missing_proportions_long %>% 
  summarise(mean_NA_proportion = mean(value)) 
# on average, 10% of values are missing -> number of recommended imputed datasets will be 10

nrow(na.omit(UiN_data_analysis))





###########################
### MULTIPLE IMPUTATION ###
###########################

library(mice)
# dry run for specifying predictors and imputation methods
init_imp <- mice(NorPD_UiO_agg, maxit = 0)

# predictor matrix
predictor_matrix <- init_imp$predictorMatrix
predictor_matrix[,"id"] <- 0
predictor_matrix["id",] <- 0
# we are including all variables in imputation models

# imputation methods
imp_method <- init_imp$method
# predictive mean matching (pmm) for continuous data,
# logistic regression (logreg) for binary, and
# polytomous regression (polyreg) for un-ordered categorical, 

# MULTIPLE IMPUTATION
set.seed(12345)
Imp1 <- mice(NorPD_UiO_agg, 
             m = 10, # number of imputed datasets
             maxit = 5, # number of iterations
             method = imp_method, 
             predictorMatrix = predictor_matrix)

# convergence diagnostics
plot(Imp1) 

NorPD_UiN_analysis_imp <- complete(Imp1, action = "long", include = FALSE) %>% 
  left_join(., NorPD_meds, by = "id") 
  
# saving data: long-format
write.csv(NorPD_UiN_analysis_imp, "UiN_analysis_imp")




# SCALE RELIABILITY

outcome_4i <- ' outcome =~ item_1 + item_2 + item_3 + item_4 '

loneli_reliability <- NorPD_UiO_agg %>% 
  select(id, starts_with("UCLA") & !ends_with("5") & !ends_with(c("n", "m"))) %>%  # UCLA loneliness scale
  pivot_longer(cols = -id,
               names_to = c("time_wave", "item"),
               names_sep = "_",
               values_to = "symptoms") %>% 
  pivot_wider(names_from = "item",
              values_from = "symptoms") %>% 
  filter(time_wave %in% c("UCLA1", "UCLA2", "UCLA3", "UCLA4")) %>% 
  mutate(item_1 = `1`,
         item_2 = `2`,
         item_3 = `3`,
         item_4 = `4`) %>% 
  group_by(time_wave) %>% 
  nest() %>% 
  mutate(model = map(.x = data,
                     ~ lavaan::cfa(model = outcome_4i,
                                   data = .x,
                                   estimator = "MLR", # estimator: robust maximum likelihood
                                   missing = "ML")), # full-information likelihood for missing values
         omega = map_dbl(.x = model,
                         ~ semTools::reliability(.x)["omega",])) %>%  # McDonald omega: scale reliability
  select(omega) %>% 
  unnest()
loneli_reliability
