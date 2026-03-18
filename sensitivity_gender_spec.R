#############################################################
### EXPLORING HETEROGENEITY: GENDER-SPECIFIC ASSOCIATIONS ### 
#############################################################


# loading packages
library(readr) # data loading
library(tidyverse) # data manipulation
library(survival) # Cox proportional hazards models
library(WeightIt) # inverse-probability-of-treatment-weighting (IPTW)
library(cobalt) # covariate balance assessments
library(marginaleffects) # computing marginal contrasts


# loading data
NorPD_UiN_analysis <- read_csv("/Data/NorPD_UiN_analysis")
NorPD_UiN_analysis_imp <- read_csv("/Data/NorPD_UiN_analysis_imp")

# pre-processing
NorPD_UiN_analysis_imp <- NorPD_UiN_analysis_imp %>% 
  group_by(.imp) %>% 
  mutate(across(starts_with("loneli_direct_"),
                ~ (.x - mean(.x))/sd(.x))) %>% 
  ungroup() %>% 
  select(id, .imp, age:bmi_3, loneliness_1:loneliness_4, loneli_direct_1:loneli_direct_4, 
         starts_with(c("event_", "time_")), ends_with(c("_DDD", "_nr"))) %>% 
  pivot_longer(names_to = c("measure", "timepoint"),
               names_pattern = "(.*)_(.*)",
               values_to = "zscore",
               cols = `loneliness_1`:`loneli_direct_4`) %>% 
  pivot_wider(names_from = timepoint,
              values_from = zscore) %>% 
  mutate(loneliness_1 = `1`,
         loneliness_2 = `2`,
         loneliness_3 = `3`,
         loneliness_4 = `4`,
         measure = ifelse(measure == "loneli_direct", "direct", "indirect")) %>% 
  select(id, .imp, measure, loneliness_1:loneliness_4, age:bmi_3, starts_with(c("event_", "time_")), ends_with(c("_DDD", "_nr")))




### analgesics and antibiotics use: number of prescriptions and defined daily doses (DDD)
nonadj_models_acute_gender <- NorPD_UiN_analysis_imp %>%
  select(.imp, id, gender, measure, loneliness_2:loneliness_4, starts_with(c("antibiotics", "analgesics"))) %>% 
  pivot_longer(names_to = "outcome",
               values_to = "prescr",
               cols = antibiotics_DDD:analgesics_nr) %>% 
  group_by(.imp, measure, outcome) %>% 
  nest() %>% 
  mutate(longi_model = map(.x = data,
                           ~ glm(prescr ~ gender * (loneliness_2 + loneliness_3 + loneliness_4), data = .x, family = poisson)), 
         avg_comp_all = map(.x = longi_model, 
                            ~ avg_predictions(.x,
                                              type = "response",
                                              variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                              by = c("loneliness_2", "loneliness_3","loneliness_4", "gender"))),
         avg_comp_all_male = map(.x = avg_comp_all,
                                ~ hypotheses(.x, hypothesis =  "b16 - b2 = 0", vcov = "HC3")), 
         avg_comp_all_female = map(.x = avg_comp_all,
                                 ~ hypotheses(.x, hypothesis =  "b15 - b1 = 0", vcov = "HC3")),
         avg_comp_all_int = map(.x = avg_comp_all,
                                ~ hypotheses(.x, hypothesis =  "(b16 - b2) - (b15 - b1) = 0", vcov = "HC3"))) # ref.: female

nonadj_results_acute_gender <- select(nonadj_models_acute_gender, ends_with(c("male", "female", "int"))) %>% 
  pivot_longer(cols = avg_comp_all_male:avg_comp_all_int, names_to = "effect", values_to = "contrasts") %>% 
  group_by(measure, outcome, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)


table_nonadj_results_acute_gender <- nonadj_results_acute_gender %>% 
  mutate(gender = case_when(str_detect(effect, "_male") ~ "male", 
                            str_detect(effect, "_female") ~ "female",
                            str_detect(effect, "_int") ~ "interaction"),
         effect = case_when(str_detect(effect, "t1") ~ "adolescence (17y)",
                            str_detect(effect, "t2") ~ "emerging adulthod (22y)",
                            str_detect(effect, "t3") ~ "young adulthood (28y)",
                            str_detect(effect, "all") ~ "cumulative (17-28y)")) %>% 
  mutate(estimate = paste0(formatC(round(estimate, 2), format = "f", digits = 2), " (", 
                           formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                           formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")")) %>% 
  select(measure, outcome, gender, effect, lonely_measure, estimate) %>% 
  arrange(measure, outcome, gender, effect) %>% 
  filter(outcome %in% c("analgesics_DDD", "antibiotics_DDD")) %>% 
  pivot_wider(names_from = "gender",
              values_from = "estimate") %>% 
  select(measure, outcome, male, female, interaction)  
  
write.csv(table_nonadj_results_acute_gender, "/Results/Tables/table_nonadj_results_acute_gender")




### risk of prescriptions for chronic physical health-related conditions: cardiovascular, metabolic, gastrointestinal, autoimmune and migraine 

### gender-specific outcome counts: in women all >= 30, in men only C and H02
NorPD_UiN_analysis_imp %>% 
  select(id, .imp, gender, loneliness_2:loneliness_4, starts_with(c("event", "time")), -ends_with(c("N05B"))) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_C`:`time_N02CC`) %>% 
  pivot_wider(names_from = time_event,
              values_from = prescr) %>% 
  group_by(gender, outcome) %>% 
  summarise(sum(event)/10) 

nonadj_models_chronic_gender <- NorPD_UiN_analysis_imp %>%
  select(id, measure, gender, .imp, loneliness_2:loneliness_4, starts_with(c("event", "time")), -ends_with(c("N05B"))) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_C`:`time_N02CC`) %>% 
  pivot_wider(names_from = time_event,
              values_from = prescr) %>% 
  group_by(.imp, measure, outcome) %>% 
  nest() %>% 
  mutate(cox_longi_model = map(.x = data,
                               ~ coxph(Surv(time, event) ~ gender * (loneliness_2 + loneliness_3 + loneliness_4), robust = T, data = .x)), 
         avg_comp_all = map(.x = cox_longi_model,
                               ~ avg_predictions(.x, 
                                                 type = "risk", 
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4", "gender"))),
         avg_comp_all_male = map(.x = avg_comp_all,
                                 ~ hypotheses(.x, hypothesis =  "b16/b2 = 0", type = "risk", comparison = "ratio")), 
         avg_comp_all_female = map(.x = avg_comp_all,
                                   ~ hypotheses(.x, hypothesis =  "b15/b1 = 0", type = "risk", comparison = "ratio")),
         avg_comp_all_int = map(.x = avg_comp_all,
                                ~ hypotheses(.x, hypothesis =  "(b16/b2) / (b15/b1) = 0", type = "risk", comparison = "ratio")))

nonadj_results_chronic_gender <- select(nonadj_models_chronic_gender, ends_with(c("male", "female", "int"))) %>% 
  pivot_longer(cols = avg_comp_all_male:avg_comp_all_int, names_to = "effect", values_to = "contrasts") %>% 
  group_by(measure, outcome, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)

table_chronic_gender <- nonadj_results_chronic_gender %>% 
  mutate(gender = case_when(str_detect(effect, "_male") ~ "male", 
                            str_detect(effect, "_female") ~ "female",
                            str_detect(effect, "_int") ~ "interaction"),
         outcome = factor(outcome, levels = c("H03A", "C", "A07", "H02", "L04A", "N02CC", "N06A"), ordered = TRUE),
         effect = factor(case_when(str_detect(effect, "t1") ~ "adolescence (17y)",
                                   str_detect(effect, "t2") ~ "emerging adulthod (22y)",
                                   str_detect(effect, "t3") ~ "young adulthood (28y)",
                                   str_detect(effect, "all") ~ "cumulative (17-28y)"),
                         levels = c("adolescence (17y)", "emerging adulthod (22y)", "young adulthood (28y)", "cumulative (17-28y)"),
                         ordered = T)) %>% 
  mutate(estimate = paste0(formatC(round(estimate, 2), format = "f", digits = 2), " (", 
                           formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                           formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")"),
         estimate = ifelse(gender %in% c("male", "interaction") & outcome %in% c("H03A", "A07", "L04A", "N02CC"), NA, estimate)) %>% 
  filter(effect == "cumulative (17-28y)") %>% 
  select(measure, outcome, lonely_measure, gender, effect, estimate) %>% 
  arrange(measure, outcome, gender, effect) %>% 
  pivot_wider(names_from = "gender",
              values_from = "estimate") %>% 
  select(measure, outcome, lonely_measure, male, female, interaction)

write.csv(table_chronic_gender, "/Results/Tables/table_chronic_gender")




