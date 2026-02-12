#############################################################################
### LONELINESS AND PRESCRIPTIONS FOR PHYSICAL HEALTH-RELATED MEDICATION ####
############################################################################



##############################
### STATISTICAL MODELLING ###
#############################

# loading packages
library(readr) # data loading
library(tidyverse) # data manipulation
library(survival) # Cox proportional hazards models
library(WeightIt) # inverse-probability-of-treatment-weighting (IPTW)
library(cobalt) # covariate balance assessments
library(marginaleffects) # computing marginal contrasts


# loading data
NorPD_UiN_analysis <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_UiN_analysis")
NorPD_UiN_analysis_imp <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_UiN_analysis_imp")
NorPD_nonresp_weights <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_nonresp_weights") %>% 
  mutate(id = norpd_id, nonresp_weights = weights) %>% select(id, nonresp_weights)


# pre-processing
NorPD_UiN_analysis_imp <- merge(NorPD_UiN_analysis_imp, NorPD_nonresp_weights, by = "id") %>% 
  group_by(.imp) %>% 
  mutate(across(starts_with("loneli_direct_"),
                ~ (.x - mean(.x))/sd(.x))) %>% 
  ungroup() %>% 
  select(id, nonresp_weights, .imp, age:bmi_3, loneliness_1:loneliness_4, loneli_direct_1:loneli_direct_4, 
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
  select(id, nonresp_weights, .imp, measure, loneliness_1:loneliness_4, age:bmi_3, starts_with(c("event_", "time_")), ends_with(c("_DDD", "_nr")))


#########################
### UNADJUSTED MODELS ###
#########################

# analgesics and antibiotics use: number of prescriptions and defined daily doses (DDD)
nonadj_models_acute <- NorPD_UiN_analysis_imp %>%
  select(.imp, id, nonresp_weights, measure, loneliness_2:loneliness_4, starts_with(c("antibiotics", "analgesics"))) %>% 
  pivot_longer(names_to = "outcome",
               values_to = "prescr",
               cols = antibiotics_DDD:analgesics_nr) %>% 
  group_by(.imp, measure, outcome) %>% 
  nest() %>% 
  mutate(longi_model = map(.x = data,
                           ~ glm(prescr ~ loneliness_2 + loneliness_3 + loneliness_4, data = .x, weights = .x$nonresp_weights, family = poisson)), 
         avg_predictions = map2(.x = longi_model, .y = data,
                               ~ avg_predictions(.x,
                                                 wts = .y$nonresp_weights,
                                                 type = "response",
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         model_t1 = map(.x = data,
                        ~ glm(prescr ~ loneliness_2, data = .x, weights = .x$nonresp_weights, family = poisson)),
         model_t2 = map(.x = data,
                        ~ glm(prescr ~ loneliness_3, data = .x, weights = .x$nonresp_weights, family = poisson)),
         model_t3 = map(.x = data,
                        ~ glm(prescr ~ loneliness_4, data = .x, weights = .x$nonresp_weights, family = poisson)),
         avg_comp_all = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b8 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)), # 1,1,1 vs 0,0,0
         avg_comp_te_t1 = map2(.x = model_t1, .y = data,
                            ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_2 = c(0,1)))), # 1 vs 0
         avg_comp_te_t2 = map2(.x = model_t2, .y = data,
                            ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_3 = c(0,1)))), # 1 vs 0
         avg_comp_te_t3 = map2(.x = model_t3, .y = data,
                            ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_4 = c(0,1)))), # 1 vs 0
         avg_comp_cde_t3 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b2 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)), # 0,0,1 vs 0,0,0
         avg_comp_cde_t2 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b3 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)), # 0,1,0 vs 0,0,0
         avg_comp_cde_t1 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b5 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights))) # 1,0,0 vs 0,0,0 


nonadj_results_acute <- select(nonadj_models_acute, starts_with("avg_comp")) %>% 
  pivot_longer(cols = avg_comp_all:avg_comp_cde_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)
write.csv(nonadj_results_acute, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/raw_nonadj_results_acute")


# risk of prescriptions for chronic physical health-related conditions: cardiovascular, metabolic, gastrointestinal, autoimmune and migraine 
nonadj_models_chronic <- NorPD_UiN_analysis_imp %>%
  select(id, nonresp_weights, .imp, measure, loneliness_2:loneliness_4, starts_with(c("event", "time")), -ends_with(c("N05B", "N06A"))) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_C`:`time_N02CC`) %>% 
  pivot_wider(names_from = time_event,
              values_from = prescr) %>% 
  group_by(.imp, measure, outcome) %>% 
  nest() %>% 
  mutate(cox_longi_model = map(.x = data,
                              ~ coxph(Surv(time, event) ~ loneliness_2 + loneliness_3 + loneliness_4, robust = T, data = .x, weights = .x$nonresp_weights)), 
         avg_predictions = map2(.x = cox_longi_model, .y = data,
                               ~ avg_predictions(.x, 
                                                 wts = .y$nonresp_weights,
                                                 type = "risk", # hazard rate
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         t1_model = map(.x = data,
                        ~ coxph(Surv(time, event) ~ loneliness_2, robust = T, data = .x, weights = .x$nonresp_weights)), 
         t2_model = map(.x = data,
                        ~ coxph(Surv(time, event) ~ loneliness_3, robust = T, data = .x, weights = .x$nonresp_weights)), 
         t3_model = map(.x = data,
                        ~ coxph(Surv(time, event) ~ loneliness_4, robust = T, data = .x, weights = .x$nonresp_weights)), 
         avg_comp_all = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b8/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  # 1,1,1 vs 0,0,0
         avg_comp_te_t1 = map(.x = t1_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_2 = c(0,1)))), # 1 vs 0
         avg_comp_te_t2 = map(.x = t2_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_3 = c(0,1)))), # 1 vs 0
         avg_comp_te_t3 = map(.x = t3_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_4 = c(0,1)))), # 1 vs 0
         avg_comp_cde_t3 = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b2/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  # 0,0,1 vs 0,0,0
         avg_comp_cde_t2 = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b3/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  # 0,1,0 vs 0,0,0
         avg_comp_cde_t1 = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b5/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)))  # 1,0,0 vs 0,0,0

nonadj_results_chronic <- select(nonadj_models_chronic, starts_with("avg_comp")) %>% 
  pivot_longer(cols = avg_comp_all:avg_comp_cde_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)
write.csv(nonadj_results_chronic, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/raw_nonadj_results_chronic")
  



########################
### EXPOSURE MODELS ###
#######################

# specifying longitudinal exposure (treatment assignment) models
exposure_models <- list(loneliness_2 ~ gender + age + ethnicity + parental_education + urbanity + 
                                       warm_parenting + parental_alcoholuse + parental_smoking + 
                                       asthma_allergy + phys_disability +
                                       social_support_1 + depression_1 + 
                                       living_situation_1 + relationship_1 + friends_1 + employment_1 +
                                       smoking_1 + alcohol_use_1 + illicit_drug_1 + phys_exer_1 + bmi_1 +
                                       loneliness_1, 
                        loneliness_3 ~ gender + age + ethnicity + parental_education + urbanity + 
                                       warm_parenting + parental_alcoholuse + parental_smoking + 
                                       asthma_allergy + phys_disability +
                                       social_support_2 + depression_2 + 
                                       living_situation_2 + relationship_2 + friends_2 + employment_2 +
                                       smoking_2 + alcohol_use_2 + illicit_drug_2 + phys_exer_2 + bmi_2 +
                                       loneliness_2,
                        loneliness_4 ~ gender + age + ethnicity + parental_education + urbanity + 
                                       warm_parenting + parental_alcoholuse + parental_smoking + 
                                       asthma_allergy + phys_disability +
                                       social_support_3 + depression_3 + 
                                       living_situation_3 + relationship_3 + friends_3 + employment_3 +
                                       smoking_3 + alcohol_use_3 + illicit_drug_3 + phys_exer_3 + bmi_3 +
                                       loneliness_3)

# Assessing covariate balance before weighting
# (mean differences for binary, and correlations for continuous covariates (default))
balance_bef_indi <- bal.tab(exposure_models,
                       data = NorPD_UiN_analysis_imp[NorPD_UiN_analysis_imp$measure == "indirect",], 
                       s.weights = NorPD_UiN_analysis_imp[NorPD_UiN_analysis_imp$measure == "indirect",]$nonresp_weights,
                       thresholds = 0.10, 
                       which.time = TRUE) # balance summary across all timepoints

balance_bef_dir <- bal.tab(exposure_models,
                       data = NorPD_UiN_analysis_imp[NorPD_UiN_analysis_imp$measure == "direct",], 
                       s.weights = NorPD_UiN_analysis_imp[NorPD_UiN_analysis_imp$measure == "direct",]$nonresp_weights,
                       thresholds = 0.10, 
                       which.time = TRUE) 

balance_bef_indi$Balanced.correlations["Not Balanced, >0.1", "count"] # number of non-balanced covariates (> 0.10) before weighting
balance_bef_dir$Balanced.correlations["Not Balanced, >0.1", "count"] 

# Constructing Inverse-Probability-of-Treatment-Weights (IPTWs) for Marginal Structural Models (MSM)
weighted_dataset <- NorPD_UiN_analysis_imp %>% 
  group_by(measure, .imp) %>% 
  nest() %>% 
  mutate(weighted_data = map(.x = data,
                                   ~ weightitMSM(exposure_models,
                                                 data = .x,
                                                 s.weights = .x$nonresp_weights,
                                                 method = "cbps", 
                                                 stabilize = TRUE)), 
         balance_after = map(.x = weighted_data,
                              ~ bal.tab(.x, thresholds = 0.10, which.time = TRUE)),
         non_balanced = map_dbl(.x = balance_after,
                                ~ .x$Balanced.correlations["Not Balanced, >0.1", "count"]),
         ESS = map_dbl(.x = balance_after,
                       ~ .x$Observations$loneliness_2["Adjusted",]),
         t1_weighted_data = map(.x = data,
                                ~ weightit(exposure_models[[1]], # loneliness T1 (adolescence)
                                           data = .x,
                                           s.weights = .x$nonresp_weights,
                                           method = "cbps", 
                                           stabilize = TRUE)),
         t2_weighted_data = map(.x = data,
                                ~ weightit(exposure_models[[2]], # loneliness T2 (emerging adulthood)
                                           data = .x,
                                           s.weights = .x$nonresp_weights,
                                           method = "cbps", 
                                           stabilize = TRUE)),
         t3_weighted_data = map(.x = data,
                                ~ weightit(exposure_models[[3]], # loneliness T3 (young adulthood)
                                           data = .x,
                                           s.weights = .x$nonresp_weights,
                                           method = "cbps", 
                                           stabilize = TRUE)))

weighted_dataset <- weighted_dataset %>% 
  mutate(weighted_data_trim = map(.x = weighted_data, 
                                  ~ trim(.x, at = .99)),
         balance_after_trim = map(.x = weighted_data_trim,
                                  ~ bal.tab(.x, thresholds = 0.10, which.time = TRUE)),
         non_balanced_trim = map_dbl(.x = balance_after_trim,
                                    ~ .x$Balanced.correlations["Not Balanced, >0.1", "count"]),
         ESS_trim = map_dbl(.x = balance_after_trim,
                          ~ .x$Observations$loneliness_2["Adjusted",]),
         t1_weighted_data_trim = map(.x = t1_weighted_data, 
                                  ~ trim(.x, at = .99)),
         t2_weighted_data_trim = map(.x = t2_weighted_data, 
                                  ~ trim(.x, at = .99)),
         t3_weighted_data_trim = map(.x = t3_weighted_data, 
                                  ~ trim(.x, at = .99)))
         

# assessing balance after weighting
max(weighted_dataset$non_balanced) # number of non-balanced covariates (> 0.10) = 0
mean(weighted_dataset$ESS) # effective sample size (ESS) low -> trimming at 99th percentile
summary(weighted_dataset$weighted_data[[2]]$weights) # extreme weights (> 20) -> trimming at 99th percentile

max(weighted_dataset$non_balanced_trim)  # number of non-balanced covariates still 0
mean(weighted_dataset$ESS_trim) # improved ESS
summary(weighted_dataset$weighted_data_trim[[2]]$weights) # no extreme weights (> 10)


weighted_data_analysis <- select(weighted_dataset, .imp, measure, weighted_data, weighted_data_trim, t1_weighted_data:t3_weighted_data, t1_weighted_data_trim:t3_weighted_data_trim)
save(weighted_data_analysis, file = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/weighted_data_analysis.RData")




######################
### OUTCOME MODELS ###
######################

# risk of prescriptions for chronic physical health-related conditions: cardiovascular, metabolic, gastrointestinal, autoimmune and migraine 
MSM_chronic_models <- unnest(weighted_dataset, data) %>% 
  select(id, nonresp_weights, loneliness_2:loneliness_4, starts_with(c("event", "time")),-ends_with(c("N05B", "N06A"))) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_C`:`time_N02CC`) %>% 
  pivot_wider(names_from = time_event,
              values_from = prescr) %>% 
  group_by(.imp, measure, outcome) %>% 
  nest() %>% 
  left_join(., weighted_data_analysis, by = c(".imp", "measure")) %>% 
   mutate(longi_model = map2(.x = data, .y = weighted_data_trim, 
                           ~ coxph_weightit(survival::Surv(time, event) ~ loneliness_2 + loneliness_3 + loneliness_4, 
                                            vcov = "HC0",
                                            data = .x, weightit = .y, weights = .x$nonresp_weights)),
         avg_predictions = map(.x = longi_model, .y = data,
                               ~ avg_predictions(.x, 
                                                 type = "risk",
                                                 vcov = "HC0", 
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         t1_model = map2(.x = data, .y = t1_weighted_data_trim,
                         ~ coxph_weightit(survival::Surv(time, event) ~ loneliness_2, 
                                          vcov = "HC0",
                                          data = .x, weightit = .y, weights = .x$nonresp_weights)),
         t2_model = map2(.x = data, .y = t2_weighted_data_trim,
                         ~ coxph_weightit(survival::Surv(time, event) ~ loneliness_3, 
                                          vcov = "HC0",
                                          data = .x, weightit = .y, weights = .x$nonresp_weights)),
         t3_model = map2(.x = data, .y = t3_weighted_data_trim,
                         ~ coxph_weightit(survival::Surv(time, event) ~ loneliness_4, 
                                          vcov = "HC0",
                                          data = .x, weightit = .y, weights = .x$nonresp_weights)),
         avg_comp_all = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b8/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),
         avg_comp_te_t1 = map(.x = t1_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_2 = c(0,1)))),
         avg_comp_te_t2 = map(.x = t2_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_3 = c(0,1)))),
         avg_comp_te_t3 = map(.x = t3_model, .y = data,
                           ~ avg_comparisons(.x, comparison = "ratio", wts = .y$nonresp_weights, variables = list(loneliness_4 = c(0,1)))),
         avg_comp_cde_t3 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b2/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_cde_t2 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b3/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_cde_t1 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b5/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)))

MSM_chronic_results <- select(MSM_chronic_models, starts_with("avg_comp")) %>% 
  pivot_longer(cols = avg_comp_all:avg_comp_cde_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)
write.csv(MSM_chronic_results, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/raw_MSM_chronic_results")

# analgesics and antibiotics use: number of prescriptions and defined daily doses (DDD)
MSM_acute_models <- unnest(weighted_dataset, data) %>% 
  select(id, nonresp_weights, measure, loneliness_2:loneliness_4, starts_with(c("antibiotics", "analgesics"))) %>% 
    pivot_longer(names_to = "outcome",
                 values_to = "prescr",
                 cols = starts_with(c("antibiotics", "analgesics"))) %>% 
    group_by(.imp, measure, outcome) %>%   
    nest() %>% 
  left_join(., weighted_data_analysis, by = c(".imp", "measure")) %>% 
  mutate(regres_models = map2(.x = data, .y = weighted_data_trim,
                              ~ glm_weightit(prescr ~ loneliness_2 + loneliness_3 + loneliness_4,
                                             family = poisson, weights = .x$nonresp_weights,
                                             data = .x, weightit = .y)),
         avg_predictions = map(.x = regres_models, .y = data,
                               ~ avg_predictions(.x,
                                                 wts = .y$nonresp_weights,
                                                 type = "response",
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         model_t1 = map2(.x = data, .y = t1_weighted_data_trim,
                        ~ glm_weightit(prescr ~ loneliness_2, data = .x, weightit = .y, weights = .x$nonresp_weights, family = poisson)),
         model_t2 = map2(.x = data, .y = t2_weighted_data_trim,
                        ~ glm_weightit(prescr ~ loneliness_3, data = .x, weightit = .y, weights = .x$nonresp_weights, family = poisson)),
         model_t3 = map2(.x = data, .y = t3_weighted_data_trim,
                        ~ glm_weightit(prescr ~ loneliness_4, data = .x, weightit = .y, weights = .x$nonresp_weights, family = poisson)),
         avg_comp_all = map2(.x = avg_predictions,  .y = data,
                             ~ hypotheses(.x, hypothesis =  "b8 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)),
         avg_comp_te_t1 = map2(.x = model_t1, .y = data,
                               ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_2 = c(0,1)))),
         avg_comp_te_t2 = map2(.x = model_t2, .y = data,
                               ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_3 = c(0,1)))),
         avg_comp_te_t3 = map2(.x = model_t3, .y = data,
                               ~ avg_comparisons(.x, type = "response", vcov = "HC3", wts = .y$nonresp_weights, variables = list(loneliness_4 = c(0,1)))), 
         avg_comp_cde_t3 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b2 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)),
         avg_comp_cde_t2 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b3 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)),
         avg_comp_cde_t1 = map2(.x = avg_predictions,  .y = data,
                            ~ hypotheses(.x, hypothesis =  "b5 - b1 = 0", vcov = "HC3", wts = .y$nonresp_weights)))
  
MSM_acute_results <- select(MSM_acute_models, starts_with("avg_comp")) %>% 
  pivot_longer(cols = avg_comp_all:avg_comp_cde_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)
write.csv(MSM_acute_results, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/raw_MSM_acute_results")





###############
### RESULTS ###
###############

# chronic health conditions
combined_results_chronic <- mutate(nonadj_results_chronic, estimation = "unadjusted") %>% 
  bind_rows(mutate(MSM_chronic_results, estimation = "IPTW-MSM")) %>% 
  mutate(outcome = factor(case_when(outcome == "C" ~ "cardiovascular medication",
                                    outcome == "H03A" ~ "thyroid preparations", 
                                    outcome == "A07" ~ "gastrointestinal medication",
                                    outcome == "H02" ~ "systemic corticosteroids", 
                                    outcome == "L04A" ~ "immunosuppressants",
                                    outcome == "N02CC" ~ "antimigraine preparations"),
                          levels = c("antimigraine preparations", "immunosuppressants", "systemic corticosteroids",
                                     "gastrointestinal medication", "thyroid preparations", "cardiovascular medication"),
                          ordered = T),
         estimation_detail = case_when(estimation == "unadjusted" & str_detect(effect, "_te_|_de_") ~ "unadjusted (TE)",
                                       estimation == "unadjusted" & str_detect(effect, "_cde_|_all") ~ "unadjusted (CDE)",
                                       estimation == "IPTW-MSM" & str_detect(effect, "_cde_|_all") ~ "IPTW-MSM (CDE)",
                                       estimation == "IPTW-MSM" & str_detect(effect, "_te_|_de_") ~ "IPTW (TE)"),
         measure = ifelse(measure == "direct", "Direct/explicit measure: 'I feel lonely'", "Indirect/implicit measure: UCLA-4"),
         effect = factor(case_when(str_detect(effect, "t1") ~ "adolescence (17y)",
                                   str_detect(effect, "t2") ~ "emerging adulthod (22y)",
                                   str_detect(effect, "t3") ~ "young adulthood (28y)",
                                   str_detect(effect, "all") ~ "cumulative (17-28y)"),
                         levels = c("adolescence (17y)", "emerging adulthod (22y)", "young adulthood (28y)", "cumulative (17-28y)"),
                         ordered = T)) 

plot_combined_chronic <- combined_results_chronic %>% 
  filter(effect == "cumulative (17-28y)") %>% 
  arrange(effect) %>% 
  ggplot(aes(x = estimate, y = outcome, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  lemon::facet_rep_wrap(~ measure, 
                        ncol = 2,
                        repeat.tick.labels = FALSE) +
  labs(title = "Cumulative effect of loneliness on risk of prescriptions for physical health-related medications in young adulthood",
       x = "Hazard ratio (95% confidence interval)",
       color = "Estimation",
       y = "") +
  theme_minimal() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(strip.text = element_text(face = "bold", size = 14, family = "Times"),
        axis.text.y = element_text(size = 14, family = "Times", color = "black"),
        legend.text = element_text(size = 14, family = "Times", color = "black"),
        text = element_text(size = 14, family = "Times"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 20)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_combined_chronic.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 13.5, 
       height = 8,  
       bg="white",
       dpi=900)
  
plot_combined_chronic_timespec <- combined_results_chronic %>% 
  arrange(effect) %>% 
  ggplot(aes(x = estimate, y = outcome, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.3, size = 1) +
  lemon::facet_rep_wrap(~ measure + effect, 
                        ncol = 4, nrow = 2,
                        repeat.tick.labels = FALSE) +
  labs(title = "Loneliness and subsequent risk of prescriptions for physical health-related medications in young adulthood",
       x = "Hazard ratio (95% confidence interval)",
       color = "Estimation",
       y = "") +
  theme_minimal() +
  geom_vline(xintercept = 1, linetype = "dashed") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(text = element_text(size = 12, family = "Times"),
        plot.title = element_text(hjust = 0.5, margin = margin(b = 20)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_combined_chronic_timespec.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 13.5, 
       height = 8,  
       bg="white",
       dpi=900)

table_combined_chronic <- combined_results_chronic %>% 
  mutate(estimate = paste0(formatC(round(estimate, 2), format = "f", digits = 2), " (", 
                           formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                           formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")"),
         measure = ifelse(measure == "Direct/explicit measure: 'I feel lonely'", "direct", "indirect")) %>% 
  select(-c(`2.5 %`, `97.5 %`)) %>% 
  pivot_wider(names_from = "estimation",
              values_from = "estimate") %>% 
  pivot_wider(names_from = "measure",
              values_from = c("unadjusted", `IPTW-MSM`)) %>% 
  arrange(desc(outcome), effect) %>% 
  select(outcome, effect, ends_with(c("direct", "indirect")))
write.csv(table_combined_chronic, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/table_combined_chronic")


# antibiotics and analgesics
combined_results_acute <- mutate(nonadj_results_acute, estimation = "unadjusted") %>% 
  bind_rows(mutate(MSM_acute_results, estimation = "IPTW-MSM")) %>% 
  separate(outcome, into = c("outcome", "measure_prescr"), sep = "_") %>% 
  mutate(measure_prescr = case_when(measure_prescr == "nr" ~ "number of prescriptions",
                                    measure_prescr == "DDD" ~ "defined daily doses"),
         measure = ifelse(measure == "direct", "Direct/explicit measure: 'I feel lonely'", "Indirect/implicit measure: UCLA-4"),
         effect = factor(case_when(str_detect(effect, "t1") ~ "adolescence (17y)",
                                   str_detect(effect, "t2") ~ "emerging adulthod (22y)",
                                   str_detect(effect, "t3") ~ "young adulthood (28y)",
                                   str_detect(effect, "all") ~ "cumulative (17-28y)"),
                         levels = c("adolescence (17y)", "emerging adulthod (22y)", "young adulthood (28y)", "cumulative (17-28y)"),
                         ordered = T))

plot_combined_acute <- combined_results_acute %>% 
  filter(effect == "cumulative (17-28y)",
         measure_prescr == "defined daily doses") %>% 
  arrange(effect) %>% 
  ggplot(aes(x = estimate, y = outcome, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.15, size = 1) +
  lemon::facet_rep_wrap(~ measure, 
                        ncol = 2,
                        scales = "free_x",
                        repeat.tick.labels = FALSE) +
  labs(title = "Cumulative effect of loneliness on antibiotics and analgesics use in young adulthood",
       color = "Estimation",
       x = "Mean difference, DDD (95% confidence interval)", 
       y = "") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(strip.text = element_text(face = "bold", size = 14, family = "Times"),
        axis.text.y = element_text(size = 14, family = "Times", color = "black"),
        legend.text = element_text(size = 14, family = "Times", color = "black"),
        text = element_text(size = 14, family = "Times"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 20)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_combined_acute.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 10, 
       height = 5,  
       bg="white",
       dpi=700)


plot_combined_acute_prescr_nr <- combined_results_acute %>% 
  filter(effect == "cumulative (17-28y)",
         measure_prescr == "number of prescriptions") %>% 
  arrange(effect) %>% 
  ggplot(aes(x = estimate, y = outcome, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.15, size = 1) +
  lemon::facet_rep_wrap(~ measure, 
                        ncol = 2,
                        scales = "free_x",
                        repeat.tick.labels = FALSE) +
  labs(title = "Cumulative effect of loneliness on antibiotics and analgesics use in young adulthood",
       color = "Estimation",
       x = "Mean difference, number of prescriptions (95% confidence interval)", 
       y = "") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(strip.text = element_text(face = "bold", size = 14, family = "Times"),
        axis.text.y = element_text(size = 14, family = "Times", color = "black"),
        legend.text = element_text(size = 14, family = "Times", color = "black"),
        text = element_text(size = 14, family = "Times"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 20)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_combined_acute_prescr_nr.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 10, 
       height = 5,  
       bg="white",
       dpi=700)


plot_combined_acute_timespec <- combined_results_acute %>% 
  filter(measure_prescr == "defined daily doses") %>%
  arrange(effect) %>% 
  ggplot(aes(x = estimate, y = outcome, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.15, size = 1) +
  lemon::facet_rep_wrap(~ measure + effect, 
                        ncol = 4,
                        scales = "free_x",
                        repeat.tick.labels = FALSE) +
  labs(title = "Loneliness and antibiotics and analgesics use in young adulthood",
       x = "Mean difference, DDD (95% confidence interval)", 
       y = "") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(color = guide_legend(reverse = TRUE)) +
  theme(strip.text = element_text(size = 12, family = "Times"),
        axis.text.y = element_text(size = 12, family = "Times", color = "black"),
        legend.text = element_text(size = 12, family = "Times", color = "black"),
        text = element_text(size = 12, family = "Times"),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 20)),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),  
        axis.line = element_line(color = "black"))

ggsave(filename = "plot_combined_acute_timespec.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 13.8, 
       height = 8,  
       bg="white",
       dpi=700)


table_combined_acute <- combined_results_acute %>% 
  filter(measure_prescr == "defined daily doses") %>% 
  mutate(estimate = paste0(formatC(round(estimate, 2), format = "f", digits = 2), " (", 
                           formatC(round(`2.5 %`, 2), format = "f", digits = 2), ", ", 
                           formatC(round(`97.5 %`, 2), format = "f", digits = 2), ")"),
          measure = ifelse(measure == "Direct/explicit measure: 'I feel lonely'", "direct", "indirect")) %>% 
  select(-c(`2.5 %`, `97.5 %`)) %>% 
  pivot_wider(names_from = "estimation",
              values_from = "estimate") %>% 
  pivot_wider(names_from = "measure",
              values_from = c("unadjusted", `IPTW-MSM`)) %>% 
  arrange(desc(outcome), effect) %>% 
  select(outcome, effect, ends_with(c("direct", "indirect")))
write.csv(table_combined_acute, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/table_combined_acute")


