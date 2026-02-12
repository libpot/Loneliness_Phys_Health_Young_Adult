##################################################
### POSITIVE CONTROL ANALYSIS: ANTIDEPRESSANTS ###
##################################################


# loading data
NorPD_UiN_analysis <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_UiN_analysis")
NorPD_UiN_analysis_imp <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_UiN_analysis_imp")
NorPD_antidep_washout <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_antidep_washout")
NorPD_nonresp_weights <- read_csv("N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Data/NorPD_nonresp_weights") %>% 
  mutate(id = norpd_id, nonresp_weights = weights) %>% select(id, nonresp_weights)

library(contsurvplot)
library(survminer)
library(cowplot)

# pre-processing
NorPD_UiN_analysis_dep_imp <- merge(NorPD_UiN_analysis_imp, NorPD_nonresp_weights, by = "id") %>% 
  group_by(.imp) %>% 
  mutate(across(starts_with("loneli_direct_"),
                ~ (.x - mean(.x))/sd(.x))) %>% 
  ungroup() %>% 
  select(id, nonresp_weights, .imp, age:bmi_3, loneliness_1:loneliness_4, loneli_direct_1:loneli_direct_4, ends_with(c("N06A"))) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_N06A`:`time_N06A`) %>% 
  pivot_wider(names_from = time_event,
              values_from = prescr) %>% 
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
  select(id, nonresp_weights, .imp, age:bmi_3, outcome, event, time, measure, loneliness_1:loneliness_4)

# fitting unadjusted Cox regression models
nonadj_models_psychotropic <- NorPD_UiN_analysis_dep_imp %>%
  group_by(.imp, outcome, measure) %>% 
  nest() %>% 
  mutate(cox_longi_model = map(.x = data,
                               ~ coxph(Surv(time, event) ~ loneliness_2 + loneliness_3 + loneliness_4, robust = T, data = .x, weights = .x$nonresp_weights)), 
         avg_predictions = map(.x = cox_longi_model, .y = data,
                               ~ avg_predictions(.x, 
                                                 wts = .y$nonresp_weights,
                                                 type = "risk", # hazard rate
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         avg_comp_all = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b8/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t3 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b2/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t2 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b3/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t1 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b5/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights))) 

nonadj_results_psychotropic <- select(nonadj_models_psychotropic, starts_with("avg_comp")) %>% 
  pivot_longer(cols = avg_comp_all:avg_comp_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(outcome, measure, estimate, `2.5 %`, `97.5 %`) 




#######################
### ADJUSTED MODELS ###
######################

adj_models_psychotropic <- unnest(weighted_dataset, data) %>% 
  select(id, nonresp_weights, loneliness_2:loneliness_4, ends_with("N06A")) %>% 
  pivot_longer(names_to = c("time_event", "outcome"),
               names_pattern = "(.*)_(.*)",
               values_to = "prescr",
               cols = `event_N06A`:`time_N06A`) %>% 
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
         avg_comp_all = map2(.x = avg_predictions, .y = data,
                             ~ hypotheses(.x, hypothesis =  "b8/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t3 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b2/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t2 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b3/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)),  
         avg_comp_t1 = map2(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b5/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)))

MSM_results_psychotropic <- select(adj_models_psychotropic, starts_with("avg_comp")) %>%
  pivot_longer(cols = avg_comp_all:avg_comp_t1, names_to = "effect", values_to = "contrasts") %>% 
  group_by(outcome, measure, effect) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$contrasts), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, outcome, estimate, `2.5 %`, `97.5 %`)



### RESULTS ###

results_combined_antidepressants <- mutate(nonadj_results_psychotropic, estimation = "unadjusted") %>% 
  bind_rows(mutate(MSM_results_psychotropic, estimation = "IPTW-MSM")) %>% 
  mutate(measure = ifelse(measure == "direct", "Direct: 'I feel lonely'", "Indirect: UCLA-4"),
         effect = factor(case_when(str_detect(effect, "t1") ~ "adolescence (17y)",
                                   str_detect(effect, "t2") ~ "emerging adulthod (22y)",
                                   str_detect(effect, "t3") ~ "young adulthood (28y)",
                                   str_detect(effect, "all") ~ "cumulative (17-28y)"),
                         levels = c("adolescence (17y)", "emerging adulthod (22y)", "young adulthood (28y)", "cumulative (17-28y)"),
                         ordered = T))


plot_combined_antidepress <- results_combined_antidepressants %>% 
  ggplot(aes(x = estimate, y = measure, color = estimation)) +
  geom_point(size = 4, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(xmin = `2.5 %`, xmax = `97.5 %`), 
                position = position_dodge(width = 0.3),
                width = 0.15, size = 1) +
  lemon::facet_rep_wrap(~ effect, 
                        ncol = 4,
                        repeat.tick.labels = FALSE) +
  labs(title = "Loneliness and risk of prescriptions for antidepressants in young adulthood",
       color = "Estimation",
       x = "Hazard ratio (95% confidence interval)",
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

ggsave(filename = "plot_combined_antidepress.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 13.8, 
       height = 5,  
       bg="white",
       dpi=700)

results_combined_antidepressants <- results_combined_antidepressants %>% 
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
write.csv(results_combined_antidepressants, "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Tables/results_combined_antidepressants")



#########################################################
# Estimating and fitting survival curves (unadjusted) ###
#########################################################

Cox_direct_antidepress <- coxph(Surv(time, event) ~ loneliness_4, 
                                data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "direct",]$data[[1]], x = TRUE)
Cox_indirect_antidepress <- coxph(Surv(time, event) ~ loneliness_4, 
                                  data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "indirect",]$data[[1]], x = TRUE)

Cox_direct_antidepress_adj <- coxph(Surv(time, event) ~ loneliness_4 + gender + age + ethnicity + parental_education + urbanity + 
                                  warm_parenting + parental_alcoholuse + parental_smoking + 
                                  asthma_allergy + phys_disability +
                                  social_support_3 + depression_3 + 
                                  living_situation_3 + relationship_3 + friends_3 + employment_3 +
                                  smoking_3 + alcohol_use_3 + illicit_drug_3 + phys_exer_3 + bmi_3 +
                                  loneliness_3, 
                                data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "direct",]$data[[1]], x = TRUE)
Cox_indirect_antidepress_adj <- coxph(Surv(time, event) ~ loneliness_4 + gender + age + ethnicity + parental_education + urbanity + 
                                    warm_parenting + parental_alcoholuse + parental_smoking + 
                                    asthma_allergy + phys_disability +
                                    social_support_3 + depression_3 + 
                                    living_situation_3 + relationship_3 + friends_3 + employment_3 +
                                    smoking_3 + alcohol_use_3 + illicit_drug_3 + phys_exer_3 + bmi_3 +
                                    loneliness_3, 
                                  data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "indirect",]$data[[1]], x = TRUE)
# without adjustment and non-response weighting
windowsFonts(Times = windowsFont("Times New Roman"))

# Survival curves
surv_curve_direct <- plot_surv_lines(time ="time", status = "event", 
                                     variable = "loneliness_4", 
                                     data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "direct",]$data[[1]],
                                     model = Cox_direct_antidepress_adj, horizon=c(0, 1)) +
  labs(x = "Time (in years)",
       y = "Probability of remaining anti-depressant free",
       color = "Loneliness \n (std)",
       title = "Direct measure: 'I feel lonely'") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14,
                                  face = "bold",
                                  margin = margin(b = 20, t = 30)),
        text = element_text(size = 14, family = "Times"))

surv_curve_indirect <-plot_surv_lines(time="time", status="event", 
                                      variable= "loneliness_4", 
                                      data= nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "indirect",]$data[[1]], 
                                      model= Cox_indirect_antidepress_adj, horizon=c(0, 1)) +
  labs(x = "Time (in years)",
       y = "",
       color = "Loneliness \n (std)",
       title = "Indirect measure: UCLA-4") +
  theme(plot.title = element_text(hjust = 0.5, 
                                  size = 14,
                                  face = "bold",
                                  margin = margin(b = 20, t = 30)),
        text = element_text(size = 14, family = "Times"))

# Survival probability at specific time-points of follow-up
surv_curve_t_direct <- plot_surv_at_t(time = "time",
                                      status = "event",
                                      variable = "loneliness_4",
                                      data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "direct",]$data[[1]],
                                      model = Cox_direct_antidepress_adj,
                                      t = 1:10) +
  scale_x_continuous(limits = c(-1.1, 2.1)) +
  scale_y_continuous(limits = c(0.80, 1.00)) +
  labs(x = "Loneliness (std)",
       y = "Probability of remaining anti-depressant free", # survival probability
       color = "Time \n (in years)",
       title = "") +
  theme(text = element_text(size = 14, family = "Times"))

surv_curve_t_indirect <- plot_surv_at_t(time = "time",
                                        status = "event",
                                        variable = "loneliness_4",
                                        data = nonadj_models_psychotropic[nonadj_models_psychotropic$measure == "indirect",]$data[[1]],
                                        model = Cox_indirect_antidepress_adj,
                                        t = 1:10) +
  scale_x_continuous(limits = c(-1.1, 2.1)) +
  scale_y_continuous(limits = c(0.80, 1.00)) +
  labs(x = "Loneliness (std)",
       y = "", # survival probability
       color = "Time \n (in years)",
       title = "") +
  theme(text = element_text(size = 14, family = "Times"))


combined_surv_curves <- plot_grid(surv_curve_direct, surv_curve_indirect, 
                                  surv_curve_t_direct, surv_curve_t_indirect,
                                  ncol = 2) +
  draw_label("Adjusted survival curves for antidepressant use by loneliness levels in young adulthood (28y)", 
             fontface = "bold",
             x = 0.5, hjust = 0.5,   # horizontal centering
             y = 0.995, vjust = 0.995,
             size = 16, fontfamily = "Times")

ggsave(filename = "combined_surv_curves.jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 18, 
       height = 10,  
       bg="white",
       dpi=700)

combined_antidep_MSM_surv_curves <- plot_grid(plot_combined_antidepress, NULL, combined_surv_curves, 
                                              ncol = 1,
                                              rel_heights = c(1, 0.1, 1.5))

ggsave(filename = "combined_antidep_MSM_surv_curves .jpeg",
       path = "N:/durable/Data_analyses/Libor/Chronic_Loneliness_Health_Adulthood/Results/Graphs", 
       width = 15, 
       height = 15,  
       bg="white",
       dpi=900)








##############################
### 1-YEAR WASHOUT PERIOD ###
##############################
# excluding individuals with antidepressant use in 2004

NorPD_UiN_antidep_washout <- select(NorPD_UiN_analysis_dep_imp, -starts_with(c("event", "time"))) %>% 
  left_join(., NorPD_antidep_washout, by = "id") %>% 
  mutate(event = ifelse(is.na(event), 0, event),
         time = ifelse(is.na(time), 9.33, time)) %>% 
  filter(no_preexp_med == TRUE|is.na(no_preexp_med)) # 216 individuals with antidepressants (after excluding 42 with pre-exposure use)


weighted_dataset_washout <- NorPD_UiN_antidep_washout %>% 
  group_by(measure, .imp) %>% 
  nest() %>% 
  mutate(weighted_data = map(.x = data,
                             ~ weightitMSM(exposure_models,
                                           data = .x,
                                           s.weights = .x$nonresp_weights,
                                           method = "cbps", 
                                           stabilize = TRUE)))

weighted_data_analysis_washout <- select(weighted_dataset_washout, .imp, weighted_data) 


adj_models_psychotropic_washout <- unnest(weighted_dataset_washout, data) %>% 
  select(id, nonresp_weights, loneliness_2:loneliness_4, event, time) %>% 
  group_by(.imp, measure) %>% 
  nest() %>% 
  left_join(., weighted_data_analysis_washout, by = c(".imp", "measure")) %>% 
  mutate(longi_model = map2(.x = data, .y = weighted_data, 
                            ~ coxph_weightit(survival::Surv(time, event) ~ loneliness_2 + loneliness_3 + loneliness_4, 
                                             vcov = "HC0",
                                             data = .x, weightit = .y, weights = .x$nonresp_weights)),
         avg_predictions = map(.x = longi_model, .y = data,
                               ~ avg_predictions(.x, 
                                                 type = "risk",
                                                 vcov = "HC0", 
                                                 variables = list(loneliness_2 = c(0,1), loneliness_3 = c(0,1), loneliness_4 = c(0,1)), 
                                                 by = c("loneliness_2", "loneliness_3","loneliness_4"))),
         avg_comp_all = map(.x = avg_predictions, .y = data,
                            ~ hypotheses(.x, hypothesis =  "b8/b1 = 0", type = "risk", comparison = "ratio", wts = .y$nonresp_weights)))  

MSM_results_psychotropic_washout <- select(adj_models_psychotropic_washout, starts_with("avg_comp")) %>%
  group_by(measure) %>% 
  nest() %>% 
  mutate(contrasts = map(.x = data, 
                         ~ summary(mice::pool(.x$avg_comp_all), conf.int = TRUE))) %>% # pooling according to Rubin's rules
  select(contrasts) %>% unnest() %>% 
  select(measure, estimate, `2.5 %`, `97.5 %`)

