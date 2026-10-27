#### Surgery Receipt Likelihood for Pre-/Post-Expansion ####
# whole population # 

seer_nonmet <- seer_nonmet%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )


surg_data <- seer_nonmet %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    # Time_Period  = factor(Time_Period, levels = c("Pre","Post")),
    surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0)
  )

surgery_did_model <- glm(
  surgery_binary ~ Treatment_ME * Time_Period,
  data = surg_data,
  family = binomial()
)

# Model summaries
coeftest(surgery_did_model, vcov = sandwich)
summary(surgery_did_model)

# Export tidy output
write.csv(
  tidy(surgery_did_model, conf.int = TRUE, exponentiate = TRUE),
  "nonmet/surgery/glm_model_surgery.csv",
  row.names = FALSE
)


#### Survival  Likelihood for Pre-/Post-Expansion ####
# 1-1 analysis, do expansion patients survive longer? 
# can be an preceding analysis to the group-specific exp data
# matched population # 

combined_matched_data <- combined_matched_data%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )

surv_1 <- combined_matched_data %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    # Time_Period  = factor(Time_Period, levels = c("Pre","Post")),
    surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0)
  )


surv_1$capped_time <- pmin(surv_1$Survival.months, 36)
surv_1$capped_event <- ifelse(surv_1$Survival.months > 36, 0, surv_1$Event)



capped_overall_did <- coxph(Surv(capped_time, capped_event) ~ Time_Period * Treatment_ME, # ,
                            data = surv_1, robust = TRUE)


summary(capped_overall_did)

capped_overall_did_tidy <- tidy(capped_overall_did, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

write.csv(capped_overall_did_tidy, "nonmet/DiD_survival_capped.csv")

##
overall_did <- coxph(Surv(Survival.months, Event) ~ Time_Period * Treatment_ME, # ,
                     data = surv_1, robust = TRUE)


summary(overall_did)

overall_did_tidy <- tidy(overall_did, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

write.csv(overall_did_tidy, "nonmet/DiD_survival_full.csv")


cox.zph(overall_did)
cox.zph(capped_overall_did)

#### Kaplan-Meier Curves ####

seer_nonmet <- seer_nonmet%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )


seer_nonmet <- seer_nonmet %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    # Time_Period  = factor(Time_Period, levels = c("Pre","Post")),
    surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0)
  )


seer_nonmet$capped_time <- pmin(seer_nonmet$Survival.months, 36)
seer_nonmet$capped_event <- ifelse(seer_nonmet$Survival.months > 36, 0, seer_nonmet$Event)


seer_nonmet <- seer_nonmet %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    # Time_Period  = factor(Time_Period, levels = c("Pre","Post")),
    surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0)
  )

# Whole population expansion vs non-expansion patients

# Log Rank Testing
logrank_test <- survdiff(Surv(capped_time, capped_event) ~ Treatment_ME, data = seer_nonmet)

# Extract key components from the survdiff object
logrank_df <- data.frame(
  Group = names(logrank_test$n),           # Group names
  N = as.numeric(logrank_test$n),          # Number of patients in each group
  Observed = logrank_test$obs,             # Observed events
  Expected = logrank_test$exp,             # Expected events
  `O_minus_E^2_over_E` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$exp,
  `O_minus_E^2_over_V` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$var
)

# Add the Chi-squared statistic and p-value as separate columns
logrank_df$Chisq <- logrank_test$chisq
logrank_df$df <- length(logrank_test$n) - 1
logrank_df$p.value <- 1 - pchisq(logrank_test$chisq, df = logrank_df$df[1])

# Export to CSV
write.csv(logrank_df, "nonmet/logrank_allpts_ME.csv", row.names = FALSE)

# KM Fit
fit_ME <- survfit(Surv(capped_time, capped_event) ~ Treatment_ME, data = seer_nonmet)

# To get approximate p-value
p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
p_value

fit_ME_plot <- ggsurvplot(fit_ME, conf.int = TRUE,
                          risk.table = TRUE,
                          pval = TRUE,
                          xlim = c(0, 36),
                          ylim = c(0.4, 1),
                          break.time.by = 6,
                          xlab = "Months",
                          ylab = "Percent Alive",
                          title = "Kaplan-Meier: Patients Exposed to Medicaid Expansion vs Non-Expansion",
                          palette = "Dark2"
)
# If you want the full ggsurvplot (survival + risk table), use cowplot::ggsave or patchwork
png("nonmet/exp_v_nonexp_KM.png", width = 3000, height = 2000, res = 300)
print(fit_ME_plot)   # must explicitly print for base graphics device
dev.off()



fit_groups <- survfit(Surv(Survival.months, capped_event) ~ Medicaid_Expansion_Status, data = seer_nonmet)
km_ME <- ggsurvplot(fit_groups, conf.int = TRUE,
                    risk.table = TRUE,
                    pval = TRUE,
                    xlim = c(0, 36),
                    ylim = c(0.3, 1),
                    break.time.by = 6,
                    xlab = "Months",
                    ylab = "Percent Alive",
                    title = "Kaplan-Meier: Expansion Timing")


png("nonmet/km_all_groups.png", width = 2100, height = 2000, res = 300)
print(km_ME)   # must explicitly print for base graphics device
dev.off()

## Group Site.recode into HN, Lung, Esophagus
cleaned_data$Site.recode.ICD.O.3.2023.Revision <- factor(cleaned_data$Site.recode.ICD.O.3.2023.Revision)  
cleaned_data <- cleaned_data %>%
  mutate(Site = case_when(
    (Site.recode.ICD.O.3.2023.Revision == "Esophagus") ~ "Esophagus",
    (Site.recode.ICD.O.3.2023.Revision == "Floor of Mouth") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Gum") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Hypopharynx") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Larynx") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Lip") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Nasopharynx") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Oropharynx") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Tongue Anterior") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Palate Excluding Soft And Uvula") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Pharynx And Oral Cavity Other") ~ "H&N",
    (Site.recode.ICD.O.3.2023.Revision == "Lung And Bronchus") ~ "Lung",
    (Site.recode.ICD.O.3.2023.Revision == "Trachea, And Respiratory Other") ~ "Lung")
  )



# cleaned_data <- cleaned_data %>%
#   mutate(
#     Site = case_when(
#       Site.recode.ICD.O.3.2023.Revision == "Esophagus" ~ "Esophagus",
#       Site.recode.ICD.O.3.2023.Revision %in% c(
#         "Floor of Mouth", "Gum", "Hypopharynx", "Larynx",
#         "Lip", "Nasopharynx", "Oropharynx", "Tongue Anterior",
#         "Palate Excluding Soft And Uvula", "Pharynx And Oral Cavity Other"
#       ) ~ "H&N",
#       Site.recode.ICD.O.3.2023.Revision %in% c(
#         "Lung And Bronchus", "Trachea, And Respiratory Other"
#       ) ~ "Lung",
#       TRUE ~ NA_character_   # everything else gets NA
#     )
#   ) 


# Disease-specific KM (lung, esophageal, H&N)
fit_disease <- survfit(Surv(Survival.months, Event) ~ Site, data = cleaned_data)

km_dztype <- ggsurvplot(fit_disease, conf.int = TRUE,
                        risk.table = TRUE,
                        pval = TRUE,
                        xlim = c(0, 36),
                        ylim = c(0.2, 1),
                        break.time.by = 6,
                        xlab = "Months",
                        ylab = "Percent Alive",
                        title = "Kaplan-Meier: Expansion Timing",
                        palette = "Dark2"
)

# If you want the full ggsurvplot (survival + risk table), use cowplot::ggsave or patchwork
png("nonmet/km_dztype.png", width = 2100, height = 2000, res = 300)
print(km_dztype)   # must explicitly print for base graphics device
dev.off()


## Lung KM (Non-Exp vs Exp) 

# recreate the true expansion criteria (intersection between diagnosis before/after ME and state-identity)
seer_nonmet_lung <- seer_nonmet_lung%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )
seer_nonmet_lung <- seer_nonmet_lung %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    Time_Period  = factor(Time_Period, levels = c("Pre","Post"))
  )

fit_lung <- survfit(Surv(Survival.months, Event) ~ Time_Period, data = seer_nonmet_lung)

# Log Rank Testing
logrank_test <- survdiff(Surv(capped_time, capped_event) ~ Treatment_ME, data = seer_nonmet_lung)

# Extract key components from the survdiff object
logrank_df <- data.frame(
  Group = names(logrank_test$n),           # Group names
  N = as.numeric(logrank_test$n),          # Number of patients in each group
  Observed = logrank_test$obs,             # Observed events
  Expected = logrank_test$exp,             # Expected events
  `O_minus_E^2_over_E` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$exp,
  `O_minus_E^2_over_V` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$var
)

# Add the Chi-squared statistic and p-value as separate columns
logrank_df$Chisq <- logrank_test$chisq
logrank_df$df <- length(logrank_test$n) - 1
logrank_df$p.value <- 1 - pchisq(logrank_test$chisq, df = logrank_df$df[1])

# Export to CSV
write.csv(logrank_df, "nonmet/lung/logrank_allpts_ME.csv", row.names = FALSE)

km_lung <- ggsurvplot(fit_lung, conf.int = TRUE,
                      risk.table = TRUE,
                      pval = TRUE,
                      xlim = c(0, 36),
                      ylim = c(0.2, 1),
                      break.time.by = 6,
                      xlab = "Months",
                      ylab = "Percent Alive",
                      title = "Lung Cancer Survival",
                      palette = "Dark2"
)
# Relabel x-axis from months to years
km_surg$plot <- km_surg$plot + scale_x_continuous(breaks = seq(0, 36, 12),
                                                  labels = 0:3)

km_surg$table <- km_surg$table + scale_x_continuous(breaks = seq(0, 36, 12),
                                                    labels = 0:3)

# If you want the full ggsurvplot (survival + risk table), use cowplot::ggsave or patchwork
png("nonmet/lung/km_survival_ME.png", width = 2800, height = 2000, res = 300)
print(km_lung)   # must explicitly print for base graphics device
dev.off()


## Head and Neck KM (Non-Exp vs Exp) 

# recreate the true expansion criteria (intersection between diagnosis before/after ME and state-identity)
seer_nonmet_HN <- seer_nonmet_HN%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )
seer_nonmet_HN <- seer_nonmet_HN %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    Time_Period  = factor(Time_Period, levels = c("Pre","Post"))
  )

fit_hn <- survfit(Surv(Survival.months, Event) ~ Time_Period, data = seer_nonmet_HN)

# Log Rank Testing
logrank_test <- survdiff(Surv(capped_time, capped_event) ~ Treatment_ME, data = seer_nonmet_HN)

# Extract key components from the survdiff object
logrank_df <- data.frame(
  Group = names(logrank_test$n),           # Group names
  N = as.numeric(logrank_test$n),          # Number of patients in each group
  Observed = logrank_test$obs,             # Observed events
  Expected = logrank_test$exp,             # Expected events
  `O_minus_E^2_over_E` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$exp,
  `O_minus_E^2_over_V` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$var
)

# Add the Chi-squared statistic and p-value as separate columns
logrank_df$Chisq <- logrank_test$chisq
logrank_df$df <- length(logrank_test$n) - 1
logrank_df$p.value <- 1 - pchisq(logrank_test$chisq, df = logrank_df$df[1])

# Export to CSV
write.csv(logrank_df, "nonmet/hn/logrank_allpts_ME.csv", row.names = FALSE)

km_hn <- ggsurvplot(fit_hn, conf.int = TRUE,
                    risk.table = TRUE,
                    pval = TRUE,
                    xlim = c(0, 36),
                    ylim = c(0.7, 1),
                    break.time.by = 6,
                    xlab = "Months",
                    ylab = "Percent Alive",
                    title = "Head and Neck SCC Survival",
                    palette = "Dark2"
)

# If you want the full ggsurvplot (survival + risk table), use cowplot::ggsave or patchwork
png("nonmet/hn/km_survival_ME.png", width = 2800, height = 2000, res = 300)
print(km_hn)   # must explicitly print for base graphics device
dev.off()



## Lung KM (Non-Exp vs Exp) 

# recreate the true expansion criteria (intersection between diagnosis before/after ME and state-identity)
seer_nonmet_esoph <- seer_nonmet_esoph%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )
seer_nonmet_esoph <- seer_nonmet_esoph %>%
  mutate(
    # Treatment is time-invariant: did the state ever expand?
    Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
    
    # Time is relative to each state's expansion year
    Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
  ) %>%
  mutate(
    Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    Time_Period  = factor(Time_Period, levels = c("Pre","Post"))
  )

fit_esophagus <- survfit(Surv(Survival.months, Event) ~ Time_Period, data = seer_nonmet_esoph)

# Log Rank Testing
logrank_test <- survdiff(Surv(capped_time, capped_event) ~ Treatment_ME, data = seer_nonmet_esoph)

# Extract key components from the survdiff object
logrank_df <- data.frame(
  Group = names(logrank_test$n),           # Group names
  N = as.numeric(logrank_test$n),          # Number of patients in each group
  Observed = logrank_test$obs,             # Observed events
  Expected = logrank_test$exp,             # Expected events
  `O_minus_E^2_over_E` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$exp,
  `O_minus_E^2_over_V` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$var
)

# Add the Chi-squared statistic and p-value as separate columns
logrank_df$Chisq <- logrank_test$chisq
logrank_df$df <- length(logrank_test$n) - 1
logrank_df$p.value <- 1 - pchisq(logrank_test$chisq, df = logrank_df$df[1])

# Export to CSV
write.csv(logrank_df, "nonmet/esophagus/logrank_allpts_ME.csv", row.names = FALSE)

km_esophagus <- ggsurvplot(fit_esophagus, conf.int = TRUE,
                           risk.table = TRUE,
                           pval = TRUE,
                           xlim = c(0, 36),
                           ylim = c(0.2, 1),
                           break.time.by = 6,
                           xlab = "Months",
                           ylab = "Percent Alive",
                           title = "Esophagus SCC Survival",
                           palette = "Dark2"
)

# If you want the full ggsurvplot (survival + risk table), use cowplot::ggsave or patchwork
png("nonmet/esophagus/km_survival_ME.png", width = 2800, height = 2000, res = 300)
print(km_esophagus)   # must explicitly print for base graphics device
dev.off()

#### Race.Ethnicity Absolute Survival Differences ####

# Restrict to patients with at least 3 years follow-up
seer_3yr <- seer_nonmet %>% mutate(capped_event = ifelse(Survival.months > 36, 0, Event),
                                  capped_time = pmin(Survival.months, 36))

# Function to compute 2-year survival by group
library(dplyr)
library(survminer)

race_groups <- seer_3yr %>% group_by(Race.Ethnicity, Time_Period) %>%
  summarize(
    surv_fit = list(survfit(Surv(capped_time, capped_event) ~ 1)),
    .groups = "drop"
  )

# Extract 36-month survival
race_surv <- race_groups %>%
  mutate(
    surv_36mo = map_dbl(surv_fit, ~ {
      s <- summary(.x, times = 36)
      if ("surv" %in% names(s)) s$surv else NA_real_
    })
  )

# Compute absolute differences between White and other races
race_wide <- race_surv %>%
  dplyr::select(Race.Ethnicity, Time_Period, surv_36mo) %>%
  tidyr::pivot_wider(names_from = Race.Ethnicity, values_from = surv_36mo)

race_wide <- race_wide %>% mutate(
  Black_vs_White = `Non-Hispanic Black` - `Non-Hispanic White`,
  Asian_vs_White = `Non-Hispanic Asian or Pacific Islander` - `Non-Hispanic White`,
  AIAN_vs_White = `Non-Hispanic American Indian/Alaska Native` - `Non-Hispanic White`,
  Hispanic_vs_White = `Hispanic (All Races)` - `Non-Hispanic White`
)

print(race_wide)

write.csv(race_wide, "nonmet/3yr/race-wide_comparisons.csv", row.names = FALSE)

#### Survival for Post-Operative Patients ####
# Surgical patients only
surg_patients <- seer_nonmet %>%
  filter(surgery_status == "Recommended + performed")

# Log Rank Testing
logrank_test <- survdiff(Surv(capped_time, capped_event) ~ Treatment_ME, data = surg_patients)

# Extract key components from the survdiff object
logrank_df <- data.frame(
  Group = names(logrank_test$n),           # Group names
  N = as.numeric(logrank_test$n),          # Number of patients in each group
  Observed = logrank_test$obs,             # Observed events
  Expected = logrank_test$exp,             # Expected events
  `O_minus_E^2_over_E` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$exp,
  `O_minus_E^2_over_V` = (logrank_test$obs - logrank_test$exp)^2 / logrank_test$var
)

# Add the Chi-squared statistic and p-value as separate columns
logrank_df$Chisq <- logrank_test$chisq
logrank_df$df <- length(logrank_test$n) - 1
logrank_df$p.value <- 1 - pchisq(logrank_test$chisq, df = logrank_df$df[1])

# Export to CSV
write.csv(logrank_df, "nonmet/logrank_surgpts_ME.csv", row.names = FALSE)

# Kaplan-Meier by Time_Period
fit_surg <- survfit(Surv(Survival.months, Event) ~ Time_Period, data = surg_patients)

km_surg <- ggsurvplot(fit_surg, conf.int = TRUE,
                      risk.table = TRUE,
                      pval = TRUE,
                      xlim = c(0, 36),
                      ylim = c(0.7, 1),
                      break.time.by = 6,
                      xlab = "Months",
                      ylab = "Overall Survival",
                      title = "Survival Pre- vs Post-Expansion for Surgical Patients",
                      palette = c("blue", "red"))

png("nonmet/surgery_patients_KM.png", width = 2100, height = 2000, res = 300)
print(km_surg)   # must explicitly print for base graphics device
dev.off()

## should this be a cox DiD model..?

# Check reference levels
contrasts(surg_patients$Time_Period)
contrasts(surg_patients$Treatment_ME)

# Step 4: Run Cox regression model
surg_survival <- coxph(Surv(Survival.months, Event) ~ Time_Period * Treatment_ME ,
                       data = surg_patients, robust = TRUE)

# Step 5: Summary of the model
summary(surg_survival)

surg_survival_tidy <- tidy(surg_survival, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(surg_survival_tidy)

write.csv(surg_survival_tidy, "nonmet/did_surgicalpt_survival.csv")

cox.zph(surg_survival)


### Three-Year Survival Cox Endpoint ### 

#  Run Cox regression model
surg_patients$capped_time <- pmin(surg_patients$Survival.months, 40)
surg_patients$capped_event <- ifelse(surg_patients$Survival.months > 40, 0, surg_patients$Event)

surg_survival_capped <- coxph(Surv(capped_time, capped_event) ~ Time_Period * Treatment_ME ,
                              data = surg_patients, robust = TRUE)

# Step 5: Summary of the model
summary(surg_survival_capped)

surg_survival_capped_tidy <- tidy(surg_survival_capped, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(surg_survival_capped_tidy)

write.csv(surg_survival_capped_tidy, "nonmet/did_surgicalpt_survival_capped36.csv")

cox.zph(surg_survival_capped)



#### Placebo Falsification Models ####

# Subset pre-expansion years only
placebo_data_1 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2006-2010"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 1)) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0)
  )

# Run placebo DiD model
placebo_model <- lm(capped_event ~ treat * fake_time, data = placebo_data_1)
summary(placebo_model)

library(modelsummary)
library(pandoc)
modelsummary(placebo_model, output = "nonmet/placebo_survival_model_output.docx")

# Output
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_survival_1.csv", row.names = FALSE)

## survival for group 2 ##
# Subset pre-expansion years only
placebo_data_2 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2006-2010","2011-2013"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 2)) %>%
  mutate(
    fake_time = ifelse(dx.yr_grouped == "2011-2013", 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 2, 1, 0)
  )

# Run placebo DiD model
placebo_model <- lm(capped_event ~ treat * fake_time, data = placebo_data_2)
summary(placebo_model)

modelsummary(placebo_model, output = "nonmet/placebo_survival_model_output_2.docx")

# Output
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_survival_2.csv", row.names = FALSE)


## survival for group 3 ##
# Subset pre-expansion years only #
placebo_data_3 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2014-2016"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 3)) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2016, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 3, 1, 0)
  )

# Run placebo DiD model
placebo_model <- lm(capped_event ~ treat * fake_time, data = placebo_data_3)
summary(placebo_model)

modelsummary(placebo_model, output = "nonmet/placebo_survival_model_output_3.docx")

# Output #
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_survival_3.csv", row.names = FALSE)


## Surgery Reciept ##

# Subset pre-expansion years only
surgery_data_1 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2006-2010"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 1)) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2010, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0)
  )

# Run placebo DiD model
surg_placebo <- glm(surgery_binary ~ treat * fake_time, data = surgery_data_1)
summary(surg_placebo)

modelsummary(surg_placebo, output = "nonmet/placebo_surgery_model_output_1.docx")

# Output
tidy(placebo_model, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(surg_placebo, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_surgery_1.csv", row.names = FALSE)


## 2014 Expansion Surgery ##
# Subset pre-expansion years only
surgery_data_2 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2011-2013"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 2)) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2013, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 2, 1, 0)
  )

# Run placebo DiD model
surg_placebo <- glm(surgery_binary ~ treat * fake_time, data = surgery_data_2)
summary(surg_placebo)

modelsummary(surg_placebo, output = "nonmet/placebo_surgery_model_output_2.docx")

# Output
tidy(surg_placebo, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(surg_placebo, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_surgery_2.csv", row.names = FALSE)

## Late Expansion Surgery ##
# Subset pre-expansion years only
surgery_data_3 <- seer_nonmet %>%
  filter(dx.yr_grouped %in% c("2014-2016"),   # Ensure dx.yr_grouped is treated as a string
         Medicaid_Expansion_Status %in% c(0, 3)) %>%
  mutate(
    fake_time = ifelse(dx.yr == 2016, 1, 0),  # Pretend 2010 is "post"
    treat = ifelse(Medicaid_Expansion_Status == 3, 1, 0)
  )

# Run placebo DiD model
surg_placebo <- glm(surgery_binary ~ treat * fake_time, data = surgery_data_3)
summary(surg_placebo)

modelsummary(surg_placebo, output = "nonmet/placebo_surgery_model_output_2.docx")

# Output
tidy(surg_placebo, conf.int = TRUE, exponentiate = TRUE)


# Assuming you have already run tidy() on the placebo model
placebo_tidy <- tidy(surg_placebo, conf.int = TRUE, exponentiate = TRUE)

# Save the tidy result as a CSV
write.csv(placebo_tidy, "nonmet/placebofalse_surgery_3.csv", row.names = FALSE)







