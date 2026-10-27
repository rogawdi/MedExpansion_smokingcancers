#### Cox Proportional Hazards Model for Full Follow-up Survival ####
cleaned_data <- cleaned_data %>%
  filter(!grepl("Unknown Race.Ethnicity", trimws(Race.Ethnicity)))

#Combine previous 3 Matched Data tables into 1 group
# Combine the matched datasets and add a group identifier
matched_data_1$group <- "Group 1"
matched_data_2$group <- "Group 2"
matched_data_3$group <- "Group 3"

# # Combine all matched datasets into one
# combined_matched_data <- bind_rows(matched_data_1, matched_data_2, matched_data_3)
# 
# 
# combined_matched_data$capped_time <- pmin(combined_matched_data$Survival.months, 36)
# combined_matched_data$capped_event <- ifelse(combined_matched_data$Survival.months > 36, 0, combined_matched_data$Event)
# combined_matched_data <- combined_matched_data[combined_matched_data$capped_time > 0, ]
# 
# cox_psm_final <- coxph(
#   Surv(Survival.months, Event) ~ age_category + Sex + Marital_Grouped +
#     Rural.Urban.Continuum.Code + Race.Ethnicity + Medicaid_Expansion_Status,
#   # cluster(matched_id), 
#   data = combined_matched_data,
#   robust = TRUE  # Adds robust standard errors
# )
# 
# capped_fin_cox_psm_final_sum <- tidy(cox_psm_final, conf.int = TRUE) %>%
#   mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
#          CI_lower = exp(conf.low),
#          CI_upper = exp(conf.high)) %>%
#   dplyr::select(term, HR, CI_lower, CI_upper, p.value)
# 
# # # View(capped_fin_cox_psm_final_sum)
# 
# write.csv(capped_fin_cox_psm_final_sum, "nonmet/coxdata_full/psm_covariates.csv")

# Group by Group Difference in Differences analysis
# 0 vs 1

# Create expansion1_time variable (Pre = 2006-2010, Post = 2011-2019)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_time = case_when(
    dx.yr_grouped %in% c("2006-2010") ~ "Pre",
    dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion1 states, 0 = non-expansion states)
matched_data_1 <- matched_data_1 %>%
  mutate(expansion1_group = case_when(
    Medicaid_Expansion_Status == 1 ~ 1,  # Expansion 1 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))

# Subset to only include subjects from Expansion Group 0 and 1
matched_data_1_did <- matched_data_1 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 1)) %>%
  # Create a binary treatment variable:
  # did_treat = 1 if Expansion Group 1, 0 if Expansion Group 0
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0),
         # Create a time variable: "Pre" for 2006-2010, "Post" for later years.
         did_time = ifelse(dx.yr_grouped == "2006-2010", "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_1_did <- matched_data_1_did %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
  )

matched_data_1_did$capped_time <- pmin(matched_data_1_did$Survival.months, 36)
matched_data_1_did$capped_event <- ifelse(matched_data_1_did$Survival.months > 36, 0, matched_data_1_did$Event)

n_pairs <- nrow(matched_data_1_did) / 2
matched_data_1_did$matched_id <- rep(1:n_pairs, each = 2)
matched_data_1_did$matched_id <- matched_data_1$subclass


capped_fin_cox_psm_final_did_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                                        data = matched_data_1_did, robust = TRUE)


summary(capped_fin_cox_psm_final_did_1)

capped_fin_cox_psm_final_did_1 <- tidy(capped_fin_cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_1

write.csv(capped_fin_cox_psm_final_did_1, "nonmet/coxdata_full/psm_0v1.csv")



# EARLY - consider early months (post-implementation periods)
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_early_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_early_1 <- matched_data_1_did_early_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_early_1$did_time)
contrasts(matched_data_1_did_early_1$did_treat)

matched_data_1_did_early_1$matched_id <- as.numeric(as.factor(
  matched_data_1_did_early_1$subclass)
)

n_pairs <- nrow(matched_data_1_did_early_1) / 2
matched_data_1_did_early_1$matched_id <- rep(1:n_pairs, each = 2)
matched_data_1_did_early_1$matched_id <- matched_data_1_did_early_1$subclass

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat,
                                          data = matched_data_1_did_early_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_1)

capped_fin_cox_psm_final_did_early_1 <- tidy(capped_fin_cox_psm_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_early_1

write.csv(capped_fin_cox_psm_final_did_early_1, "nonmet/coxdata_full/psm_0v1_early.csv")


# LATER
# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016", "2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

matched_data_1_did_later_1$matched_id <- as.numeric(as.factor(
  matched_data_1_did_later_1$subclass)
)

n_pairs <- nrow(matched_data_1_did_later_1) / 2
matched_data_1_did_later_1$matched_id <- rep(1:n_pairs, each = 2)
matched_data_1_did_later_1$matched_id <- matched_data_1_did_later_1$subclass

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_later_1

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_full/psm_0v1_later.csv")

### LAST ###

# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_full/psm_0v1_LAST.csv")


### Mid Group ###

# Step 1: Subset data to focus on 2014-2019
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_full/psm_0v1_MID.csv")



##### LATE + MID #
# Step 1: Subset data to focus on 2011-2016
matched_data_1_did_later_1 <- matched_data_1 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_time = ifelse(dx.yr_grouped %in% c("2011-2013", "2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(did_treat = ifelse(Medicaid_Expansion_Status == 1, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_1_did_later_1 <- matched_data_1_did_later_1 %>%
  mutate(
    did_time = factor(did_time, levels = c("Pre", "Post")),
    did_treat = factor(did_treat, levels = c(0, 1))
    
  )

# Check reference levels
contrasts(matched_data_1_did_later_1$did_time)
contrasts(matched_data_1_did_later_1$did_treat)

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_full/psm_0v1_2011-2016.csv")



# 0 vs 2
# Create expansion1_time variable (Pre = 2006-2013, Post = 2014-2019)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_time = case_when(
    dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
    dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
    TRUE ~ NA_character_  # Assign NA for years outside the study period
  ))

# Create expansion1_group variable (1 = expansion2 states, 0 = non-expansion states)
matched_data_2 <- matched_data_2 %>%
  mutate(expansion2_group = case_when(
    Medicaid_Expansion_Status == 2 ~ 1,  # Expansion 2 (early)
    Medicaid_Expansion_Status == 0 ~ 0,  # Non-expansion (control)
    TRUE ~ NA_real_  # Assign NA for other groups
  ))



# Subset to only include subjects from Expansion Group 0 and 2
matched_data_2_did <- matched_data_2 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 2)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0),
         # Define Pre (2006-2013) and Post (2014+)
         did_time_2 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_2_did <- matched_data_2_did %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

matched_data_2_did <- matched_data_2_did %>%
  filter(!grepl("Unknown Race.Ethnicity", trimws(Race.Ethnicity)))

matched_data_2_did$matched_id <- as.numeric(as.factor(
  matched_data_2_did$subclass)
)

n_pairs <- nrow(matched_data_2_did) / 2
matched_data_2_did$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2_did$matched_id <- matched_data_2_did$subclass




capped_fin_cox_psm_final_did_2 <- coxph(Surv(Survival.months, Event) ~ did_time_2 * did_treat_2 ,
                                        data = matched_data_2_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_2)

capped_fin_cox_psm_final_did_2 <- tidy(capped_fin_cox_psm_final_did_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_2

write.csv(capped_fin_cox_psm_final_did_2, "nonmet/coxdata_full/psm_0v2.csv")

# EARLY - consider early months (post-implementation periods)
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_early_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2014-2016"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_early_2 <- matched_data_2_did_early_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_early_2$did_time_2)
contrasts(matched_data_2_did_early_2$did_treat_2)

n_pairs <- nrow(matched_data_2_did_early_2) / 2
matched_data_2_did_early_2$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2_did_early_2$matched_id <- matched_data_2_did_early_2$subclass


# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_2 <- coxph(Surv(Survival.months, Event) ~ did_time_2 * did_treat_2 ,
                                          data = matched_data_2_did_early_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_2)

capped_fin_cox_psm_final_did_early_2 <- tidy(capped_fin_cox_psm_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_2)

write.csv(capped_fin_cox_psm_final_did_early_2, "nonmet/coxdata_full/psm_0v2_early.csv")




# LATER
# Step1: Subset data to focus on 2017-2019
matched_data_2_did_later_2 <- matched_data_2 %>%
  filter(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2017-2019"))

# Step 2: Update the time variable for 2014-2016 and 2017-2019 periods
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_time_2 = ifelse(dx.yr_grouped %in% c("2017-2019"), "Post", "Pre"))

# Create treatment variable (same as before)
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(did_treat_2 = ifelse(Medicaid_Expansion_Status == 2, 1, 0))

# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_later_2$did_time_2)
contrasts(matched_data_2_did_later_2$did_treat_2)

n_pairs <- nrow(matched_data_2_did_later_2) / 2
matched_data_2_did_later_2$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2_did_later_2$matched_id <- matched_data_2_did_later_2$subclass


# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_2 <- coxph(Surv(Survival.months, Event) ~ did_time_2 * did_treat_2 ,
                                          data = matched_data_2_did_later_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_2)

capped_fin_cox_psm_final_did_later_2 <- tidy(capped_fin_cox_psm_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_2)

write.csv(capped_fin_cox_psm_final_did_later_2, "nonmet/coxdata_full/psm_0v2_later.csv")



# 0 vs 3
# Subset to only include subjects from Expansion Group 0 and 3
matched_data_3_did <- matched_data_3 %>%
  filter(Medicaid_Expansion_Status %in% c(0, 3)) %>%
  # Assign treatment variable: 1 for expansion (2), 0 for non-expansion (0)
  mutate(did_treat_3 = ifelse(Medicaid_Expansion_Status == 3, 1, 0),
         # Define Pre (2006-2016) and Post (2017+)
         did_time_3 = ifelse(dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016"), "Pre", "Post"))

# Convert did_time and did_treat to factors with explicit ordering:
matched_data_3_did <- matched_data_3_did %>%
  mutate(
    did_time_3 = factor(did_time_3, levels = c("Pre", "Post")),
    did_treat_3 = factor(did_treat_3, levels = c(0, 1))
  )

matched_data_3_did <- matched_data_3_did %>%
  filter(!grepl("Unknown Race.Ethnicity", trimws(Race.Ethnicity)))

n_pairs <- nrow(matched_data_3_did) / 2
matched_data_3_did$matched_id <- rep(1:n_pairs, each = 2)
matched_data_3_did$matched_id <- matched_data_3_did$subclass

capped_fin_cox_psm_final_did_3 <- coxph(Surv(Survival.months, Event) ~ did_time_3 * did_treat_3 +
                                          cluster(matched_id),
                                        data = matched_data_3_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_3)

capped_fin_cox_psm_final_did_3 <- tidy(capped_fin_cox_psm_final_did_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_3

write.csv(capped_fin_cox_psm_final_did_3, "nonmet/coxdata_full/psm_0v3.csv")




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



#### Table 2 ####
library(dplyr)
library(survival)
library(rstpm2)

# Use seer_smoke to reproduce the original all-stage Table 2.
# Later, use seer_nonmet for the stage I-III replacement.
table2_data <- seer_nonmet %>%
  mutate(
    exp_yr = case_when(
      as.character(Medicaid_Expansion_Status) == "1" ~ 2011,
      as.character(Medicaid_Expansion_Status) == "2" ~ 2014,
      as.character(Medicaid_Expansion_Status) == "3" ~ 2017,
      TRUE ~ NA_real_
    ),
    
    # Expansion × Post indicator
    did_interaction = as.integer(
      !is.na(exp_yr) & dx.yr >= exp_yr
    ),
    
    Survival.months = ifelse(
      Survival.months <= 0,
      0.01,
      Survival.months
    ),
    
    age_category = relevel(
      factor(age_category),
      ref = "Middle-Aged"
    ),
    
    Sex = relevel(
      factor(Sex),
      ref = "Male"
    ),
    
    Race.Ethnicity = relevel(
      factor(Race.Ethnicity),
      ref = "Non-Hispanic White"
    ),
    
    Marital_Grouped = relevel(
      factor(Marital_Grouped),
      ref = "Married"
    ),
    
    Rural.Urban.Continuum.Code = relevel(
      factor(Rural.Urban.Continuum.Code),
      ref = "Counties in metropolitan areas ge 1 million pop"
    ),
    
    stage_group = relevel(
      factor(stage_group),
      ref = "I"
    )
  ) %>%
  filter(
    complete.cases(
      Survival.months,
      Event,
      did_interaction,
      dx.yr,
      age_category,
      Sex,
      Race.Ethnicity,
      Marital_Grouped,
      Rural.Urban.Continuum.Code,
      stage_group,
      Medicaid.Expansion.Status
    )
  ) %>%
  droplevels()


# Exact original Table 2 Royston-Parmar model
table2_model <- stpm2(
  Surv(Survival.months, Event) ~
    did_interaction +
    dx.yr +
    age_category +
    Sex +
    Race.Ethnicity +
    Marital_Grouped +
    Rural.Urban.Continuum.Code +
    stage_group,
  data = table2_data,
  df = 5,
  robust = TRUE,
  cluster = table2_data$Medicaid.Expansion.Status
)

summary(table2_model)

# Extract robust model results
b  <- coef(table2_model)
se <- sqrt(diag(vcov(table2_model)))

rp_tidy <- data.frame(
  term = names(b),
  HR = exp(b),
  CI_lower = exp(b - 1.96 * se),
  CI_upper = exp(b + 1.96 * se),
  p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE),
  row.names = NULL
)


# Formatting helpers
header_row <- function(label) {
  data.frame(
    Variable = label,
    HR = "",
    `95% CI` = "",
    `p-value` = "",
    check.names = FALSE
  )
}

result_row <- function(term_name, label) {
  
  x <- rp_tidy[rp_tidy$term == term_name, ]
  
  if (nrow(x) != 1) {
    stop(paste("Could not find model term:", term_name))
  }
  
  p_display <- if (
    x$p_value < 0.001
  ) {
    "<0.001"
  } else if (
    x$p_value < 0.10
  ) {
    sprintf("%.3f", x$p_value)
  } else {
    sprintf("%.2f", x$p_value)
  }
  
  data.frame(
    Variable = label,
    HR = sprintf("%.2f", x$HR),
    `95% CI` = sprintf(
      "%.2f – %.2f",
      x$CI_lower,
      x$CI_upper
    ),
    `p-value` = p_display,
    check.names = FALSE
  )
}


# Assemble Main Table 2
table2_final <- bind_rows(
  header_row("Medicaid Expansion"),
  result_row(
    "did_interaction",
    "Expansion × Post"
  ),
  
  result_row(
    "dx.yr",
    "Year of diagnosis (per year)"
  ),
  
  header_row("Age (ref: 50–59)"),
  result_row(
    "age_categoryYoung Adult",
    "Young Adult (<50)"
  ),
  result_row(
    "age_categoryOlder Adult",
    "Older Adult (60–64)"
  ),
  
  header_row("Sex (ref: Male)"),
  result_row(
    "SexFemale",
    "Female"
  ),
  
  header_row("Race/Ethnicity (ref: NH White)"),
  result_row(
    "Race.EthnicityHispanic (All Races)",
    "Hispanic (All Races)"
  ),
  result_row(
    "Race.EthnicityNon-Hispanic Black",
    "Black"
  ),
  result_row(
    "Race.EthnicityNon-Hispanic Asian or Pacific Islander",
    "Asian / Pacific Islander"
  ),
  result_row(
    "Race.EthnicityNon-Hispanic American Indian/Alaska Native",
    "AI/AN"
  ),
  
  header_row("Marital Status (ref: Married)"),
  result_row(
    "Marital_GroupedPreviously Married",
    "Previously Married"
  ),
  result_row(
    "Marital_GroupedSingle/Unknown",
    "Single / Unknown"
  ),
  
  header_row("Rurality (ref: Metro >1M)"),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Counties in metropolitan areas of 250,000 to 1 million pop"
    ),
    "Metro 250k–1M"
  ),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Counties in metropolitan areas of lt 250 thousand pop"
    ),
    "Metro <250k"
  ),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Nonmetropolitan counties adjacent to a metropolitan area"
    ),
    "Nonmetro, metro-adjacent"
  ),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Nonmetropolitan counties not adjacent to a metropolitan area"
    ),
    "Nonmetro, nonmetro-adjacent"
  ),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Unknown/missing/no match (Alaska or Hawaii - Entire State)"
    ),
    "Unknown / AK–HI"
  ),
  result_row(
    paste0(
      "Rural.Urban.Continuum.Code",
      "Unknown/missing/no match/Not 1990-2023"
    ),
    "Other Unknown"
  ),
  
  header_row("Stage (ref: I)"),
  result_row(
    "stage_groupII",
    "II"
  ),
  result_row(
    "stage_groupIII",
    "III"
  )
)


print(table2_final, row.names = FALSE)

write.csv(
  table2_final,
  "nonmet/table2_mortality_recreated.csv",
  row.names = FALSE,
  na = ""
)


