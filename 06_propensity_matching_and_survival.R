#### Propensity Score Matching ####

cleaned_data <- seer_nonmet

cleaned_data <- cleaned_data %>%
  mutate(Match_Group = sample(c("match_1", "match_2", "match_3"), n(), replace = TRUE))
# List to store results
cox_results <- list()


cleaned_data$binary_group_1 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0,
                                      ifelse(cleaned_data$Medicaid_Expansion_Status == 1, 1, NA))
cleaned_data$binary_group_2 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0,
                                      ifelse(cleaned_data$Medicaid_Expansion_Status == 2, 1, NA))
cleaned_data$binary_group_3 <- ifelse(cleaned_data$Medicaid_Expansion_Status == 0, 0,
                                      ifelse(cleaned_data$Medicaid_Expansion_Status == 3, 1, NA))


# 
# cleaned_data$binary_group_1 <- ifelse(cleaned_data$Medicaid_Expansion_Status == "Non-Expansion", 0, 
#                                       ifelse(cleaned_data$Medicaid_Expansion_Status == "Early Expansion", 1, NA))
# cleaned_data$binary_group_2 <- ifelse(cleaned_data$Medicaid_Expansion_Status == "Non-Expansion", 0, 
#                                       ifelse(cleaned_data$Medicaid_Expansion_Status == "2014 Expansion", 1, NA))
# cleaned_data$binary_group_3 <- ifelse(cleaned_data$Medicaid_Expansion_Status == "Non-Expansion", 0, 
#                                       ifelse(cleaned_data$Medicaid_Expansion_Status == "Late Expansion", 1, NA))
# 


cleaned_data_cleaned_1 <- cleaned_data[!is.na(cleaned_data$binary_group_1), ]

cleaned_data_cleaned_2 <- cleaned_data[!is.na(cleaned_data$binary_group_2), ]
cleaned_data_cleaned_3 <- cleaned_data[!is.na(cleaned_data$binary_group_3), ]


# Create propensity score model using logistic regression



# Group 1 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_1_complete <- cleaned_data_cleaned_1[complete.cases(
  cleaned_data_cleaned_1[, c("age_category", "Sex", 
                             "Race.Ethnicity",
                             "rad_group", "stage_group", "surgery_status", 
                             "income_quintile_equal", "Rural.Urban.Continuum.Code", 
                             "Marital_Grouped")]), ]


# Step 1: Fit Propensity Score Model
ps_model_1 <- glm(binary_group_1 ~ age_category + Sex + Combined.Stage +
                    Race.Ethnicity + 
                    Rural.Urban.Continuum.Code, 
                  data = cleaned_data_cleaned_1_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_1_complete$pscore <- predict(ps_model_1, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_1_complete, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_before_1.png", width = 8, height = 6, dpi = 300)


cleaned_data_cleaned_1_complete$logit_pscore <- log(cleaned_data_cleaned_1_complete$pscore / (1 - cleaned_data_cleaned_1_complete$pscore))

match_1 <- matchit(
  binary_group_1 ~ 1,
  data = cleaned_data_cleaned_1_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_1_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 1,
  replace = FALSE
)

# Step 5: Extract Matched Data
matched_data_1 <- match.data(match_1)

n_pairs <- nrow(matched_data_1) / 2
matched_data_1$matched_id <- rep(1:n_pairs, each = 2)
matched_data_1$matched_id <- matched_data_1$subclass

# Step 6: Plot After Matching
ggplot(matched_data_1, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_after_1.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_1 <- love.plot(match_1, 
                    binary = "raw", 
                    threshold = 0.1, 
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    theme(text = element_text(family = "Times"))
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_1)
ggsave("nonmet/psm/loveplot_1.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_1)

# Step 5: Perform Matching with all covariates in the model
match_1 <- matchit(
  binary_group_1 ~ age_category +
    Sex + 
    Combined.Stage +
    Race.Ethnicity +
    Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_1_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_1_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 1,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_1 <- match.data(match_1)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_1, aes(x = pscore, fill = as.factor(binary_group_1))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 1", 
       x = "Propensity Score", y = "Density") +
  theme_minimal()
# ggsave("nonmet/psm/propensity_scores_after_1.png", width = 8, height = 6, dpi = 300)


# Check balance again
bal.tab(match_1)
# Step 8: Covariate Balance Plot
plot_1 <- love.plot(match_1, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")


#Print the love plot
print(plot_1)
ggsave("nonmet/psm/full_loveplot_1.png", width = 16, height = 4, dpi = 300)

bal.tab(match_1)


# Group 2 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_2_complete <- cleaned_data_cleaned_2[complete.cases(
  cleaned_data_cleaned_2[, c("age_category", "Sex", 
                             "Race.Ethnicity",
                             "rad_group", "surgery_status", 
                             "income_quintile_equal", "Rural.Urban.Continuum.Code", 
                             "Marital_Grouped")]), ]

# Step 1: Fit Propensity Score Model
ps_model_2 <- glm(binary_group_2 ~ age_category + Sex + Combined.Stage +
                    Race.Ethnicity + 
                    Rural.Urban.Continuum.Code
                  , 
                  data = cleaned_data_cleaned_2_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_2_complete$pscore <- predict(ps_model_2, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_2_complete, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_before_2.png", width = 8, height = 6, dpi = 300)

# Step 4: Perform Matching
cleaned_data_cleaned_2_complete$logit_pscore <- log(cleaned_data_cleaned_2_complete$pscore / (1 - cleaned_data_cleaned_2_complete$pscore))

match_2 <- matchit(
  binary_group_2 ~ 1,
  data = cleaned_data_cleaned_2_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_2_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 1,
  replace = FALSE
)

# match_2 <- matchit(binary_group_2 ~ 1, 
#                    data = cleaned_data_cleaned_2_complete, 
#                    method = "nearest", 
#                    distance = cleaned_data_cleaned_2_complete$pscore, 
#                    ratio = 1, 
#                    caliper = 0.2, 
#                    replace = FALSE)

# Step 5: Extract Matched Data
matched_data_2 <- match.data(match_2)

n_pairs <- nrow(matched_data_2) / 2
matched_data_2$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2$matched_id <- matched_data_2$subclass

# Step 6: Plot After Matching
ggplot(matched_data_2, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_after_2.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_2 <- love.plot(match_2, 
                    binary = "raw", 
                    threshold = 0.1, 
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_2)
ggsave("nonmet/psm/loveplot_2.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_2)


# Step 5: Perform Matching with all covariates in the model
match_2 <- matchit(
  binary_group_2 ~ age_category + Sex + Combined.Stage +
    Race.Ethnicity +  Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_2_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_2_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 1,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_2 <- match.data(match_2)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_2, aes(x = pscore, fill = as.factor(binary_group_2))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 2", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_after_2.png", width = 8, height = 6, dpi = 300)

# Step 8: Covariate Balance Plot
plot_2 <- love.plot(match_2, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group 2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

#Print the love plot
print(plot_2)
ggsave("nonmet/psm/full_loveplot_2.png", width = 16, height = 4, dpi = 300)

# Step 8: Show Table of Balance
bal.tab(match_2)



# Group 3 vs 0
# Step 0: Filter data to complete cases for modeling
cleaned_data_cleaned_3_complete <- cleaned_data_cleaned_3[complete.cases(
  cleaned_data_cleaned_3[, c("age_category", "Sex", 
                             "Race.Ethnicity",
                             "rad_group", "stage_group", "surgery_status", 
                             "income_quintile_equal", "Rural.Urban.Continuum.Code", 
                             "Marital_Grouped")]), ]

# Step 1: Fit Propensity Score Model
ps_model_3 <- glm(binary_group_3 ~ age_category + Sex + Combined.Stage +
                    Race.Ethnicity + 
                    Rural.Urban.Continuum.Code
                  , 
                  data = cleaned_data_cleaned_3_complete, 
                  family = "binomial")

# Step 2: Attach predicted propensity scores to this filtered dataset
cleaned_data_cleaned_3_complete$pscore <- predict(ps_model_3, type = "response")

# Step 3: Plot Before Matching
ggplot(cleaned_data_cleaned_3_complete, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions Before Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_before_3.png", width = 8, height = 6, dpi = 300)

# Step 4: Perform Matching
cleaned_data_cleaned_3_complete$logit_pscore <- log(cleaned_data_cleaned_3_complete$pscore / (1 - cleaned_data_cleaned_3_complete$pscore))

match_3 <- matchit(
  binary_group_3 ~ 1,
  data = cleaned_data_cleaned_3_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_3_complete$logit_pscore,
  caliper = 0.2,  # on logit scale
  ratio = 1,
  replace = FALSE
)
# 
# 
# match_3 <- matchit(binary_group_3 ~ 1, 
#                    data = cleaned_data_cleaned_3_complete, 
#                    method = "nearest", 
#                    distance = cleaned_data_cleaned_3_complete$pscore, 
#                    ratio = 1, 
#                    caliper = 0.2, 
#                    replace = FALSE)

# Step 5: Extract Matched Data
matched_data_3 <- match.data(match_3)

n_pairs <- nrow(matched_data_3) / 2
matched_data_3$matched_id <- rep(1:n_pairs, each = 2)
matched_data_3$matched_id <- matched_data_3$subclass

# Step 6: Plot After Matching
ggplot(matched_data_3, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal() + theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_after_3.png", width = 8, height = 6, dpi = 300)

# Step 7: Covariate Balance Plot
plot_3 <- love.plot(match_3, 
                    binary = "raw", 
                    threshold = 0.1, 
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group 1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

print(plot_3)


print(plot_3)
ggsave("nonmet/psm/loveplot_3.png", width = 16, height = 4, dpi = 300)


# Step 8: Show Table of Balance
bal.tab(match_3)


# Step 5: Perform Matching with all covariates in the model
match_3 <- matchit(
  binary_group_3 ~ age_category + Sex + Combined.Stage +
    Race.Ethnicity + Rural.Urban.Continuum.Code,  # Include all covariates used in the model
  data = cleaned_data_cleaned_3_complete,
  method = "nearest",
  distance = cleaned_data_cleaned_3_complete$logit_pscore,
  caliper = 0.2,  # On logit scale
  ratio = 1,
  replace = FALSE
)

# Step 6: Extract Matched Data
matched_data_3 <- match.data(match_3)

# Step 7: Plot Propensity Score Distributions After Matching
ggplot(matched_data_3, aes(x = pscore, fill = as.factor(binary_group_3))) +
  geom_density(alpha = 0.4) +
  labs(title = "Propensity Score Distributions After Matching for Group 3", 
       x = "Propensity Score", y = "Density") +
  theme_minimal()+ theme(text = element_text(family = "Times"))
ggsave("nonmet/psm/propensity_scores_after_3.png", width = 8, height = 6, dpi = 300)

# Step 8: Covariate Balance Plot
plot_3 <- love.plot(match_3, 
                    binary = "raw", 
                    threshold = 0.1,  # Set threshold for standardized mean differences
                    sample.names = c("Before PSM", "After PSM"),                      title = "Covariate Balance: Absolute Standardized Mean Differences for Group ") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Covariates", y = "Absolute Standardized Mean Difference")

#Print the love plot
print(plot_3)
ggsave("nonmet/psm/full_loveplot_3.png", width = 16, height = 4, dpi = 300)

# Step 8: Show Table of Balance
bal.tab(match_3)


bal.tab(match_1)  # For Group 1 vs 0
bal.tab(match_2)  # For Group 2 vs 0
bal.tab(match_3)  # For Group 3 vs 0


#### Cox Proportional Hazards Ratio for 2-year Survival ####
cleaned_data <- cleaned_data %>%
  filter(!grepl("Unknown Race.Ethnicity", trimws(Race.Ethnicity)))

#Combine previous 3 Matched Data tables into 1 group
# Combine the matched datasets and add a group identifier
matched_data_1$group <- "Group 1"
matched_data_2$group <- "Group 2"
matched_data_3$group <- "Group 3"

# Combine all matched datasets into one

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


n_pairs <- nrow(matched_data_1_did) / 2
matched_data_1_did$matched_id <- rep(1:n_pairs, each = 2)
matched_data_1_did$matched_id <- matched_data_1$subclass

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

capped_fin_cox_psm_final_did_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat ,
                                        data = matched_data_1_did, robust = TRUE)


summary(capped_fin_cox_psm_final_did_1)

capped_fin_cox_psm_final_did_1 <- tidy(capped_fin_cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_1

write.csv(capped_fin_cox_psm_final_did_1, "nonmet/coxdata_3/psm_0v1.csv")

cox_psm_final_did_1 <- coxph(Surv(Survival.months, Event) ~ did_time * did_treat ,
                             data = matched_data_1_did, robust = TRUE)
summary(cox_psm_final_did_1)

cox_psm_final_did_1 <- tidy(cox_psm_final_did_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

cox_psm_final_did_1



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

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat, #,
                                          data = matched_data_1_did_early_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_1)

capped_fin_cox_psm_final_did_early_1 <- tidy(capped_fin_cox_psm_final_early_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

cox.zph(capped_fin_cox_psm_final_early_1)
# # # # View(capped_fin_cox_psm_final_did_early_1)

write.csv(capped_fin_cox_psm_final_did_early_1, "nonmet/coxdata_3/psm_0v1_early.csv")


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

# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat, #,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_3/psm_0v1_later.csv")

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
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_3/psm_0v1_LAST.csv")


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
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_3/psm_0v1_MID.csv")



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
capped_fin_cox_psm_final_later_1 <- coxph(Surv(capped_time, capped_event) ~ did_time * did_treat ,
                                          data = matched_data_1_did_later_1, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_1)

capped_fin_cox_psm_final_did_later_1 <- tidy(capped_fin_cox_psm_final_later_1, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

# # # # View(capped_fin_cox_psm_final_did_later_1)

write.csv(capped_fin_cox_psm_final_did_later_1, "nonmet/coxdata_3/psm_0v1_2011-2016.csv")



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


n_pairs <- nrow(matched_data_2_did) / 2
matched_data_2_did$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2_did$matched_id <- matched_data_1$subclass


matched_data_2_did <- matched_data_2_did %>%
  filter(!grepl("Unknown Race.Ethnicity", trimws(Race.Ethnicity)))


capped_fin_cox_psm_final_did_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 ,
                                        data = matched_data_2_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_2)

capped_fin_cox_psm_final_did_2 <- tidy(capped_fin_cox_psm_final_did_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_2

write.csv(capped_fin_cox_psm_final_did_2, "nonmet/coxdata_3/psm_0v2.csv")

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



# Step 4: Run Cox regression model
capped_fin_cox_psm_final_early_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 ,
                                          data = matched_data_2_did_early_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_early_2)

capped_fin_cox_psm_final_did_early_2 <- tidy(capped_fin_cox_psm_final_early_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_early_2)

write.csv(capped_fin_cox_psm_final_did_early_2, "nonmet/coxdata_3/psm_0v2_early.csv")




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


n_pairs <- nrow(matched_data_2_did_later_2) / 2
matched_data_2_did_later_2$matched_id <- rep(1:n_pairs, each = 2)
matched_data_2_did_later_2$matched_id <- matched_data_2_did_later_2$subclass


# Step 3: Ensure 'did_time' is a factor with 'Pre' as reference
matched_data_2_did_later_2 <- matched_data_2_did_later_2 %>%
  mutate(
    did_time_2 = factor(did_time_2, levels = c("Pre", "Post")),
    did_treat_2 = factor(did_treat_2, levels = c(0, 1))
  )

# Check reference levels
contrasts(matched_data_2_did_later_2$did_time_2)
contrasts(matched_data_2_did_later_2$did_treat_2)




# Step 4: Run Cox regression model
capped_fin_cox_psm_final_later_2 <- coxph(Surv(capped_time, capped_event) ~ did_time_2 * did_treat_2 ,
                                          data = matched_data_2_did_later_2, robust = TRUE)

# Step 5: Summary of the model
summary(capped_fin_cox_psm_final_later_2)

capped_fin_cox_psm_final_did_later_2 <- tidy(capped_fin_cox_psm_final_later_2, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(capped_fin_cox_psm_final_did_later_2)

write.csv(capped_fin_cox_psm_final_did_later_2, "nonmet/coxdata_3/psm_0v2_later.csv")



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


capped_fin_cox_psm_final_did_3 <- coxph(Surv(capped_time, capped_event) ~ did_time_3 * did_treat_3 ,
                                        data = matched_data_3_did, robust = TRUE)

summary(capped_fin_cox_psm_final_did_3)

capped_fin_cox_psm_final_did_3 <- tidy(capped_fin_cox_psm_final_did_3, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

capped_fin_cox_psm_final_did_3

write.csv(capped_fin_cox_psm_final_did_3, "nonmet/coxdata_3/psm_0v3.csv")


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




