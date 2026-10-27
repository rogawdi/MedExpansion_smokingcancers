## Healthcare Disparities for total Pre/Post populations # 

seer_nonmet_disp <- seer_nonmet%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )

# ### Trying this again ###
# 
# # Step 0: Data prep
# 
# seer_nonmet_disp <- seer_nonmet_disp %>%
#   mutate(
#     Treatment_ME = factor(ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0), levels = c(0,1)),
#     Time_Period  = factor(ifelse(dx.yr < exp_yr, "Pre","Post"), levels = c("Pre","Post")),
#     surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0),
#     income_quintile_equal = relevel(factor(income_quintile_equal), ref = "Q3"),
#     Race.Ethnicity = relevel(factor(Race.Ethnicity), ref = "Non-Hispanic White"),
#     stage_group = relevel(factor(stage_group), ref = "III"),
#     age_category = relevel(factor(age_category), ref = "Middle-Aged"),
#     Sex = relevel(factor(Sex), ref = "Male"),
#     Marital_Grouped = relevel(factor(Marital_Grouped), ref = "Single/Unknown"),
#     Rural.Urban.Continuum.Code = factor(
#       case_when(
#         Rural.Urban.Continuum.Code == "Unknown/missing/no match (Alaska or Hawaii - Entire State)" ~ 
#           "Nonmetropolitan counties adjacent to a metropolitan area",
#         TRUE ~ as.character(Rural.Urban.Continuum.Code)
#       )
#     )
#   ) %>%
#   filter(!grepl("Non-Hispanic Unknown Race", trimws(Race.Ethnicity)))
# 
# # Step 1: Time stratification
# 
# cut_times <- seq(from = 12, to = 180, by = 60) # 1-15 years
# 
# stratified_data <- survSplit(
#   Surv(Survival.months, Event) ~ .,
#   data = seer_nonmet_disp,
#   cut = cut_times,
#   event = "Event",
#   start = "tstart",
#   end = "tstop",
#   zero = -1   # ensures it's < any observed survival time
# )
# stratified_data$time_strata <- cut(
#   stratified_data$tstop,
#   breaks = c(0, cut_times),
#   labels = c("0-24m", "24-60m", "60-180m"),
#   include.lowest = TRUE, right = TRUE
# )
# 
# 
# # Step 2: Final Cox model with time-varying covariates
# 
# final_cox_model <- coxph(
#   Surv(tstart, tstop, Event) ~
#     Treatment_ME + Time_Period + Treatment_ME:Time_Period +
#     # The tt() function correctly models time-varying effects for your covariates.
#     tt(age_category) + tt(Marital_Grouped) + tt(Sex) + tt(Race.Ethnicity) +
#     tt(income_quintile_equal) + tt(Rural.Urban.Continuum.Code) +
#     stage_group + surgery_status + rad_group +
#     strata(time_strata),
#   data = stratified_data,
#   robust = TRUE,
#   cluster = Medicaid.Expansion.Status
# )
# 
# 
# # Step 3: Extract and tidy ALL results for manuscript
# 
# final_results <- tidy(final_cox_model, conf.int = TRUE) %>%
#   mutate(
#     HR = exp(estimate),
#     CI_lower = exp(conf.low),
#     CI_upper = exp(conf.high),
#     sig = ifelse(p.value < 0.05, "Yes", "No")
#   )
# 
# # Add a 'type' column to distinguish between fixed and time-varying effects
# final_results <- final_results %>%
#   mutate(
#     type = case_when(
#       grepl("tt\\(age_category\\)", term) ~ "Time-Varying: Age",
#       grepl("tt\\(Marital_Grouped\\)", term) ~ "Time-Varying: Marital Status",
#       grepl("tt\\(Sex\\)", term) ~ "Time-Varying: Sex",
#       grepl("tt\\(Race.Ethnicity\\)", term) ~ "Time-Varying: Race",
#       grepl("tt\\(income_quintile_equal\\)", term) ~ "Time-Varying: Income",
#       grepl("tt\\(Rural.Urban.Continuum.Code\\)", term) ~ "Time-Varying: Rurality",
#       grepl("Time_PeriodPost", term) ~ "Interaction",
#       TRUE ~ "Fixed"
#     )
#   )
# 
# # View the full results.
# print(final_results, n=Inf)
# write.csv(final_results, "stageIV/supplement_tt()_validation.csv", row.names = FALSE)
# 
# final_disparities_model <- coxph(
#   Surv(tstart, tstop, Event) ~
#     Treatment_ME + Time_Period + Treatment_ME:Time_Period +
#     # We use explicit interactions to get the hazard ratios for each group
#     Time_Period:age_category +
#     Time_Period:Marital_Grouped +
#     Time_Period:Sex +
#     Time_Period:Race.Ethnicity +
#     Time_Period:income_quintile_equal +
#     Time_Period:Rural.Urban.Continuum.Code +
#     stage_group + surgery_status + rad_group +
#     strata(time_strata),
#   data = stratified_data,
#   robust = TRUE,
#   cluster = Medicaid.Expansion.Status
# )
# 
# summary(final_disparities_model)
# 
# 
# # Step 3: Extract and tidy results for manuscript
# 
# disparities_results <- tidy(final_disparities_model, conf.int = TRUE) %>%
#   filter(grepl("Time_PeriodPost", term)) %>%       # keep only Post interactions
#   mutate(
#     # Assign group labels for faceting
#     group = case_when(
#       grepl("age_category", term) ~ "Age",
#       grepl("Marital_Grouped", term) ~ "Marital Status",
#       grepl("Sex", term) ~ "Sex",
#       grepl("Race.Ethnicity", term) ~ "Race/Ethnicity",
#       grepl("income_quintile_equal", term) ~ "Income Quintile",
#       grepl("Rural\\.Urban\\.Continuum\\.Code", term) ~ "Rurality",
#       TRUE ~ NA_character_  # mark anything not a subgroup as NA
#     )
#   ) %>%
#   filter(!is.na(group)) %>%  # remove Overall rows
#   mutate(
#     HR = exp(estimate),
#     CI_lower = exp(conf.low),
#     CI_upper = exp(conf.high),
#     sig = ifelse(p.value < 0.05, "Yes", "No"),
#     
#     # Clean term names for plotting: keep only subgroup identity
#     term_clean = term %>%
#       gsub("^Time_PeriodPost:", "", .) %>%
#       gsub("age_category", "", .) %>%
#       gsub("Marital_Grouped", "", .) %>%
#       gsub("Sex", "", .) %>%
#       gsub("Race.Ethnicity", "", .) %>%
#       gsub("income_quintile_equal", "", .) %>%
#       gsub("Rural\\.Urban\\.Continuum\\.Code", "", .) %>%
#       trimws()
#   ) %>%
#   arrange(group, term_clean)
# 
# 
# 
# # Forest plot
# 
# ggplot(disparities_results, aes(x = HR, y = term_clean)) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
#   geom_point(aes(color = sig), size = 3) +
#   geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = sig), height = 0.2) +
#   scale_color_manual(values = c("Yes" = "#1b7837", "No" = "#762a83")) +
#   labs(
#     title = "Changes in Disparities from Pre- to Post-Medicaid Expansion",
#     x = "Interaction HR (Post vs Pre, relative to Reference Group)",
#     y = NULL,
#     color = "Statistically Significant"
#   ) +
#   facet_grid(group ~ ., scales = "free_y", space = "free", switch = "y") +
#   theme_minimal(base_size = 14) +
#   theme(
#     strip.text.y.left = element_text(angle = 0, size = 13, hjust = 0.5),
#     strip.placement = "outside",
#     strip.background = element_blank(),
#     axis.text.y = element_text(size = 10),
#     legend.position = "bottom",
#     panel.grid.major.y = element_blank(),
#     panel.spacing = unit(1.2, "lines")
#   )
# 
# ggsave("stageIV/stratified_forestplot_disparityCLN.png", width = 12, height = 8, dpi = 300)
# 
# 
# # Print and save the final results
# print(disparities_results, n = Inf)
# write.csv(disparities_results, "stageIV/timestratified_interaction_disparities.csv", row.names = FALSE)
# 



#### Demographics ####

# 
# # Define variables
# vars <- c("Sex",
#           "Race.Ethnicity",
#           "age_category",
#           "stage_group",
#           "income_quintile_equal",
#           "Rural.Urban.Continuum.Code",
#           "Marital_Grouped",
#           "rad_group",
#           "surgery_status")
# 
# # Create TableOne object
# table1_expansion <- CreateTableOne(
#   vars = vars,
#   strata = "Medicaid_Expansion_Status",
#   data = seer_nonmet,
#   addOverall = TRUE
# )
# 
# # Ensure the n (%) format is clean
# # Extract summary with n (%) formatting
# tab1_mat <- print(table1_expansion,
#                   showAllLevels = TRUE,
#                   test = TRUE,
#                   quote = FALSE,
#                   noSpaces = TRUE,
#                   format = "fp",        # <-- gives counts + percents
#                   digits = 2,
#                   printToggle = FALSE)
# 
# # Convert to data frame
# tab1_fmt <- as.data.frame(tab1_mat)
# 
# # Fix the ugly X, X.1, etc. by moving rownames into a proper column
# tab1_fmt$Variable <- rownames(tab1_fmt)
# rownames(tab1_fmt) <- NULL
# 
# # Reorder so Variable column is first
# tab1_fmt <- tab1_fmt[, c(ncol(tab1_fmt), 1:(ncol(tab1_fmt)-1))]
# 
# # Now it should look like: "150860 (61.9)" instead of raw counts
# head(tab1_fmt)
# 
# write.csv(tab1_fmt, "nonmet/demographics_MEgroups.csv", row.names = FALSE)
# 
# 
# # Create TableOne object
# seer_nonmet <- seer_nonmet%>%
#   mutate(
#     exp_yr = case_when(
#       (Medicaid_Expansion_Status == "0") ~ "2011",
#       (Medicaid_Expansion_Status == "1") ~ "2014",
#       (Medicaid_Expansion_Status == "2") ~ "2017",
#       (Medicaid_Expansion_Status == "3") ~ "2020"
#     )
#   )
# 
# 
# seer_nonmet <- seer_nonmet %>%
#   mutate(
#     # Treatment is time-invariant: did the state ever expand?
#     Treatment_ME = ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0),
#     
#     # Time is relative to each state's expansion year
#     Time_Period = factor(case_when(dx.yr <= 2010 ~ "2006-2010", dx.yr <= 2013 ~ "2011-2013", dx.yr <= 2016 ~ "2014-2016", dx.yr <= 2019 ~ "2017-2019"))
#   ) %>%
#   mutate(
#     Treatment_ME = factor(Treatment_ME, levels = c(0,1)),
    # Time_Period  = factor(Time_Period, levels = c("Pre","Post")),
#     surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0)
#   )
# 
# table2_expansion <- CreateTableOne(
#   vars = vars,
#   strata = "Treatment_ME",
#   data = seer_nonmet,
#   addOverall = TRUE
# )
# 
# # Ensure the n (%) format is clean
# # Extract summary with n (%) formatting
# tab2_mat <- print(table2_expansion,
#                   showAllLevels = TRUE,
#                   test = TRUE,
#                   quote = FALSE,
#                   noSpaces = TRUE,
#                   format = "fp",        # <-- gives counts + percents
#                   digits = 2,
#                   printToggle = FALSE)
# 
# # Convert to data frame
# tab2_fmt <- as.data.frame(tab2_mat)
# 
# # Fix the ugly X, X.1, etc. by moving rownames into a proper column
# tab2_fmt$Variable <- rownames(tab2_fmt)
# rownames(tab2_fmt) <- NULL
# 
# # Reorder so Variable column is first
# tab2_fmt <- tab2_fmt[, c(ncol(tab2_fmt), 1:(ncol(tab2_fmt)-1))]
# 
# # Now it should look like: "150860 (61.9)" instead of raw counts
# head(tab2_fmt)
# 
# write.csv(tab2_fmt, "overall/demographics_exp_binary.csv", row.names = FALSE)


#### Censoring at 3 years #### 

### Healthcare Disparities for total Pre/Post populations ### 
seer_censored <- seer_nonmet

seer_censored$capped_time <- pmin(seer_censored$Survival.months, 36)
seer_censored$capped_event <- ifelse(seer_censored$Survival.months > 36, 0, seer_censored$Event)


# Step 1: Define expansion timing for each group
seer_censored <- seer_censored %>%
  mutate(expansion1_time = if_else(
    Medicaid_Expansion_Status == 1,
    case_when(
      dx.yr_grouped %in% c("2006-2010") ~ "Pre",
      dx.yr_grouped %in% c("2011-2013", "2014-2016", "2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

seer_censored <- seer_censored %>%
  mutate(expansion2_time = if_else(
    Medicaid_Expansion_Status == 2,
    case_when(
      dx.yr_grouped %in% c("2006-2010", "2011-2013") ~ "Pre",
      dx.yr_grouped %in% c("2014-2016", "2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

seer_censored <- seer_censored %>%
  mutate(expansion3_time = if_else(
    Medicaid_Expansion_Status == 3,
    case_when(
      dx.yr_grouped %in% c("2006-2010", "2011-2013", "2014-2016") ~ "Pre",
      dx.yr_grouped %in% c("2017-2019") ~ "Post",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

seer_censored <- seer_censored %>%
  mutate(expansion0_time = if_else(
    Medicaid_Expansion_Status == 0,
    case_when(
      dx.yr_grouped %in% c("2006-2010","2011-2013", "2014-2016", "2017-2019") ~ "Pre",
      TRUE ~ NA_character_
    ),
    NA_character_
  ))

# Step 2: Define Pre/Post periods across all expansion groups
seer_disp_censored <- seer_censored %>%
  mutate(
    Time_Period = case_when(
      (Medicaid_Expansion_Status == "1" & expansion1_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "1" & expansion1_time == "Post") ~ "Post",
      (Medicaid_Expansion_Status == "2" & expansion2_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "2" & expansion2_time == "Post") ~ "Post",
      (Medicaid_Expansion_Status == "3" & expansion3_time == "Pre") ~ "Pre",
      (Medicaid_Expansion_Status == "3" & expansion3_time == "Post") ~ "Post",
      (Medicaid_Expansion_Status == "0" & expansion0_time == "Pre") ~ "Pre",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Time_Period))

seer_disp_censored$Time_Period <- factor(seer_disp_censored$Time_Period)

seer_disp_censored$Time_Period <- relevel(seer_disp_censored$Time_Period, ref = "Pre")


# Set references
seer_disp_censored <- seer_disp_censored %>%
  mutate(
    income_quintile_equal = relevel(factor(income_quintile_equal), ref = "Q3"),
    Race.Ethnicity             = relevel(factor(Race.Ethnicity),             ref = "Non-Hispanic White"),
    stage_group      = relevel(factor(stage_group),      ref = "III"),
    age_category     = relevel(factor(age_category),     ref = "Middle-Aged"),
    Sex = relevel(factor(Sex), ref = "Male"),
    Marital_Grouped  = relevel(factor(Marital_Grouped),  ref = "Single/Unknown"),
    Rural.Urban.Continuum.Code = relevel(factor(Rural.Urban.Continuum.Code), ref = "Counties in metropolitan areas ge 1 million pop")
  )

seer_disp_censored <- seer_disp_censored %>%
  mutate(
    Rural.Urban.Continuum.Code = as.character(Rural.Urban.Continuum.Code),
    Rural.Urban.Continuum.Code = case_when(
      Rural.Urban.Continuum.Code == "Unknown/missing/no match (Alaska or Hawaii - Entire State)" ~ 
        "Nonmetropolitan counties adjacent to a metropolitan area",
      TRUE ~ Rural.Urban.Continuum.Code
    ),
    Rural.Urban.Continuum.Code = factor(Rural.Urban.Continuum.Code)
  )

# Step 3: Create a Cox model with interaction terms to test for disparity change
cox_model <- coxph(
  Surv(capped_time, capped_event) ~ 
    (age_category + Marital_Grouped + Sex + 
       + Race.Ethnicity + income_quintile_equal + surgery_status + stage_group + rad_group + Rural.Urban.Continuum.Code) * Time_Period,
  data = seer_disp_censored,
  robust = TRUE,
  cluster = Medicaid.Expansion.Status
)

print(cox_model)


# Step 4: Extract interaction terms (i.e., change in disparities pre/post)
interaction_results <- tidy(cox_model, conf.int = TRUE) %>%
  filter(grepl(":Time_PeriodPost", term)) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No")
  ) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value, sig)

# # View interaction effects
write.csv(file ="nonmet/3yr/interaction.effects.csv", interaction_results)

library(forcats)

# Step 4: Extract and clean interaction terms
interaction_results <- tidy(cox_model, conf.int = TRUE) %>%
  filter(grepl(":Time_PeriodPost", term)) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No"),
    
    # Clean term names for display
    term_clean = term %>%
      gsub(":Time_PeriodPost", "", .) %>%
      gsub("age_category", "", .) %>%
      gsub("Sex", "", .) %>%
      gsub("Marital_Grouped", "", .) %>%
      gsub("Race.Ethnicity", "", .) %>%
      gsub("income_quintile_equal", "", .) %>%
      gsub("stage_group", "", .) %>%
      gsub("surgery_status", "", .) %>%
      gsub("rad_group", "", .) %>%
      gsub("Rural\\.Urban\\.Continuum\\.Code", "", .) %>%
      gsub("Site\\.recode\\.ICD\\.O\\.3\\.2023\\.Revision", "", .) %>%
      gsub("Unknown/missing/no match \\(Alaska or Hawaii - Entire State\\)", 
           "Nonmetropolitan counties adjacent to a metropolitan area", .) %>%
      gsub("Young Adult", "<50 years old", .) %>%
      gsub("Older Adult", "60–64 years old", .) %>%
      trimws(),
    
    # Assign group labels
    group = case_when(
      grepl("Race.Ethnicity", term) ~ "Race.Ethnicity",
      grepl("age_category", term) ~ "Age",
      grepl("Sex", term) ~ "Sex",
      grepl("Marital_Grouped", term) ~ "Marital Status",
      grepl("income_quintile_equal", term) ~ "Income Quintile",
      grepl("stage_group", term) ~ "Stage",
      grepl("surgery_status", term) ~ "Surgery",
      grepl("rad_group", term) ~ "Radiation",
      grepl("Rural\\.Urban\\.Continuum\\.Code", term) ~ "Rurality",
      grepl("Site\\.recode", term) ~ "Primary Tumor Site",
      TRUE ~ "Other"
    ),
    
    # Combine label for plotting
    term_label = paste0(group, ": ", term_clean)
  ) %>%
  filter(!grepl("Race: Unknown|Rurality:.*Not.*1990.*2023", term_label)) %>%
  mutate(
    term_label = fct_inorder(term_label),
    group = factor(group, levels = unique(group))  # lock in group order
  ) %>%
  arrange(group, term_label)

# Step 5: Forest plot with facet strip headers
ggplot(interaction_results, aes(x = HR, y = term_label)) +
  geom_point(aes(color = sig), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = sig), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Yes" = "#1b7837", "No" = "#762a83")) +
  labs(
    title = "Interaction Effects for Smoking-Related Aerodigestive Cancers (36 Month Cap): Pre- vs Post-Expansion",
    x = "Interaction HR (Post vs Pre, vs Reference Group)",
    y = NULL,
    color = "Statistically Significant"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free", switch = "y") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text.y.left = element_text(angle = 0, size = 13, hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

ggsave("nonmet/3yr/forest_plot_interactions.png", width = 14, height = 10, dpi = 300)

interaction_results %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No")
  ) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value, sig) %>%
  print(n = Inf)
# View(interaction_results)
write.csv(file ="nonmet/3yr/full.interaction.effects.csv", interaction_results)



