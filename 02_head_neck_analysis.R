#### Head and Neck Analysis ####
# Define Medicaid expansion categories
seer_HN$Medicaid.Expansion.Status <- tolower(seer_HN$Medicaid.Expansion.Status)

seer_HN$Medicaid_Expansion_Status <- as.factor(seer_HN$Medicaid.Expansion.Status)

seer_HN <- seer_HN %>%
  mutate(dx.yr_grouped = case_when(
    dx.yr >= 2006 & dx.yr <= 2010 ~ "2006-2010",
    dx.yr >= 2011 & dx.yr <= 2013 ~ "2011-2013",
    dx.yr >= 2014 & dx.yr <= 2016 ~ "2014-2016",
    dx.yr >= 2017 & dx.yr <= 2019 ~ "2017-2019",
    TRUE ~ NA_character_
  ))
seer_HN$dx.yr_grouped <- factor(seer_HN$dx.yr_grouped)

seer_HN <- seer_HN %>%
  mutate(
    Medicaid_Expansion_Status = case_when(
      # 2006-2010: Assign pre-expansion groups based on future expansion
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  # future EARLY expander
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new mexico") ~ 2,  # future ON TIME expander
      dx.yr_grouped == "2006-2010" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future LATE expander
      dx.yr_grouped == "2006-2010" ~ 0,  # Never expanders stay 0
      
      # 2010-2013: First expansion
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  # NEW expander
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new mexico") ~ 2,  # future expander
      dx.yr_grouped == "2011-2013" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future expander
      dx.yr_grouped == "2011-2013" ~ 0,  
      
      # 2014-2016: Second expansion
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new mexico") ~ 2,  #new expander
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  #already expanded
      dx.yr_grouped == "2014-2016" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # future expander
      dx.yr_grouped == "2014-2016" ~ 0,  
      
      # 2017-2019: Final expansion
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("louisiana", "alaska natives") ~ 3,  # new expander
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("iowa", "hawaii", "kentucky", "new mexico") ~ 2,  #already expanded
      dx.yr_grouped == "2017-2019" & Medicaid.Expansion.Status %in% c("california", "connecticut", "new jersey", "seattle (puget sound)") ~ 1,  #already expanded
      dx.yr_grouped == "2017-2019" ~ 0,  
      
      TRUE ~ 0  # Default: No expansion
    )
  )

# seer_HN$Medicaid.Expansion.Status <- relevel(seer_HN$Medicaid.Expansion.Status, ref = "0")



# Check if mapping was successful:
table(seer_HN$Medicaid_Expansion_Status)

# Convert categorical variables into factors
seer_HN$Age <- factor(seer_HN$Age)  # Age grouped by 5 years
seer_HN$Combined.Stage <- factor(seer_HN$Combined.Stage)

stage_map <- data.frame(
  stage_range = c("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA","Stage IIIB+C", "Stage IV"),
  stage_group = c("I", "II", "II", "III","III", "IV")
) 
# Merge SEER dataset with age categories
seer_HN <- seer_HN %>%
  left_join(stage_map, by = c("Combined.Stage" = "stage_range"))

# # # View new age categories
table(seer_HN$stage_group, useNA = "ifany")

seer_HN$Sex <- factor(seer_HN$Sex, levels = c("Male", "Female"))
seer_HN$Median.household.income.inflation.adj.to.2023 <- factor(seer_HN$Median.household.income.inflation.adj.to.2023)  # Median household 
seer_HN$Rural.Urban.Continuum.Code <- factor(seer_HN$Rural.Urban.Continuum.Code)  # Rural/urban continuum
#seer_HN <- seer_HN %>%
#  filter(!grepl("Unknown/missing/no match", trimws(Rural.Urban.Continuum.Code)))
seer_HN$Marital.status.at.diagnosis <- factor(seer_HN$Marital.status.at.diagnosis)  # Marital status
# seer_HN$Marital.status.at.diagnosis <- relevel(seer_HN$Marital.status.at.diagnosis, ref = "Married (including common law)") # make Married the reference

seer_HN$Race.Ethnicity <- factor(seer_HN$Race.Ethnicity)  # Race/ethnicity
seer_HN <- seer_HN %>%
  filter(!grepl("Non-Hispanic Unknown Race", trimws(Race.Ethnicity)))

seer_HN$Race.Ethnicity <- relevel(seer_HN$Race.Ethnicity, ref = "Non-Hispanic White") # make NH White the reference


# Group Radiation recodes
seer_HN$Radiation.recode <- factor(seer_HN$Radiation.recode)  # Radiation
rad_map <- data.frame(
  rad_range = c("Beam radiation", 
                "Combination of beam with implants or isotopes", 
                "None/Unknown",
                "Radiation, NOS  method or source not specified",
                "Radioactive implants (includes brachytherapy) (1988+)",
                "Radioisotopes (1988+)","Recommended, unknown if administered",
                "Refused (1988+)"),
  rad_group = c("Radiation recieved", "Radiation recieved","none/unknown/refused",
                "Radiation recieved",
                "Radiation recieved", "Radiation recieved",
                "none/unknown/refused","none/unknown/refused")
)
# Merge SEER dataset with radiation categories
seer_HN <- seer_HN %>%
  left_join(rad_map, by = c("Radiation.recode" = "rad_range"))

# # # View new rad categories
table(seer_HN$rad_group, useNA = "ifany")

# seer_HN$rad_group <- relevel(seer_HN$rad_group, ref = "none/unknown/refused") # should radiation performed be baseline? mmm actually no i think should not be is baseline since binary. 


# Group Surgery recodes
seer_HN$Surgery.Performed <- factor(seer_HN$Surgery.Performed)  # Surgery Status

seer_HN$Surgery.Performed <- factor(seer_HN$Surgery.Performed)  # Surgery

surg_map <- data.frame(
  surg_range = c("Surgery performed", 
                 "Not performed, patient died prior to recommended surgery", 
                 "Not recommended",
                 "Not recommended, contraindicated due to other cond; autopsy only (1973-2002)",
                 "Recommended but not performed, patient refused",
                 "Recommended but not performed, unknown reason",
                 "Recommended, unknown if performed",
                 "Unknown; death certificate; or autopsy only (2003+)"),
  surgery_status = c("Recommended + performed", 
                     "Recommended, not performed",
                     "Not recommended/unknown",
                     "Not recommended/unknown",
                     "Recommended, not performed",
                     "Recommended, not performed",
                     "Recommended, not performed",
                     "Not recommended/unknown")
)

# Merge SEER dataset with surgical categories
seer_HN <- seer_HN %>%
  left_join(surg_map, by = c("Surgery.Performed" = "surg_range"))

# Convert surgery_status to a factor if it's not already
seer_HN$surgery_status <- factor(seer_HN$surgery_status)

# Now, relevel the factor
# seer_HN$surgery_status <- relevel(seer_HN$surgery_status, ref = "Recommended + performed")

# # # View new age categories
table(seer_HN$surgery_status, useNA = "ifany")



# Create a binary event indicator (1 = death, 0 = alive)
seer_HN$Event <- ifelse(seer_HN$Vital.status.recode..study.cutoff.used. == "Dead", 1, 0)

# Convert Time into numeric variable
seer_HN$Survival.months <- as.numeric(as.character(seer_HN$Survival.months))
seer_HN <- seer_HN[!is.na(seer_HN$Survival.months), ]

# INCOME QUINTILE
# Define income ranges and assign midpoints
income_map <- data.frame(
  income_range = c("< $40,000", "$40,000 - $44,999", "$45,000 - $49,999",
                   "$50,000 - $54,999", "$55,000 - $59,999", "$60,000 - $64,999",
                   "$65,000 - $69,999", "$70,000 - $74,999", "$75,000 - $79,999",
                   "$80,000 - $84,999", "$85,000 - $89,999", "$90,000 - $94,999",
                   "$95,000 - $99,999", "$100,000 - $109,999", "$110,000 - $119,999",
                   "$120,000+"),
  midpoint = c(35000, 42500, 47500, 52500, 57500, 62500, 67500, 72500, 77500, 
               82500, 87500, 92500, 97500, 105000, 115000, 125000)
)

seer_HN <- seer_HN %>%
  left_join(income_map, by = c("Median.household.income.inflation.adj.to.2023" = "income_range"))

quintile_cutoffs <- quantile(seer_HN$midpoint, probs = seq(0, 1, 0.2), na.rm = TRUE)

seer_HN <- seer_HN %>%
  mutate(income_quintile_equal = ntile(midpoint, 5)) %>%
  mutate(income_quintile_equal = factor(income_quintile_equal,
                                        levels = 1:5,
                                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# seer_HN$income_quintile <- relevel(seer_HN$income_quintile, ref = "Q3")

seer_HN <- seer_HN %>%
  filter(!grepl("<NA>", trimws(income_quintile_equal)))


table(seer_HN$income_quintile_equal, useNA = "ifany")

### State by State population income quintiles ### 
# 
# # 1. Create a new variable for state-specific income quintiles
# seer_HN <- seer_HN %>%
#   group_by(Medicaid.Expansion.Status) %>%
#   mutate(
#     state_income_quintile = ntile(midpoint, 5)  # Split income into quintiles within each state
#   ) %>%
#   ungroup()
# 
# # 2. Create a 'low_income' flag for the lowest 2 quintiles within each state
# seer_HN <- seer_HN %>%
#   mutate(
#     low_income_flag = ifelse(state_income_quintile %in% c(1, 2), 1, 0)
#   )
# 
# # 3. Optional: If you want to see the quintile ranges by state
# state_income_ranges <- seer_HN %>%
#   filter(!is.na(midpoint)) %>%  # Exclude NA values
#   group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
#   summarize(
#     min_income = min(midpoint, na.rm = TRUE),
#     max_income = max(midpoint, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# # 5. # View or export the quintile cutoffs table
# print(state_quintile_cutoffs)
# # write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)


# AGE GROUPS
# Create a mapping of age groups
age_map <- data.frame(
  age_range = c("20-24 years", "25-29 years", "30-34 years", "35-39 years",
                "40-44 years", "45-49 years", "50-54 years", "55-59 years", "60-64 years"),
  age_category = c("Young Adult", "Young Adult", "Young Adult", "Young Adult",
                   "Young Adult", "Young Adult", "Young Adult", "Middle-Aged",
                   "Older Adult")
)

# Merge SEER dataset with age categories
seer_HN <- seer_HN %>%
  left_join(age_map, by = c("Age" = "age_range"))

# View new age categories
table(seer_HN$age_category, useNA = "ifany")

#Truncate survival at 36 months
seer_HN <- seer_HN %>%
  mutate(capped_time = pmin(Survival.months, 36))


seer_HN$Marital_Grouped <- factor(ifelse(
  seer_HN$Marital.status.at.diagnosis %in% c("Married (including common law)", "Unmarried or Domestic Partner"), "Married",
  ifelse(seer_HN$Marital.status.at.diagnosis %in% c("Divorced", "Separated", "Widowed"), "Previously Married", "Single/Unknown")
))

table(seer_HN$Marital_Grouped)

# seer_HN$Marital_Grouped <- relevel(seer_HN$Marital_Grouped, ref = "Single/Unknown")

seer_HN$capped_time <- pmin(seer_HN$Survival.months, 36)
seer_HN$capped_event <- ifelse(seer_HN$Survival.months > 36, 0, seer_HN$Event)

#### Stage IV Excluded ####
seer_nonmet_HN <- seer_HN

seer_nonmet_HN <- seer_nonmet_HN %>%
  filter(stage_group != "IV")

### Healthcare Disparities for total Pre/Post populations ### 

seer_nonmet_HN_disp <- seer_nonmet_HN%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )

### Trying this again ###

# Step 0: Data prep

seer_nonmet_HN_disp <- seer_nonmet_HN_disp %>%
  mutate(
    Treatment_ME = factor(ifelse(Medicaid_Expansion_Status %in% c("1","2","3"), 1, 0), levels = c(0,1)),
    Time_Period  = factor(ifelse(dx.yr < exp_yr, "Pre","Post"), levels = c("Pre","Post")),
    surgery_binary = ifelse(surgery_status == "Recommended + performed", 1, 0),
    income_quintile_equal = relevel(factor(income_quintile_equal), ref = "Q3"),
    Race.Ethnicity = relevel(factor(Race.Ethnicity), ref = "Non-Hispanic White"),
    stage_group = relevel(factor(stage_group), ref = "III"),
    age_category = relevel(factor(age_category), ref = "Middle-Aged"),
    Sex = relevel(factor(Sex), ref = "Male"),
    Marital_Grouped = relevel(factor(Marital_Grouped), ref = "Single/Unknown"),
    Rural.Urban.Continuum.Code = factor(
      case_when(
        Rural.Urban.Continuum.Code == "Unknown/missing/no match (Alaska or Hawaii - Entire State)" ~ 
          "Nonmetropolitan counties adjacent to a metropolitan area",
        TRUE ~ as.character(Rural.Urban.Continuum.Code)
      )
    )
  ) %>%
  filter(!grepl("Non-Hispanic Unknown Race", trimws(Race.Ethnicity)))

# Step 1: Time stratification

cut_times <- seq(from = 12, to = 180, by = 60) # 1-15 years

stratified_data <- survSplit(
  Surv(Survival.months, Event) ~ .,
  data = seer_nonmet_HN_disp,
  cut = cut_times,
  event = "Event",
  start = "tstart",
  end = "tstop",
  zero = -1   # ensures it's < any observed survival time
)
stratified_data$time_strata <- cut(
  stratified_data$tstop,
  breaks = c(0, cut_times),
  labels = c("0-24m", "24-60m", "60-180m"),
  include.lowest = TRUE, right = TRUE
)


# Step 2: Final Cox model with time-varying covariates

final_cox_model <- coxph(
  Surv(tstart, tstop, Event) ~
    Treatment_ME + Time_Period + Treatment_ME:Time_Period +
    # The tt() function correctly models time-varying effects for your covariates.
    tt(age_category) + tt(Marital_Grouped) + tt(Sex) + tt(Race.Ethnicity) +
    tt(income_quintile_equal) + tt(Rural.Urban.Continuum.Code) +
    stage_group + surgery_status + rad_group +
    strata(time_strata),
  data = stratified_data,
  robust = TRUE,
  cluster = Medicaid.Expansion.Status
)


# Step 3: Extract and tidy ALL results for manuscript

final_results <- tidy(final_cox_model, conf.int = TRUE) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No")
  )

# Add a 'type' column to distinguish between fixed and time-varying effects
final_results <- final_results %>%
  mutate(
    type = case_when(
      grepl("tt\\(age_category\\)", term) ~ "Time-Varying: Age",
      grepl("tt\\(Marital_Grouped\\)", term) ~ "Time-Varying: Marital Status",
      grepl("tt\\(Sex\\)", term) ~ "Time-Varying: Sex",
      grepl("tt\\(Race.Ethnicity\\)", term) ~ "Time-Varying: Race",
      grepl("tt\\(income_quintile_equal\\)", term) ~ "Time-Varying: Income",
      grepl("tt\\(Rural.Urban.Continuum.Code\\)", term) ~ "Time-Varying: Rurality",
      grepl("Time_PeriodPost", term) ~ "Interaction",
      TRUE ~ "Fixed"
    )
  )

# View the full results.
print(final_results, n=Inf)
write.csv(final_results, "nonmet/hn/supplement_tt()_validation.csv", row.names = FALSE)

final_disparities_model <- coxph(
  Surv(tstart, tstop, Event) ~
    Treatment_ME + Time_Period + Treatment_ME:Time_Period +
    # We use explicit interactions to get the hazard ratios for each group
    Time_Period:age_category +
    Time_Period:Marital_Grouped +
    Time_Period:Sex +
    Time_Period:Race.Ethnicity +
    Time_Period:income_quintile_equal +
    Time_Period:Rural.Urban.Continuum.Code +
    stage_group + surgery_status + rad_group +
    strata(time_strata),
  data = stratified_data,
  robust = TRUE,
  cluster = Medicaid.Expansion.Status
)

summary(final_disparities_model)


# Step 3: Extract and tidy results for manuscript

disparities_results <- tidy(final_disparities_model, conf.int = TRUE) %>%
  filter(grepl("Time_PeriodPost", term)) %>%       # keep only Post interactions
  mutate(
    # Assign group labels for faceting
    group = case_when(
      grepl("age_category", term) ~ "Age",
      grepl("Marital_Grouped", term) ~ "Marital Status",
      grepl("Sex", term) ~ "Sex",
      grepl("Race.Ethnicity", term) ~ "Race/Ethnicity",
      grepl("income_quintile_equal", term) ~ "Income Quintile",
      grepl("Rural\\.Urban\\.Continuum\\.Code", term) ~ "Rurality",
      TRUE ~ NA_character_  # mark anything not a subgroup as NA
    )
  ) %>%
  filter(!is.na(group)) %>%  # remove Overall rows
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No"),
    
    # Clean term names for plotting: keep only subgroup identity
    term_clean = term %>%
      gsub("^Time_PeriodPost:", "", .) %>%
      gsub("age_category", "", .) %>%
      gsub("Marital_Grouped", "", .) %>%
      gsub("Sex", "", .) %>%
      gsub("Race.Ethnicity", "", .) %>%
      gsub("income_quintile_equal", "", .) %>%
      gsub("Rural\\.Urban\\.Continuum\\.Code", "", .) %>%
      trimws()
  ) %>%
  arrange(group, term_clean)



# Forest plot

ggplot(disparities_results, aes(x = HR, y = term_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = sig), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = sig), height = 0.2) +
  scale_color_manual(values = c("Yes" = "#1b7837", "No" = "#762a83")) +
  labs(
    title = "Changes in Disparities from Pre- to Post-Medicaid Expansion",
    x = "Interaction HR (Post vs Pre, relative to Reference Group)",
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

ggsave("nonmet/hn/stratified_forestplot_disparityCLN.png", width = 12, height = 8, dpi = 300)


# Print and save the final results
print(disparities_results, n = Inf)
write.csv(disparities_results, "nonmet/hn/timestratified_interaction_disparities.csv", row.names = FALSE)




