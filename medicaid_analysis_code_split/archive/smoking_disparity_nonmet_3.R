## Analysis of Pre-Post Disparities Analaysis ##


#### Load necessary libraries ####

library(survival)
library(survminer)
library(dplyr)
library(broom)
library(Matching)
library(MatchIt)
library(nnet)
library(cobalt)
library(data.table)
library(rgenoud)
library(emmeans)
library(tidyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(broom)
library(scales)
library(sandwich)
library(lmtest)
library(purrr)
library(ggpubr)
library(usmap)
library(rstpm2)
library(dplyr)
library(forcats)
library(timereg)
library(tableone)
library(bacondecomp)
library(did)

setwd("/Users/rgawdi/Documents/LungCA_ACA/SMOKING_DISPARITY")

# final
# Read the dataset
seer_esoph <- read.csv("inputs/esophagus.csv")

seer_lung <- read.csv("inputs/lung.csv")

seer_HN <- read.csv("inputs/H&N.csv")

# Combine into one dataframe
seer_smoke <- bind_rows(seer_esoph, seer_lung, seer_HN)

# Save to a new CSV (optional)
write.csv(seer_smoke, "inputs/smoking_all.csv")



#### Variable Coding ####
seer_lung <- read.csv("inputs/lung.csv")
# Define Medicaid expansion categories
seer_smoke$Medicaid.Expansion.Status <- tolower(seer_smoke$Medicaid.Expansion.Status)

seer_smoke$Medicaid_Expansion_Status <- as.factor(seer_smoke$Medicaid.Expansion.Status)

seer_smoke <- seer_smoke %>%
  mutate(dx.yr_grouped = case_when(
    dx.yr >= 2006 & dx.yr <= 2010 ~ "2006-2010",
    dx.yr >= 2011 & dx.yr <= 2013 ~ "2011-2013",
    dx.yr >= 2014 & dx.yr <= 2016 ~ "2014-2016",
    dx.yr >= 2017 & dx.yr <= 2019 ~ "2017-2019",
    TRUE ~ NA_character_
  ))
seer_smoke$dx.yr_grouped <- factor(seer_smoke$dx.yr_grouped)

seer_smoke <- seer_smoke %>%
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

# seer_smoke$Medicaid.Expansion.Status <- relevel(seer_smoke$Medicaid.Expansion.Status, ref = "0")



# Check if mapping was successful:
table(seer_smoke$Medicaid_Expansion_Status)

# Convert categorical variables into factors
seer_smoke$Age <- factor(seer_smoke$Age)  # Age grouped by 5 years
seer_smoke$Combined.Stage <- factor(seer_smoke$Combined.Stage)

stage_map <- data.frame(
  stage_range = c("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA","Stage IIIB+C", "Stage IV"),
  stage_group = c("I", "II", "II", "III","III", "IV")
) 
# Merge SEER dataset with age categories
seer_smoke <- seer_smoke %>%
  left_join(stage_map, by = c("Combined.Stage" = "stage_range"))

# # # View new age categories
table(seer_smoke$stage_group, useNA = "ifany")

seer_smoke$Sex <- factor(seer_smoke$Sex, levels = c("Male", "Female"))
seer_smoke$Median.household.income.inflation.adj.to.2023 <- factor(seer_smoke$Median.household.income.inflation.adj.to.2023)  # Median household 
seer_smoke$Rural.Urban.Continuum.Code <- factor(seer_smoke$Rural.Urban.Continuum.Code)  # Rural/urban continuum
#seer_smoke <- seer_smoke %>%
#  filter(!grepl("Unknown/missing/no match", trimws(Rural.Urban.Continuum.Code)))
seer_smoke$Marital.status.at.diagnosis <- factor(seer_smoke$Marital.status.at.diagnosis)  # Marital status
# seer_smoke$Marital.status.at.diagnosis <- relevel(seer_smoke$Marital.status.at.diagnosis, ref = "Married (including common law)") # make Married the reference

seer_smoke$Race.Ethnicity <- factor(seer_smoke$Race.Ethnicity)  # Race/ethnicity

seer_smoke$Race.Ethnicity <- relevel(
  seer_smoke$Race.Ethnicity, 
  ref = "Non-Hispanic White") # make NH White the reference
seer_smoke <- seer_smoke %>%
  filter(!grepl("Non-Hispanic Unknown Race.Ethnicity", trimws(Race.Ethnicity)))

# Group Radiation recodes
seer_smoke$Radiation.recode <- factor(seer_smoke$Radiation.recode)  # Radiation
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
seer_smoke <- seer_smoke %>%
  left_join(rad_map, by = c("Radiation.recode" = "rad_range"))

# # # View new rad categories
table(seer_smoke$rad_group, useNA = "ifany")

# seer_smoke$rad_group <- relevel(seer_smoke$rad_group, ref = "none/unknown/refused") # should radiation performed be baseline? mmm actually no i think should not be is baseline since binary. 


# Group Surgery recodes
seer_smoke$Surgery.Performed <- factor(seer_smoke$Surgery.Performed)  # Surgery Status

seer_smoke$Surgery.Performed <- factor(seer_smoke$Surgery.Performed)  # Surgery

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
seer_smoke <- seer_smoke %>%
  left_join(surg_map, by = c("Surgery.Performed" = "surg_range"))

# Convert surgery_status to a factor if it's not already
seer_smoke$surgery_status <- factor(seer_smoke$surgery_status)

# Now, relevel the factor
# seer_smoke$surgery_status <- relevel(seer_smoke$surgery_status, ref = "Recommended + performed")

# # # View new age categories
table(seer_smoke$surgery_status, useNA = "ifany")



# Create a binary event indicator (1 = death, 0 = alive)
seer_smoke$Event <- ifelse(seer_smoke$Vital.status.recode..study.cutoff.used. == "Dead", 1, 0)

# Convert Time into numeric variable
seer_smoke$Survival.months <- as.numeric(as.character(seer_smoke$Survival.months))
seer_smoke <- seer_smoke[!is.na(seer_smoke$Survival.months), ]

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

seer_smoke <- seer_smoke %>%
  left_join(income_map, by = c("Median.household.income.inflation.adj.to.2023" = "income_range"))

quintile_cutoffs <- quantile(seer_smoke$midpoint, probs = seq(0, 1, 0.2), na.rm = TRUE)

seer_smoke <- seer_smoke %>%
  mutate(income_quintile_equal = ntile(midpoint, 5)) %>%
  mutate(income_quintile_equal = factor(income_quintile_equal,
                                        levels = 1:5,
                                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# seer_smoke$income_quintile <- relevel(seer_smoke$income_quintile, ref = "Q3")

seer_smoke <- seer_smoke %>%
  filter(!grepl("<NA>", trimws(income_quintile_equal)))

# # # View result
table(seer_smoke$income_quintile_equal, useNA = "ifany")

income_ranges <- seer_smoke %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(income_quintile_equal) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )
print(income_ranges)
write.csv(income_ranges, file = "nonmet/income_quintiles.csv")

### State by State population income quintiles ### 

# 1. Create a new variable for state-specific income quintiles
seer_smoke <- seer_smoke %>%
  group_by(Medicaid.Expansion.Status) %>%
  mutate(
    state_income_quintile = ntile(midpoint, 5)  # Split income into quintiles within each state
  ) %>%
  ungroup()

state_income_ranges <- seer_smoke %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Create a 'low_income' flag for the lowest 2 quintiles within each state
seer_smoke <- seer_smoke %>%
  mutate(
    low_income_flag = ifelse(state_income_quintile %in% c(1, 2), 1, 0)
  )

# 3. Optional: If you want to see the quintile ranges by state
state_income_ranges <- seer_smoke %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )

# 5. # View or export the quintile cutoffs table
print(state_quintile_cutoffs)
# write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)


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
seer_smoke <- seer_smoke %>%
  left_join(age_map, by = c("Age" = "age_range"))

# # # View new age categories
table(seer_smoke$age_category, useNA = "ifany")

#Truncate survival at 36 months
seer_smoke <- seer_smoke %>%
  mutate(capped_time = pmin(Survival.months, 36))


seer_smoke$Marital_Grouped <- factor(ifelse(
  seer_smoke$Marital.status.at.diagnosis %in% c("Married (including common law)", "Unmarried or Domestic Partner"), "Married",
  ifelse(seer_smoke$Marital.status.at.diagnosis %in% c("Divorced", "Separated", "Widowed"), "Previously Married", "Single/Unknown")
))

table(seer_smoke$Marital_Grouped)

# seer_smoke$Marital_Grouped <- relevel(seer_smoke$Marital_Grouped, ref = "Single/Unknown")

seer_smoke$capped_time <- pmin(seer_smoke$Survival.months, 36)
seer_smoke$capped_event <- ifelse(seer_smoke$Survival.months > 36, 0, seer_smoke$Event)

#### Stage IV Excluded ####
seer_nonmet <- seer_smoke

seer_nonmet <- seer_nonmet %>%
  filter(stage_group != "IV")


#
#### Existing Mortality Hazards (non-matched) ####

seer_nonmet <- seer_nonmet %>%
  mutate(
    income_quintile_equal = relevel(factor(income_quintile_equal), ref = "Q3"),
    Race.Ethnicity = relevel(factor(
      Race.Ethnicity), ref = "Non-Hispanic White"),
    stage_group      = relevel(factor(stage_group),      ref = "III"),
    age_category     = relevel(factor(age_category),     ref = "Middle-Aged"),
    Sex = relevel(factor(Sex), ref = "Male"),
    Marital_Grouped  = relevel(factor(Marital_Grouped),  ref = "Single/Unknown"),
    Rural.Urban.Continuum.Code = relevel(factor(Rural.Urban.Continuum.Code), ref = "Counties in metropolitan areas ge 1 million pop")
  )

cox_unadj <- coxph(
  Surv(Survival.months, Event) ~ age_category + Sex + Marital_Grouped + income_quintile_equal +
    Rural.Urban.Continuum.Code + Race.Ethnicity + surgery_status +
    stage_group + rad_group + Medicaid_Expansion_Status,
  # cluster(matched_id), 
  data = seer_nonmet,
  robust = TRUE  
)

cox_unadj_tidy <- tidy(cox_unadj, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(cox_unadj_tidy)

write.csv(cox_unadj_tidy, "nonmet/mortalityHR_fullcohort.csv")


## capped at 3 years##

seer_nonmet$capped_time <- pmin(seer_nonmet$Survival.months, 36)
seer_nonmet$capped_event <- ifelse(seer_nonmet$Survival.months > 36, 0, seer_nonmet$Event)
seer_nonmet <- seer_nonmet[seer_nonmet$capped_time > 0, ]

cox_unadj_capped <- coxph(
  Surv(capped_time, capped_event) ~ age_category + Sex + Marital_Grouped + income_quintile_equal +
    Rural.Urban.Continuum.Code + Race.Ethnicity + surgery_status +
    stage_group + rad_group + Medicaid_Expansion_Status,
  # cluster(matched_id), 
  data = seer_nonmet,
  robust = TRUE  # Adds robust standard errors
)

cox_unadj_capped_tidy <- tidy(cox_unadj_capped, conf.int = TRUE) %>%
  mutate(HR = exp(estimate),  # Convert coefficients to hazard ratios
         CI_lower = exp(conf.low),
         CI_upper = exp(conf.high)) %>%
  dplyr::select(term, HR, CI_lower, CI_upper, p.value)

print(cox_unadj_capped_tidy)

write.csv(cox_unadj_capped_tidy, "nonmet/mortalityHR_fullcohort_3yrCapped.csv")


## Univariate Model ## 
univar_results <- lapply(
  c("age_category", "Sex", "income_quintile_equal", "stage_group", 
    "Marital_Grouped", "Rural.Urban.Continuum.Code", "Race.Ethnicity"),
  function(var) {
    fit <- coxph(as.formula(paste0("Surv(Survival.months, Event) ~ ", var)), data = seer_nonmet)
    tidy(fit, conf.int = TRUE) %>%
      mutate(
        Variable_Group = var,
        HR = exp(estimate),
        CI_lower = exp(conf.low),
        CI_upper = exp(conf.high),
        CI = paste0(round(CI_lower, 2), "-", round(CI_upper, 2)),
        p.value = signif(p.value, 3)
      ) %>%
      dplyr::select(Variable_Group, term, HR, CI, p.value)
  }
)

univar_results <- do.call(rbind, univar_results)

write.csv(univar_results, "nonmet/univariable_results.csv", row.names = FALSE)


## Multivariable analysis ##
multi_fit <- coxph(Surv(capped_time, capped_event) ~ age_category + Sex + Marital_Grouped + income_quintile_equal +
                     Rural.Urban.Continuum.Code + Race.Ethnicity + surgery_status +
                     stage_group + rad_group + Medicaid_Expansion_Status,
                   data = seer_nonmet)

# Tidy and compute HR + 95% CI

multivar_results <- tidy(multi_fit, conf.int = TRUE) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    CI = paste0(round(CI_lower,2), "-", round(CI_upper,2)),
    p.value = signif(p.value, 3)
  ) %>%
  dplyr::select(term, HR, CI, p.value)

# Optional: rename term to include variable group
multivar_results <- multivar_results %>%
  mutate(Variable_Group = sub("(\\w+)_.*", "\\1", term)) %>%
  dplyr::select(Variable_Group, term, HR, CI, p.value)


# View results
print(multivar_results)

write.csv(multivar_results, "nonmet/multivariable_results.csv", row.names = FALSE)

#### Healthcare Disparities Alleviation for Total Population #### 

seer_disp <- seer_nonmet%>%
  mutate(
    exp_yr = case_when(
      Medicaid_Expansion_Status == "0" ~ NA_real_,
      Medicaid_Expansion_Status == "1" ~ 2011,
      Medicaid_Expansion_Status == "2" ~ 2014,
      Medicaid_Expansion_Status == "3" ~ 2017
    )
  )

seer_disp$income_quintile_equal <- relevel(seer_disp$income_quintile_equal, ref = "Q5")


seer_disp <- seer_disp %>%
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

library(car)

# Just test Treatment_ME and Time_Period
col_test <- lm(as.numeric(Treatment_ME) ~ Time_Period, data = seer_disp)
summary(col_test)

# Variance inflation factors (VIF)
vif(lm(as.numeric(Treatment_ME) ~ Time_Period, data = seer_disp))
vif(lm(as.numeric(Time_Period) ~ Treatment_ME, data = seer_disp))

# Step 1: Time stratification
# Recode 0 months to 0.5 months
seer_disp$Survival.months <- ifelse(seer_disp$Survival.months == 0,
                                    0.5,
                                    seer_disp$Survival.months)

# Now survSplit works
cut_times <- seq(from = 12, to = 180, by = 60)
stratified_data <- survSplit(
  Surv(Survival.months, Event) ~ .,
  data = seer_disp,
  cut = cut_times,
  event = "Event",
  start = "tstart",
  end = "tstop",
  id = "id",
  zero = 0
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
    age_category + Marital_Grouped + Sex + Race.Ethnicity +
    tt(income_quintile_equal) + tt(Rural.Urban.Continuum.Code) +
    tt(stage_group) + tt(surgery_status) + tt(rad_group) +
    strata(time_strata),
  data = stratified_data,
  robust = TRUE,
  cluster = Medicaid.Expansion.Status
)


final_cox_model_clean <- coxph(
  Surv(tstart, tstop, Event) ~
    Treatment_ME * Time_Period +
    age_category + Marital_Grouped + Sex + Race.Ethnicity +
    income_quintile_equal + Rural.Urban.Continuum.Code +
    strata(stage_group, surgery_status, rad_group, time_strata),
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
      grepl("tt\\(income_quintile_equal\\)", term) ~ "Time-Varying: Income",
      grepl("tt\\(Rural.Urban.Continuum.Code\\)", term) ~ "Time-Varying: Rurality",
      grepl("tt\\(stage_group\\)", term) ~ "Time-Varying: Stage",
      grepl("tt\\(surgery_status\\)", term) ~ "Time-Varying: Surgery Status",
      grepl("tt\\(rad_group\\)", term) ~ "Time-Varying: Radiation",
      grepl("Time_PeriodPost", term) ~ "Interaction",
      TRUE ~ "Fixed"
    )
  )

# View the full results.
print(final_results, n=Inf)
write.csv(final_results, "nonmet/supplement_tt()_validation.csv", row.names = FALSE)

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
) ### Replace Time_Period with Medicaid_Expansion_Status ### 

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

ggsave("nonmet/stratified_forestplot_disparityCLN.png", width = 12, height = 8, dpi = 300)




#### disparities part 2 ####
#### Healthcare Disparities: Triple-Interaction DiD ####

library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(stringr)

# -------------------------------------------------------------------
# 1. Prepare common stage I-III cohort
# -------------------------------------------------------------------

seer_disp_base <- seer_nonmet %>%
  mutate(
    patient_id = row_number(),
    
    income_quintile_equal = relevel(
      factor(income_quintile_equal),
      ref = "Q5"
    ),
    
    Race.Ethnicity = relevel(
      factor(Race.Ethnicity),
      ref = "Non-Hispanic White"
    ),
    
    stage_group = relevel(
      factor(stage_group),
      ref = "III"
    ),
    
    age_category = relevel(
      factor(age_category),
      ref = "Middle-Aged"
    ),
    
    Sex = relevel(
      factor(Sex),
      ref = "Male"
    ),
    
    Marital_Grouped = relevel(
      factor(Marital_Grouped),
      ref = "Married"
    ),
    
    Rural.Urban.Continuum.Code = factor(
      case_when(
        Rural.Urban.Continuum.Code ==
          "Unknown/missing/no match (Alaska or Hawaii - Entire State)" ~
          "Nonmetropolitan counties adjacent to a metropolitan area",
        TRUE ~ as.character(Rural.Urban.Continuum.Code)
      )
    ),
    
    Survival.months = ifelse(
      Survival.months <= 0,
      0.5,
      Survival.months
    )
  ) %>%
  filter(
    !grepl(
      "Non-Hispanic Unknown Race",
      trimws(Race.Ethnicity)
    )
  )


# -------------------------------------------------------------------
# 2. Function to construct each valid DiD comparison
# -------------------------------------------------------------------

make_did_stack <- function(data, treated_group, expansion_year, stack_name) {
  
  data %>%
    filter(
      as.character(Medicaid_Expansion_Status) %in%
        c("0", treated_group)
    ) %>%
    
    mutate(
      stack_id = stack_name,
      
      Treatment_ME = factor(
        ifelse(
          as.character(Medicaid_Expansion_Status) == treated_group,
          "Expansion",
          "Control"
        ),
        levels = c("Control", "Expansion")
      ),
      
      # IMPORTANT:
      # BOTH expansion and control patients use the same cutoff
      Time_Period = factor(
        ifelse(
          dx.yr < expansion_year,
          "Pre",
          "Post"
        ),
        levels = c("Pre", "Post")
      )
    )
}


# -------------------------------------------------------------------
# 3. Construct the 3 staggered DiD comparisons
# -------------------------------------------------------------------

disp_early <- make_did_stack(
  seer_disp_base,
  treated_group = "1",
  expansion_year = 2011,
  stack_name = "Early Expansion"
)

disp_2014 <- make_did_stack(
  seer_disp_base,
  treated_group = "2",
  expansion_year = 2014,
  stack_name = "2014 Expansion"
)

disp_late <- make_did_stack(
  seer_disp_base,
  treated_group = "3",
  expansion_year = 2017,
  stack_name = "Late Expansion"
)


# Stack all three comparisons
stratified_data <- bind_rows(
  disp_early,
  disp_2014,
  disp_late
)


# Sanity check: THERE SHOULD BE NO NAs HERE
print(
  table(
    stratified_data$Treatment_ME,
    stratified_data$Time_Period,
    useNA = "ifany"
  )
)


# -------------------------------------------------------------------
# 4. Triple-interaction Cox model
# -------------------------------------------------------------------

final_disparities_model <- coxph(
  Surv(Survival.months, Event) ~
    
    Treatment_ME * Time_Period * (
      age_category +
        Marital_Grouped +
        Sex +
        Race.Ethnicity +
        income_quintile_equal +
        Rural.Urban.Continuum.Code
    ) +
    
    stage_group +
    dx.yr +
    
    # Separate baseline hazard for each expansion comparison
    strata(stack_id),
  
  data = stratified_data,
  
  # Controls appear in multiple stacks and patients within states
  # are correlated; retain state-level robust clustering
  robust = TRUE,
  cluster = Medicaid.Expansion.Status
)

summary(final_disparities_model)


# -------------------------------------------------------------------
# 5. Extract ONLY actual 3-way interaction terms
# -------------------------------------------------------------------

all_disparity_terms <- tidy(
  final_disparities_model,
  conf.int = TRUE
)


disparities_results <- all_disparity_terms %>%
  
  # A true 3-way term must contain BOTH policy components
  # AND one of the subgroup variables
  filter(
    str_detect(term, "Treatment_MEExpansion"),
    str_detect(term, "Time_PeriodPost"),
    str_detect(
      term,
      paste(
        "age_category",
        "Marital_Grouped",
        "Sex",
        "Race.Ethnicity",
        "income_quintile_equal",
        "Rural\\.Urban\\.Continuum\\.Code",
        sep = "|"
      )
    )
  ) %>%
  
  mutate(
    group = case_when(
      str_detect(term, "age_category") ~ "Age",
      str_detect(term, "Marital_Grouped") ~ "Marital Status",
      str_detect(term, "Sex") ~ "Sex",
      str_detect(term, "Race.Ethnicity") ~ "Race/Ethnicity",
      str_detect(term, "income_quintile_equal") ~ "Income Quintile",
      str_detect(term, "Rural\\.Urban\\.Continuum\\.Code") ~ "Rurality",
      TRUE ~ NA_character_
    ),
    
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    
    sig = ifelse(
      p.value < 0.05,
      "Yes",
      "No"
    ),
    
    term_clean = term %>%
      str_replace_all("Treatment_MEExpansion", "") %>%
      str_replace_all("Time_PeriodPost", "") %>%
      str_replace_all("age_category", "") %>%
      str_replace_all("Marital_Grouped", "") %>%
      str_replace_all("Sex", "") %>%
      str_replace_all("Race.Ethnicity", "") %>%
      str_replace_all("income_quintile_equal", "") %>%
      str_replace_all(
        "Rural\\.Urban\\.Continuum\\.Code",
        ""
      ) %>%
      str_replace_all(":", "") %>%
      str_trim()
  ) %>%
  
  filter(!is.na(group)) %>%
  arrange(group, term_clean)

income_quintile_equal = relevel(
  factor(income_quintile_equal),
  ref = "Q5"
)

Marital_Grouped = relevel(
  factor(Marital_Grouped),
  ref = "Married"
)

triple_check_clean <- triple_check %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high)
  ) %>%
  dplyr::select(
    term,
    HR,
    CI_lower,
    CI_upper,
    p.value
  )

print(
  triple_check_clean,
  n = Inf,
  width = Inf
)

View(disparities_results)

write.csv(
  disparities_results,
  "nonmet/figure2_triple_interaction_results.csv",
  row.names = FALSE
)


triple_check_clean <- triple_check %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high)
  ) %>%
  dplyr::select(
    term,
    HR,
    CI_lower,
    CI_upper,
    p.value
  )

print(
  triple_check_clean,
  n = Inf,
  width = Inf
)
write.csv(
  triple_check_clean,
  "nonmet/TRUE_triple_interaction_raw_results.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------
# 6. Forest plot
# -------------------------------------------------------------------

ggplot(
  disparities_results,
  aes(
    x = HR,
    y = term_clean
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey50"
  ) +
  
  geom_point(
    aes(color = sig),
    size = 3
  ) +
  
  geom_errorbarh(
    aes(
      xmin = CI_lower,
      xmax = CI_upper,
      color = sig
    ),
    height = 0.2
  ) +
  
  scale_color_manual(
    values = c(
      "Yes" = "#1b7837",
      "No" = "#762a83"
    )
  ) +
  
  labs(
    title = "Differential Association of Medicaid Expansion With Survival Disparities",
    x = "Three-Way Interaction HR",
    y = NULL,
    color = "Statistically Significant"
  ) +
  
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free",
    switch = "y"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    strip.text.y.left =
      element_text(
        angle = 0,
        size = 13,
        hjust = 0.5
      ),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )


ggsave(
  "nonmet/stratified_forestplot_disparityCLN_triple.png",
  width = 12,
  height = 8,
  dpi = 300
)

#### RP model####
library(rstpm2)

# Step 1: Royston-Parmar flexible parametric model
# Use restricted cubic splines (df = 5 is typical, adjust for smoothness)
rp_model <- stpm2(
  Surv(Survival.months, Event) ~ 
    Treatment_ME * Time_Period +
    (age_category + Marital_Grouped + Sex + Race.Ethnicity +
       income_quintile_equal + Rural.Urban.Continuum.Code) * Time_Period +
    stage_group + surgery_status + rad_group,
  data = seer_disp,
  df = 5,              # degrees of freedom for baseline spline
  cluster = seer_disp$Medicaid.Expansion.Status, 
  robust = TRUE
)

summary(rp_model)

write.csv(rp_model, "nonmet/roypalm_validation_disparities.csv")

# Step 2: Extract interaction terms
rp_interactions <- tidy(rp_model, conf.int = TRUE) %>%
  filter(grepl("Time_PeriodPost", term)) %>%
  mutate(
    HR = exp(estimate),
    CI_lower = exp(conf.low),
    CI_upper = exp(conf.high),
    sig = ifelse(p.value < 0.05, "Yes", "No")
  )

# Step 3: Forest plot
ggplot(rp_interactions, aes(x = HR, y = term)) +
  geom_point(aes(color = sig), size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, color = sig), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("Yes" = "hotpink", "No" = "orange3")) +
  labs(
    title = "Royston–Parmar Model: Pre vs Post Expansion",
    x = "Interaction HR",
    y = NULL
  ) +
  theme_minimal()
 ggsave("nonmet/forest_plot_RoystonParmar_interactions.png", width = 14, height = 10, dpi = 300)

capture.output(summary(rp_model), file = "nonmet/roypalm_validation_disparities.txt")

# install.packages("timereg")
# 
# 
# # Step 1: Aalen's additive hazard model
# aalen_model <- aalen(
#   Surv(Survival.months, Event) ~ 
#     Treatment_ME * Time_Period +
#     (age_category + Marital_Grouped + Sex + Race.Ethnicity +
#        income_quintile_equal + Rural.Urban.Continuum.Code) * Time_Period +
#     stage_group + surgery_status + rad_group,
#   data = seer_disp,
#   robust = TRUE,
#   n.sim = 1000   # optional: to get resampled CI bands
# )
# 
# summary(aalen_model)
# 
# # Step 2: Extract interaction terms
# aalen_df <- as.data.frame(aalen_model$cum) %>%
#   filter(grepl("Time_PeriodPost", names(.)))
# 
# # Step 3: Plot cumulative effects (more standard than HRs)
# plot(aalen_model, sim.ci = TRUE, col = 1:3)




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




#### Lung Analysis ####
# Define Medicaid expansion categories
seer_lung$Medicaid.Expansion.Status <- tolower(seer_lung$Medicaid.Expansion.Status)

seer_lung$Medicaid_Expansion_Status <- as.factor(seer_lung$Medicaid.Expansion.Status)

seer_lung <- seer_lung %>%
  mutate(dx.yr_grouped = case_when(
    dx.yr >= 2006 & dx.yr <= 2010 ~ "2006-2010",
    dx.yr >= 2011 & dx.yr <= 2013 ~ "2011-2013",
    dx.yr >= 2014 & dx.yr <= 2016 ~ "2014-2016",
    dx.yr >= 2017 & dx.yr <= 2019 ~ "2017-2019",
    TRUE ~ NA_character_
  ))
seer_lung$dx.yr_grouped <- factor(seer_lung$dx.yr_grouped)

seer_lung <- seer_lung %>%
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

# seer_lung$Medicaid.Expansion.Status <- relevel(seer_lung$Medicaid.Expansion.Status, ref = "0")



# Check if mapping was successful:
table(seer_lung$Medicaid_Expansion_Status)

# Convert categorical variables into factors
seer_lung$Age <- factor(seer_lung$Age)  # Age grouped by 5 years
seer_lung$Combined.Stage <- factor(seer_lung$Combined.Stage)

stage_map <- data.frame(
  stage_range = c("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA","Stage IIIB+C", "Stage IV"),
  stage_group = c("I", "II", "II", "III","III", "IV")
) 
# Merge SEER dataset with age categories
seer_lung <- seer_lung %>%
  left_join(stage_map, by = c("Combined.Stage" = "stage_range"))

# # # View new age categories
table(seer_lung$stage_group, useNA = "ifany")

seer_lung$Sex <- factor(seer_lung$Sex, levels = c("Male", "Female"))
seer_lung$Median.household.income.inflation.adj.to.2023 <- factor(seer_lung$Median.household.income.inflation.adj.to.2023)  # Median household 
seer_lung$Rural.Urban.Continuum.Code <- factor(seer_lung$Rural.Urban.Continuum.Code)  # Rural/urban continuum
#seer_lung <- seer_lung %>%
#  filter(!grepl("Unknown/missing/no match", trimws(Rural.Urban.Continuum.Code)))
seer_lung$Marital.status.at.diagnosis <- factor(seer_lung$Marital.status.at.diagnosis)  # Marital status
# seer_lung$Marital.status.at.diagnosis <- relevel(seer_lung$Marital.status.at.diagnosis, ref = "Married (including common law)") # make Married the reference

seer_lung$Race.Ethnicity <- factor(seer_lung$Race.Ethnicity)  # Race/ethnicity
# seer_lung$Race.Ethnicity <- relevel(seer_lung$Race.Ethnicity, ref = "Non-Hispanic White") # make NH White the reference

seer_lung <- seer_lung %>%
  filter(!grepl("Non-Hispanic Unknown Race", trimws(Race.Ethnicity)))

# Group Radiation recodes
seer_lung$Radiation.recode <- factor(seer_lung$Radiation.recode)  # Radiation
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
seer_lung <- seer_lung %>%
  left_join(rad_map, by = c("Radiation.recode" = "rad_range"))

# # # View new rad categories
table(seer_lung$rad_group, useNA = "ifany")

# seer_lung$rad_group <- relevel(seer_lung$rad_group, ref = "none/unknown/refused") # should radiation performed be baseline? mmm actually no i think should not be is baseline since binary. 


# Group Surgery recodes
seer_lung$Surgery.Performed <- factor(seer_lung$Surgery.Performed)  # Surgery Status

seer_lung$Surgery.Performed <- factor(seer_lung$Surgery.Performed)  # Surgery

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
seer_lung <- seer_lung %>%
  left_join(surg_map, by = c("Surgery.Performed" = "surg_range"))

# Convert surgery_status to a factor if it's not already
seer_lung$surgery_status <- factor(seer_lung$surgery_status)

# Now, relevel the factor
# seer_lung$surgery_status <- relevel(seer_lung$surgery_status, ref = "Recommended + performed")

# # # View new age categories
table(seer_lung$surgery_status, useNA = "ifany")



# Create a binary event indicator (1 = death, 0 = alive)
seer_lung$Event <- ifelse(seer_lung$Vital.status.recode..study.cutoff.used. == "Dead", 1, 0)

# Convert Time into numeric variable
seer_lung$Survival.months <- as.numeric(as.character(seer_lung$Survival.months))
seer_lung <- seer_lung[!is.na(seer_lung$Survival.months), ]

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

seer_lung <- seer_lung %>%
  left_join(income_map, by = c("Median.household.income.inflation.adj.to.2023" = "income_range"))

quintile_cutoffs <- quantile(seer_lung$midpoint, probs = seq(0, 1, 0.2), na.rm = TRUE)

seer_lung <- seer_lung %>%
  mutate(income_quintile_equal = ntile(midpoint, 5)) %>%
  mutate(income_quintile_equal = factor(income_quintile_equal,
                                        levels = 1:5,
                                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# seer_lung$income_quintile <- relevel(seer_lung$income_quintile, ref = "Q3")

seer_lung <- seer_lung %>%
  filter(!grepl("<NA>", trimws(income_quintile_equal)))

# # # View result
table(seer_lung$income_quintile_equal, useNA = "ifany")

### State by State population income quintiles ### 

# 1. Create a new variable for state-specific income quintiles
seer_lung <- seer_lung %>%
  group_by(Medicaid.Expansion.Status) %>%
  mutate(
    state_income_quintile = ntile(midpoint, 5)  # Split income into quintiles within each state
  ) %>%
  ungroup()

# 2. Create a 'low_income' flag for the lowest 2 quintiles within each state
seer_lung <- seer_lung %>%
  mutate(
    low_income_flag = ifelse(state_income_quintile %in% c(1, 2), 1, 0)
  )

# 3. Optional: If you want to see the quintile ranges by state
state_income_ranges <- seer_lung %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )

# 5. # View or export the quintile cutoffs table
print(state_quintile_cutoffs)
# write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)


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
seer_lung <- seer_lung %>%
  left_join(age_map, by = c("Age" = "age_range"))

# # # View new age categories
table(seer_lung$age_category, useNA = "ifany")

#Truncate survival at 36 months
seer_lung <- seer_lung %>%
  mutate(capped_time = pmin(Survival.months, 36))


seer_lung$Marital_Grouped <- factor(ifelse(
  seer_lung$Marital.status.at.diagnosis %in% c("Married (including common law)", "Unmarried or Domestic Partner"), "Married",
  ifelse(seer_lung$Marital.status.at.diagnosis %in% c("Divorced", "Separated", "Widowed"), "Previously Married", "Single/Unknown")
))

table(seer_lung$Marital_Grouped)

# seer_lung$Marital_Grouped <- relevel(seer_lung$Marital_Grouped, ref = "Single/Unknown")

seer_lung$capped_time <- pmin(seer_lung$Survival.months, 36)
seer_lung$capped_event <- ifelse(seer_lung$Survival.months > 36, 0, seer_lung$Event)


#### Stage IV Excluded ####
seer_nonmet_lung <- seer_lung

seer_nonmet_lung <- seer_nonmet_lung %>%
  filter(stage_group != "IV")

### Healthcare Disparities for total Pre/Post populations ### 


seer_nonmet_lung_disp <- seer_nonmet_lung%>%
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

seer_nonmet_lung_disp <- seer_nonmet_lung_disp %>%
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
  data = seer_nonmet_lung_disp,
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
write.csv(final_results, "nonmet/lung/supplement_tt()_validation.csv", row.names = FALSE)

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

ggsave("nonmet/lung/stratified_forestplot_disparityCLN.png", width = 12, height = 8, dpi = 300)


# Print and save the final results
print(disparities_results, n = Inf)
write.csv(disparities_results, "nonmet/lung/timestratified_interaction_disparities.csv", row.names = FALSE)



#### Esophagus Analysis ####
# Define Medicaid expansion categories
seer_esoph$Medicaid.Expansion.Status <- tolower(seer_esoph$Medicaid.Expansion.Status)

seer_esoph$Medicaid_Expansion_Status <- as.factor(seer_esoph$Medicaid.Expansion.Status)

seer_esoph <- seer_esoph %>%
  mutate(dx.yr_grouped = case_when(
    dx.yr >= 2006 & dx.yr <= 2010 ~ "2006-2010",
    dx.yr >= 2011 & dx.yr <= 2013 ~ "2011-2013",
    dx.yr >= 2014 & dx.yr <= 2016 ~ "2014-2016",
    dx.yr >= 2017 & dx.yr <= 2019 ~ "2017-2019",
    TRUE ~ NA_character_
  ))
seer_esoph$dx.yr_grouped <- factor(seer_esoph$dx.yr_grouped)

seer_esoph <- seer_esoph %>%
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

# seer_esoph$Medicaid.Expansion.Status <- relevel(seer_esoph$Medicaid.Expansion.Status, ref = "0")



# Check if mapping was successful:
table(seer_esoph$Medicaid_Expansion_Status)

# Convert categorical variables into factors
seer_esoph$Age <- factor(seer_esoph$Age)  # Age grouped by 5 years
seer_esoph$Combined.Stage <- factor(seer_esoph$Combined.Stage)

stage_map <- data.frame(
  stage_range = c("Stage I", "Stage IIA", "Stage IIB", "Stage IIIA","Stage IIIB+C", "Stage IV"),
  stage_group = c("I", "II", "II", "III","III", "IV")
) 
# Merge SEER dataset with age categories
seer_esoph <- seer_esoph %>%
  left_join(stage_map, by = c("Combined.Stage" = "stage_range"))

# # # View new age categories
table(seer_esoph$stage_group, useNA = "ifany")

seer_esoph$Sex <- factor(seer_esoph$Sex, levels = c("Male", "Female"))
seer_esoph$Median.household.income.inflation.adj.to.2023 <- factor(seer_esoph$Median.household.income.inflation.adj.to.2023)  # Median household 
seer_esoph$Rural.Urban.Continuum.Code <- factor(seer_esoph$Rural.Urban.Continuum.Code)  # Rural/urban continuum
#seer_esoph <- seer_esoph %>%
#  filter(!grepl("Unknown/missing/no match", trimws(Rural.Urban.Continuum.Code)))
seer_esoph$Marital.status.at.diagnosis <- factor(seer_esoph$Marital.status.at.diagnosis)  # Marital status
# seer_esoph$Marital.status.at.diagnosis <- relevel(seer_esoph$Marital.status.at.diagnosis, ref = "Married (including common law)") # make Married the reference

seer_esoph$Race.Ethnicity <- factor(seer_esoph$Race.Ethnicity)  # Race/ethnicity
# seer_esoph$Race.Ethnicity <- relevel(seer_esoph$Race.Ethnicity, ref = "Non-Hispanic White") # make NH White the reference
seer_esoph <- seer_esoph %>%
  filter(!grepl("Non-Hispanic Unknown Race", trimws(Race.Ethnicity)))

# Group Radiation recodes
seer_esoph$Radiation.recode <- factor(seer_esoph$Radiation.recode)  # Radiation
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
seer_esoph <- seer_esoph %>%
  left_join(rad_map, by = c("Radiation.recode" = "rad_range"))

# # # View new rad categories
table(seer_esoph$rad_group, useNA = "ifany")

# seer_esoph$rad_group <- relevel(seer_esoph$rad_group, ref = "none/unknown/refused") # should radiation performed be baseline? mmm actually no i think should not be is baseline since binary. 


# Group Surgery recodes
seer_esoph$Surgery.Performed <- factor(seer_esoph$Surgery.Performed)  # Surgery Status

seer_esoph$Surgery.Performed <- factor(seer_esoph$Surgery.Performed)  # Surgery

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
seer_esoph <- seer_esoph %>%
  left_join(surg_map, by = c("Surgery.Performed" = "surg_range"))

# Convert surgery_status to a factor if it's not already
seer_esoph$surgery_status <- factor(seer_esoph$surgery_status)

# Now, relevel the factor
# seer_esoph$surgery_status <- relevel(seer_esoph$surgery_status, ref = "Recommended + performed")

# # # View new age categories
table(seer_esoph$surgery_status, useNA = "ifany")



# Create a binary event indicator (1 = death, 0 = alive)
seer_esoph$Event <- ifelse(seer_esoph$Vital.status.recode..study.cutoff.used. == "Dead", 1, 0)

# Convert Time into numeric variable
seer_esoph$Survival.months <- as.numeric(as.character(seer_esoph$Survival.months))
seer_esoph <- seer_esoph[!is.na(seer_esoph$Survival.months), ]

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

seer_esoph <- seer_esoph %>%
  left_join(income_map, by = c("Median.household.income.inflation.adj.to.2023" = "income_range"))

quintile_cutoffs <- quantile(seer_esoph$midpoint, probs = seq(0, 1, 0.2), na.rm = TRUE)

seer_esoph <- seer_esoph %>%
  mutate(income_quintile_equal = ntile(midpoint, 5)) %>%
  mutate(income_quintile_equal = factor(income_quintile_equal,
                                        levels = 1:5,
                                        labels = c("Q1", "Q2", "Q3", "Q4", "Q5")))

# seer_esoph$income_quintile <- relevel(seer_esoph$income_quintile, ref = "Q3")

seer_esoph <- seer_esoph %>%
  filter(!grepl("<NA>", trimws(income_quintile_equal)))

# # # View result
table(seer_esoph$income_quintile_equal, useNA = "ifany")

### State by State population income quintiles ### 

# 1. Create a new variable for state-specific income quintiles
seer_esoph <- seer_esoph %>%
  group_by(Medicaid.Expansion.Status) %>%
  mutate(
    state_income_quintile = ntile(midpoint, 5)  # Split income into quintiles within each state
  ) %>%
  ungroup()


# 2. Create a 'low_income' flag for the lowest 2 quintiles within each state
seer_esoph <- seer_esoph %>%
  mutate(
    low_income_flag = ifelse(state_income_quintile %in% c(1, 2), 1, 0)
  )

# 3. Optional: If you want to see the quintile ranges by state
state_income_ranges <- seer_esoph %>%
  filter(!is.na(midpoint)) %>%  # Exclude NA values
  group_by(Medicaid.Expansion.Status, state_income_quintile) %>%
  summarize(
    min_income = min(midpoint, na.rm = TRUE),
    max_income = max(midpoint, na.rm = TRUE),
    .groups = "drop"
  )

# 5. # View or export the quintile cutoffs table
print(state_quintile_cutoffs)
# write.csv(state_quintile_cutoffs, "state_income_quintile_cutoffs.csv", row.names = FALSE)


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
seer_esoph <- seer_esoph %>%
  left_join(age_map, by = c("Age" = "age_range"))

# # # View new age categories
table(seer_esoph$age_category, useNA = "ifany")

#Truncate survival at 36 months
seer_esoph <- seer_esoph %>%
  mutate(capped_time = pmin(Survival.months, 36))


seer_esoph$Marital_Grouped <- factor(ifelse(
  seer_esoph$Marital.status.at.diagnosis %in% c("Married (including common law)", "Unmarried or Domestic Partner"), "Married",
  ifelse(seer_esoph$Marital.status.at.diagnosis %in% c("Divorced", "Separated", "Widowed"), "Previously Married", "Single/Unknown")
))

table(seer_esoph$Marital_Grouped)

# seer_esoph$Marital_Grouped <- relevel(seer_esoph$Marital_Grouped, ref = "Single/Unknown")

seer_esoph$capped_time <- pmin(seer_esoph$Survival.months, 36)
seer_esoph$capped_event <- ifelse(seer_esoph$Survival.months > 36, 0, seer_esoph$Event)



#### Stage IV Excluded ####
seer_nonmet_esoph <- seer_esoph

seer_nonmet_esoph <- seer_nonmet_esoph %>%
  filter(stage_group != "IV")

### Healthcare Disparities for total Pre/Post populations ### 

seer_nonmet_esophagus_disp <- seer_nonmet_esoph%>%
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

seer_nonmet_esophagus_disp <- seer_nonmet_esophagus_disp %>%
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
  data = seer_nonmet_esophagus_disp,
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
write.csv(final_results, "nonmet/esophagus/supplement_tt()_validation.csv", row.names = FALSE)

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

ggsave("nonmet/esophagus/stratified_forestplot_disparityCLN.png", width = 12, height = 8, dpi = 300)


# Print and save the final results
print(disparities_results, n = Inf)
write.csv(disparities_results, "nonmet/esophagus/timestratified_interaction_disparities.csv", row.names = FALSE)




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







#### Figures ####
# Required packages
if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

# Panel A: Conceptual Timeline
#################################

panelA <- ggplot() +
  # Full study period line (now drawn between y=1 and y=1.1 so it has "height")
  geom_segment(aes(x = 2006, xend = 2019, y = 1, yend = 1), size = 1, color = "black") +
  
  # Vertical ACA expansion marker
  geom_vline(xintercept = 2014, linetype = "dashed", color = "gray40", size = 1) +
  
  # Labels for eras
  annotate("text", x = 2010, y = 1.2, label = "Pre-Expansion", size = 4, fontface = "bold") +
  annotate("text", x = 2016.5, y = 1.2, label = "Post-Expansion", size = 4, fontface = "bold") +
  
  # ACA expansion label
  annotate("text", x = 2014, y = 1.4, label = "ACA Expansion\n(Jan 1, 2014)",
           hjust = 0.5, size = 4) +
  
  # Force y-axis to include the line & labels
  scale_y_continuous(limits = c(0.8, 1.5), breaks = NULL) +
  
  # X-axis for years
  scale_x_continuous(breaks = seq(2006, 2019, 2), limits = c(2006, 2019)) +
  
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Panel A. Conceptual Study Timeline and Policy Event")

print(panelA)

#################################
# Panel B: Cohort-Specific Timelines
#################################

# Define bars for each cohort (correct ordering, no swaps)
panelB_data <- data.frame(
  cohort = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion"),
  pre_start = c(2006, 2006, 2006, 2006),
  pre_end   = c(2020, 2011, 2014, 2017),
  post_start = c(NA, 2011, 2014, 2017),
  post_end   = c(NA, 2020, 2020, 2020)
)

# Preserve cohort order explicitly
panelB_data$cohort <- factor(panelB_data$cohort,
                             levels = c("Non-Expansion", "Early Expansion", "2014 Expansion", "Late Expansion"))

panelB <- ggplot() +
  # Pre-expansion bars
  geom_rect(data = panelB_data,
            aes(xmin = pre_start, xmax = pre_end,
                ymin = as.numeric(cohort) - 0.4,
                ymax = as.numeric(cohort) + 0.4),
            fill = "#E31A1C") +
  
  # Post-expansion bars
  geom_rect(data = subset(panelB_data, !is.na(post_start)),
            aes(xmin = post_start, xmax = post_end,
                ymin = as.numeric(cohort) - 0.4,
                ymax = as.numeric(cohort) + 0.4),
            fill = "#1F78B4") +
  
  # ACA expansion marker
  geom_vline(xintercept = 2014, linetype = "dashed", color = "gray40", size = 1) +
  
  scale_x_continuous(breaks = seq(2006, 2020, 2), limits = c(2006, 2020)) +
  scale_y_continuous(breaks = 1:4, labels = levels(panelB_data$cohort), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Panel B. Cohort-Specific Time Periods")

#################################
# Print panels
#################################
print(panelA)
print(panelB)

# Save as needed:
ggsave("nonmet/panelA_timeline.png", panelA, width = 7, height = 2.5, dpi = 300)
ggsave("nonmet/panelB_cohorts.png", panelB, width = 7, height = 4, dpi = 300)



library(ggplot2)
library(patchwork)


# ---------- Panel A: Conceptual Study Timeline ----------
pA <- ggplot() +
  geom_segment(aes(x = 2006, xend = 2020, y = 0, yend = 0), size = 0.8) +
  geom_vline(xintercept = 2014, linetype = "dashed", size = 0.6) +
  annotate("text", x = 2014, y = 0.04, label = "ACA Expansion\n(Jan 1, 2014)", 
           hjust = -.2, size = 2) +
  geom_vline(xintercept = 2011, linetype = "dotted", size = 0.5) +
  annotate("text", x = 2011, y = 0.04, label = "Early Expansion\n(Section 1115 Waivers)", 
           hjust = -.1, size = 2) +
  geom_vline(xintercept = 2017, linetype = "dotted", size = 0.5) +
  annotate("text", x = 2017, y = 0.04, label = "Late Expansion\n(State Legislation)", 
           hjust = -.1, size = 2) +
  annotate("text", x = 2008, y = 0.25, label = "Pre-Expansion", 
           hjust = 0.3, fontface = "bold", size = 3) +
  annotate("text", x = 2019, y = 0.25, label = "Post-Expansion", 
           hjust = 0.7, fontface = "bold", size = 3) +
  scale_x_continuous(breaks = seq(2006, 2020, 2), limits = c(2006, 2020)) +
  scale_y_continuous(limits = c(-0.05, 0.25)) +
  theme_minimal() +
  labs(title = "Conceptual Study Timeline and Policy Event") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.margin = margin(2, 2, 2, 2)
  )

pA
# ---------- Panel B: Cohort-Specific Periods ----------
cohort_data <- data.frame(
  cohort = c("Non-Expansion", "Early Expansion", "On-Time Expansion (2014)", "Late Expansion"),
  start = c(2006, 2006, 2006, 2006),
  end   = c(2020, 2011, 2014, 2017),
  status = "Pre-Expansion"
)

expanded_data <- data.frame(
  cohort = c("Early Expansion", "On-Time Expansion (2014)", "Late Expansion"),
  start = c(2011, 2014, 2017),
  end   = c(2020, 2020, 2020),
  status = "Post-Expansion"
)

plot_data <- rbind(cohort_data, expanded_data)
plot_data$cohort <- factor(plot_data$cohort,
                           levels = c("Non-Expansion", "Early Expansion", "On-Time Expansion (2014)", "Late Expansion"))

pB <- ggplot(plot_data, aes(xmin = start, xmax = end, ymin = as.numeric(cohort) - 0.4,
                            ymax = as.numeric(cohort) + 0.4, fill = status)) +
  geom_rect(color = "black") +
  geom_vline(xintercept = 2014, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 2011, linetype = "dotted", color = "gray40") +
  geom_vline(xintercept = 2017, linetype = "dotted", color = "gray40") +
  scale_fill_manual(values = c("Pre-Expansion" = "#E31A1C", "Post-Expansion" = "#1F78B4")) +
  scale_y_continuous(breaks = 1:4, labels = levels(plot_data$cohort)) +
  scale_x_continuous(breaks = seq(2006, 2020, 2), limits = c(2006, 2020)) +
  labs(title = "Cohort-Specific Time Periods", x = "Year", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
pB
# ---------- Combine Panels ----------
final_plot <- pA / pB + plot_annotation(tag_levels = "A")
final_plot

# Save
ggsave("nonmet/supplemental_timeline_figure.png", final_plot, width = 9, height = 6, dpi = 300)
ggsave("nonmet/timeline_only_Central.png", pB, width = 7, height = 4, dpi = 300)


#### STATE MAP ####
install.packages("usmap")

install.packages("ggplot2")


# Your Medicaid expansion classification
state_expansion <- data.frame(
  state = c(state.name),  # built-in state names
  status = NA_character_
)

# Define statuses manually based on your DiD design
early_states <- c("California", "Connecticut", "New Jersey", "Washington")
mid_states <- c("Kentucky", "new mexico", "Hawaii", "Iowa")
late_states <- c("Louisiana", "Alaska")
non_expanded <- c("Utah", "Georgia")
not_considered <- setdiff(state.name, c(early_states, mid_states, late_states, non_expanded))

# Assign categories
state_expansion$status[state_expansion$state %in% early_states] <- "Early Expansion (2011)"
state_expansion$status[state_expansion$state %in% mid_states] <- "On-Time Expansion (2014)"
state_expansion$status[state_expansion$state %in% late_states] <- "Late Expansion (2017)"
state_expansion$status[state_expansion$state %in% non_expanded] <- "Non-Expansion"
state_expansion$status[state_expansion$state %in% not_considered] <- "N/A"

# Plot
state_map <- usmap::plot_usmap(data = state_expansion, values = "status") +
  scale_fill_manual(
    values = c(
      "Non-Expansion" = "azure4",          
      "Early Expansion (2011)" = "blue4",    
      "On-Time Expansion (2014)" = "yellow3",      
      "Late Expansion (2017)" = "red3"
    ),
    name = "Medicaid Expansion",
    na.value = "floralwhite"
  ) +
  labs(title = "Medicaid Expansion Status by State") +
  theme(
    legend.position = "left",          # move legend to right side
    legend.direction = "vertical",      # stack items vertically
    legend.title = element_text(size = 14, family = "Times"),
    legend.text  = element_text(size = 12, family = "Times"),
    text = element_text(family = "Times")
  )

print(state_map)

ggsave("nonmet/CRC_StatePlot.pdf", 
       plot = state_map, 
       width = 10, height = 6)



# plot if nonempty
if (nrow(plot_36) > 0) {
  ggplot(plot_36, aes(x = event_bin, y = HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray30") + 
    geom_point(size = 3, color = "navy") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, color = "navy") +
    labs(
      x = "Patients' Year Diagnosed Relative to State's Medicaid Expansion",
      y = "Hazard Ratio for Mortality (All Cause)",
      title = "Event-Study for All-Cause 3-Year Mortality Hazard",
      subtitle = "Reference = year −1 (pre-expansion)"
    ) +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
} else {
  message("No event-time rows matched. Check the exact term labels in tidy_cox_es_36$term.")
}
# optional save
ggsave("binned/cox_36_survivalbins.png", width = 8, height = 6, dpi = 300)


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


