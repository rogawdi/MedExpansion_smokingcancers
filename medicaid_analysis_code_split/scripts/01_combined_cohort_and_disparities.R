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




