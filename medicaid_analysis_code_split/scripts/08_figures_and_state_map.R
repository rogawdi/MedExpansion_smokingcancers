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


