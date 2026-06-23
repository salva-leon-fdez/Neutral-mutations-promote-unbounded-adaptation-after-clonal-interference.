# Reproducible Code for: Neutral mutations promote unbounded adaptation after clonal interference.

### Salvador León Fernández, Miguel Ángel Fortuna Alcolado.

# We are going to reproduce the figures and statistical analyses from the Results and Discussion section: **Strictly neutral epistasis unlocks inaccessible adaptive pathways.**

# ---

## Setup & Environment

library("tidyverse")

## Data Loading & Data Preparation

db_local = read_csv(file = "../../DATA/db_local_peaks.csv")
db_local %>% head(3)
db_local %>% tail(3)

# ---

### Fitness peaks analysis

db_high_mut <- db_local %>% filter(mut_rate == 0.01)
db_low_mut  <- db_local %>% filter(mut_rate == 0.001)

low_mut_peaks <- db_low_mut %>%
    summarise(pct_trapped = mean(p_beneficial_mut_previous_non_neutral_sequence == 0) * 100)
low_mut_peaks

high_mut_lineages_strict <- db_high_mut %>%
  summarise(pct_trapped = mean(p_beneficial_mut_previous_non_neutral_sequence == 0) * 100)
high_mut_lineages_strict

# > **We identified local fitness peaks under the high mutation rate (21 %) but not under the low mutation rate**

high_mut_experiments <- db_high_mut %>%
  group_by(rep_id) %>%
  summarise(has_peak = any(p_beneficial_mut_previous_non_neutral_sequence == 0)) %>%
  ungroup() %>%
  summarise(pct_experiments_trapped = mean(has_peak) * 100)
high_mut_experiments

# > **Representing 23% of high-mutation-rate experiments**

high_mut_lineages_relaxed <- db_high_mut %>%
  summarise(pct_trapped = mean(p_beneficial_mut_previous_non_neutral_sequence < 0.001) * 100)
high_mut_lineages_relaxed

# > **Relaxing the peak definition reveals that 63% of lineages are functionally trapped due to extreme scarcity of beneficial neighbors (<0.1%).**

### Epistatic effects on ancestral background

# Classify fitness effects by comparing gestation times (lower time = higher fitness)
db_classified <- db_local %>%
  mutate(
    effect_on_ancestor = case_when(
      gest_time_conditionally_neutral_sequence == 0 | is.na(gest_time_conditionally_neutral_sequence) ~ "Lethal",      # Fails to replicate
      gest_time_conditionally_neutral_sequence < gest_time_previous_non_neutral_sequence ~ "Beneficial",  # Faster replication
      gest_time_conditionally_neutral_sequence == gest_time_previous_non_neutral_sequence ~ "Neutral",     # Equal replication
      gest_time_conditionally_neutral_sequence > gest_time_previous_non_neutral_sequence ~ "Deleterious"  # Slower replication
    )
  )

# Calculate outcome percentages for the high mutation rate (0.01)
epistasis_high <- db_classified %>%
  filter(dplyr::near(mut_rate, 0.01)) %>%
  group_by(effect_on_ancestor) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2))

cat("\n--- High Mutation Rate (0.01) ---\n")
print(epistasis_high)

# Calculate outcome percentages for the low mutation rate (0.001)
epistasis_low <- db_classified %>%
  filter(dplyr::near(mut_rate, 0.001)) %>%
  group_by(effect_on_ancestor) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = round((count / sum(count)) * 100, 2))

cat("\n--- Low Mutation Rate (0.001) ---\n")
print(epistasis_low)

# > **the fitness effect of the final beneficial mutation changed in 27% of cases (becoming lethal in 11% of cases, neutral in 10%, and deleterious in 6%).**

# ---

## Reproducibility Session Info

sessionInfo()

# ---
# ---
