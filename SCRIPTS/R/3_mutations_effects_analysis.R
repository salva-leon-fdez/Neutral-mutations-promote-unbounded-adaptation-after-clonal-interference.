# Reproducible Code for: Neutral mutations promote unbounded adaptation after clonal interference.

### Salvador León Fernández, Miguel Ángel Fortuna Alcolado.

# We are going to reproduce the figures and statistical analyses from the Results and Discussion section: **Beneficial exhaustion drives a shift toward neutral evolution**.

# ---

## Setup & Environment

library("tidyverse")
library("scales") 
library("patchwork")
library("data.table")

## Data Loading

#### Read mutations fixed by selection:

db_fixed <- read_csv("../../DATA/db_fixed_mutations.csv") %>% mutate(rep_id_paired = ifelse(mut_rate == 0.01, rep_id - 200, rep_id))
db_fixed %>% head(3)
db_fixed %>% tail(3)

#### Read lineages fixed by selection from files:

db_lineages <- read_csv("../../DATA/db_lineages.csv") %>% mutate(rep_id = ifelse(mut_rate == 0.01, rep_id - 200, rep_id))
    db_lineages %>% nrow
    db_lineages %>% head(3)
    db_lineages %>% tail(3)

db_lineages_unique <- db_lineages %>%
    group_by(rep_id, mut_rate, lineage_id) %>%
    filter(update != max(update)) %>%
    ungroup() %>%
    group_by(rep_id, mut_rate, sequence) %>%
    slice_min(update, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(rep_id, mut_rate, update) %>%
    select(-lineage_id, -ancestor_dist)
db_lineages_unique %>% head(5)

#### Read muutants alive every 1000 updates:

db_load_lineages <- fread("../../DATA/db_mut_load.csv", 
                 colClasses = c("integer", "integer", "character", "numeric", 
                                "integer", "integer", "character", "integer")) %>%
    filter(mu > 0.0005) %>%
    mutate(rep_id = ifelse(rep_id > 200, rep_id - 200, rep_id)) %>%
    arrange(rep_id, update)
    
db_load_lineages %>% head
db_load_lineages %>% tail

#### Read mutational supply data

# ==============================================================================
# DATA GENERATION PIPELINE (For Reproducibility Only)
# ==============================================================================
# The following block of code was used to process the massive raw dataset 
# ("DATA/db_mut_supply.csv") and generate the summary file used in this script. 
# It is commented out to save execution time and memory during standard runs, 
# as the final processed file is already provided in the repository.
#
# db_available <- fread("../../DATA/db_mut_supply.csv", colClasses = cls) %>%
#  filter(mut_rate > 0.0005) %>%
#  mutate(relative_fitness = ifelse(update == 0, NA, relative_fitness)) %>%
#  arrange(rep_id, update)
#
# db_available_summary <- db_available %>%
#  filter(
#    !is.na(relative_fitness), 
#    !is.na(mut_rate),
#    resource == "NOT",
#    update != 0                
#  ) %>%
#  mutate(
#    update_bin = ceiling(update / 1000) * 1000, 
#    effect_type = case_when(
#      relative_fitness > 1  ~ "Beneficial",
#      relative_fitness == 1 ~ "Neutral",
#      relative_fitness < 1 & relative_fitness > 0  ~ "Deleterious",
#      relative_fitness == 0  ~ "Lethal"
#    ),
#    effect_type = factor(effect_type, levels = c("Lethal", "Beneficial", "Neutral", "Deleterious")),
#    mut_rate    = factor(mut_rate)
#  ) %>%
#  count(mut_rate, update_bin, effect_type) %>%
#  group_by(mut_rate, update_bin) %>%
#  mutate(percentage = (n / sum(n)) * 100) %>%
#  ungroup()
#
# write_csv(db_available_summary, file = "../../DATA/db_mut_supply_summary.csv")
# ==============================================================================

db_available_summary = read_csv(file = "../../DATA/db_mut_supply_summary.csv")
db_available_summary %>% head
db_available_summary %>% tail

# ---

## Statistical Analysis Function:

perform_paired_analysis <- function(data, rep_col, cond_col, val_col, group1, group2) {
  
  # Pivot data to align paired observations
  df_paired <- data %>%
    select(all_of(c(rep_col, cond_col, val_col))) %>%
    pivot_wider(names_from = all_of(cond_col), values_from = all_of(val_col)) %>%
    drop_na() # Drop incomplete pairs
  
  vec_g1 <- df_paired[[group1]]
  vec_g2 <- df_paired[[group2]]
  
  differences <- vec_g1 - vec_g2
  
  # Check normality of differences
  shapiro_res <- shapiro.test(differences)
  
  cat("\n--------------------------------------------------\n")
  cat("1. NORMALITY TEST RESULTS (Shapiro-Wilk)\n")
  cat("--------------------------------------------------\n")
  cat("P-value:", round(shapiro_res$p.value, 5), "\n\n")
  
  # Select and run the appropriate test
  if (shapiro_res$p.value > 0.05) {
    
    cat("=> Differences FOLLOW a normal distribution (p > 0.05).\n")
    cat("=> Running Paired T-Test (Parametric)...\n\n")
    
    stat_test_res <- t.test(vec_g1, vec_g2, paired = TRUE)
    
  } else {
    
    cat("=> Differences DO NOT follow a normal distribution (p < 0.05).\n")
    cat("=> Running Paired Wilcoxon Signed-Rank Test (Non-parametric)...\n\n")
    
    # Handle zero-difference ties smoothly
    stat_test_res <- wilcox.test(vec_g1, vec_g2, paired = TRUE, exact = FALSE) 
  }
  
  cat("--------------------------------------------------\n")
  cat("2. STATISTICAL TEST RESULTS\n")
  cat("--------------------------------------------------\n")
  print(stat_test_res)
  
  # Plot the distribution of differences
  p_diff <- ggplot(data.frame(diff = differences), aes(x = diff)) +
    geom_density(fill = "gray80", alpha = 0.5, color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "Distribution of Paired Differences",
         subtitle = paste("Selected test:", stat_test_res$method),
         x = paste("Difference (", group1, "-", group2, ")"),
         y = "Density") +
    theme_minimal()
  
  print(p_diff)
  
  # Return objects silently for optional assignment
  invisible(list(
    shapiro = shapiro_res,
    stat_test = stat_test_res,
    plot = p_diff
  ))
}

# ---

## Statistical Analysis & Figure Generation

# ==============================================================================
# FIGURE EXPORT TEMPLATE (Optional)
# ==============================================================================
# If you wish to extract and save any of the generated plots in high resolution,
# you can uncomment and use the following base template. 
# Just replace 'your_plot_variable' with the actual name of the plot object 
# you want to save (e.g., relative_fitness_ci, final_plot, etc.).

# dir.create("Figures", showWarnings = FALSE)

# ggsave(
#   filename = "Figures/your_figure_name.png", 
#   plot     = your_plot_variable, 
#   width    = 5, 
#   height   = 6, 
#   units    = "in", 
#   dpi      = 600, 
#   bg       = "white", 
#   scale    = 1.2
# )

### Figure 3A. Effects of fixed mutations over time

#### Relative frequency of fixed mutations over time

db_fixed_per100_NOT <- db_fixed %>%
  filter(
    !is.na(relative_fitness),
    mut_rate > 0.0005,
    !is.na(mut_rate),
    resource == "NOT",
    update != 0               
  ) %>%
  mutate(
    # Binning: Rounds updates up to the nearest 1000-step interval
    update_bin = ceiling(update / 1000) * 1000, 
    
    # Classification of mutational effects based on relative fitness
    effect_type = case_when(
      relative_fitness > 1  ~ "Beneficial",
      relative_fitness == 1 ~ "Neutral",
      relative_fitness < 1  ~ "Deleterious"
    ),
    
    # Define factor levels for consistent plotting order
    effect_type = factor(effect_type, levels = c("Beneficial", "Neutral", "Deleterious")),
    mut_rate    = factor(mut_rate)
  ) %>%
  
  # Calculate frequency and relative proportions (percentage)
  count(mut_rate, update_bin, effect_type) %>%
  group_by(mut_rate, update_bin) %>%
  mutate(percentage = (n / sum(n)) * 100) %>%
  ungroup()

db_fixed_per100_NOT %>% head()

options(repr.plot.width = 12, repr.plot.height = 6)

global_distribution = ggplot(db_fixed_per100_NOT %>%
       mutate(
            mut_rate = factor(mut_rate, 
                              levels = c(0.001, 0.01), 
                              labels = c("Low Mutation Rate", "High Mutation Rate"))), 
       aes(x = update_bin, y = percentage, color = effect_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5, alpha = 0.7) +
  
  # Stratify by mutation rate
  facet_wrap(~ mut_rate) +
  scale_x_continuous(labels = function(x) x / 100000) +
                     
  # Color palette for fitness effects
  scale_color_manual(values = c(
    "Beneficial"  = "#2ca02c",  # Green
    "Neutral"     = "#7f7f7f",  # Gray
    "Deleterious" = "#d62728"   # Red
  )) +
  
labs(
    x = expression(paste("Time (Updates [", "" %*% 10^5, "])")),    
    y     = "Percentage of Mutations (%)",
    color = "Fitness effect",
    fill  = "Fitness effect" 
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
  theme_bw() +
  theme(
      strip.text = element_text(size = 18), 
      strip.background = element_blank(),
      
      axis.title.x = element_text(size = 18, margin = margin(t = 5), color = "black"),
      axis.title.y = element_text(size = 20, margin = margin(r = 15), color = "black"),
      
      axis.text.x = element_text(color = "black", size = 16),
      axis.text.y = element_text(color = "black", size = 16),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid.major = element_line(color = "gray90", linetype = "solid"),
      panel.grid.minor = element_blank(),
      
      legend.position = "top",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
  )
global_distribution

### Statistical analysis Figure 3A: Beneficial mutations dominate the entire adaptive process at a low mutation rate

db_fixed_per100_NOT %>% head(3)
db_fixed_per100_NOT %>% tail(3)

db_fixed_per100_NOT <- db_fixed_per100_NOT %>% filter(mut_rate == "0.001") %>% complete(update_bin, effect_type, fill = list(percentage = 0))

results_dominance <- perform_paired_analysis(
  data     = db_fixed_per100_NOT,
  rep_col  = "update_bin", 
  cond_col = "effect_type", 
  val_col  = "percentage",   
  group1   = "Beneficial",
  group2   = "Neutral"
)

### Figure 3B-E. Combined plot: results at the end of adaptation

#### Mutation effects dristribution at the end of experiment:

db_rel_fitness <- db_lineages_unique %>%
    filter(!is.na(mut_type)) %>%
    group_by(mut_rate, rep_id) %>%
    summarise(
        n_total_mut = n(),
        n_mut_beneficial = sum(relative_fitness > 1),
        n_mut_detrimental = sum(relative_fitness < 1),
        n_mut_neutral = sum(relative_fitness == 1), .groups = "drop")
db_rel_fitness %>% nrow
db_rel_fitness %>% head(3)
db_rel_fitness %>% tail(3)

db_rel_fitness_long <- db_rel_fitness %>%
    pivot_longer(
        cols = c(n_mut_beneficial, n_mut_detrimental, n_mut_neutral),
        names_to = "effect_mut",
        names_prefix = "n_mut_",
        values_to = "n_mut") %>%
    mutate(p_effect_mut = n_mut/n_total_mut)
db_rel_fitness_long %>% nrow
db_rel_fitness_long %>% head(3)
db_rel_fitness_long %>% tail(3)

#### Pairwise genome distance at the end of adaptation:

LEVENSHTEIN_BIN <- "levenshtein/./levenshtein"
 
run_levenshtein <- function(sequences) {
    seqs <- unique(sequences)
    if (length(seqs) < 2) {
        return(tibble(n=length(seqs), pairs=0L, mean=NA_real_, variance=NA_real_))
    }
    tmp <- tempfile()
    write.table(seqs, tmp, quote=FALSE, row.names=FALSE, col.names=FALSE)
    out <- system(paste(LEVENSHTEIN_BIN, tmp, "2>&1"), intern=TRUE)
    file.remove(tmp)
    read.csv(text=tail(out, 1), header=FALSE,
             col.names=c("n", "pairs", "mean", "variance"))
}

# ==============================================================================
# 1. GET SEQUENCES AT THE PENULTIMATE update_sample PER REPLICATE
# ==============================================================================
seqs_final <- db_load_lineages %>%
    group_by(rep_id, mu) %>%
    mutate(penultimate = sort(unique(update), decreasing=TRUE)[2]) %>%
    filter(update == penultimate) %>%
    ungroup()
 
cat("Penultimate update_sample per group:\n")
seqs_final %>%
    distinct(rep_id, mu, update) %>%
    group_by(mu) %>%
    summarise(update = first(update), n_reps = n(), .groups="drop") %>%
    print()
 
cat("\nSequences per replicate:\n")
seqs_final %>%
    group_by(mu, rep_id) %>%
    summarise(n_seqs = n(), .groups="drop") %>%
    group_by(mu) %>%
    summarise(median_n = median(n_seqs),
              min_n    = min(n_seqs),
              max_n    = max(n_seqs), .groups="drop") %>%
    print()

# ==============================================================================
# 2. LEVENSHTEIN WITHIN (rep_id x mut_rate)
# ==============================================================================
lev_final <- seqs_final %>%
    group_by(rep_id, mu) %>%
    group_modify(~ run_levenshtein(.x$sequence)) %>%
    ungroup() %>%
    filter(!is.na(mean)) 
 
cat("\nValid replicates:\n")
lev_final %>% count(mu) %>% print()
 
cat("\nSummary:\n")
lev_final %>%
    group_by(mu) %>%
    summarise(n=n(), median=median(mean), sd=sd(mean), .groups="drop") %>%
    print()

### Combined plot. Figure 3B-E

# Define shared color scale to ensure consistency across all panels
shared_fill <- scale_fill_manual(
  name   = "Mutation Rate",
  values = c("0.001" = "#F8766D", "0.01" = "#00BFC4"),
  labels = c("0.001" = "0.001", "0.01" = "0.01")
)
 
# Main panel: Proportion of fixed mutations by effect type
# (This plot holds the master legend for the final combined figure)
p_rel_fitness_NOT <- ggplot(
    db_rel_fitness_long %>%
        mutate(
          effect_mut = factor(effect_mut, levels = c("beneficial", "neutral", "detrimental")),
          mut_rate   = factor(mut_rate)
        ),
    aes(x = effect_mut, y = p_effect_mut, fill = mut_rate)
  ) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  shared_fill +
  scale_y_continuous(labels = percent) +
  scale_x_discrete(
    labels = c("beneficial" = "+", "neutral" = "Neutral", "detrimental" = "-")
  ) +
  labs(x = NULL, y = "Fixed Mutations", fill = "Mutation Rate") +
  theme_bw() +
  theme(
    plot.title        = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 22, margin = margin(r = 10)),
    axis.text.x       = element_text(color = "black", size = 20),
    axis.text.y       = element_text(color = "black", size = 20),
    axis.ticks        = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major  = element_line(color = "gray90", linetype = "solid"),
    legend.position   = "top",
    legend.title      = element_text(size = 22),
    legend.text       = element_text(size = 20),
    legend.key.size   = unit(1.0, "cm")
  )
 
# Subpanel: Total count of fixed mutations prior to the final update
df_fixed_end <- db_lineages_unique %>%
  filter(mut_rate %in% c(0.001, 0.01)) %>%
  filter(update < 499000) %>%
  group_by(mut_rate, rep_id) %>%
  summarise(total_fixed = n(), .groups = "drop") %>%
  mutate(mut_rate = factor(mut_rate))
 
p_fixed_count <- ggplot(
    df_fixed_end,
    aes(x = mut_rate, y = total_fixed, fill = mut_rate)
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, outlier.size = 1.5, width = 0.5) +
  shared_fill +
  scale_x_discrete(labels = c("0.001" = "Low", "0.01" = "High")) +
  labs(x = NULL, y = "Fixed Mutations (Count)") +
  theme_bw() +
  theme(
    plot.title        = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 18, margin = margin(r = 5)),
    axis.text.x       = element_text(color = "black", size = 20),
    axis.text.y       = element_text(color = "black", size = 18),
    axis.ticks        = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major  = element_line(color = "gray90", linetype = "solid"),
    legend.position   = "none"
  )
 
# Subpanel: Number of surviving genotypes at the end of the experiment
db_alive <- db_load_lineages %>%
  filter(update == 499000) %>%
  group_by(mu, rep_id) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(mu = factor(mu))
 
p_living_genotypes <- ggplot(
    db_alive,
    aes(x = mu, y = n, fill = mu)
  ) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, outlier.size = 1.5, width = 0.5) +
  shared_fill +
  scale_x_discrete(labels = c("0.001" = "Low", "0.01" = "High")) +
  labs(x = NULL, y = "Number of Surviving Genotypes") +
  theme_bw() +
  theme(
    plot.title        = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 18, margin = margin(r = 5)),
    axis.text.x       = element_text(color = "black", size = 20),
    axis.text.y       = element_text(color = "black", size = 18),
    axis.ticks        = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major  = element_line(color = "gray90", linetype = "solid"),
    legend.position   = "none"
  )
 
# Subpanel: Mean pairwise genotypic distance at the final update
p_lev_final <- ggplot(
    lev_final %>% mutate(mu = factor(mu)),
    aes(x = mu, y = mean, fill = mu)
  ) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  shared_fill +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0.0, 0.15)) 
  ) +
  coord_cartesian(ylim = c(0,0.3)) +
  scale_x_discrete(labels = c("0.001" = "Low", "0.01" = "High")) +
  labs(x = NULL, y = "Mean Pairwise Genotype Distance") +
  theme_bw() +
  theme(
    plot.title        = element_blank(),
    axis.title.x      = element_blank(),
    axis.title.y      = element_text(size = 20, margin = margin(r = 10)),
    axis.text.x       = element_text(color = "black", size = 18),
    axis.text.y       = element_text(color = "black", size = 18),
    axis.ticks        = element_line(linewidth = 0.5),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major  = element_line(color = "gray90", linetype = "solid"),
    legend.position   = "none"
  )
 
# Assemble final figure layout
options(repr.plot.width = 12, repr.plot.height = 6)
 
combined_plot <- p_rel_fitness_NOT + p_fixed_count + p_living_genotypes + p_lev_final +
  plot_layout(widths = c(3, 1, 1, 1)) 
 
combined_plot

### Statistical analysis. Figure 3B-E

#### Figure 3B. Mutation effects dristribution at the end of experiment

# DESCRIPTIVE STATISTICS
stats_panelB <- db_rel_fitness_long %>%
  group_by(mut_rate, effect_mut) %>%
  summarise(
    median_val = median(p_effect_mut, na.rm = TRUE),
    q1_val     = quantile(p_effect_mut, 0.25, na.rm = TRUE),
    q3_val     = quantile(p_effect_mut, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

print("--- PANEL B: FIXED MUTATIONS BY EFFECT TYPE ---")
stats_panelB


db_rel_fitness_long %>% head(3)

interaction_model <- lm(p_effect_mut ~ as.factor(mut_rate) * as.factor(effect_mut), data = db_rel_fitness_long)

summary(interaction_model)

anova(interaction_model)

#### Figure 3C. Length of the adaptive walk

stats_panelC <- df_fixed_end %>%
  group_by(mut_rate) %>%
  summarise(
    median_val = median(total_fixed, na.rm = TRUE),
    q1_val     = quantile(total_fixed, 0.25, na.rm = TRUE),
    q3_val     = quantile(total_fixed, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

print("--- PANEL C: TOTAL FIXED MUTATIONS ---")
stats_panelC


df_fixed_end %>% head(3)

results_panelC <- perform_paired_analysis(
  data     = df_fixed_end,     
  rep_col  = "rep_id",         
  cond_col = "mut_rate",       
  val_col  = "total_fixed",   
  group1   = "0.001",
  group2   = "0.01"
)

#### Figure 3D. Genotypes alive at the end of experiments

stats_panelD <- db_alive %>%
  group_by(mu) %>%
  summarise(
    median_val = median(n, na.rm = TRUE),
    q1_val     = quantile(n, 0.25, na.rm = TRUE),
    q3_val     = quantile(n, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

print("--- PANEL D: SURVIVING GENOTYPES ---")
stats_panelD

db_alive %>% head(3)

results_panelD <- perform_paired_analysis(
  data     = db_alive,         
  rep_col  = "rep_id",         
  cond_col = "mu",            
  val_col  = "n",              
  group1   = "0.001",
  group2   = "0.01"
)

#### Figure 3E. Mean pairwise Levenshtein distance 

stats_panelE <- lev_final %>%
  group_by(mu) %>%
  summarise(
    median_val = median(mean, na.rm = TRUE),
    q1_val     = quantile(mean, 0.25, na.rm = TRUE),
    q3_val     = quantile(mean, 0.75, na.rm = TRUE),
    .groups    = "drop"
  )

print("--- PANEL E: LEVENSHTEIN DISTANCE ---")
stats_panelE

lev_final %>% head(3)

results_panelE <- perform_paired_analysis(
  data     = lev_final,        
  rep_col  = "rep_id",        
  cond_col = "mu",             
  val_col  = "mean",           
  group1   = "0.001",
  group2   = "0.01"
)

# ---

## Reproducibility Session Info

sessionInfo()

# ---
# ---
