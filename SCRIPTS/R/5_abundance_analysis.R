# Reproducible Code for: Neutral mutations promote unbounded adaptation after clonal interference.

### Salvador León Fernández, Miguel Ángel Fortuna Alcolado.

# We are going to reproduce the figures and statistical analyses from the Results and Discussion section: **Clonal monopolization restricts the accumulation of neutral mutations.**

# ---

## Setup & Environment

library("tidyverse")
library("scales") 
library("patchwork")
library("data.table")

## Data Loading & Data Preparation

db_load <- fread("../../DATA/db_mut_load.csv", 
                 colClasses = c("integer", "integer", "character", "numeric", 
                                "integer", "integer", "character", "integer")) %>%
  arrange(desc(resource), rep_id, update)
db_load %>% head(3)
db_load %>% tail(3)

db_load_lineages = read_csv("../../DATA/db_joined_mut_load_lineages.csv",
col_types = cols("i", "i", "d", "c", "d", "i", "l", "d", "d", "i", "i"))
db_load_lineages %>% head(3)
db_load_lineages %>% tail(3)

## Mapping seq_ids to genomes

seq_dict = db_load %>%
  distinct(founder_id, mu, rep_id, sequence) %>%
  group_by(founder_id, mu, rep_id) %>%
  mutate(seq_id = row_number()) %>%
  ungroup()

db_load_id = db_load %>%
  left_join(seq_dict, by = c("founder_id", "mu", "rep_id", "sequence"))

db_load_id %>% head(3)

## Extracting the most abundant organism at each update

max_abundance_update = db_load_id %>% 
  filter(abundance == max(abundance), .by = c(founder_id, mu, rep_id, update))

## Extracting main lineage distribution

# The extracted organisms represent the most abundant genotypes of their respective lineages at each sampled update.

db_main = db_load_lineages %>% 
  group_by(founder_id, rep_id, lineage_id, update_sample) %>% 
  slice_max(order_by = abundance, n = 1, with_ties = FALSE) %>% arrange(founder_id, rep_id, update_sample)
db_main %>% head(3)

# --- 

## Functions:

#### Abundance of all organisms and main lineage distribution. SUPPLEMENTARY FIGURE 1

abundance_time_course_effects <- function(founder_id_number, mu_val, rep_id_number = NULL) {
  
  # Filter and prepare background abundance data for all sequences
  bg_data <- db_load_id %>% 
    filter(founder_id == founder_id_number, dplyr::near(mu, mu_val))
    
  if(!is.null(rep_id_number)) {
    bg_data <- bg_data %>% filter(rep_id == rep_id_number)
  }
  
  if(nrow(bg_data) == 0) {
    stop("No data found for the specified parameters in db_load_id.")
  }

  # Shuffle factor levels to assign highly distinct colors to adjacent sequences
  bg_data <- bg_data %>%
  mutate(seq_id_factor = factor(seq_id, levels = sample(unique(seq_id))))
    
  num_sequences <- length(unique(bg_data$seq_id))
  distinct_palette <- colorRampPalette(palette.colors(n = 36, palette = "Polychrome 36"))(num_sequences)    
    
  # Filter and prepare trajectory data for the main evolutionary lineage
  main_data <- db_main %>% 
    filter(founder_id == founder_id_number, 
           dplyr::near(mut_rate, mu_val),
           is_main_lineage == TRUE)
    
  if(!is.null(rep_id_number)) {
    main_data <- main_data %>% filter(rep_id == rep_id_number)
  }

  # Construct the visualization overlaying main lineage and mutational events
  p <- ggplot() +
    
    # Background: Individual sequence trajectories
    geom_line(data = bg_data, 
              aes(x = update, y = abundance, color = seq_id_factor), 
              alpha = 0.7, linewidth = 0.5) +
              
    # Foreground: Dominant lineage trajectory
    geom_line(data = main_data, 
              aes(x = update_sample, y = abundance, group = lineage_id), 
              color = "black", linewidth = 1.5) +
              
    # Highlight: First appearance of deleterious mutations (red)
    geom_point(data = main_data %>% 
                 filter(relative_fitness < 1, relative_fitness > 0) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#F8766D", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of neutral mutations (blue)
    geom_point(data = main_data %>% 
                 filter(relative_fitness == 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#619CFF", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of beneficial mutations (green)
    geom_point(data = main_data %>% 
                 filter(relative_fitness > 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#00BA38", color = "black", size = 3.5, stroke = 1) +
    
     scale_color_manual(values = distinct_palette) +

    xlim(0, 500000) +  
    ylim(0, 10000) +
    
    labs(
      title = paste0("Evolutionary Dynamics (Founder: ", founder_id_number, " | mu: ", mu_val,
                     ifelse(is.null(rep_id_number), " | All Replicates", paste0(" | Replicate: ", rep_id_number)), ")"),
      subtitle = "Colors: All seq_id | Black: Main Lineage\nDots: Green = Beneficial (>1) | Blue = Neutral (==1) | Red = Deleterious (<1)",
      x = "Update",
      y = "Abundance"
    ) +
    theme_minimal() +
    theme(legend.position = "none") 
  
  # Facet by replicate if querying the entire population
  if(is.null(rep_id_number)) {
    p <- p + facet_wrap(~rep_id) 
  }
  
  return(p)
}

#### Abundance of all organisms and main lineage distribution. FIGURE 4A-B

abundance_main_image <- function(founder_id_number, mu_val, rep_id_number = NULL) {
  
  # Filter and prepare background abundance data for all sequences
  bg_data <- db_load_id %>% 
    filter(founder_id == founder_id_number, dplyr::near(mu, mu_val))
    
  if(!is.null(rep_id_number)) {
    bg_data <- bg_data %>% filter(rep_id == rep_id_number)
  }
  
  if(nrow(bg_data) == 0) {
    stop("No data found for the specified parameters in db_load_id.")
  }

  # Shuffle factor levels to assign highly distinct colors to adjacent sequences
  bg_data <- bg_data %>%
  mutate(seq_id_factor = factor(seq_id, levels = sample(unique(seq_id))))
    
  num_sequences <- length(unique(bg_data$seq_id))
  distinct_palette <- colorRampPalette(palette.colors(n = 36, palette = "Polychrome 36"))(num_sequences)    
    
  # Filter and prepare trajectory data for the main evolutionary lineage
  main_data <- db_main %>% 
    filter(founder_id == founder_id_number, 
           dplyr::near(mut_rate, mu_val),
           is_main_lineage == TRUE)
    
  if(!is.null(rep_id_number)) {
    main_data <- main_data %>% filter(rep_id == rep_id_number)
  }

  # Construct the main visualization
  p <- ggplot() +
    
    # Background: Individual sequence trajectories
    geom_line(data = bg_data, 
              aes(x = update, y = abundance, color = seq_id_factor), 
              alpha = 0.7, linewidth = 0.5) +
              
    # Foreground: Dominant lineage trajectory
    geom_line(data = main_data, 
              aes(x = update_sample, y = abundance, group = lineage_id), 
              color = "black", linewidth = 1.5) +
              
    # Highlight: First appearance of deleterious mutations (red)
    geom_point(data = main_data %>% 
                 filter(relative_fitness < 1, relative_fitness > 0) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#F8766D", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of neutral mutations (blue)
    geom_point(data = main_data %>% 
                 filter(relative_fitness == 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#619CFF", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of beneficial mutations (green)
    geom_point(data = main_data %>% 
                 filter(relative_fitness > 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#00BA38", color = "black", size = 3.5, stroke = 1) +
    
     scale_color_manual(values = distinct_palette) +

     scale_x_continuous(labels = function(x) x / 100000, expand = expansion(mult = c(0, 0.05))) +
     scale_y_continuous(limits = c(0, 10000), expand = expansion(mult = c(0, 0.05))) +    
    labs(
      x = expression(paste("Time (Updates [", "" %*% 10^5, "])")),    
      y = "Abundance"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14),
      
      strip.background = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # Facet by replicate if querying the entire population
  if(is.null(rep_id_number)) {
    p <- p + facet_wrap(~rep_id) 
  }
  
  return(p)
}

#### Visualizing lineage competition at the end of the adaptive period

abundance_time_course_lineages <- function(founder_id_number, mu_val, rep_id_number = NULL, min_update = 0, max_update = 500000) {
  
  # Load and filter primary lineage dataset
  plot_data <- db_main %>% 
    filter(founder_id == founder_id_number, dplyr::near(mut_rate, mu_val))
    
  if(!is.null(rep_id_number)) {
    plot_data <- plot_data %>% filter(rep_id == rep_id_number)
  }
  
  if(nrow(plot_data) == 0) {
    stop("No data found for the specified parameters.")
  }
  
  # Calculate dynamic y-axis maximum within the specified zoom window
  max_y_zoom <- plot_data %>%
    filter(update_sample >= min_update, update_sample <= max_update) %>%
    pull(abundance) %>%
    max(na.rm = TRUE)
    
  if(is.infinite(max_y_zoom) || is.na(max_y_zoom)) max_y_zoom <- 10000 
  
  # Segregate highly successful lineages (abundance > 500) from background noise
  successful_ids <- plot_data %>%
    group_by(lineage_id) %>%
    filter(max(abundance, na.rm = TRUE) > 500) %>% 
    pull(lineage_id) %>%
    unique()
    
  bg_data <- plot_data %>% filter(!lineage_id %in% successful_ids)
  color_data <- plot_data %>% filter(lineage_id %in% successful_ids)
  
  # Construct the lineage competence visualization
  p <- ggplot() +
    
    # Background: Low-abundance or unsuccessful lineages
    geom_line(data = bg_data, 
              aes(x = update_sample, y = abundance, group = lineage_id), 
              color = "gray80", alpha = 0.4, linewidth = 0.4) +
    
    # Foreground: Trajectories of highly successful lineages
    geom_line(data = color_data, 
              aes(x = update_sample, y = abundance, color = as.factor(lineage_id), group = lineage_id), 
              linewidth = 1.2, alpha = 0.9) +
              
    geom_point(data = color_data, 
               aes(x = update_sample, y = abundance, color = as.factor(lineage_id)), 
               size = 1.5, alpha = 0.7) +
               
    # Highlight: First appearance of deleterious mutations (red)
    geom_point(data = color_data %>% 
                 filter(relative_fitness < 1, relative_fitness > 0) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#F8766D", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of neutral mutations (blue)
    geom_point(data = color_data %>% 
                 filter(relative_fitness == 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#619CFF", color = "black", size = 2.5, stroke = 1) +
               
    # Highlight: First appearance of beneficial mutations (green)
    geom_point(data = color_data %>% 
                 filter(relative_fitness > 1) %>% 
                 group_by(lineage_id, update_born) %>% 
                 slice_min(order_by = update_sample, n = 1, with_ties = FALSE), 
               aes(x = update_sample, y = abundance), 
               shape = 21, fill = "#00BA38", color = "black", size = 3.5, stroke = 1) +

    # Apply Cartesian zoom to preserve underlying data structure for faceting
    coord_cartesian(xlim = c(min_update, max_update), ylim = c(0, max_y_zoom * 1.05)) +
    
    scale_x_continuous(labels = function(x) x / 100000, expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(limits = c(0, 10000), expand = expansion(mult = c(0, 0.05))) +    
    
    labs(
      x = expression(paste("Time (Updates [", "" %*% 10^5, "])")),    
      y = "Abundance"
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      
      axis.text = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size=14, margin=margin(t=1)),
      axis.title.y = element_text(size=14, margin=margin(r=10)),    
      
      strip.background = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # Facet by replicate with free y-axes to accommodate scale variance
  if(is.null(rep_id_number)) {
    p <- p + facet_wrap(~rep_id, scales = "free_y") 
  }
  
  return(p)
}

## Statistical Analysis & Figure Generation

# ==============================================================================
# FIGURE EXPORT TEMPLATE (Optional)
# ==============================================================================
# Uncomment the lines below to save a high-resolution copy of the figure.
# NOTE: By default, ggsave() captures the last plot displayed. Ensure you run 
# this command immediately after generating/viewing the desired plot.

# dir.create("Figures", showWarnings = FALSE)
# ggsave("Figures/your_figure_name.png", dpi = 600, width = X, height = Y)

### Figure 4A

options(repr.plot.width = 8, repr.plot.height = 5)

abundance_main_image(1271, 0.001, 6)

### Figure 4B

options(repr.plot.width = 8, repr.plot.height = 5)

abundance_main_image(1271, 0.01, 206) # Note: rep_ids > 200 are paired counterparts (e.g., 206 corresponds to rep_id 6)

### Figure 4C

plot_data <- db_load_id %>%
  
  filter(mu > 0.0005) %>% 
  group_by(mu, founder_id, rep_id, update) %>%
  summarise(
    n_seqs = n_distinct(seq_id), 
    .groups = "drop"
  ) %>%
  
  group_by(mu, update) %>%
  summarise(
    global_mean = mean(n_seqs),
    global_sd = sd(n_seqs),
    
    sd_lower = global_mean - global_sd,
    sd_upper = global_mean + global_sd,
    
    .groups = "drop"
  ) %>%
  
  mutate(mu = as.factor(mu))

options(repr.plot.width = 8, repr.plot.height = 5)

genome_diversification = ggplot(plot_data, aes(x = update, y = global_mean, color = mu, fill = mu, group = mu)) +
  
  geom_ribbon(aes(ymin = sd_lower, ymax = sd_upper), alpha = 0.4, color = NA) +
  
  geom_line(linewidth = 1.2) +
  
  scale_x_continuous(labels = function(x) x / 100000) +
  
  labs(
    x = expression(paste("Time (Updates [", "" %*% 10^5, "])")),    
    y = "Mean Count of Unique Genotypes",
    color = "Mutation Rate:",
    fill = "Mutation Rate:" 
  ) +
  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20, margin = margin(t = 15)),
    axis.title.y = element_text(size = 18, margin = margin(r = 12)),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(size = 16, color = "black"),    
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray90", linetype = "solid"),
    
    legend.position = "top",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.7, "cm"),
    legend.margin = margin(b = 10) 
  )
genome_diversification

### Figure 4D

options(repr.plot.width = 8, repr.plot.height = 5)

# Isolate and categorize first-appearance mutation events across all replicates
mutation_events_global <- db_main %>% 
  filter(is_main_lineage == TRUE) %>%
  mutate(mut_type = case_when(
    relative_fitness < 1 ~ "Deleterious",
    relative_fitness == 1 ~ "Neutral",
    relative_fitness > 1 ~ "Beneficial"
  )) %>%
  # Group by lineage parameters to ensure mutational events are counted uniquely
  group_by(founder_id, mut_rate, rep_id, lineage_id, update_born) %>% 
  # Extract the exact update when the mutation is first sampled
  slice_min(order_by = update_sample, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Exclude deleterious mutations to focus exclusively on neutral and beneficial dynamics
  filter(mut_type != "Deleterious")

# Construct overlapping histogram faceted by mutation rate
ggplot(mutation_events_global, aes(x = abundance, fill = mut_type)) +
  
  geom_histogram(position = "identity", alpha = 0.6, bins = 50, color = "black", linewidth = 0.3) +
  
  scale_fill_manual(values = c("Neutral" = "#619CFF", "Beneficial" = "#00BA38")) +
  
  scale_x_continuous(labels = function(x) x / 1000, expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  
  facet_wrap(~ mut_rate, scales = "free_y", 
             labeller = labeller(mut_rate = c("0.001" = "Low Mutation Rate", "0.01" = "High Mutation Rate"))) +  
  
  labs(
    x = expression(paste("Lineage Abundance [", "" %*% 10^3, "]")),    
    y = "Frequency (Count)",
    fill = "Mutation Effect:"
  ) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 16, color = "black"),
    
    axis.text.x = element_text(color = "black", size = 15.5), 
    axis.text.y = element_text(color = "black", size = 15.5),
    
    axis.ticks.x = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.y = element_line(linewidth = 0.8, color = "black"),
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "gray90", linetype = "solid"),
    panel.grid.minor = element_blank(),
    
    strip.background = element_blank(), 
    strip.text = element_text(size = 14, color = "black", face = "bold"),
    
    legend.position = "top",
    legend.title = element_text(size = 17), 
    legend.text = element_text(size = 15), 
    legend.key.size = unit(0.55, "cm") 
  )

### Statistical analysis. Figure 4D: Neutral mutations can only accumulate when there is no genome monopolization

#### Unpaired analysis function

perform_unpaired_analysis <- function(data, group_col, val_col, group1, group2) {
  
  # Isolate complete cases and extract numerical vectors per group
  df_clean <- data %>% drop_na(all_of(c(group_col, val_col)))
  
  vec_g1 <- df_clean %>% filter(!!sym(group_col) == group1) %>% pull(!!sym(val_col))
  vec_g2 <- df_clean %>% filter(!!sym(group_col) == group2) %>% pull(!!sym(val_col))
  
  # Validate presence of observations in both comparison groups
  if(length(vec_g1) == 0 | length(vec_g2) == 0) {
    stop("One of the groups has 0 observations. Check group names.")
  }
  
  cat("\n--------------------------------------------------\n")
  cat("1. DESCRIPTIVE STATISTICS\n")
  cat("--------------------------------------------------\n")
  cat(group1, "-> N:", length(vec_g1), "| Median:", median(vec_g1), "\n")
  cat(group2, "-> N:", length(vec_g2), "| Median:", median(vec_g2), "\n\n")
  
  # Default to non-parametric Wilcoxon Rank-Sum Test to accommodate 
  # large sample sizes (N > 5000 limitations in Shapiro-Wilk) and highly skewed distributions
  cat("=> Running Unpaired Wilcoxon Rank-Sum Test (Non-parametric)...\n\n")
  
  stat_test_res <- wilcox.test(vec_g1, vec_g2, paired = FALSE, exact = FALSE) 
  
  cat("--------------------------------------------------\n")
  cat("2. STATISTICAL TEST RESULTS\n")
  cat("--------------------------------------------------\n")
  print(stat_test_res)
  
  # Plot overlapping density distributions for visual comparison
  p_dist <- ggplot(df_clean %>% filter(!!sym(group_col) %in% c(group1, group2)), 
                   aes(x = !!sym(val_col), fill = !!sym(group_col))) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Distribution of", val_col),
         subtitle = paste("Test:", stat_test_res$method),
         x = val_col,
         y = "Density") +
    theme_minimal()
  
  print(p_dist)
  
  # Return objects silently for optional assignment
  invisible(list(
    stat_test = stat_test_res,
    plot = p_dist
  ))
}

mutation_events_global %>% head(3)

high_mut_rate<- mutation_events_global %>% filter(mut_rate == 0.01)

resultados <- perform_unpaired_analysis(
  data = high_mut_rate, 
  group_col = "mut_type", 
  val_col = "abundance", 
  group1 = "Neutral", 
  group2 = "Beneficial"
)

low_mut_rate<- mutation_events_global %>% filter(mut_rate == 0.001)

resultados <- perform_unpaired_analysis(
  data = low_mut_rate, 
  group_col = "mut_type", 
  val_col = "abundance", 
  group1 = "Neutral", 
  group2 = "Beneficial"
)

### Figure 4E

#### Mean pairwise Levenshtein distance over time

db_load_id_paired =  db_load_id %>% 
    filter(mu > 0.0005, resource == "NOT") %>%
    mutate(rep_id = ifelse(rep_id > 200, rep_id - 200, rep_id)) %>%
    arrange(rep_id, update)
    
db_load_id_paired %>% head
db_load_id_paired %>% tail

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
# MAIN FUNCTION
# Usage examples:
#   plot_diversity_trajectory(rep_ids = 1)          # single replicate
#   plot_diversity_trajectory(rep_ids = c(1, 2, 3)) # several replicates
#   plot_diversity_trajectory(rep_ids = 1:10)        # range
# ==============================================================================
plot_diversity_trajectory <- function(rep_ids) {
 
    cat("Computing Levenshtein trajectories for rep_ids:",
        paste(rep_ids, collapse=", "), "\n")
 
    # Filter to selected replicates
    db_subset <- db_load_id_paired %>%
        filter(rep_id %in% rep_ids)
 
    if (nrow(db_subset) == 0) {
        stop("No data found for the requested rep_ids. Check your ids.")
    }
 
    # Compute within-group distance per (rep_id x mu x update)
    lev_traj <- db_subset %>%
        group_by(rep_id, mu, update) %>%
        group_modify(~ run_levenshtein(.x$sequence)) %>%
        ungroup() %>%
        filter(!is.na(mean)) %>%
        mutate(
            mu_label = ifelse(mu == 0.01, "mu = 0.01", "mu = 0.001"),
            mu_label = factor(mu_label,
                                    levels = c("mu = 0.001", "mu = 0.01"))
        )
 
    # X axis: normalise update to [0, 1] so all replicates are comparable
    # regardless of total length (start = 0, end = 1)
    lev_traj <- lev_traj %>%
        group_by(rep_id, mu) %>%
        mutate(
            upd_min  = min(update),
            upd_max  = max(update),
            upd_norm = (update - upd_min) / (upd_max - upd_min)
        ) %>%
        ungroup()
 
    # Build plot
    n_reps <- length(unique(lev_traj$rep_id))
 
    # Alpha and linewidth scale with number of replicates
    line_alpha <- ifelse(n_reps == 1, 1, max(0.15, 1 / sqrt(n_reps)))
    lw_rep     <- ifelse(n_reps == 1, 1.2, 0.5)
 
    p <- ggplot(lev_traj,
                aes(x     = upd_norm,
                    y     = mean,
                    color = mu_label,
                    group = interaction(rep_id, mu)))
 
    # If more than one replicate: draw individual lines faintly + median on top
    if (n_reps > 1) {
        # Compute median trajectory across selected replicates
        median_traj <- lev_traj %>%
            # Bin the normalised x into 100 equal bins for smooth median
            mutate(x_bin = round(upd_norm, 2)) %>%
            group_by(mu_label, x_bin) %>%
            summarise(median_mean = median(mean),
                      q25         = quantile(mean, 0.25),
                      q75         = quantile(mean, 0.75),
                      .groups     = "drop")
 
        p <- p +
            geom_line(alpha = line_alpha, linewidth = lw_rep) +
            geom_ribbon(data = median_traj,
                        aes(x = x_bin, ymin = q25, ymax = q75,
                            fill = mu_label, group = mu_label),
                        alpha = 0.15, color = NA, inherit.aes = FALSE) +
            geom_line(data = median_traj,
                      aes(x = x_bin, y = median_mean,
                          color = mu_label, group = mu_label),
                      linewidth = 1.4, inherit.aes = FALSE)
    } else {
        p <- p + geom_line(linewidth = lw_rep)
    }
 
    title_str <- if (n_reps == 1) {
        paste0("Genotypic diversity trajectory - replicate ", rep_ids)
    } else {
        paste0("Genotypic diversity trajectory - ",
               n_reps, " replicates (bold = median, ribbon = IQR)")
    }
 
    p <- p +
        scale_color_discrete(labels = c("mu = 0.001" = "Low Mutation Rate",
                                        "mu = 0.01"  = "High Mutation Rate")) +
        scale_x_continuous(breaks = c(0, 1),
                           labels = c("start", "end"),
                           expand = expansion(mult = 0.02)) +
        scale_y_continuous(labels = scales::label_percent()) +
        labs(title  = title_str,
             x      = "Adaptation Process",
             y      = "Pairwise Genotype Distance",
             color  = NULL,
             fill   = NULL) +
        theme_bw() +
        theme(
            plot.title        = element_blank(),
            axis.title.x      = element_text(size=20, margin=margin(t=1)),

            axis.title.y      = element_text(size=20, margin=margin(r=10)),
            axis.text.x       = element_text(color="black", size=20),
            axis.text.y       = element_text(color="black", size=18),
            axis.ticks        = element_line(linewidth=0.5),
            axis.ticks.length = unit(0.25, "cm"),
            panel.border      = element_rect(color="black", fill=NA, linewidth=1),
            panel.grid.major  = element_line(color="gray90", linetype="solid"),
            legend.position   = "top",
            legend.title      = element_text(size=20),
            legend.text       = element_text(size=18),
            legend.key.size   = unit(1.0, "cm")
        )
 
    print(p)
    invisible(lev_traj)   # return data silently for further use
}

options(repr.plot.width = 14, repr.plot.height = 5)

plot_diversity_trajectory(rep_ids = 6)          

### Figure 5A. Low mutation rate 

abundance_time_course_lineages(1405, 0.001, 97, min_update = 300000, max_update = 500000)

### Figure 5B. High mutation rate

options(repr.plot.width = 10, repr.plot.height = 6)

abundance_time_course_lineages(1405, 0.01, 297, min_update = 300000, max_update = 500000)

# ---

### Supplementary Figures

#### Supplementary Figure 2

options(repr.plot.width = 16, repr.plot.height = 12)

abundance_time_course_effects(1776,0.001) / abundance_time_course_effects(1776,0.01)
abundance_time_course_effects(1405,0.001) / abundance_time_course_effects(1405,0.01)

#### Supplementary Figure 3

options(repr.plot.width = 14, repr.plot.height = 5)
plot_diversity_trajectory(rep_ids = 1)          

plot_diversity_trajectory(rep_ids = 97)          

plot_diversity_trajectory(rep_ids = 188)          

#### Supplementary Figure 4

options(repr.plot.width = 14, repr.plot.height = 5)

abundance_time_course_lineages(1221, 0.001, 75, min_update = 370000, max_update = 500000)

options(repr.plot.width = 14, repr.plot.height = 5)

abundance_time_course_lineages(1405, 0.001, 37, min_update = 300000, max_update = 500000)

options(repr.plot.width = 14, repr.plot.height = 5)

abundance_time_course_lineages(1271, 0.01, 346, min_update = 350000, max_update = 500000)

options(repr.plot.width = 14, repr.plot.height = 5)

abundance_time_course_lineages(1776, 0.01, 250, min_update = 370000, max_update = 500000)

# ---

## Reproducibility Session Info

sessionInfo()

# ---
# ---
