# Reproducible Code for: Neutral mutations promote unbounded adaptation after clonal interference.

### Salvador León Fernández, Miguel Ángel Fortuna Alcolado.

# Here, we will merge the data of the most abundant organism every 1000 updates from `db_mut_load.csv` with `db_mut_supply.csv` to extract its fitness measurement. We also included an ancestor file because `db_mut_supply` only contains mutants and not the original sequences.

# ---

## Setup & Environment

library("data.table")
library("tidyverse")
library("scales")

## Ancestor file

db_fixed <- read_csv("../../DATA/db_fixed_mutations.csv",
                     col_types = cols("i", "i", "i", "i", "c", "d", "c", "i", "d", "d", "c", "i")) %>%
    mutate(sequence = str_remove_all(sequence, "\\*")) %>%
    arrange(desc(resource), rep_id, update) %>%
    mutate(length = nchar(sequence))
db_fixed %>% nrow
db_fixed %>% head(5)
db_fixed %>% tail(5)

ancestors = db_fixed %>% filter(update == 0, mut_rate > 0.0005) %>% group_by(founder_id, rep_id, mut_rate, update, fitness, relative_fitness, sequence) %>% summarise(.groups = "drop")
ancestors %>% head
ancestors %>% tail

write_csv(ancestors, file = "../../DATA/db_all_ancestors.csv")

## Extract most abundant organism's fitness values

cat("1. Loading mutational_load.csv and filtering dominant lineages...\n")
db_load <- fread("../../DATA/db_mut_load.csv")

dominant_lineages <- db_load %>%
  filter(abundance == max(abundance), .by = c(update, mu, rep_id, founder_id)) %>%
  rename(sample_update = update) %>%
  as.data.table()

rm(db_load)
gc()

cat("2. Loading the unified ancestor file...\n")
ancestors_db <- fread("../../DATA/db_all_ancestors.csv")

cat("3. Loading fitness_all_mutations.csv (essential columns only)...\n")
all_mutations_db <- fread(
  "../../DATA/db_mut_supply.csv",
  select = c("founder_id", "rep_id", "mut_rate", "update", "fitness_org_id", "relative_fitness", "seq")
)

cat("4. Merging data and filling missing values with ancestors...\n")

final_dataset <- dominant_lineages %>%
  left_join(
    all_mutations_db,
    by = c("founder_id" = "founder_id", 
           "rep_id" = "rep_id", 
           "mu" = "mut_rate", 
           "sequence" = "seq")
  ) %>%
  rename(origin_update_mut = update) %>%
  left_join(
    ancestors_db,
    by = c("founder_id" = "founder_id", 
           "rep_id" = "rep_id", 
           "mu" = "mut_rate", 
           "sequence" = "sequence"),
    suffix = c("_mut", "_anc")
  ) %>%
  mutate(
    fitness_org_id = coalesce(fitness_org_id, fitness),
    relative_fitness = coalesce(relative_fitness_mut, relative_fitness_anc),
    origin_update = coalesce(as.numeric(origin_update_mut), 0)
  ) %>%
  select(-fitness, -relative_fitness_mut, -relative_fitness_anc, -origin_update_mut) %>%
  filter(!is.na(fitness_org_id)) %>%
  arrange(founder_id, mu, rep_id, sample_update)

cat("5. Process complete. Releasing RAM...\n")
rm(all_mutations_db, ancestors_db, dominant_lineages)
gc()

head(final_dataset)

final_dataset = final_dataset %>% filter(mu > 0.0005)

fwrite(final_dataset, "../../DATA/db_joined_most_abundant_organism.csv")
