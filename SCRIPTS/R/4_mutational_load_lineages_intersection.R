# Reproducible Code for: Neutral mutations promote unbounded adaptation after clonal interference.

### Salvador León Fernández, Miguel Ángel Fortuna Alcolado.

# Here, we will merge `db_mut_load.csv` with `db_lineages.csv` to extract the main lineage of organisms belonging to the evolutionary lines that survived until the end of the experiment, thus obtaining their abundance data.

# ---

## Setup & Environment

library("tidyverse")
library("data.table")

db_abundance = fread("../../DATA/db_mut_load.csv")
db_abundance %>% head
db_abundance %>% tail
db_abundance %>% nrow

db_lineages = fread("../../DATA/db_lineages.csv")
db_lineages %>% head
db_lineages %>% tail
db_lineages %>% nrow

## Extract lineage abundance

# Deduplicate lineages, keeping only the earliest occurrence per sequence 
# to prevent many-to-many joins
db_lineages_prep <- db_lineages %>%
  select(
    founder_id, 
    rep_id, 
    mut_rate,             
    sequence,             
    update_born = update, 
    lineage_id, 
    is_main_lineage,      
    fitness, 
    relative_fitness      
  ) %>%
  arrange(update_born) %>% 
  distinct(founder_id, rep_id, sequence, mut_rate, .keep_all = TRUE)

# Subset abundance data to essential columns for memory efficiency
db_abundance_prep <- db_abundance %>%
  select(
    founder_id, 
    rep_id, 
    mu,                   
    sequence,             
    update_sample = update, 
    abundance
  )

# Merge datasets via primary keys and order chronologically by sampling time
db_final <- db_lineages_prep %>%
  inner_join(
    db_abundance_prep, 
    by = c("founder_id", "rep_id", "sequence", "mut_rate" = "mu")
  ) %>%
  arrange(update_sample)

db_final <- db_final %>% arrange(founder_id, rep_id, update_sample)
db_final %>% head
db_final %>% tail
db_final %>% nrow

write_csv(db_final, file = "../../DATA/db_joined_mut_load_lineages.csv")
