library(dplyr)
library(readr)
library(tidyr)
library(purrr)

# Load guide metadata
guide_file <- "/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/harmonized_guide_file_poolabcd.tsv"
guides <- read_tsv(guide_file)

# 1. Check for duplicated spacer sequences
spacer_dups <- guides %>%
  group_by(spacer) %>%
  filter(n() > 1) %>%
  ungroup()

cat("=== Total duplicated spacers: ", n_distinct(spacer_dups$spacer), " ===\n\n")

compare_within_spacer <- spacer_dups %>%
  mutate(across(everything(), as.character)) %>%
  group_by(spacer) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-c(spacer, row_id)) %>%
  group_by(spacer, name) %>%
  summarise(
    n_unique = n_distinct(value),
    values = paste(unique(value), collapse = " | "),
    .groups = "drop"
  ) %>%
  filter(n_unique > 1) %>%
  arrange(spacer, desc(n_unique))

cat("=== Columns with differing values across duplicated spacers ===\n")
print(compare_within_spacer)

# 3. Optional: Show raw rows for a few example duplicated spacers
example_spacers <- spacer_dups %>%
  count(spacer, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(spacer)

cat("\n=== Example raw rows for duplicated spacers (differing columns only) ===\n")

for (s in example_spacers) {
  sub <- spacer_dups %>% filter(spacer == s)
  differing_cols <- sub %>%
    mutate(across(everything(), as.character)) %>%
    summarise(across(everything(), ~ n_distinct(.) > 1)) %>%
    pivot_longer(everything(), names_to = "column", values_to = "differs") %>%
    filter(differs) %>%
    pull(column)
  
  cat("\n--- Spacer:", s, "---\n")
  print(sub[, unique(c("spacer", differing_cols))], n = Inf)
}

# 2. Check for spacer length inconsistencies
guide_lengths <- guides %>%
  mutate(spacer_length = nchar(spacer)) %>%
  group_by(spacer_length) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cat("\n=== Spacer length distribution ===\n")
print(guide_lengths)

# List all non-modal length spacers
modal_length <- guide_lengths$spacer_length[1]

off_length <- guides %>%
  filter(nchar(spacer) != modal_length)

cat("\n=== Spacers with non-standard length ===\n")
print(off_length)
