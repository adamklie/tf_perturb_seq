############################################################################
### CREATE SPEC FILE FOR HELEN'S TF PERTURB-SEQ DATA
### Sara Geraghty, June 2025
############################################################################

library(dplyr)
library(stringr)
library(tibble)
library(readr)
library(tidyr)

PATH <- '/hpc/group/gersbachlab/'
PATH_OUT <- paste0(PATH, 'seg95/crispr-pipeline-personal/example-data')

# Data paths (note: no hash files)
data_crispr_p1 <- paste0(PATH, 'hls/IGVF_TF_Perturbseq/Complete_Screen/3-23-25_TF_perturbseq_full_sequencing/analysis/mkfastq/CRISPR_novaseq_1/fastqs')
data_gex_p1 <- paste0(PATH, 'hls/IGVF_TF_Perturbseq/Complete_Screen/3-23-25_TF_perturbseq_full_sequencing/analysis/mkfastq/GEX_novaseq_1/fastqs')
data_crispr_p2 <- paste0(PATH, 'hls/IGVF_TF_Perturbseq/Complete_Screen/3-23-25_TF_perturbseq_full_sequencing/analysis/mkfastq/CRISPR_novaseq_2/fastqs')
data_gex_p2 <- paste0(PATH, 'hls/IGVF_TF_Perturbseq/Complete_Screen/3-23-25_TF_perturbseq_full_sequencing/analysis/mkfastq/GEX_novaseq_2/fastqs')

# Static paths
barcode_onlist_path <- '/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/737K-august-2016.txt'
guide_metadata_path <- '/hpc/group/gersbachlab/seg95/crispr-pipeline-personal/example-data/harmonized_guide_file_poolabcdf.tsv'

# Required format: tsv file with the following columns: R1_path,R2_path,file_modality,measurement_sets,sequencing_run,lane,seqspec,barcode_onlist,guide_design,barcode_hashtag_map

# Build a dataframe from each set of files
gather_seqspec_entries <- function(data_path, modality, sequencing_run_label) {
  files <- list.files(data_path, pattern = "_R[12]_001.fastq.gz$", full.names = TRUE)
  if (length(files) == 0) {
    warning(paste("No FASTQ files found in:", data_path))
    return(tibble())
  }

  # Extract sample ID (e.g., A1) and read type (R1 or R2)
  info <- tibble(
    path = files,
    basename = basename(files),
    sample = str_extract(basename, "^[^_]+"),  # e.g., A1, B2
    read = str_extract(basename, "_R[12]_") %>% str_remove_all("_"),
    unique_id = str_extract(basename, "^[^_]+_S\\d+_L\\d{3}"),
    lane = str_extract(basename, "L\\d{3}") %>% str_remove("L") %>% as.integer()
  )

  # Spread into R1 and R2 columns
  paired <- info %>%
    select(unique_id, read, path, lane) %>%
    pivot_wider(names_from = read, values_from = path, values_fn = first) %>%
    filter(!is.na(R1) & !is.na(R2)) %>%
    mutate(
      file_modality = modality,
      measurement_sets = paste0(unique_id, "_r", sequencing_run_label),
      sequencing_run = as.numeric(sequencing_run_label),
      seqspec = case_when(
        modality == "scRNA" ~ "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/example_data/yaml_files/rna_seqspec.yml",
        modality == "gRNA" ~ "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/example_data/yaml_files/guide_seqspec.yml",
        TRUE ~ ""
      ),
      barcode_onlist = barcode_onlist_path,
      guide_design = guide_metadata_path,
      barcode_hashtag_map = ""    
  ) %>%
    select(
      R1_path = R1, 
      R2_path = R2, 
      file_modality,
      measurement_sets,
      sequencing_run,
      lane,
      seqspec,
      barcode_onlist,
      guide_design,
      barcode_hashtag_map
    )
  
  return(paired)
}

# Generate spec entries for all 4 data sources
spec_crispr_1 <- gather_seqspec_entries(data_crispr_p1, "gRNA", "1")
spec_crispr_2 <- gather_seqspec_entries(data_crispr_p2, "gRNA", "2")
spec_gex_1    <- gather_seqspec_entries(data_gex_p1,    "scRNA", "1")
spec_gex_2    <- gather_seqspec_entries(data_gex_p2,    "scRNA", "2")

# Combine into one data frame
samplesheet <- bind_rows(spec_crispr_1, spec_crispr_2, spec_gex_1, spec_gex_2) %>%
  arrange(measurement_sets, file_modality, sequencing_run, lane)

# Write to output TSV, comma-separated (weird)
write.table(samplesheet, file = file.path(PATH_OUT, "helen_samplesheet_seqspec_poolabcdf.tsv"), sep = ",", quote = F, row.names = F)

# Optional: split samplesheet by a grouping variable to make smaller versions
split_by <- "sequencing_run"  
samplesheet_split <- split(samplesheet, samplesheet[[split_by]])
# Write each group to a separate file
for (grp in names(samplesheet_split)) {
  outfile <- file.path(PATH_OUT, paste0("helen_samplesheet_", split_by, "_", grp, ".tsv"))
  write.table(samplesheet_split[[grp]], file = outfile, sep = ",", quote = FALSE, row.names = FALSE)
  message("Wrote: ", outfile)
}
