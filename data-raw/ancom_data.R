## ECAM DATASET
ecam_table_taxa <-
  "data-raw/ecam-table-taxa.tsv" %>%
  readr::read_tsv(skip = 1) %>%
  dplyr::rename(otu_id = `feature-id`) %>%
  data.frame(check.names = FALSE)

rownames(ecam_table_taxa) <- ecam_table_taxa$otu_id

ecam_sample_metadata <-
  "data-raw/ecam-sample-metadata.tsv" %>%
  readr::read_tsv() %>%
  .[-1, ] %>%
  dplyr::rename(Sample.ID = `#SampleID`)

ecam_prepro <- feature_table_pre_process(
  feature_table = ecam_table_taxa,
  meta_data = ecam_sample_metadata,
  sample_var = "Sample.ID",
  group_var = "delivery",
  out_cut = 0.05,
  zero_cut = 0.90,
  lib_cut = 0,
  neg_lb = TRUE
)

## MOVING PICS DATASET
moving_pics_table_taxa <-
  "data-raw/moving-pics-table.tsv" %>%
  readr::read_tsv(skip = 1) %>%
  dplyr::rename(otu_id = `feature-id`) %>%
  data.frame(row.names = 1, check.names = FALSE)

# rownames(moving_pics_table_taxa) <- moving_pics_table_taxa$otu_id

moving_pics_sample_metadata <-
  "data-raw/moving-pics-sample-metadata.tsv" %>%
  readr::read_tsv() %>%
  .[-1, ] %>%
  dplyr::rename(Sample.ID = SampleID)

moving_pics_prepro <- feature_table_pre_process(
  feature_table = moving_pics_table_taxa,
  meta_data = moving_pics_sample_metadata,
  sample_var = "Sample.ID",
  group_var = NULL,
  out_cut = 0.05,
  zero_cut = 0.90,
  lib_cut = 1000,
  neg_lb = TRUE
)

usethis::use_data(
  ecam_sample_metadata,
  ecam_table_taxa,
  ecam_prepro,
  moving_pics_sample_metadata,
  moving_pics_table_taxa,
  moving_pics_prepro,
  internal = TRUE,
  overwrite = TRUE
)
