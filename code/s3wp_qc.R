# Mai Nguyen 11/11/2024
# S3WP analysis (Olink HT and Olink Target)
# QC code 

# import require lib
library(dplyr)
library(impute)
library(openxlsx)
library(ggplot2)
library(umap)
library(tidyr)

# Functions for QC
# Function to filter a data frame by a specified column and value
filter_by_column_value = function(data, column, value) {
  filtered_data = data[data[[column]] == value, ]
  return(filtered_data)}

# Function to divide a df into sub-df based on a column
divide_df_by_column = function(data, column) {
  unique_values = unique(data[[column]])
  for (value in unique_values) {
    new_var_name = paste0(deparse(substitute(data)), '_', column, value)
    sub_data = data[data[[column]] == value, ]
    assign(new_var_name, sub_data, envir = .GlobalEnv)}}

# Function to calculate warning counts for samples
calculate_sample_warnings = function(data) {
  data %>%
    group_by(DAid) %>%
    summarize(
      sample_warn_count = sum(SampleQC == 'FAIL', na.rm = TRUE),
      sample_total_count = n(),
      sample_warn_percentage = (sample_warn_count / sample_total_count) * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(sample_warn_count))}

# Function to calculate warning counts for proteins
calculate_protein_warnings = function(data) {
  data %>%
    group_by(UniProt) %>%
    summarize(
      protein_warn_count = sum(AssayQC == 'WARN', na.rm = TRUE),
      protein_total_count = n(),
      protein_warn_percentage = (protein_warn_count / protein_total_count) * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(protein_warn_count))}

# Function to filter samples with more than 50% warnings
filter_samples = function(data, sample_warning_counts) {
  filtered_data = data %>%
    inner_join(sample_warning_counts, by = 'DAid') %>%
    {
      removed_samples = filter(., sample_warn_percentage > 50) %>%
        group_by(DAid) %>%
        summarize(records_removed = n(), .groups = 'drop')
      cat('Samples removed:', nrow(removed_samples), '(', sum(removed_samples$records_removed), 'records)\n')
      filter(., sample_warn_percentage <= 50)
    } %>%
    select(-sample_warn_count, -sample_total_count, -sample_warn_percentage)
  return(filtered_data)}

# Function to filter proteins with more than 50% warnings
filter_proteins = function(data, protein_warning_counts) {
  filtered_data = data %>%
    inner_join(protein_warning_counts, by = 'UniProt') %>%
    {
      removed_proteins = filter(., protein_warn_percentage > 50) %>%
        group_by(UniProt) %>%
        summarize(records_removed = n(), .groups = 'drop')
      cat('Proteins removed:', nrow(removed_proteins), '(', sum(removed_proteins$records_removed), 'records)\n')
      filter(., protein_warn_percentage <= 50)
    } %>%
    select(-protein_warn_count, -protein_total_count, -protein_warn_percentage)
  return(filtered_data)}

# Function to plot LOD distribution
plot_lod_distribution = function(data, lod_col, lod_type) {
  ggplot(data, aes(x = .data[[lod_col]])) +
    geom_histogram(binwidth = 0.2, fill = 'skyblue', color = 'black', alpha = 0.7) +
    labs(
      title = paste0('LOD distribution of ', lod_type),
      x = 'LOD',
      y = 'Frequency'
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))}

# Function to plot QC UMAP
plot_qc_umap = function(x, y, method, num) {
  title = paste0(method, '\nNumber of samples=', num)
  plot(x, y, pch = 19, xlab = 'UMAP1', ylab = 'UMAP2', main = title)}


# Run code
# Read in data
olink_ht = read.csv('C:/Users/maihu/wellness/data/onlink_ht/olink_ht.csv', header = T)
olink_target = read.delim('C:/Users/maihu/wellness/data/olink_target/olink_target.txt', sep = '\t', header = T)

# QC for Olink HT
# Calculate warning counts
sample_warning_counts = calculate_sample_warnings(olink_ht)
protein_warning_counts = calculate_protein_warnings(olink_ht)

# Filter samples with > 50% warnings
olink_ht_fil1 = filter_samples(olink_ht, sample_warning_counts)

# Filter proteins with > 50% warnings
olink_ht_fil2 = filter_proteins(olink_ht_fil1, protein_warning_counts)

# Filter measurements with warnings 
# MISSING: information how to assess measurements to filter warnings
# (no other columns having WARN/FAIL)

# Select one of the replicate assays
# Proteins with multiple measurements: GBP1 (UniProt: P32455), MAP2K1 (UniProt: Q02750)
# Filtered by UniProt 
gbp1_measurements = filter_by_column_value(olink_ht_fil2, 'UniProt', 'P32455')
map2k1_measurements = filter_by_column_value(olink_ht_fil2, 'UniProt', 'Q02750')
# There are multiple blocks for these proteins
# Separate by blocks to see LOD distribution
divide_df_by_column(gbp1_measurements, 'Block')
divide_df_by_column(map2k1_measurements, 'Block')
# LOD distribution 
plot_lod_distribution(data = gbp1_measurements_Block3, lod_col = 'LOD', lod_type = 'GBP1 Block 3') # highest
plot_lod_distribution(data = gbp1_measurements_Block4, lod_col = 'LOD', lod_type = 'GBP1 Block 4')
plot_lod_distribution(data = gbp1_measurements_Block5, lod_col = 'LOD', lod_type = 'GBP1 Block 5')
plot_lod_distribution(data = map2k1_measurements_Block3, lod_col = 'LOD', lod_type = 'MAP2K1 Block 3')
plot_lod_distribution(data = map2k1_measurements_Block4, lod_col = 'LOD', lod_type = 'MAP2K1 Block 4') # highest
plot_lod_distribution(data = map2k1_measurements_Block5, lod_col = 'LOD', lod_type = 'MAP2K1 Block 5')

# LOD distribution plot to assess quality
plot_lod_distribution(data = olink_ht_fil2, lod_col = 'LOD', lod_type = 'Olink HT')

# Plot QC UMAP
# Divide data into meta data and count data 
metadata_ht = olink_ht_fil2 %>%
  select(DAid, OlinkID, visit, barcode) %>%
  distinct()  

ht_count = olink_ht_fil2 %>%
  select(UniProt, NPX, barcode)

npx_ht_wide = reshape(data = ht_count, timevar = 'UniProt', idvar = 'barcode', direction = 'wide')

npx_ht = npx_ht_wide %>%
  select(-barcode)

npx_ht = npx_ht %>%
  mutate(across(everything(), ~ as.numeric(.)))

# Remove NA values
na_npx_ht = sum(is.na(npx_ht))
cat('Total number of missing values (NA):', na_npx_ht, '\n')
clean_ht_wide = na.omit(npx_ht)

umapxy_ht = umap(clean_ht_wide, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')

umap_ht = data.frame(x = umapxy_ht$layout[, 1], y = umapxy_ht$layout[, 2])

plot_qc_umap(x = umap_ht$x, y = umap_ht$y, method = 'UMAP QC Olink HT', num = nrow(umap_ht))


# QC for Olink target
# Assumption: done filter warnings (no warning data)
# NA values (remove, can change to impute if wanted)
na_target = sum(is.na(olink_target))
cat('Total number of missing values (NA):', na_target, '\n')
clean_olink_target = na.omit(olink_target)

# LOD distribution plot 
# lod_table_long = calculate_lod_wide(data = clena_olink_target,
#   id_col = 'sampleID',
#   blank_ids = # list of blank ids,
#   exclude_cols = c('subject_id', 'visit'))

# Plot LOD distribution
# plot_lod_distribution(data = lod_table_long, lod_col = 'lod', protein_id_col = 'protein')

# Plot QC UMAP
metadata_target = clean_olink_target %>%
  select(sampleID, subject_id, visit)

npx_target = clean_olink_target %>%
  select(-sampleID, -subject_id, -visit)

umapxy_target = umap(npx_target, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')

umap_target = data.frame(x = umapxy_target$layout[, 1], y = umapxy_target$layout[, 2])

plot_qc_umap(x = umap_target$x, y = umap_target$y, method = 'UMAP QC Olink Target', num = nrow(umap_target))


