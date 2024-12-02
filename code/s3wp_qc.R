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
library(stringr)
library(effectsize)

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

# Function to extract information from a specific column
extract_info_from_column = function(df, source_column, pattern, new_column_name) {
  df %>%
    mutate(!!new_column_name := str_extract(.data[[source_column]], pattern) %>%
        str_replace('.*?: ', ''))}

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

# Function to analyze protein detection with LOD
analyze_proteins = function(data, detection_percentage = 100) {
  sample_count = length(unique(data$barcode))
  detection_threshold = ceiling((detection_percentage / 100) * sample_count)
  protein_summary = data %>%
    group_by(UniProt) %>%
    summarise(total_measurements = n(),
      detected_in_samples = sum(!is.na(NPX) & !is.na(LOD) & NPX >= LOD),
      detected = ifelse(detected_in_samples >= detection_threshold, 'Yes', 'No'))
  higher_than_lod_count = sum(protein_summary$detected == 'Yes')
  lower_than_lod_count = sum(protein_summary$detected == 'No')
  total_proteins = nrow(protein_summary)
  higher_than_lod_percentage = (higher_than_lod_count / total_proteins) * 100
  lower_than_lod_percentage = (lower_than_lod_count / total_proteins) * 100
  cat('Number of proteins detected in at least', detection_percentage, '% of samples:', higher_than_lod_count, '(', round(higher_than_lod_percentage, 1), '%)', '\n')
  cat('Number of proteins detected in less than', detection_percentage, '% of samples:', lower_than_lod_count, '(', round(lower_than_lod_percentage, 1), '%)', '\n')
  return(protein_summary)}

# Run code
# Read in data
olink_ht_raw = read.csv('C:/Users/maihu/wellness/data/onlink_ht/olink_ht.csv', header = T)
olink_target = read.delim('C:/Users/maihu/wellness/data/olink_target/olink_target.txt', sep = '\t', header = T)
subject_info = read.xlsx('C:/Users/maihu/wellness/data/onlink_ht/subjects_2024-11-21.xlsx', sheet = 1)
target_id = read.csv('C:/Users/maihu/wellness/data/Annoteringsfil_Olink.csv', header = T, sep = ';')

# QC for Olink HT
# Calculate warning counts
sample_warning_counts = calculate_sample_warnings(olink_ht_raw)
protein_warning_counts = calculate_protein_warnings(olink_ht_raw)

# Filter samples with > 50% warnings
olink_ht_fil1 = filter_samples(olink_ht_raw, sample_warning_counts)

# Filter proteins with > 50% warnings
olink_ht_fil2 = filter_proteins(olink_ht_fil1, protein_warning_counts)

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

# Filter low LOD assays
olink_ht = olink_ht_fil2 %>%  filter(!(UniProt == 'P32455' & (Block == 4 | Block == 5)))
olink_ht = olink_ht %>%  filter(!(UniProt == 'Q02750' & (Block == 3 | Block == 5)))

# LOD distribution plot to assess quality
plot_lod_distribution(data = olink_ht_fil2, lod_col = 'LOD', lod_type = 'Olink HT')

# Test how many protein detected more than LOD
olink_ht_fil2 = olink_ht_fil2 %>% mutate(under_LOD = ifelse(NPX < LOD, 'Yes', 'No'))

# Number of protein compared with LOD
protein_detected_100 = analyze_proteins(olink_ht_fil2)
protein_detected_80 = analyze_proteins(olink_ht_fil2, detection_percentage = 80)

# Extract subject_id for meta data 
olink_ht = extract_info_from_column(olink_ht, 'Extra.data', 'subject_id: \\S+', 'subject_id')
# Add sampleID to the column
olink_ht$sampleID = paste0('WEL_', olink_ht$visit, '_', olink_ht$subject_id)
olink_ht = olink_ht %>% mutate(subject_id = str_replace(subject_id, '.*?-', ''))

# Add gender to subjects
subject_info = subject_info %>% rename(subject_id = Subject)
subject_info = subject_info %>% mutate(subject_id = as.character(subject_id))
olink_ht = olink_ht %>% left_join(subject_info %>% select(subject_id, Sex), by = 'subject_id')

# Plot QC 
# Divide data into meta data and count data 
metadata_ht = olink_ht %>% select(DAid, visit, barcode, subject_id, Sex, sampleID) %>% distinct()  

# count data go with barcode bcs there are some dup in patient-visit (sampleID)
ht_count = olink_ht %>% select(UniProt, NPX, barcode)
npx_ht_wide = ht_count %>% pivot_wider(names_from = UniProt, values_from = 'NPX')
npx_ht = npx_ht_wide %>% select(-barcode)
npx_ht_matrix = as.matrix(npx_ht)

# Run the code if data is not yet numeric
# npx_ht = npx_ht %>% mutate(across(everything(), ~ as.numeric(.)))

# Impute missing data
imputed_olink_ht_npx = impute.knn(npx_ht_matrix)
imputed_olink_ht = imputed_olink_ht_npx$data
imputed_olink_ht = as.data.frame(imputed_olink_ht)

# Plot UMAP
umapxy_ht = umap(imputed_olink_ht, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_ht = data.frame(x = umapxy_ht$layout[, 1], y = umapxy_ht$layout[, 2])
umap_ht$barcode = npx_ht_wide$barcode
umap_ht = umap_ht %>% left_join(metadata_ht %>% select(barcode, Sex, visit, subject_id), by = 'barcode')

# UMAP colored by visit
ggplot(umap_ht, aes(x = x, y = y, color = factor(visit))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink HT colored by Visit',
    x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  scale_color_manual(values = c('#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#e41a1c', '#ffff33')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    axis.title = element_text(size = 14), legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# UMAP colored by sex
ggplot(umap_ht, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink HT colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    axis.title = element_text(size = 14), legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# UMAP colored by patient (for the consistency of patient profile throughout visits)
# No legend bcs there are too many
ggplot(umap_ht, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink HT colored by Patient',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  # facet_wrap(~subject_id) + # if you want to divide into sub-plot of different factors
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# QC for Olink target
# Assumption: done filter warnings (no warning data)
# Divide into meta data and count 
metadata_target = olink_target %>% select(sampleID, subject_id, visit)
metadata_target = metadata_target %>% mutate(subject_id = str_replace(subject_id, '.*?-', ''))
metadata_target = metadata_target %>% left_join(subject_info %>% select(subject_id, Sex), by = 'subject_id')

npx_target = olink_target %>% select(-sampleID, -subject_id, -visit)

# Impute missing data
npx_target_matrix = as.matrix(npx_target)
imputed_olink_target_npx = impute.knn(npx_target_matrix)
imputed_olink_target = imputed_olink_target_npx$data
imputed_olink_target = as.data.frame(imputed_olink_target)

# Plot UMAP
umapxy_target = umap(imputed_olink_target, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_target = data.frame(x = umapxy_target$layout[, 1], y = umapxy_target$layout[, 2])
umap_target = cbind(umap_target, metadata_target)

# UMAP colored by visit
ggplot(umap_target, aes(x = x, y = y, color = factor(visit))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink Target colored by Visit',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  scale_color_manual(values = c('#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#e41a1c', '#ffff33')) +
#  facet_wrap(~visit) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# UMAP colored by sex
ggplot(umap_target, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink Target colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title = element_text(size = 14), legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# UMAP colored by patient
ggplot(umap_target, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink Target colored by Patient ID',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  # scale_color_manual(values = c('#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#e41a1c', '#ffff33')) +
  #facet_wrap(~visit) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))


# Platform comparison
# Replace OID with Uniprot ID where relevant
uniprot_target = target_id$Uniprot[match(names(imputed_olink_target), target_id$OlinkID)]
# Only replace OID with single value UniProt
new_colnames = names(imputed_olink_target)
for (i in seq_along(new_colnames)) {replacement = uniprot_target[i]
if (!is.na(replacement) && !grepl(',', replacement)) {
  new_colnames[i] = replacement}}
names(imputed_olink_target) = new_colnames
# If there are OID with no or multiple UniProt IDs, filter them out for manual evaluation 
problematic_cases = data.frame(OID = names(imputed_olink_target), Uniprot = uniprot_target
) %>% filter(is.na(Uniprot) | grepl(',', Uniprot))
# Replace OID with UniProt ID after manual inspection 
colnames(imputed_olink_target)[colnames(imputed_olink_target) == 'OID00368'] = 'P29459'
colnames(imputed_olink_target)[colnames(imputed_olink_target) == 'OID00402'] = 'Q8NEV9'
colnames(imputed_olink_target)[colnames(imputed_olink_target) == 'OID00723'] = 'Q29983'
colnames(imputed_olink_target)[colnames(imputed_olink_target) == 'OID01492'] = 'Q11128'
# Extract 1 protein due to its missing Uniprot (cannot find the correct unique Uniprot)
removed_column = imputed_olink_target[['OID01214']]
imputed_olink_target = imputed_olink_target[ , !(colnames(imputed_olink_target) %in% 'OID01214')]


# Metadata comparison
metadata_ht_subset = metadata_ht %>% select(sampleID)
metadata_target_subset = metadata_target %>% select(sampleID)

metadata_shared = inner_join(metadata_ht_subset, metadata_target_subset)
metadata_ht_unique = anti_join(metadata_ht_subset, metadata_target_subset)
metadata_target_unique = anti_join(metadata_target_subset, metadata_ht_subset)

cat('\nThere are', nrow(metadata_shared), 'shared data between Olink HT and Olink Target data sets.\n')
cat('There are', nrow(metadata_ht_unique), 'unique data in Olink HT dataset.\n')
cat('There are', nrow(metadata_target_unique), 'unique data in Olink Target dataset.\n\n')

# Protein detected comparison
protein_shared = intersect(colnames(imputed_olink_ht), colnames(imputed_olink_target))
protein_ht_unique = setdiff(colnames(imputed_olink_ht), colnames(imputed_olink_target))
protein_target_unique = setdiff(colnames(imputed_olink_target), colnames(imputed_olink_ht))

cat('\nNumber of proteins detected in both datasets:', length(protein_shared), '\n')
cat('Number of proteins only detected in Olink HT:', length(protein_ht_unique), '\n')
cat('Number of proteins only detected in Olink Target:', length(protein_target_unique), '\n\n')


# Protein correlation
# Insert barcode to map
imputed_olink_ht$barcode = npx_ht_wide$barcode
imputed_olink_ht$sampleID = metadata_ht$sampleID[match(imputed_olink_ht$barcode, metadata_ht$barcode)]
imputed_olink_ht = imputed_olink_ht %>% select(-barcode)

imputed_olink_target$sampleID = metadata_target$sampleID

# Filter for shared sampleID and proteins
imputed_olink_ht_filtered = imputed_olink_ht %>% filter(sampleID %in% metadata_shared$sampleID) %>%
  select(sampleID, all_of(protein_shared))
imputed_olink_target_filtered = imputed_olink_target %>% filter(sampleID %in% metadata_shared$sampleID) %>%
  select(sampleID, all_of(protein_shared))

# Convert to long format to merge
imputed_olink_ht_long = imputed_olink_ht_filtered %>%
  pivot_longer(cols = -sampleID, names_to = 'protein_name', values_to = 'npx_ht')
imputed_olink_target_long = imputed_olink_target_filtered %>%
  pivot_longer(cols = -sampleID, names_to = 'protein_name', values_to = 'npx_target')

# Merge all data
correlation = inner_join(imputed_olink_ht_long, imputed_olink_target_long, by = c('sampleID', 'protein_name'))

# Correlation 
# Protein wise 
correlation = correlation %>% group_by(protein_name) %>%
  mutate(correlation = cor(npx_ht, npx_target)) %>% ungroup()

# Choose protein P41439 with highest correlation
protein_data = correlation[correlation$protein_name == 'P41439', ]
# Pearson correlation coefficient
corr_coeff = cor(protein_data$npx_ht, protein_data$npx_target)
# Plot parameters
par(mar = c(4, 4, 3, 1), mgp = c(1.5, 0.3, 0), tcl = -0.2, las = 1, cex.axis = 0.8, cex.lab = 0.9)           
# Scatter plot
plot(protein_data$npx_ht, protein_data$npx_target,
     main = bquote(bold('Scatter Plot for Protein P41439')),
     xlab = 'Olink HT', ylab = 'Olink Target',
     pch = 19, col = rgb(70, 130, 180, max = 255, alpha = 100), cex = 1.2, asp = 1, bty = 'l') 
# Trendline
abline(lm(npx_target ~ npx_ht, data = protein_data), col = 'black', lwd = 2, lty = 2)
# Corr coeff
mtext(bquote(italic('Correlation Coefficient (r): ') ~ .(round(corr_coeff, 2))), side = 3, line = -0.1, cex = 0.8) 
# Histogram of corr distribution
# Correlation per protein
unique_correlations = correlation %>% distinct(protein_name, correlation)
# Histogram
hist(unique_correlations$correlation, breaks = 20, col = 'steelblue', border = 'white',
     main = 'Distribution of protein correlation', xlab = 'Correlation', ylab = 'Frequency')

# Merge data 
imputed_olink_ht_long = imputed_olink_ht_long %>% 
  left_join(select(metadata_ht, sampleID, visit, Sex), by = 'sampleID')

# Anova
anova_result = aov(npx_ht ~ Sex + visit, data = imputed_olink_ht_long)
summary(anova_result)

# ETA squared
eta_squared_result = eta_squared(anova_result)








