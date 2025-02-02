# Mai Nguyen 11/11/2024
# S3WP analysis (Olink HT and Olink Target)
# QC code 

# Import require lib
library(dplyr)
library(impute)
library(openxlsx)
library(ggplot2)
library(umap)
library(tidyr)
library(stringr)
library(effectsize)
library(FactoMineR)
library(factoextra)
library(lme4)

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

# Function to create scatter plot for a given protein
plot_correlation = function(protein_data, protein_name) {
  corr_coeff = cor(protein_data$npx_ht, protein_data$npx_target)
  par(mar = c(4, 4, 3, 1), mgp = c(1.5, 0.3, 0), tcl = -0.2, las = 1, cex.axis = 0.8, cex.lab = 0.9)           
  plot(protein_data$npx_ht, protein_data$npx_target,
       main = bquote(bold('Scatter Plot for Protein' ~ .(protein_name))),
       xlab = 'Olink HT', ylab = 'Olink Target',
       pch = 19, col = rgb(70, 130, 180, max = 255, alpha = 100), cex = 1.2, asp = 1, bty = 'l') 
  abline(lm(npx_target ~ npx_ht, data = protein_data), col = 'black', lwd = 2, lty = 2)
  mtext(bquote(italic('Correlation Coefficient (r): ') ~ .(round(corr_coeff, 2))), side = 3, line = -0.1, cex = 0.8)}

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

# Extract subject_id for meta data 
olink_ht_fil2 = extract_info_from_column(olink_ht_fil2, 'Extra.data', 'subject_id: \\S+', 'subject_id')
# Add sampleID to the column
olink_ht_fil2$sampleID = paste0('WEL_', olink_ht_fil2$visit, '_', olink_ht_fil2$subject_id)
olink_ht_fil2 = olink_ht_fil2 %>% mutate(subject_id = str_replace(subject_id, '.*?-', ''))

# Add gender to subjects
subject_info = subject_info %>% rename(subject_id = Subject)
subject_info = subject_info %>% mutate(subject_id = as.character(subject_id))
olink_ht_fil2 = olink_ht_fil2 %>% left_join(subject_info %>% select(subject_id, Sex), by = 'subject_id')

# Filter duplicated measurements
repeated_measurements = olink_ht_fil2 %>% group_by(subject_id, visit) %>%
  filter(n_distinct(DAid) > 1) %>%  ungroup() %>% distinct(subject_id, visit, DAid) %>%  
  arrange(subject_id, visit, DAid) 

# Remove first measurement of the duplicates (DAid == DA12126, DA12092, DA11971)
olink_ht_fil2 = olink_ht_fil2 %>%  filter(!(DAid == 'DA12126'| DAid == 'DA12092' | DAid == 'DA11971'))

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
olink_ht_fil2 = olink_ht_fil2 %>% mutate(under_LOD = ifelse(is.na(NPX) | NPX < LOD, 'Yes', 'No'))
olink_ht_fil2$under_LOD = as.factor(olink_ht_fil2$under_LOD)

# Number of protein compared with LOD
protein_detected_100 = analyze_proteins(olink_ht_fil2)
protein_detected_80 = analyze_proteins(olink_ht_fil2, detection_percentage = 80)

# Per protein, calculate the percentage of samples detected above LOD
protein_LOD_comparison = olink_ht_fil2 %>% group_by(UniProt) %>% summarise(
    more_than_LOD = sum(under_LOD == 'No'), total_samples = n(),
    percent_more_than_LOD = round((more_than_LOD / total_samples) * 100, 2)) %>%
  select(UniProt, more_than_LOD, percent_more_than_LOD) %>% arrange(desc(percent_more_than_LOD))

# Separate 0% into one bin
protein_LOD_comparison = protein_LOD_comparison %>%
  mutate(percent_more_than_LOD = ifelse(percent_more_than_LOD == 0, -1, percent_more_than_LOD))
# Maximum count to have adaptive y axis
max_count = max(hist(ifelse(protein_LOD_comparison$percent_more_than_LOD == 0, -1, 
                            protein_LOD_comparison$percent_more_than_LOD), 
                     breaks = seq(-5, 105, by = 5), plot = FALSE)$counts)
# Above LOD percentage histogram
ggplot(protein_LOD_comparison, aes(x = ifelse(percent_more_than_LOD == 0, -1, percent_more_than_LOD))) +
  geom_histogram(binwidth = 5, fill = 'skyblue', color = 'black', alpha = 0.7, boundary = 0, closed = 'left') +
  geom_segment(aes(x = 0, xend = 0, y = 0, yend = max_count), color = 'red', linetype = 'dashed', size = 0.6, inherit.aes = FALSE) +
  geom_segment(aes(x = 50, xend = 50, y = 0, yend = max_count), color = 'blue', linetype = 'dashed', size = 0.6, inherit.aes = FALSE) +
  labs(title = 'Distribution of percentage above LOD', x = 'Proteins detected rate compared to LOD', y = 'Frequency') +
  theme_minimal() + theme(
    plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)),
    axis.title.x = element_text(size = 10, margin = margin(t = 10)), 
    axis.title.y = element_text(size = 10, margin = margin(r = 10)), 
    axis.text.x = element_text(size = 9), 
    axis.text.y = element_text(size = 9),
    axis.line = element_line(size = 0.5, color = 'black'), 
    panel.grid = element_blank()) + scale_x_continuous(
    breaks = c(-1, seq(5, 100, by = 10)), 
    labels = c('0%', seq(10, 100, by = 10))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + coord_cartesian(xlim = c(-5, 105)) 


# Plot QC 
# Divide data into meta data and count data 
metadata_ht = olink_ht %>% select(DAid, visit, barcode, subject_id, Sex, sampleID) %>% distinct() 
# Add age column into meta data 
subject_info$Birth.date = as.Date(subject_info$Birth.date)
metadata_ht$age = 2015 - as.numeric(format(subject_info$Birth.date[match(metadata_ht$subject_id, subject_info$subject_id)], '%Y'))

# count data go with barcode bcs there are some dup in patient-visit (sampleID)
ht_count = olink_ht %>% select(UniProt, NPX, barcode)
npx_ht_wide = ht_count %>% pivot_wider(names_from = UniProt, values_from = 'NPX')
npx_ht = npx_ht_wide %>% select(-barcode)
npx_ht_matrix = as.matrix(npx_ht)

# Run the code if data is not yet numeric
# npx_ht = npx_ht %>% mutate(across(everything(), ~ as.numeric(.)))

# Impute missing data
imputed_olink_ht_npx = impute.knn(npx_ht_matrix)
imputed_olink_ht_matrix = imputed_olink_ht_npx$data
imputed_olink_ht = as.data.frame(imputed_olink_ht_matrix)

# Plot UMAP
umapxy_ht = umap(imputed_olink_ht, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_ht = data.frame(x = umapxy_ht$layout[, 1], y = umapxy_ht$layout[, 2])
umap_ht$barcode = npx_ht_wide$barcode
umap_ht = umap_ht %>% left_join(metadata_ht %>% select(barcode, Sex, visit, subject_id), by = 'barcode')

# UMAP colored by sex
ggplot(umap_ht, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink HT colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
    axis.title = element_text(size = 14), legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

# Filter outliers from UMAP
# print(umap_ht[which(umap_ht$x < -30), ])
# print(umap_ht[which(umap_ht$y < -20), ])
# Outliers detected: samples belong to subject_id: 3283, 3217, 3119, 3293
# Remove outlier in all data: npx_ht_wide, imputed_olink_ht, metadata_ht
npx_ht_wide = npx_ht_wide[-which(umap_ht$x < -30 | umap_ht$y < -20), ]
imputed_olink_ht = imputed_olink_ht[-which(umap_ht$x < -30 | umap_ht$y < -20), ]
imputed_olink_ht_matrix = as.matrix(imputed_olink_ht)
metadata_ht = metadata_ht[-which(umap_ht$x < -30 | umap_ht$y < -20), ]

umapxy_ht = umap(imputed_olink_ht, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_ht = data.frame(x = umapxy_ht$layout[, 1], y = umapxy_ht$layout[, 2])
umap_ht$barcode = npx_ht_wide$barcode
umap_ht = umap_ht %>% left_join(metadata_ht %>% select(barcode, Sex, visit, subject_id), by = 'barcode')

# Re-run UMAP colored by sex
ggplot(umap_ht, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink HT colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title = element_text(size = 14), legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# UMAP colored by visit
ggplot(umap_ht, aes(x = x, y = y, color = factor(visit))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink HT colored by Visit',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  scale_color_manual(values = c('#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#e41a1c', '#ffff33')) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# UMAP colored by patient 
# No legend bcs there are too many
ggplot(umap_ht, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink HT colored by Patient',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# UMAP colored by patient separately - assess profile consistency over visits
ggplot(umap_ht, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink HT colored by Patient',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  facet_wrap(~subject_id) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# Plot PCA
olink_ht_pca = PCA(imputed_olink_ht_matrix, graph = FALSE)
fviz_pca_ind(olink_ht_pca, geom.ind = 'point', col.ind = recode(metadata_ht$Sex, 'f' = 'female', 'm' = 'male'), 
             legend.title = 'Gender', title = 'Olink HT PCA, colored by Gender')
# There was a significant outlier to check on
pca_coordinates = as.data.frame(olink_ht_pca$ind$coord)
print(metadata_ht[which(pca_coordinates$Dim.2 < -100), ])

# Meta data of outlier (might not be a human sample, like contamination)        
# DAid visit barcode subject_id Sex     sampleID age
# 72 DA12242     6  138004       3829   m WEL_6_1-3829  56

# Re-run PCA
imputed_olink_ht = imputed_olink_ht_matrix[-which(pca_coordinates$Dim.2 < -100), ]
metadata_ht = metadata_ht[-which(pca_coordinates$Dim.2 < -100), ]
npx_ht_wide = npx_ht_wide[-which(pca_coordinates$Dim.2 < -100), ]
olink_ht_pca_removed = PCA(imputed_olink_ht, graph = FALSE)
fviz_pca_ind(olink_ht_pca_removed, geom.ind = 'point', 
             col.ind = recode(metadata_ht$Sex, 'f' = 'female', 'm' = 'male'),
             legend.title = 'Gender', title = 'Olink HT PCA (removed outlier)')
imputed_olink_ht = data.frame(imputed_olink_ht)


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
imputed_olink_target_matrix = imputed_olink_target_npx$data
imputed_olink_target = as.data.frame(imputed_olink_target_matrix)

# Plot UMAP
umapxy_target = umap(imputed_olink_target, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_target = data.frame(x = umapxy_target$layout[, 1], y = umapxy_target$layout[, 2])
umap_target = cbind(umap_target, metadata_target)

# UMAP colored by sex
ggplot(umap_target, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink Target colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title = element_text(size = 14), legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Filter outliers from UMAP
print(umap_target[which(umap_target$y < -20), ])
# Outliers detected: samples belong to subject_id: 3293
# Remove outlier in all data: npx_target, imputed_olink_target, metadata_target
npx_target = npx_target[-which(umap_target$y < -20), ]
imputed_olink_target = imputed_olink_target[-which(umap_target$y < -20), ]
imputed_olink_target_matrix = as.matrix(imputed_olink_target)
metadata_target = metadata_target[-which(umap_target$y < -20), ]

umapxy_target = umap(imputed_olink_target, n_neighbors = 15, min_dist = 0.1, metric = 'euclidean')
umap_target = data.frame(x = umapxy_target$layout[, 1], y = umapxy_target$layout[, 2])
umap_target = cbind(umap_target, metadata_target)

# Re-run UMAP colored by sex
ggplot(umap_target, aes(x = x, y = y, color = Sex)) + geom_point(size = 3, alpha = 0.8) +
  labs(title = 'UMAP Olink Target colored by Sex', x = 'UMAP1', y = 'UMAP2', color = 'Sex') +
  scale_color_manual(values = c('f' = '#F8766D', 'm' = '#00BFC4')) + theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        axis.title = element_text(size = 14), legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

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

# UMAP colored by patient
ggplot(umap_target, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink Target colored by Patient ID',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  #facet_wrap(~visit) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# UMAP colored by patient separately - assess profile consistency over visits
ggplot(umap_target, aes(x = x, y = y, color = factor(subject_id))) +
  geom_point(size = 3, alpha = 0.8, show.legend = F) +
  labs(title = 'UMAP Olink Target colored by Patient ID',
       x = 'UMAP1', y = 'UMAP2', color = 'Visit') +
  facet_wrap(~subject_id) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
                          axis.title = element_text(size = 14), legend.title = element_text(size = 12),
                          legend.text = element_text(size = 10))

# Plot PCA
olink_target_pca = PCA(imputed_olink_target_matrix, graph = FALSE)
fviz_pca_ind(olink_target_pca, geom.ind = 'point', 
             col.ind = recode(metadata_target$Sex, 'f' = 'female', 'm' = 'male'), 
             legend.title = 'Gender', title = 'Olink Target PCA, colored by Gender')

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

# Correlation per protein
unique_correlations = correlation %>% distinct(protein_name, correlation) %>%
  mutate(corr_category = case_when(correlation <= 0.4 ~ 'Low', correlation > 0.4 & correlation <= 0.7 ~ 'Mid',
    correlation > 0.7 ~ 'High')) %>% arrange(correlation)

# Count proteins per corr level
category_counts = unique_correlations %>% group_by(corr_category) %>% summarize(count = n())
print(paste('Number of proteins with low correlation (0-0.4):', category_counts$count[category_counts$corr_category == 'Low']))
print(paste('Number of proteins with mid correlation (0.4-0.7):', category_counts$count[category_counts$corr_category == 'Mid']))
print(paste('Number of proteins with high correlation (0.7-1):', category_counts$count[category_counts$corr_category == 'High']))

# Histogram of corr distribution
hist(unique_correlations$correlation, breaks = 20, col = '#377eb8', border = 'white',
     main = 'Distribution of Protein Correlation', xlab = 'Correlation', ylab = 'Frequency')
# Line separated regions
abline(v = 0.4, col = 'red', lty = 2, lwd = 2) 
abline(v = 0.7, col = 'red', lty = 2, lwd = 2) 

# 5 most and least correlated proteins
lowest_proteins = unique_correlations %>% slice(1:5)
highest_proteins = unique_correlations %>% slice((n() - 4):n())

# Correlation plot
for (protein in lowest_proteins$protein_name) {
  protein_data = correlation[correlation$protein_name == protein, ]
  plot_correlation(protein_data, protein)}
for (protein in highest_proteins$protein_name) {
  protein_data = correlation[correlation$protein_name == protein, ]
  plot_correlation(protein_data, protein)}

# ANOVA
# Merge data 
imputed_olink_ht_long = imputed_olink_ht_long %>% 
  left_join(select(metadata_ht, sampleID, visit, Sex, age), by = 'sampleID')

# Anova
anova_result = aov(npx_ht ~ Sex + visit + age, data = imputed_olink_ht_long)
summary(anova_result)

# Summary: 
# Sex: anova p-value < 2e-16, eta = 2.9e-4 --> impactful but on a small part
# Age: anova p-value < 2e-16, eta = 4.6e-4 --> impactful but on a small part
# Visit: anova p-value = 0.9, eta = 5.9e-8 --> no impact

# Anova per protein (extract to columns of a df for easy read)
anova_result_protein = imputed_olink_ht_long %>% group_by(protein_name) %>%
  summarise(residual_df = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]]['Residuals', 'Df'],
    Sex_F = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['F value']][1],
    Sex_p = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['Pr(>F)']][1],
    visit_F = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['F value']][2],
    visit_p = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['Pr(>F)']][2],
    age_F = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['F value']][3],
    age_p = summary(aov(npx_ht ~ Sex + visit + age, data = cur_data()))[[1]][['Pr(>F)']][3])

# ETA squared
eta_squared_result = eta_squared(anova_result)
# Per protein
anova_result_protein = anova_result_protein %>% mutate(
    Sex_eta2 = (Sex_F * 1) / ((Sex_F * 1) + residual_df),
    visit_eta2 = (visit_F * 1) / ((visit_F * 1) + residual_df),
    age_eta2 = (age_F * 1) / ((age_F * 1) + residual_df))


# Longitudinal analysis with mixed-effect model 
imputed_olink_ht_long = imputed_olink_ht_long %>% 
  left_join(metadata_ht %>% select(sampleID, subject_id), by = 'sampleID')

# Change dtype to factors
imputed_olink_ht_long$visit = as.factor(imputed_olink_ht_long$visit)
imputed_olink_ht_long$Sex = as.factor(imputed_olink_ht_long$Sex)
imputed_olink_ht_long$subject_id = as.factor(imputed_olink_ht_long$subject_id)

# Model for general effect
model = lmer(npx_ht ~ age + Sex + visit + (1 | subject_id), data = imputed_olink_ht_long)
summary(model)

# Model for each protein (not sure what this one is)
# I completely don't know how to do this so I used ChatGPT :<
# Split data by protein_name
protein_models = imputed_olink_ht_long %>% group_by(protein_name) %>% group_split()

# Fit the original mixed-effects model (with age, Sex, and visit) for each protein
model_with_covariates = lapply(protein_models, function(data) {
  lmer(npx_ht ~ age + Sex + visit + (1 | subject_id), data = data)})

# Fit the simplified mixed-effects model (only visit) for each protein
model_only_visit = lapply(protein_models, function(data) {
  lmer(npx_ht ~ visit + (1 | subject_id), data = data)})

# Extract results from both models using broom.mixed::tidy()
library(broom.mixed)

results_with_covariates = lapply(model_with_covariates, tidy)
results_only_visit = lapply(model_only_visit, tidy)

# Combine results into dataframes for each model
results_with_covariates_df = do.call(rbind, results_with_covariates)
results_only_visit_df = do.call(rbind, results_only_visit)

# Add protein names to the results (ensure correct alignment)
results_with_covariates_df$protein_name = rep(
  sapply(protein_models, function(df) unique(df$protein_name)), 
  sapply(model_with_covariates, function(model) nrow(tidy(model))))

results_only_visit_df$protein_name = rep(
  sapply(protein_models, function(df) unique(df$protein_name)), 
  sapply(model_only_visit, function(model) nrow(tidy(model))))

















# Pipeline
# - general qc
# - platform comparison - how do assays in target and ht correlate?
# - umap - how stable is the plasma proteome across visits for target and ht?
# - ht - how many of the 5.400 proteins do we expect to detect above LOD in a healthy population? if time allows, biological differences
# - anova - how much of the variation in the plasma proteome is explained by age, sex and visit?




# Focus now: mixed-effect model for each protein --> pathway enrichment






# Later (if time allows)
# Anova plot: stacked bar plot for variance explained by factors --> Zoom in plot per factor 
# Npx against factors (violin plot or smth Wen article)
# Numb of proteins related to each factor (kinda a table?)
# Shared protein between HT and Target then plot PCA to see if similar
# Anova for target and compare protein to HT






# Done
# protein wise above LOD (percentage per protein) - done
# Edit hist of corr - done
# Do a divided umap per patient (facets) - done
# PCA for both data sets to identify patterns - done
# Anova per protein - done
