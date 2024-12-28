# Mai Nguyen 11/11/2024
# S3WP analysis (Olink HT and Olink Target)
# Function file

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




