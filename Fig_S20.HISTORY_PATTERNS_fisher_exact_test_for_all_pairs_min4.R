# Step 1: Calculate the percentage of "FAILED" and "ND" rows for each region
failed_nd_counts <- table(df_filtered$pair, df_filtered$CORE_OR_PERIPHERAL %in% c("FAILED", "ND"))
total_counts <- rowSums(failed_nd_counts)

# Calculate the percentage of "FAILED" + "ND"
failed_nd_percentages <- (failed_nd_counts[, "TRUE"] / total_counts) * 100

# Step 2: Identify regions where more than 50% of the rows are "FAILED" or "ND"
regions_to_exclude <- names(failed_nd_percentages[failed_nd_percentages > 50])

# Step 3: Filter out these regions from the data
df_filtered_clean <- df_filtered[!df_filtered$pair %in% regions_to_exclude, ]

# Step 4: Perform the core vs non-core analysis as before

# Define the core regions as before
core_regions <- c("Croatia-Romania", "Lithuania-Croatia", "Lithuania-Romania")

# Function to perform pairwise chi-square or Fisher's test and collect results
perform_pairwise_chi_square_or_fisher <- function(data, region1, region2) {
  patterns <- unique(data$PATTERN)
  patterns <- patterns[!patterns %in% c("FAILED", "ND")]  # Exclude "FAILED" and "ND"
  
  data_region1 <- data[data$pair == region1, ]
  data_region2 <- data[data$pair == region2, ]
  
  ew_region1 <- unique(data_region1$E_W)
  ew_region2 <- unique(data_region2$E_W)
  
  data_combined <- rbind(data_region1, data_region2)
  
  core_status <- ifelse(region1 %in% core_regions | region2 %in% core_regions, "Core Comparison", "Non-Core Comparison")
  
  results <- data.frame(
    Pattern = character(),
    Region_Comparison = character(),
    E_W_Comparison = character(),
    Core_Comparison = character(),
    Test = character(),
    Statistic = numeric(),
    P_Value = numeric(),
    Note = character(),
    stringsAsFactors = FALSE
  )
  
  for (pattern in patterns) {
    pattern_count_region1 <- sum(data_region1$PATTERN == pattern)
    pattern_count_region2 <- sum(data_region2$PATTERN == pattern)
    
    if (pattern_count_region1 == 0 & pattern_count_region2 == 0) {
      next
    }
    
    if (pattern_count_region1 == 0 | pattern_count_region2 == 0) {
      contingency_table <- table(data_combined$PATTERN == pattern, data_combined$pair)
      fisher_test <- fisher.test(contingency_table)
      
      results <- rbind(results, data.frame(
        Pattern = pattern,
        Region_Comparison = paste(region1, "vs", region2),
        E_W_Comparison = paste(ew_region1, "vs", ew_region2),
        Core_Comparison = core_status,
        Test = "Fisher's Exact Test",
        Statistic = fisher_test$estimate,
        P_Value = fisher_test$p.value,
        Note = "Fisher's test used due to zero observations in one region"
      ))
    } else {
      contingency_table <- table(data_combined$PATTERN == pattern, data_combined$pair)
      chi_square_test <- chisq.test(contingency_table)
      
      results <- rbind(results, data.frame(
        Pattern = pattern,
        Region_Comparison = paste(region1, "vs", region2),
        E_W_Comparison = paste(ew_region1, "vs", ew_region2),
        Core_Comparison = core_status,
        Test = "Chi-Square Test",
        Statistic = chi_square_test$statistic,
        P_Value = chi_square_test$p.value,
        Note = ""
      ))
    }
  }
  
  return(results)
}

# Step 5: Identify regions with enough observations (at least 4 observations)
region_counts_clean <- table(df_filtered_clean$pair)
regions_with_enough_data_clean <- names(region_counts_clean[region_counts_clean >= 4])

# Step 6: Perform pairwise tests (Chi-square or Fisher's) and collect results
all_results_clean <- data.frame()

for (i in 1:(length(regions_with_enough_data_clean) - 1)) {
  for (j in (i + 1):length(regions_with_enough_data_clean)) {
    region1 <- regions_with_enough_data_clean[i]
    region2 <- regions_with_enough_data_clean[j]
    
    pairwise_results <- perform_pairwise_chi_square_or_fisher(df_filtered_clean, region1, region2)
    
    all_results_clean <- rbind(all_results_clean, pairwise_results)
  }
}

# Step 7: Organize and print the results
all_results_clean <- all_results_clean[order(all_results_clean$Pattern), ]
print(all_results_clean)

# Optionally, save the organized results to a CSV file
write.csv(all_results_clean, "pairwise_chi_square_or_fisher_results_filtered.csv", row.names = FALSE)
