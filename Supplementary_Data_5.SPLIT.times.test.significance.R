# Load necessary libraries
library(dplyr)
library(readr)
library(ggpubr)
library(combinat)

# Load the data
data <- read.delim2("/Users/tuomastoivainen/Documents/Dataset_S04.txt", header = TRUE, sep = "\t", dec = ",")

# Remove problematic samples
data2 <- subset(data, sample1 != "ES21")
data3 <- subset(data2, sample1 != "ES3")

# Convert split_time_ya to numeric
data3 <- data3 %>%
  mutate(split_time_ya = as.numeric(gsub(",", ".", split_time_ya)))

# Get all unique pairs
unique_pairs <- unique(data3$pair)

# Generate all unique pairwise combinations
pair_combinations <- t(combn(unique_pairs, 2))

# Initialize list to store test results
results <- list()

# Loop over each pairwise comparison
for (i in 1:nrow(pair_combinations)) {
  pair1 <- pair_combinations[i, 1]
  pair2 <- pair_combinations[i, 2]
  
  # Extract split times for both region pairs
  group1 <- data3 %>% filter(pair == pair1) %>% pull(split_time_ya) %>% na.omit()
  group2 <- data3 %>% filter(pair == pair2) %>% pull(split_time_ya) %>% na.omit()
  
  # Perform test only if both groups have ≥3 values
  if (length(group1) >= 3 & length(group2) >= 3) {
    
    # Normality tests (only apply if n between 3 and 5000)
    normality_group1 <- if (length(group1) <= 5000) shapiro.test(group1)$p.value else NA
    normality_group2 <- if (length(group2) <= 5000) shapiro.test(group2)$p.value else NA
    
    # Choose test based on normality
    if (!is.na(normality_group1) & normality_group1 > 0.05 & 
        !is.na(normality_group2) & normality_group2 > 0.05) {
      test_result <- t.test(group1, group2, var.equal = TRUE)
      test_type <- "t-test"
    } else {
      test_result <- wilcox.test(group1, group2)
      test_type <- "Wilcoxon"
    }
    
    # Store result
    results[[i]] <- list(
      pair1 = pair1,
      pair2 = pair2,
      test_type = test_type,
      p_value = test_result$p.value
    )
  }
}

# Remove any NULL results
results <- results[!sapply(results, is.null)]

# Convert to data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Convert p_value column to numeric
results_df$p_value <- as.numeric(results_df$p_value)

# Apply multiple testing corrections
results_df$p_value_FDR <- p.adjust(results_df$p_value, method = "BH")
results_df$p_value_Bonferroni <- p.adjust(results_df$p_value, method = "bonferroni")

# View significant comparisons (FDR < 0.05)
sig_fdr <- results_df %>% filter(p_value_FDR < 0.05)
print(sig_fdr)

# Save full results
write.table(results_df, "pair_within-regions_comparisons_results_split_time_all_data_multiple_corrected.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
