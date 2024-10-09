# Load necessary libraries

library(dplyr)
library(readr)
library(ggpubr)
library(combinat)

# Load the data
data <- read.delim2("C:/Users/ttoivain/Downloads/ttoivain/Downloads/table12.tsv", head=T, sep="\t",dec=",")
data2<-subset(data,data$sample!="ES21")
data3<-subset(data2,data2$sample!="ES3")



data3 <- data3 %>%
  mutate(split_time_ya2 = as.numeric(gsub(",", ".", split_time_ya2)))
# Get unique pairs
unique_pairs <- unique(data3$pair)

# Generate all unique pair combinations
pair_combinations <- t(combn(unique_pairs, 2))

# Initialize a list to store results
results <- list()

# Loop through each combination of pairs
for (i in 1:nrow(pair_combinations)) {
  pair1 <- pair_combinations[i, 1]
  pair2 <- pair_combinations[i, 2]
  
  # Extract split times for the two pairs
  group1 <- data3 %>% filter(pair == pair1) %>% pull(split_time_ya2) %>% na.omit()
  group2 <- data3 %>% filter(pair == pair2) %>% pull(split_time_ya2) %>% na.omit()
  
  # Check if both groups have at least 3 samples
  if (length(group1) >= 3 & length(group2) >= 3) {
    # Check normality only if sample size is between 3 and 5000
    normality_group1 <- if (length(group1) >= 3 & length(group1) <= 5000) shapiro.test(group1)$p.value else NA
    normality_group2 <- if (length(group2) >= 3 & length(group2) <= 5000) shapiro.test(group2)$p.value else NA
    
    # Determine if both groups are normally distributed
    if (!is.na(normality_group1) & normality_group1 > 0.05 & !is.na(normality_group2) & normality_group2 > 0.05) {
      test_result <- t.test(group1, group2, var.equal = TRUE)
      test_type <- "t-test"
    } else {
      test_result <- wilcox.test(group1, group2)
      test_type <- "Wilcoxon"
    }
    
    # Store the result
    results[[i]] <- list(
      pair1 = pair1,
      pair2 = pair2,
      test_type = test_type,
      p_value = test_result$p.value
    )
  }
}

# Remove NULL entries from the results list
results <- results[!sapply(results, is.null)]

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Print the results
print(results_df)
setwd("C:/Users/ttoivain/Downloads/ttoivain/Downloads/msmc-im-home/figures/Fig4")
# Save the results to a CSV file
write.table(results_df, "pair_within-regions_comparisons_results_split_time_all_data.txt", row.names = FALSE, col.names=T,quote=F, sep="\t")