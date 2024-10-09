df<-read.delim2("/Users/tuomastoivainen/Downloads/Table.S3.history.processed_runs.txt", head=T, sep="\t",dec=",")

# Clean the PATTERN column: replace spaces, remove extra spaces, and standardize case
df$PATTERN <- gsub("\\s+", "_", df$PATTERN)  # Replace any spaces with underscores
df$PATTERN <- trimws(df$PATTERN)             # Remove leading/trailing spaces
df$PATTERN <- toupper(df$PATTERN) 

df <- df[!is.na(df$PATTERN2), ]

# Step 1: Add the "Core_Status" column based on the predefined core regions
core_regions <- c("Croatia-Romania", "Lithuania-Croatia", "Lithuania-Romania")

# Assign "Core" status to the pairs that are in the core regions and "Non-Core" to others
df$Core_Status <- ifelse(df$pair %in% core_regions, "Core", "Non-Core")

# Step 2: Combine CORE_1 and CORE_2 as "CORE", and the four peripheral patterns as "PERIPHERAL"
df$CORE_PERIPHERAL <- df$PATTERN
df$CORE_PERIPHERAL[df$PATTERN %in% c("CORE_1", "CORE_2")] <- "CORE"
df$CORE_PERIPHERAL[df$PATTERN %in% c("PERIPHERAL_1", "PERIPHERAL_2", "PERIPHERAL_3", "PERIPHERAL_4")] <- "PERIPHERAL"

# Step 3: Exclude "FAILED" and "ND" patterns from the data
df_filtered_clean <- df[!df$PATTERN %in% c("FAILED", "ND"), ]

# Step 4: Check if the combined patterns exist in the Core and Non-Core groups
core_peripheral_counts <- table(df_filtered_clean$CORE_PERIPHERAL, df_filtered_clean$Core_Status)
print("Combined CORE and PERIPHERAL counts by Core_Status (excluding FAILED and ND):")
print(core_peripheral_counts)

# Step 5: If valid data exists, perform the chi-square or Fisher's exact test
if (sum(core_peripheral_counts) == 0) {
  print("No data available for CORE and PERIPHERAL patterns after combining.")
} else {
  # Perform the test if data exists
  if (min(core_peripheral_counts) < 5) {
    # Use Fisher's Exact Test if any cell has fewer than 5 counts
    core_vs_noncore_test <- fisher.test(core_peripheral_counts)
  } else {
    # Otherwise, use a chi-square test
    core_vs_noncore_test <- chisq.test(core_peripheral_counts)
  }
  
  # Step 6: Print the test result
  print(core_vs_noncore_test)
}
