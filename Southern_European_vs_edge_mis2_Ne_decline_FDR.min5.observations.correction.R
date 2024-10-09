# Load necessary libraries
library(dplyr)
library(readr)

# Set the path for the log file
log_file_path <- "/Users/tuomastoivainen/Downloads/msmc-im-home/debug_log.txt"

# Create a file connection for message logging
log_connection <- file(log_file_path, open = "wt")

# Start logging standard output to the file
sink(log_file_path, append = TRUE)

# Start logging messages to the same file connection
sink(log_connection, type = "message")

# Function to validate if a file should be included based on im_N1 values across all time points
validate_file <- function(file_path, upper_limit=200000) {
  # Suppressing the read_tsv() indexing message
  df <- suppressMessages(read_tsv(file_path, show_col_types = FALSE))
  
  # Check if required columns exist
  if (!all(c("left_time_boundary", "im_N1") %in% colnames(df))) {
    message(paste("Skipping file (missing columns):", file_path))
    return(FALSE)  # Skip if required columns are missing
  }
  
  # Calculate the max im_N1 value across **all time points** in the file
  max_im_N1_all_time <- max(df$im_N1, na.rm = TRUE)
  message(paste("Max im_N1 for file", file_path, "across all time points is:", max_im_N1_all_time))
  
  # Check if the maximum im_N1 exceeds the upper limit across all time points
  if (max_im_N1_all_time > upper_limit) {
    message(paste("File excluded (max im_N1 exceeds threshold across all time points):", file_path))
    return(FALSE)  # Return FALSE if the file exceeds the threshold across any time point
  }
  
  return(TRUE)  # Return TRUE if the file passes the check
}

# Function to process valid files and extract im_N1 values during MIS2
process_all_files <- function(directory, min_generations=7000, max_generations=14500, upper_limit=200000) {
  file_paths <- list.files(directory, pattern = "estimates", full.names = TRUE)
  
  # Validate files and keep only valid ones based on the max im_N1 across all time points
  valid_files <- file_paths[vapply(file_paths, validate_file, logical(1), upper_limit)]
  
  # If no valid files, return empty result
  if (length(valid_files) == 0) {
    message("No valid files found in directory:", directory)
    return(list(im_N1_values = numeric(), file_count = 0))
  }
  
  used_file_count <- length(valid_files)  # Count of valid files
  all_im_N1_values <- numeric()  # Initialize to collect all im_N1 values
  
  # Process valid files during the MIS2 time range
  for (file_path in valid_files) {
    df <- suppressMessages(read_tsv(file_path, show_col_types = FALSE))
    
    # Filter for rows within the MIS2 time range
    im_N1_values <- df %>%
      filter(left_time_boundary >= min_generations & left_time_boundary <= max_generations) %>%
      pull(im_N1)
    
    all_im_N1_values <- c(all_im_N1_values, im_N1_values)  # Collect values for analysis
    message(paste("File included in the analysis:", file_path))
  }
  
  return(list(im_N1_values = all_im_N1_values, file_count = used_file_count))
}

# Main analysis function
run_analysis <- function(core_directories, main_directory, output_path, upper_limit=200000) {
  # Process core population data from multiple directories
  core_results <- lapply(core_directories, process_all_files, upper_limit = upper_limit)
  im_N1_values_core <- unlist(lapply(core_results, function(res) res$im_N1_values))
  mean_core <- mean(im_N1_values_core, na.rm = TRUE)
  
  # Get all edge folders and initialize the results table
  edge_folders <- list.dirs(main_directory, recursive = FALSE)
  results_table <- data.frame(Edge_Folder = character(), Test_Type = character(), P_Value = numeric(), FDR_P_Value = numeric(), Percentage_Difference = numeric(), File_Count = integer(), stringsAsFactors = FALSE)
  
  all_p_values <- c()
  
  for (edge_folder in edge_folders) {
    # Skip core directories
    if (edge_folder %in% core_directories) next
    
    # Process edge population data
    edge_result <- process_all_files(edge_folder, upper_limit = upper_limit)
    im_N1_values_edge <- edge_result$im_N1_values
    file_count <- edge_result$file_count
    
    # Exclude edge folders with fewer than 4 estimate files used
    if (file_count < 5) {
      message(paste("Skipping folder due to fewer than 4 files:", edge_folder))
      next
    }
    
    if (length(im_N1_values_edge) == 0) {
      results_table <- rbind(results_table, data.frame(Edge_Folder = basename(edge_folder), Test_Type = "No data", P_Value = NA, FDR_P_Value = NA, Percentage_Difference = NA, File_Count = file_count))
      next
    }
    
    # Perform the comparison test
    comparison_result <- perform_comparison_test(im_N1_values_edge, im_N1_values_core)
    p_value <- ifelse(is.null(comparison_result$result), NA, comparison_result$result$p.value)
    
    # Calculate percentage difference
    mean_edge <- mean(im_N1_values_edge, na.rm = TRUE)
    percentage_difference <- ifelse(!is.na(mean_edge) && mean_core != 0, ((mean_core - mean_edge) / mean_core) * 100, NA)
    
    # Collect p-values for FDR adjustment
    if (!is.na(p_value)) all_p_values <- c(all_p_values, p_value)
    
    # Add result to the table
    results_table <- rbind(results_table, data.frame(Edge_Folder = basename(edge_folder), Test_Type = comparison_result$test_type, P_Value = p_value, FDR_P_Value = NA, Percentage_Difference = percentage_difference, File_Count = file_count))
  }
  
  # Apply FDR correction using the Benjamini-Hochberg method
  results_table$FDR_P_Value[!is.na(results_table$P_Value)] <- p.adjust(all_p_values, method = "BH")
  
  # Save the results as a CSV file
  write.csv(results_table, output_path, row.names = FALSE)
  print(results_table)
}

# Run the analysis
core_directories <- c("/Users/tuomastoivainen/Downloads/msmc-im-home/croatia_romania", "/Users/tuomastoivainen/Downloads/msmc-im-home/italy_croatia", "/Users/tuomastoivainen/Downloads/msmc-im-home/germany_italy")
main_directory <- "/Users/tuomastoivainen/Downloads/msmc-im-home/"
output_path <- "/Users/tuomastoivainen/Downloads/msmc-im-home/edge_vs_core_comparison_results_fdr.min5.csv"

run_analysis(core_directories, main_directory, output_path)

# Stop logging and close the log file
sink()
sink(type = "message")
close(log_connection)

# Output message to indicate logging is complete
message("Logging complete. Check the log file:", log_file_path)
