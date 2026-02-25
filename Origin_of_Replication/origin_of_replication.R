# Load the libraries:
library(readxl)
library(dplyr)

# Set Working Directory:
setwd("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication")

###mev
## 1- Read data from the files into a data frame:
bed <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                  sheet = "bed")
mev <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                         sheet = "mev")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- mev[mev$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "mev_result.csv", row.names = FALSE)
write.csv(negative_differences, "mev_negative_differences_result.csv", row.names = FALSE)

###N2
## Read data from the files into a data frame:
N2 <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                 sheet = "N2")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- N2[N2$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame
# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "N2_result.csv", row.names = FALSE)
write.csv(negative_differences, "N2_negative_differences_result.csv", row.names = FALSE)

###PB306
## 1- Read data from the files into a data frame:
PB306 <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                    sheet = "PB306")


# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- PB306[PB306$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "PB306_result.csv", row.names = FALSE)
write.csv(negative_differences, "PB306_negative_differences_result.csv", row.names = FALSE)

###C.snp
## 1- Read data from the files into a data frame:
C.snp <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                    sheet = "C.snp")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- C.snp[C.snp$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "C.snp_result.csv", row.names = FALSE)
write.csv(negative_differences, "C.snp_negative_differences_result.csv", row.names = FALSE)

###E.snp
## 1- Read data from the files into a data frame:
E.snp <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                    sheet = "E.snp")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- E.snp[E.snp$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "E.snp_result.csv", row.names = FALSE)
write.csv(negative_differences, "E.snp_negative_differences_result.csv", row.names = FALSE)

###G.snp
## 1- Read data from the files into a data frame:
G.snp <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                    sheet = "G.snp")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- G.snp[G.snp$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "G.snp_result.csv", row.names = FALSE)
write.csv(negative_differences, "G.snp_negative_differences_result.csv", row.names = FALSE)

###J.snp
## 1- Read data from the files into a data frame:
J.snp <- read_excel("C:/Moein/projects/Mutability/Hipergator/Origin_of_Replication/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1.xlsx", 
                    sheet = "J.snp")

# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
  # Subset the data for the current "chr" value
  subset_file1 <- J.snp[J.snp$chr == chr, ]
  subset_file2 <- bed[bed$chr == chr, ]
  
  # Initialize an empty data frame to store the results for this "chr" value
  result_df <- data.frame()
  
  # Nested loops to calculate differences
  for (i in 1:nrow(subset_file1)) {
    for (j in 1:nrow(subset_file2)) {
      difference <- subset_file2$start[j] - subset_file1$POS[i]
      # Create a data frame with the difference and the original columns
      result_row <- data.frame(chr = chr, POS = subset_file1$POS[i], start = subset_file2$start[j], difference = difference)
      # Append the result row to the result data frame
      result_df <- rbind(result_df, result_row)
    }
  }
  
  # Append the result data frame to the list
  results_list[[chr]] <- result_df
}

# Combine all the results into a single data frame
final_result <- do.call(rbind, results_list)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed
write.csv(max_negative_differences, "J.snp_result.csv", row.names = FALSE)
write.csv(negative_differences, "J.snp_negative_differences_result.csv", row.names = FALSE)


#End
