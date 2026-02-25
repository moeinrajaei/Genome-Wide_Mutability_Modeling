library(readxl)



# Set Working Directory:

setwd("C:/Moein/projects/Mutability/New_9_13_2024/Origin_of_Replication")


I have two files, each has two columns of data, file1 has "chr" and "POS" columns and file2 has "chr" and "start" columns. 
how to Calculate the differences between each value in "start" and each value in "POS" using nested loops and then have the result files as a a data frame with all the input columns?


###nonmutated

### 1- Read data from the files into a data frame:

bed <- read_excel("C:/Users/mr2683/Desktop/mutability/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1 (003).xlsx", 
    sheet = "bed")
nonmutated <- read_excel("C:/Users/mr2683/Desktop/mutability/Copy of Pourkarimi...Whitehouse.eLIFE.2016_SUPPL Table 1 (003).xlsx", 
    sheet = "nonmutated")


# Create an empty list to store the results for each combination
results_list <- list()

# Iterate through unique values in the "chr" column
unique_chrs <- unique(bed$chr)  # You can use file1$chr as well if both files have the same "chr" values
for (chr in unique_chrs) {
    # Subset the data for the current "chr" value
    subset_file1 <- nonmutated[nonmutated$chr == chr, ]
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



# Load the dplyr library if not already loaded
library(dplyr)

# Assuming you have the final_result data frame

# Filter the data to keep only rows with negative differences
negative_differences <- final_result %>%
  filter(difference < 0)

# Group by the "POS" column and summarize to find the maximum negative difference
max_negative_differences <- negative_differences %>%
  group_by(POS) %>%
  summarize(Max_Negative_Difference = max(difference))

# Print or save max_negative_differences as needed

write.csv(max_negative_differences, "nonmutated_result.csv", row.names = FALSE)

write.csv(negative_differences, "nonmutated_negative_differences_result.csv", row.names = FALSE)

















