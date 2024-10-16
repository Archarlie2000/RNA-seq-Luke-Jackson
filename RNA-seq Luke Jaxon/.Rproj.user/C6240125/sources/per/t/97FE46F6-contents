# Load necessary libraries
library(biomaRt)

# Set up biomaRt to use Ensembl for rat and human
rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert Ensembl transcript IDs to Ensembl gene IDs
convert_to_gene_ids <- function(transcript_ids) {
  gene_mapping <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), 
                        filters = "ensembl_transcript_id", 
                        values = transcript_ids, 
                        mart = rat)
  return(gene_mapping)
}

# Set the directory containing the TSV files
directory <- "../DEG_results"

# Read all TSV files and combine them into one data frame
files <- list.files(directory, pattern = "*.tsv", full.names = TRUE)

# Initialize an empty data frame
deg_data <- data.frame()

# Loop through each file, read its contents, and add it to the main data frame
for (file in files) {
  # Read the TSV file
  temp_data <- read.delim(file, header = FALSE)
  
  # Set column names for the data frame
  colnames(temp_data) <- c("Transcript_ID", "logFC", "logCPM", "LR", "PValue", "FDR")
  
  # Add a column to indicate the source file
  temp_data$source <- basename(file)
  
  # Combine the data into the main data frame
  deg_data <- rbind(deg_data, temp_data)
}

# Extract unique Ensembl transcript IDs
transcript_ids <- unique(deg_data$Transcript_ID)

# Convert transcript IDs to gene IDs
gene_mapping <- convert_to_gene_ids(transcript_ids)

# Merge the gene mapping with the main data frame
deg_data <- merge(deg_data, gene_mapping, by.x = "Transcript_ID", by.y = "ensembl_transcript_id")

# Extract unique Ensembl gene IDs
vec_ids <- unique(deg_data$ensembl_gene_id)

# Find orthologs of rat Ensembl gene IDs in human
orthologs <- getLDS(attributes = c("ensembl_gene_id"),
                    filters = "ensembl_gene_id",
                    values = vec_ids,
                    mart = rat,
                    attributesL = c("ensembl_gene_id"),
                    martL = human)

# Rename columns for clarity
colnames(orthologs) <- c("ensembl_gene_id", "human_ensembl_gene_id")

# Merge orthologs with the main data frame
deg_data <- merge(deg_data, orthologs, by = "ensembl_gene_id", all.x = TRUE)

# Create a new directory for the updated files
output_directory <- "../Human_ortholog"
dir.create(output_directory, showWarnings = FALSE)

# Save each file with the added human ortholog column
for (file in files) {
  # Filter the data for the current file
  file_data <- deg_data[deg_data$source == basename(file), ]
  
  # Remove the source column for saving
  file_data$source <- NULL
  
  # Write the updated data to a CSV file in the new directory
  output_file <- file.path(output_directory, paste0(tools::file_path_sans_ext(basename(file)), ".csv"))
  write.csv(file_data, output_file, row.names = FALSE)
}

# Print a message indicating completion
print("All files have been saved with human orthologs in the 'Human_ortholog' folder.")