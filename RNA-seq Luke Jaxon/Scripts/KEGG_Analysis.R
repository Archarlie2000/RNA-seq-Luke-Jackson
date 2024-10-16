# Load required libraries
library(clusterProfiler)
library(org.Rn.eg.db)       # Rat-specific database
library(pathview)           # For KEGG pathway visualization
library(biomaRt)            # For ID conversion

# Directory containing the map files
input_directory <- "../DEG_results"
# Directory to save the comparison results
output_directory <- "../KEGG_Results"

# List all files in the input directory
map_files <- list.files(input_directory, full.names = TRUE)

# Initialize a data frame to store all KEGG enrichment results
all_results <- data.frame()

# Set up biomaRt to use Ensembl
ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")


# file <- "../DEG_results/DEG_results_ALL_Cell_11mM_vs_Control.tsv"
# Iterate through each map file and perform the analysis
for (file in map_files) {
  # Extract the base file name for saving purposes
  file_name <- basename(file)
  
  # Load the mapped data for KEGG analysis
  mapped_data <- read.delim(file, sep = "\t", header = TRUE)
  
  result_prefix <- sub("\\.csv$", "", file_name)  # Remove .csv extension
  
  # Create a subfolder for each comparison result
  comparison_directory <- file.path(output_directory, result_prefix)
  dir.create(comparison_directory, showWarnings = FALSE)
  
  # Convert row names to a column named "ensembl_transcript_id"
  mapped_data$ensembl_transcript_id <- rownames(mapped_data)
  rownames(mapped_data) <- NULL
  
  # Reorder columns to put "ensembl_transcript_id" first
  mapped_data <- mapped_data[, c("ensembl_transcript_id", colnames(mapped_data)[-ncol(mapped_data)])]
  
  # Extract the Ensembl Transcript IDs
  ensembl_ids <- mapped_data$ensembl_transcript_id
  
  # Convert Ensembl Transcript IDs to Entrez Gene IDs using biomaRt
  conversion_results <- getBM(attributes = c('ensembl_transcript_id', 'entrezgene_id'),
                              filters = 'ensembl_transcript_id',
                              values = ensembl_ids,
                              mart = ensembl)
  
  # Merge conversion results with the original data to keep only mapped genes
  mapped_data <- merge(mapped_data, conversion_results, by = "ensembl_transcript_id")
  
  # Extract the gene list for KEGG analysis (Entrez Gene IDs)
  gene_list <- mapped_data$entrezgene_id
  
  # Perform KEGG analysis for rats ('rno')
  keggr <- enrichKEGG(gene = gene_list,
                      organism = "rno",
                      keyType = "ncbi-geneid",
                      pvalueCutoff = 0.05)
  
  # Convert KEGG results to a readable format
  if (!is.null(keggr)) {
    keggr <- setReadable(keggr, OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
    
    # Save the KEGG pathway gene enrichment analysis
    result_prefix <- sub("\\.csv$", "", file_name)  # Remove .csv extension
    result_file <- file.path(output_directory, paste0(result_prefix, "_kegg_enrichment_results.csv"))
    write.csv(as.data.frame(keggr), result_file)
    
    # Append results to the comprehensive results data frame
    enriched_df <- as.data.frame(keggr)
    enriched_df$Comparison <- result_prefix  # Add a column for identifying the comparison
    all_results <- rbind(all_results, enriched_df)
  }
}

# Save the comprehensive results to a CSV file
write.csv(all_results, file.path(output_directory, "comprehensive_kegg_enrichment_results.csv"))