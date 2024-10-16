# Load required libraries
library(org.Rn.eg.db)       # Rat-specific database
library(clusterProfiler)    # For GO analysis
library(dplyr)              # For data manipulation
library(biomaRt)            # For converting transcript IDs to gene IDs

# Set the directory containing the DEG files and create the output directory
deg_directory <- "../DEG_Results"
output_directory <- "../GO_Results"
dir.create(output_directory, showWarnings = FALSE)

# Get the list of file names in the directory
files <- list.files(deg_directory, pattern = "*.tsv", full.names = TRUE)

# Load and combine DEG results into a list, treating the first column as row names (gene IDs)
deg_list <- lapply(files, function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", row.names = 1, fill = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  data$GeneID <- rownames(data)
  return(data)
})

# Set up biomaRt to use Ensembl for rat
ensembl <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")

# Function to convert Ensembl transcript IDs to Ensembl gene IDs
convert_to_gene_ids <- function(transcript_ids) {
  gene_mapping <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), 
                        filters = "ensembl_transcript_id", 
                        values = transcript_ids, 
                        mart = ensembl)
  return(gene_mapping)
}

# Function to map Ensembl Gene IDs to Entrez IDs using the rat database
map_to_entrez <- function(df, dataset_name) {
  transcript_ids <- df$GeneID
  gene_mapping <- convert_to_gene_ids(transcript_ids)
  unmapped_genes <- setdiff(transcript_ids, gene_mapping$ensembl_transcript_id)
  unmapped_genes_list[[dataset_name]] <<- unmapped_genes
  
  # Merge with the original data to get the Ensembl gene IDs
  df_mapped <- merge(df, gene_mapping, by.x = "GeneID", by.y = "ensembl_transcript_id")
  
  # Use bitr() to convert Ensembl gene IDs to Entrez IDs
  mapped_ids <- bitr(df_mapped$ensembl_gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Rn.eg.db)
  
  # Handle the case where some IDs fail to map
  if (nrow(mapped_ids) == 0) {
    warning(paste0("No gene IDs were mapped for dataset: ", dataset_name))
    return(NULL)
  }
  
  df_final <- merge(df_mapped, mapped_ids, by.x = "ensembl_gene_id", by.y = "ENSEMBL", all.x = TRUE)
  
  # Include logFC, logCPM, and LR in the final annotation files
  df_final <- df_final %>% dplyr::select(GeneID, ensembl_gene_id, ENTREZID, logFC, logCPM, LR, FDR)
  
  # Save the intermediate mapped data frame in the output directory
  write.csv(df_final, file.path(output_directory, paste0("mapped_data_", dataset_name, ".csv")), row.names = FALSE)
  
  return(df_final)
}

# Initialize a list to store unmapped genes for each dataset
unmapped_genes_list <- list()

# Apply the mapping function to all DEG data frames
mapped_deg_list <- lapply(seq_along(deg_list), function(i) {
  dataset_name <- basename(files[i])  # Extract the file name for use in output files
  map_to_entrez(deg_list[[i]], dataset_name)
})

print(unmapped_genes_list)

# Perform GO analysis for each dataset using the rat database
go_results_list <- lapply(seq_along(mapped_deg_list), function(i) {
  df <- mapped_deg_list[[i]]
  if (is.null(df)) return(NULL)
  
  significant_genes <- df %>% filter(FDR < 0.05) %>% pull(ENTREZID)
  
  if (length(significant_genes) == 0) {
    return(NULL)
  }
  
  go_enrich <- enrichGO(gene = significant_genes, 
                        OrgDb = org.Rn.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
  
  if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
    return(NULL)
  }
  
  go_results_df <- as.data.frame(go_enrich)
  
  if (!"p.adjust" %in% colnames(go_results_df)) {
    return(NULL)
  }
  
  go_results_df <- go_results_df %>% rename(FDR = p.adjust) %>% arrange(FDR)
  go_results_df$Source <- basename(files[i])
  return(go_results_df)
})

# Combine all GO results into a single data frame
go_combined_df <- do.call(rbind, go_results_list)

# Save the combined GO results in the output directory
write.csv(go_combined_df, file.path(output_directory, "go_combined_results_with_regulation.csv"), row.names = FALSE)