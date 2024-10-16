# Load required libraries
library(clusterProfiler)
library(org.Rn.eg.db)       # Rat-specific database
library(biomaRt)            # For ID conversion
library(orthologsBioMART)   # For finding human orthologs
library(dplyr)              # For data manipulation
library(graphite)           # For pathways
library(SPIA)               # For SPIA analysis

# Set the directory containing the DEG files and create the output directory
deg_directory <- "../DEG_Results"
output_directory <- "../KEGG_Human_Ortholog_Results"
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

# Function to map mouse Ensembl gene IDs to human orthologs
map_to_human_orthologs <- function(df, dataset_name) {
  # Extract unique Ensembl gene IDs
  vec_ids <- unique(df$ensembl_gene_id)
  
  # Map rat Ensembl gene IDs to human orthologs
  orthologs <- findOrthologsMmHs(from_filters = "ensembl_gene_id",
                                 from_values = vec_ids,
                                 to_attributes = "ensembl_gene_id")
  
  # Rename columns to avoid conflicts
  colnames(orthologs) <- c("ensembl_gene_id", "human_ensembl_gene_id")
  
  # Merge orthologs with the main data frame
  df_final <- merge(df, orthologs, by = "ensembl_gene_id", all.x = TRUE)
  
  # Save the data with human orthologs in the output directory
  write.csv(df_final, file.path(output_directory, paste0("mapped_with_human_orthologs_", dataset_name, ".csv")), row.names = FALSE)
  
  return(df_final)
}

# Apply the human ortholog mapping to all DEG data frames
human_mapped_deg_list <- lapply(seq_along(mapped_deg_list), function(i) {
  dataset_name <- basename(files[i])  # Extract the file name for use in output files
  map_to_human_orthologs(mapped_deg_list[[i]], dataset_name)
})

# Load pathways and prepare for SPIA
reactome <- pathways("hsapiens", "reactome")
reactome <- convertIdentifiers(reactome, "ENTREZID")
prepareSPIA(reactome, "reactome")

wikipathways <- pathways("hsapiens", "wikipathways")
wikipathways <- convertIdentifiers(wikipathways, "ENTREZID")
prepareSPIA(wikipathways, "wikipathways")

# Set up biomaRt to use Ensembl for human
ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Function to convert Ensembl gene IDs to Entrez IDs
convert_to_entrez_ids <- function(ensembl_ids) {
  gene_mapping <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                        filters = "ensembl_gene_id", 
                        values = ensembl_ids, 
                        mart = ensembl_human)
  return(gene_mapping)
}

# Perform SPIA analysis for each dataset
spia_output_directory <- "SPIA_Human_Orthologs_Pathway"
dir.create(spia_output_directory, showWarnings = FALSE)

for (file in list.files(output_directory, pattern = "mapped_with_human_orthologs_.*\.csv", full.names = TRUE)) {
  # Read the CSV file
  deg_data <- read.csv(file)
  
  # Extract unique Ensembl gene IDs
  ensembl_ids <- unique(deg_data$human_ensembl_gene_id)
  
  # Convert Ensembl IDs to Entrez IDs
  gene_mapping <- convert_to_entrez_ids(ensembl_ids)
  
  # Merge gene mapping with the main data frame
  deg_data <- merge(deg_data, gene_mapping, by.x = "human_ensembl_gene_id", by.y = "ensembl_gene_id", all.x = TRUE)
  
  # Remove rows without Entrez IDs
  deg_data <- deg_data[!is.na(deg_data$entrezgene_id), ]
  deg_data <- deg_data[!duplicated(deg_data$entrezgene_id), ]
  
  # Extract differentially expressed genes (DE_LIST) and all genes (ALL_Genes)
  tg11 <- subset(deg_data, FDR < 0.05)
  DE_LIST <- tg11$logFC
  names(DE_LIST) <- as.vector(tg11$entrezgene_id)
  ALL_Genes <- as.character(unique(deg_data$entrezgene_id))
  
  # Perform SPIA analysis
  results_KEGG <- spia(de = DE_LIST, all = ALL_Genes, organism = "hsa", beta = NULL, nB = 2000, plots = FALSE, verbose = TRUE, combine = "fisher")
  results_reactome <- runSPIA(de = DE_LIST, all = ALL_Genes, "reactome")
  results_wikipathways <- runSPIA(de = DE_LIST, all = ALL_Genes, "wikipathways")
  
  # Concatenate results with source information
  results_KEGG$source <- "KEGG"
  results_reactome$source <- "Reactome"
  results_wikipathways$source <- "WikiPathways"
  
  combined_results <- rbind(results_KEGG, results_reactome, results_wikipathways)
  
  # Filter significant results and save to file
  if (nrow(combined_results) > 0) {
    sig_results <- combined_results[combined_results$pGFWER < 0.05, ]
    if (nrow(sig_results) > 0) {
      output_file <- file.path(spia_output_directory, paste0(tools::file_path_sans_ext(basename(file)), "_SPIA_Results.csv"))
      write.csv(sig_results, file = output_file, row.names = FALSE)
    } else {
      print("No significant hits")
    }
  }
}

print("SPIA analysis completed for all datasets.")