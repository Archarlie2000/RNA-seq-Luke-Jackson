# Load required libraries
library(tximport)
library(edgeR)

# Read sample information
sample_info <- read.table("metadata.tsv", header = TRUE)

# Define directory and files
quant_dir <- "C:/Users/archa/OneDrive/Documents/RNA-seq Luke Jaxon/raw_reads"
files <- file.path(quant_dir, paste0(sample_info$sample, ".sf"))
names(files) <- sample_info$sample

# Import quantification files using tximport at the transcript level
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Extract transcript-level counts
counts <- txi$counts

# Strip version numbers from transcript IDs to get gene-level IDs
gene_ids <- gsub("\\..*", "", rownames(counts))  # Remove version numbers (e.g., ".2" from "ENSMUST00000205391.2")

# Aggregate counts to gene level
gene_counts <- rowsum(counts, group = gene_ids)

# Create DGEList object for edgeR analysis
dge <- DGEList(counts = gene_counts, samples = sample_info)

# Filtering low-count genes
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalization
dge <- calcNormFactors(dge)

# Design matrix
group <- factor(sample_info$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the model
fit <- glmFit(dge, design)

# Generate comparisons between each group and the control
# Adjust the contrast matrix to include all groups vs. Control
group_names <- levels(group)
control_group <- "Control"
comparisons <- group_names[group_names != control_group]

# Store results for each comparison in a list of data frames
results_df_list <- list()

for (comparison_group in comparisons) {
  # Create the contrast for the current group vs. control
  contrast <- makeContrasts(contrasts = paste(comparison_group, "-", control_group), levels = design)
  lrt <- glmLRT(fit, contrast = contrast)
  results <- topTags(lrt, n = Inf)
  
  # Add the data frame to the list with the comparison group as the name
  results_df_list[[comparison_group]] <- as.data.frame(results)
  
  # Print the top results for this comparison
  cat("Top significant genes for", comparison_group, "vs", control_group, ":\n")
  print(head(results_df_list))
  cat("\n")
}

# Go one level up and create "DEG_results" folder if it doesn't exist
if (!dir.exists("../DEG_results")) {
  dir.create("../DEG_results")
}

# (Optional) Save each comparison's data frame to a file in the "DEG_results" folder
for (comparison_group in comparisons) {
  write.table(results_df_list[[comparison_group]], 
              file = paste0("../DEG_results/DEG_results_ALL_", comparison_group, "_vs_", control_group, ".tsv"), 
              sep = "\t", row.names = TRUE, quote = FALSE)
}


