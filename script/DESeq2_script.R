setwd("~/Desktop/working_directory/AppBio_project/")

#install.packages("BiocManager")

#BiocManager::install("DESeq2") 
#BiocManager::install(c("GenomicFeatures", "txdbmaker"))


library(DESeq2)
#library(GenomicFeatures)
#library(txdbmaker)
library(ggplot2)
library(Rtsne)
library(umap)
library(matrixStats)

# List all count data files
count_file_list <- list.files(path = "count_files", pattern = "*.txt", full.names = TRUE)
count_file_list

metadata <- read.csv(file ="phenotype_metadata.csv", header = TRUE, row.names = 1)
head(metadata)

### 
# This section cannot be done on HPC as GenomicFeatures is not available 
###

# Load GTF file
# txdb <- makeTxDbFromGFF(file = "genomic.gtf", format = "gtf")
# txdb

# Extracting gene data for gene reference list
# genes <- genes(txdb)
# genes_df <- as.data.frame(genes)

# Extract only the gene IDs
# gene_ids <- rownames(genes_df)

# View the first few gene IDs
# head(gene_ids)

# Save the gene list to a file
# write.table(gene_ids, file = "gene_ids.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

### End of GenomicFeatures section
 
# Instead we have generated the file locally and sftp'd it to the HPC

reference_genes <- read.table("gene_ids.txt", header = FALSE, stringsAsFactors = FALSE)
reference_genes <- reference_genes$V1  # Extract the gene names

# Function to Align Counts to Reference Genes
read_and_align <- function(file, reference_genes) {
 
   # Read the HTSeq file
  df <- read.table(file, header = FALSE, row.names = 1, col.names = c("gene_id", "count"))
  
  # Remove summary rows (e.g., __no_feature, __ambiguous)
  df <- df[!grepl("^__", rownames(df)), , drop = FALSE]
  
  # Create a full matrix aligned to the reference genes
  aligned_df <- data.frame(count = rep(0, length(reference_genes)), row.names = reference_genes)
  
  # Fill in the counts for existing genes
  aligned_df[rownames(df), "count"] <- df$count
  
  # Return the aligned counts
  colnames(aligned_df) <- gsub(".txt$", "", basename(file))  # Use the sample name as column name
  return(aligned_df)
}

# Read and align all files to the reference gene list
aligned_counts <- lapply(count_file_list, function(file) read_and_align(file, reference_genes))

# Combine all aligned data frames into a single matrix
count_matrix <- do.call(cbind, aligned_counts)

# Check the resulting matrix
head(count_matrix)
tail(count_matrix)

# Check if metadata row names match count matrix column names
all(rownames(metadata) == colnames(count_matrix))

###

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ phenotype) 

dds <- dds[rowSums(counts(dds) >= 10) >= 2, ] 
# Remove genes with low counts

dds <- DESeq(dds)
# Performs normalization and fits the model

results <- results(dds)
#results <- results[order(results$padj), ]
head(results)

### fdrtools section ###

install.packages("fdrtool")
library(fdrtool)

pvalues <- results$pvalue
  
fdr_results <- fdrtool(pvalues, statistic = "pvalue")
summary(fdr_results)

# Add q-values to the original results table
results$qval <- fdr_results$qval

head(results)

significant_genes <- subset(results, qval < 0.1)
head(significant_genes)
significant_genes

###

write.csv(as.data.frame(results), file = "differential_expression_results.csv")
write.csv(as.data.frame(significant_genes), file = "significant_genes_results.csv")

