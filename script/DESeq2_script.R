setwd("~/Desktop/working_directory/AppBio_project/")

#install.packages("BiocManager")

#BiocManager::install("DESeq2") 
#BiocManager::install(c("GenomicFeatures", "txdbmaker"))
#install.packages("fdrtool")

library(DESeq2)
library(GenomicFeatures)
library(txdbmaker)
library(ggplot2)
library(Rtsne)
library(umap)
library(matrixStats)

# List all count data files
count_file_list <- list.files(path = "count_files", pattern = "*.txt", full.names = TRUE)
count_file_list

metadata <- read.csv(file ="phenotype_metadata.csv", header = TRUE, row.names = 1)
head(metadata)

### Creating database from GTF file to extract gene names

# Load GTF file
txdb <- makeTxDbFromGFF(file = "genomic.gtf", format = "gtf")
txdb

# Extracting gene data for gene reference list
genes <- genes(txdb)
genes_df <- as.data.frame(genes)

# Extract only the gene IDs
reference_genes <- rownames(genes_df)

# View the first few gene IDs
head(reference_genes)

###

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

# Check for NA values
sum(is.na(count_matrix)) # None

### Differential Expression according to Phenotype

dds_p <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ phenotype) 

<<<<<<< HEAD
dds <- DESeq(dds)
=======
>>>>>>> PCA
# Performs normalization and fits the model
dds_p <- DESeq(dds_p)

<<<<<<< HEAD
results <- results(dds) # Automatically performs independent filtering
head(results)

# fdrtools section - Multiple testing correction
=======
results_p <- results(dds_p) # Automatically performs independent filtering
head(results_p)

# Check for NA values
sum(is.na(results_p)) # 131 NA values
# Where are they?
sum(is.na(results_p$pvalue)) # 5 NAs in p-value

# Excluding NA values from analysis
exNA_results_p <- results_p[!is.na(results_p$pvalue), ]
head(exNA_results_p)

sum(is.na(exNA_results_p)) # 106 NAs

### Differential Expression according to Strain

dds_s <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ strain) 
# Performs normalization and fits the model
dds_s <- DESeq(dds_s) 

results_s <- results(dds_s)
head(results_s)

# Check for NA values
sum(is.na(results_s)) # 25 NA values
# Where are they?
sum(is.na(results_s$pvalue)) # 5 NAs in p-value

# Excluding NA values from analysis
exNA_results_s <- results_s[!is.na(results_s$pvalue), ]
head(exNA_results_s)

sum(is.na(exNA_results_s)) # 0 NAs ???

### fdrtools section - Multiple testing correction
>>>>>>> PCA

library(fdrtool)

# Phenotype differential expression

pvalues_p <- exNA_results_p$pvalue
  
fdr_results_p <- fdrtool(pvalues_p, statistic = "pvalue")
summary(fdr_results_p)

# Add q-values to the original results table
exNA_results_p$qval <- fdr_results_p$qval

head(exNA_results_p)

<<<<<<< HEAD
significant_genes <- subset(results, qval < 0.1)
head(significant_genes)
=======
significant_genes_p <- subset(exNA_results_p, qval < 0.1)
head(significant_genes_p)
>>>>>>> PCA

summary(significant_genes_p)

# Writing to a csv file
write.csv(as.data.frame(significant_genes_p), file = "significant_genes_phenotype.csv")

<<<<<<< HEAD
### PCA plot ###

library(ggplot2)
#install.packages("ggfortify")
=======
# Strain differential expression

pvalues_s <- exNA_results_s$pvalue

fdr_results_s <- fdrtool(pvalues_s, statistic = "pvalue")
summary(fdr_results_s)

# Add q-values to the original results table
exNA_results_s$qval <- fdr_results_s$qval

head(exNA_results_s)

significant_genes_s <- subset(exNA_results_s, qval < 0.1)
head(significant_genes_s)

summary(significant_genes_s)

# Writing to a csv file
write.csv(as.data.frame(significant_genes_s), file = "significant_genes_strain.csv")

### PCA plot ###

library(ggplot2)
>>>>>>> PCA
library(ggfortify)

# Using count_matrix for PCA data
head(count_matrix)

# PCA needs variability across columns to work so constant 0 values need to be removed
zero_var_cols <- apply(count_matrix, 1, function(x) var(x) == 0)
print(colnames(count_matrix)[zero_var_cols])

# Filtering constant 0s
filtered_matrix <- count_matrix[apply(count_matrix, 1, function(x) var(x) > 0), ]

# Applying PCA
pca <- prcomp(t(filtered_matrix), scale. = TRUE)
summary(pca)

<<<<<<< HEAD
# Showing string of phenotype metadata and strain metadata - just to check names
metadata$phenotype
metadata$strain

# Plotting PCA
PCAplot <- ggplot(pca, aes(x = PC1, y = PC2, colour = metadata$phenotype, shape = metadata$strain)) +
         geom_point(size = 3) +
         labs(title = "Transcriptomics PCA Plot", x = "PC1: 35.5% variance", y = "PC2: 28.6% variance", colour = "Phenotype", shape = "Strain Type") +
=======
# Plotting PCA
PCAplot <- ggplot(pca, aes(x = PC1, y = PC2, colour = metadata$phenotype, shape = metadata$strain)) +
         geom_point(size = 3) +
         labs(title = "Transcriptomics PCA Plot", x = "PC1: 46.07% variance", y = "PC2: 13.15% variance", colour = "Phenotype", shape = "Strain Type") +
>>>>>>> PCA
         scale_color_manual(values = c("fumarate_producer" = "indianred2", "malate_producer" = "olivedrab3", "succinate_producer" = "turquoise3", "wild_type" = "orchid2"), labels = c("Fumarate Producer", "Malate Producer", "Succinate Producer", "Wild Type")) +
         scale_shape_manual(values = c("parental" = 15, "evolved" = 16, "wild_type" = 17), labels = c("Parental", "Evolved", "Wild Type")) +
         theme_bw() + 
         theme(legend.position = "bottom")

# Saving plot as image
ggsave("PCAplot.png", plot = PCAplot, width = 12, height = 8)
<<<<<<< HEAD
=======

### This section is not included in the Rmarkdown file

# Create a loading plot
library(ggforce)

loadings <- pca$rotation
print(loadings)

# Convert loadings to a data frame for ggplot2
loadings_df <- as.data.frame(loadings)
loadings_df$Variable <- rownames(loadings_df)

# Plot loadings
ggplot(loadings_df, aes(x = PC1, y = PC2, label = Variable)) +
  geom_point() +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  labs(title = "Loading Plot", x = "PC1", y = "PC2") +
  theme_minimal() +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), linetype = "dashed", inherit.aes = FALSE)
# Plot shows minimal contribution by variables to PC1 & PC2
# This is typical for transcriptomic analysis as the data shows high-dimensionality
>>>>>>> PCA
