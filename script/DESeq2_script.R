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
library(dplyr)

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

# Grouping Phenotype & Strain together to avoid redundancy in Differential Expression
metadata$group <- factor(paste(metadata$phenotype, metadata$strain, sep = "_"))
metadata

# Temporarily reading group as characters for renaming process
metadata$group <- as.character(metadata$group)

# Renaming wild_type_wild_type to just wild_type
metadata <- metadata %>%
  mutate(group = ifelse(group == "wild_type_wild_type", "wild_type", group))

metadata

# Changing back to factors for Differential Expression
metadata$group <- as.factor(metadata$group)

### Differential Expression according to Phenotype and Strain

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ group) 

# Performs normalization and fits the model
dds <- DESeq(dds)

results <- results(dds) # Automatically performs independent filtering
head(results)

# Check for NA values
sum(is.na(results)) # 133 NA values
sum(is.na(results$pvalue)) # 6 NAs in p-value

# Excluding NA values from pvalue
exNA_results <- results[!is.na(results$pvalue), ]

### fdrtools section - Multiple testing correction

library(fdrtool)

pvalues <- exNA_results$pvalue
  
fdr_results <- fdrtool(pvalues, statistic = "pvalue")
summary(fdr_results)

# Add q-values to the original results table
exNA_results$qval <- fdr_results$qval

head(exNA_results)

significant_genes <- subset(exNA_results, qval < 0.1)
head(significant_genes)

summary(significant_genes)

# Writing to a csv file
write.csv(as.data.frame(significant_genes), file = "significant_genes.csv")

### PCA plot ###

library(ggplot2)
library(ggfortify)

### Constructing a Differential Expression matrix for PCA

# Creating normalized count matrix for PCA
norm_counts <- counts(dds, normalized = TRUE)
summary(norm_counts)

# Applying PCA
pca <- prcomp(t(norm_counts), scale. = TRUE)
summary(pca)

# Reading variables as factors due to issues in PCA legend
metadata$phenotype <- factor(metadata$phenotype, levels = c("fumarate_producer", "malate_producer", "succinate_producer", "wild_type"))
metadata$strain <- factor(metadata$strain, levels = c("parental", "evolved", "wild_type"))

# Plotting PCA
PCAplot <- ggplot(pca, aes(x = PC1, y = PC2, colour = metadata$phenotype, shape = metadata$strain)) +
         geom_point(size = 3) +
         labs(title = "Transcriptomics PCA Plot", x = "PC1: 22.74% variance", y = "PC2: 15.77% variance", colour = "Phenotype", shape = "Strain Type") +
         scale_color_manual(values = c("fumarate_producer" = "indianred2", "malate_producer" = "olivedrab3", "succinate_producer" = "turquoise3", "wild_type" = "orchid2"), labels = c("Fumarate Producer", "Malate Producer", "Succinate Producer", "Wild Type")) +
         scale_shape_manual(values = c("parental" = 15, "evolved" = 16, "wild_type" = 17), labels = c("Parental", "Evolved", "Wild Type")) +
         theme_bw() + 
         theme(legend.position = "bottom")

# Saving plot as image
ggsave("PCAplot.png", plot = PCAplot, width = 12, height = 8)

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
# Plot shows minimal contribution by any specific gene to PC1 & PC2
# This is typical for transcriptomic analysis as the data shows high-dimensionality
