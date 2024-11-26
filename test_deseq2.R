# test_deseq2.R
#testing
# Check if R is running
cat("R is working.\n")

# Test loading DESeq2
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is not installed or cannot be loaded.")
} else {
  cat("DESeq2 is successfully loaded.\n")
}

# Run a simple DESeq2 example
library(DESeq2)

# Simulate a small dataset
count_data <- matrix(
  rpois(2000, lambda = 10),
  nrow = 200,
  ncol = 10,
  dimnames = list(paste0("gene", 1:200), paste0("sample", 1:10))
)
condition <- factor(rep(c("A", "B"), each = 5)) # 5 samples per condition
coldata <- data.frame(row.names = colnames(count_data), condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)
results <- results(dds)

# Print results
print("DESeq2 analysis completed successfully!")
print(results)
