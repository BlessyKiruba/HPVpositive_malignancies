#Differential expression analysis

# Load necessary libraries
library(DESeq2)
library(readxl)
library(dplyr)
library(writexl)
library(ggplot2)

library(readxl)
library(dplyr)

# Step 1: Load the full data including gene symbols
counts_df <- read_excel("count_matrix.xlsx")
# Step 2: Convert to data frame (not tibble)
counts_df <- as.data.frame(counts_df)
# Step 3: Set gene symbols as rownames
rownames(counts_df) <- counts_df$SYMBOL
# Step 4: Remove SYMBOL column
counts_df <- counts_df[, -which(names(counts_df) == "SYMBOL")]

map_data <- read_excel("Metadata.xlsx")

# Step 2.5: Verify that count matrix columns match metadata barcodes
all(colnames(counts_df) %in% map_data$Barcode)  # Should return TRUE
all(colnames(counts_df) == map_data$Barcode)    # Should return TRUE

# Step 3: Set rownames in metadata to match count column names
rownames(map_data) <- map_data$Barcode

# Step 4: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_df)),
                              colData = map_data,
                              design = ~ tissue_type)

# Step 5: Filter low-expressed genes 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
# Step 6: Set reference level for condition (e.g., "Normal")
dds$tissue_type <- relevel(dds$tissue_type, ref = "Normal")

# Step 7: Run DESeq2
dds <- DESeq(dds)

# Step 8: Extract results
res <- results(dds)
summary(res)


# Convert to data frame and remove NAs
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)

# Add gene symbol column
res_df$GeneSymbol <- rownames(res_df)

# Add significance labels
res_df$significance <- "Not Significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Reorder columns: GeneSymbol first, then rest
res_df <- res_df[, c("GeneSymbol", setdiff(names(res_df), "GeneSymbol"))]

library(writexl)
# Save to Excel
write_xlsx(res_df, path = "differentially_expressed_gene_results.xlsx")

# Define custom palette for volcano plot
palette_colors <- c(
  "#1F77B4", "#D62728", "#2CA02C",
  "#FF7F0E", "#4E79A7", "#98DF8A",
  "#FF9896", "#A0CBE8", "#FFBB78", 'grey'
)

# Assign colors by significance class
res_df$color <- ifelse(res_df$significance == "Upregulated", palette_colors[2], 
                       ifelse(res_df$significance == "Downregulated", palette_colors[1], 
                              palette_colors[10]))  # grey for NS

# Save to high-quality TIFF
tiff("C:/Users/derri/Desktop/Internship/Cervical Cancer Project/VolcanoPlot_Uniform.tiff",
     width = 8, height = 8, units = "in", 
     res = 600, compression = "lzw")

# Volcano plot 
plot(res_df$log2FoldChange, -log10(res_df$padj),
     pch = 20,
     col = res_df$color,
     xlab = "Log2 Fold Change",
     ylab = "-Log10 Adjusted P-value",
     main = "Volcano Plot",
     cex = 0.8)

# Add thresholds
abline(h = -log10(0.05), col = "black", lty = 2, lwd = 1)  # p-adj cutoff
abline(v = c(-1, 1), col = "black", lty = 2, lwd = 1)      # log2FC cutoff

# Add legend
legend("topright",
       legend = c("Upregulated", "Downregulated", "Not Significant"),
       col = c(palette_colors[2], palette_colors[1], palette_colors[10]),
       pch = 20,
       pt.cex = 1.2,
       bty = "n")

dev.off()




