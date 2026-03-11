# Get the metadata
meta <- merged.data@meta.data

# Create a mapping from orig.ident to Type
sample_types <- unique(meta[, c("orig.ident", "Type")])

# Add a numeric ID to each sample
sample_types$sample_id <- paste0(sample_types$Type, "_", seq_len(nrow(sample_types)))

# Merge back to metadata
meta <- merge(meta, sample_types[, c("orig.ident", "sample_id")], by = "orig.ident", all.x = TRUE)

# Assign to new column
merged.data@meta.data$samples <- meta$sample_id

epi_cells <- subset(merged.data, subset = cellType_tcell == "Epithelial cells")

CC_CN <- subset(merged.data, subset = Type %in% c("CC", "CN"))
HNC_HNN <- subset(merged.data, subset = Type %in% c("HNC", "HNN"))

library(Seurat)
DefaultAssay(CC_CN)
# For CC vs CN
pseudo_CC_CN <- AggregateExpression(CC_CN, 
                                    group.by = c("cellType_tcell", "samples"), 
                                    assays = "RNA", 
                                    slot = "counts")
pseudo_CC_CN <- pseudo_CC_CN$RNA

library(Matrix)

# Convert sparse matrix to dense matrix first
dense_matrix <- as.matrix(pseudo_CC_CN)

# Then convert to data frame
pseudo_CC_CN_df <- as.data.frame(dense_matrix)

# Select only columns that contain "Epithelial cells"
counts_epi_cc <- pseudo_CC_CN_df[, grepl("Epithelial cells", colnames(pseudo_CC_CN_df))]

colData <- data.frame(samples = colnames(counts_epi_cc))
# Assuming your object is called dds (e.g., DESeqDataSet)
# And colnames(dds) look like: "Epithelial cells_CC-1", "Epithelial cells_CN-6", etc.

# Extract condition from sample names
condition <- ifelse(grepl("CC", colnames(counts_epi_cc)), "CC", "CN")

# Add to colData
colData$condition <- condition
rownames(colData) <- colData$samples
colData <- as.data.frame(colData)
colData <- colData[,-1, drop = FALSE]

# Perform DESEq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = epi_only_df,
                              colData = colData,
                              design = ~condition)
keep <- rowSums(counts(dds))>=10
dds <- dds[keep,]
colData(dds)$condition <- factor(colData(dds)$condition, levels = c("CN", "CC"))
dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds, name = "condition_CC_vs_CN")
res <- as.data.frame(res)
res$gene <- rownames(res)

# For HNC vs HNN
HNC_HNN <- subset(merged.data, subset = Type %in% c("HNC", "HNN"))

library(Seurat)
DefaultAssay(HNC_HNN)
# For CC vs CN
pseudo_HNC_HNN <- AggregateExpression(HNC_HNN, 
                                    group.by = c("cellType_tcell", "samples"), 
                                    assays = "RNA", 
                                    slot = "counts")
pseudo_HNC_HNN <- pseudo_HNC_HNN$RNA

library(Matrix)

# Convert sparse matrix to dense matrix first
dense_matrix <- as.matrix(pseudo_HNC_HNN)

# Then convert to data frame
pseudo_HNC_HNN_df <- as.data.frame(dense_matrix)

# Select only columns that contain "Epithelial cells"
counts_epi_hnc <- pseudo_HNC_HNN_df[, grepl("Epithelial cells", colnames(pseudo_HNC_HNN_df))]

colData_hnc <- data.frame(samples = colnames(counts_epi_hnc))
# Assuming your object is called dds (e.g., DESeqDataSet)
# And colnames(dds) look like: "Epithelial cells_CC-1", "Epithelial cells_CN-6", etc.

# Extract condition from sample names
condition_hnc <- ifelse(grepl("HNC", colnames(counts_epi_hnc)), "HNC", "HNN")

# Add to colData
colData_hnc$condition_hnc <- condition_hnc
rownames(colData_hnc) <- colData_hnc$samples
colData_hnc <- as.data.frame(colData_hnc)
colData_hnc <- colData_hnc[,-1, drop = FALSE]

# Perform DESEq2
library(DESeq2)
dds_hnc <- DESeqDataSetFromMatrix(countData = counts_epi_hnc,
                              colData = colData_hnc,
                              design = ~condition_hnc)
keep <- rowSums(counts(dds_hnc))>=10
dds_hnc <- dds_hnc[keep,]
colData(dds_hnc)$condition_hnc <- factor(colData(dds_hnc)$condition_hnc, levels = c("HNN", "HNC"))
dds_hnc <- DESeq(dds_hnc)

resultsNames(dds_hnc)
res_hnc <- results(dds_hnc, name = "condition_hnc_HNC_vs_HNN")
res_hnc <- as.data.frame(res_hnc)
res_hnc$gene <- rownames(res_hnc)
writexl::write_xlsx(res_hnc, "HNNvsHNC_EPI_degs.xlsx")
