# script to perform WGCNA
# setwd("~/Desktop/WGCNA")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()         

# 1. Fetch Data ------------------------------------------------
data <- as.data.frame(final_counts_tcga_gtex)
colData <- as.data.frame(MET)
rownames(data) <- data$gene
data <- data[,-352]
data <- data[,colnames(data) %in% col$barcode]
MET <- MET[MET$barcode %in% colnames(data),]

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


# pca - method 2
#data <- data[complete.cases(data), ]
pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


### NOTE: If there are batch effects observed, correct for them before moving ahead

# exclude outlier samples
samples.to.be.excluded <- c('TCGA-DS-A0VN-01A-21R-A10U-07','TCGA-EK-A2RB-01A-11R-A18M-07','Data9_8','TCGA-C5-A1BQ-01C-11R-A213-07','TCGA-EA-A3HU-01A-11R-A213-07','Data10_8','Data10_4')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization ----------------------------------------------------------------------

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))
colData <- as.data.frame(MET)
MET <- MET[MET$barcode %in% colnames(data.subset),]

#final_merged_metadata <- as.data.frame(final_merged_metadata)
rownames(colData) <- colData$barcode

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75)


# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()


# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.9, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 9
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 19000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor


# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

tiff("DENDO_CCC.tiff", width = 10, height = 7, units = "in", res = 600, compression = "lzw")
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
dev.off()

tiff("POWER_CC.tiff", 
     width = 7, height = 9, units = "in", res = 600, compression = "lzw")
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
grid.arrange(a1, a2, nrow = 2)
dev.off()

tiff("PCA_CC.tiff", 
     width = 10, height = 8, units = "in", res = 600, compression = "lzw")
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
dev.off()
moduleColors = bwnet$colors
unmergedColors = bwnet$unmergedColors



# 6A. Relate modules to traits --------------------------------------------------
cor <- WGCNA::cor
colData$tissue_numeric <- ifelse(colData$tissue_type == "Tumor", 1, 0)
colData <- colData[colnames(data.subset), ]
stopifnot(all(rownames(colData) == colnames(data.subset)))

# === Reorder eigengenes ===
MEs_col <- orderMEs(module_eigengenes)
eigengene_names <- names(MEs_col)
module_colors_clean <- gsub("^ME", "", eigengene_names)
color_df <- data.frame(Module = eigengene_names,
                       Color = module_colors_clean,
                       stringsAsFactors = FALSE)

# === Trait matrix ===
trait_data <- model.matrix(~ 0 + tissue_type, data = colData)
colnames(trait_data) <- gsub("tissue_type", "", colnames(trait_data))
colnames(trait_data)[colnames(trait_data) == "Tumor"] <- "CESC"
rownames(trait_data) <- rownames(colData)

# === Correlation and p-values ===
module_trait_cor <- cor(MEs_col, trait_data, use = "p")
module_trait_pvalue <- corPvalueStudent(module_trait_cor, nSamples = nrow(colData))

# === Melted long format ===
library(reshape2)
cor_melt <- melt(module_trait_cor)
pval_melt <- melt(module_trait_pvalue)
colnames(cor_melt) <- c("Module", "Trait", "Correlation")
colnames(pval_melt) <- c("Module", "Trait", "Pvalue")
cor_melt$Pvalue <- pval_melt$Pvalue
cor_melt <- merge(cor_melt, color_df, by = "Module")

# === Significance stars + labels ===
cor_melt$Signif <- ifelse(cor_melt$Pvalue < 0.001, "*",
                          ifelse(cor_melt$Pvalue < 0.01, "",
                                 ifelse(cor_melt$Pvalue < 0.05, "*", "")))
cor_melt$Label <- paste0(signif(cor_melt$Correlation, 2), cor_melt$Signif)

# === Reorder modules ===
module_order <- unique(cor_melt[order(abs(cor_melt$Correlation), decreasing = TRUE), "Module"])
cor_melt$Module <- factor(cor_melt$Module, levels = rev(module_order))

# === Color strip ===
color_strip_df <- unique(cor_melt[, c("Module", "Color")])
color_strip_df$Module <- factor(color_strip_df$Module, levels = levels(cor_melt$Module))

color_plot <- ggplot(color_strip_df, aes(x = 1, y = Module, fill = Module)) +
  geom_tile(width = 1) +
  scale_fill_manual(values = setNames(color_df$Color, color_df$Module), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

# === Heatmap (no space between columns) ===
heatmap_plot <- ggplot(cor_melt, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(width = 1) +  # <-- set width to 1 to remove space
  geom_text(aes(label = Label), size = 4.2) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  labs(title = "Module–Trait Relationships", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 0),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(2.5, "cm")  # Extend colorbar
  )

# === Final plot ===
library(cowplot)
final_plot <- plot_grid(color_plot, heatmap_plot, 
                        ncol = 2, rel_widths = c(0.08, 1), align = "h")
print(final_plot)

# Save as TIFF with high resolution in Cervical Cancer Project folder
tiff("module_trait_relationships_CC.tiff", 
     width = 9, height = 8, units = "in", res = 600, compression = "lzw")
print(final_plot)
dev.off()

module_assignments <- bwnet$colors


modules_to_extract <- c("turquoise", "red", "blue")
gene_modules <- data.frame(Gene = names(bwnet$colors), Module = bwnet$colors)
selected_genes <- gene_modules[gene_modules$Module %in% modules_to_extract, ]
selected_genes <- selected_genes[order(selected_genes$Module), ]

library(writexl)

write_xlsx(selected_genes,"module_names_of_genes.xlsx")
#################

# Extract genes for the 'turquoise' module
turquoise_genes <- names(module_assignments[module_assignments == "turquoise"])

# Extract genes for the 'red' module
red_genes <- names(module_assignments[module_assignments == "red"])
blue_genes <- names(module_assignments[module_assignments == "blue"])

# Extract genes for the 'brown' module
gene_lists_to_save <- list(
  "Magenta_Module" = data.frame(Genes = turquoise_genes),
  "Pink_Module" = data.frame(Genes = red_genes),
  "blue_Module" = data.frame(Genes = blue_genes)
)

# Define the name for your output Excel file
output_filename <- "HPV_CC_GENES_WGCNA.xlsx"

# Write the list to an Excel file
writexl::write_xlsx(gene_lists_to_save, path = output_filename)








