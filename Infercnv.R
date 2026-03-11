
cervical_cells <- subset(merged.data, subset = cellType_tcell == "Epithelial cells")

cervical <- subset(cervical_cells, subset = Type %in% c("CC", "CN"))
#headneck <- subset(merged.data, subset = Type %in% c("HNC", "HNN"))


library(Seurat)
#cervical <- subset(merged.data, subset = cellType_tcell == "Epithelial cells")
cervical <- JoinLayers(cervical)
counts <- LayerData(cervical, assay = "RNA", layer = "counts")
cell_annotations <- data.frame(
  cell = colnames(cervical),
  group = cervical$Type
)
#write.table(cell_annotations, "cell_annotations.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(cell_annotations, "cell_annotations.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

library(rtracklayer)
gtf <- import("gencode.v48.basic.annotation.gtf.gz")
genes <- gtf[gtf$type == "gene"]
gene_order <- data.frame(
  gene = genes$gene_name,
  chr = as.character(seqnames(genes)),
  start = start(genes),
  end = end(genes)
)
gene_order <- gene_order[gene_order$gene %in% rownames(counts), ]
gene_order <- gene_order[order(gene_order$chr, gene_order$start), ]
rownames(gene_order) <- gene_order$gene
# Keep the earliest start position per gene
gene_order_unique <- gene_order[order(gene_order$start), ]
gene_order_unique <- gene_order_unique[!duplicated(gene_order_unique$gene), ]

# Now set rownames safely
rownames(gene_order_unique) <- gene_order_unique$gene

write.table(gene_order_unique, "gene_order.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
counts <- as.matrix(counts)
#storage.mode(counts) <- "numeric"
library(Matrix)
counts <- as(counts, "dgCMatrix")  # Convert to sparse format
library(infercnv)
serialize(object, connection = NULL, ascii = FALSE)
memory.limit(size = 16000)  # Set to 16GB, adjust as needed

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = "cell_annotations.txt",
  delim = "\t",
  gene_order_file = "gene_order.txt",
  ref_group_names = c("CN")  # adjust based on your `Type` values
)
#Pre-Cancer

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,  # adjust based on sparsity
  out_dir = "infercnv_output",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 16
)

