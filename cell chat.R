# Data pre-processing
library(CellChat)
merged.data@meta.data$cellType_grouped
merged.data <- JoinLayers(merged.data)
table(merged.data$Type)
cancer_subset <- subset(merged.data, subset = Type %in% c("CC"))
# Create cellchat object
cellchat <- createCellChat(object = cancer_subset, group.by = "cellType", assay = "RNA")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
options(future.globals.maxSize = 12 * 1024^3) # Set to 2 GB. You might need to adjust this value.
future::plan("multisession", workers = 4) # do parallel

##Compute the communication probability and infer cellular communication network
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), 
                 ncol = ncol(mat), 
                 dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])}
##Part III: Visualization of cell-cell communication network
#Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
cellchat@netP$pathways
pathways.show <- c("CXCL") 
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")


par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")












