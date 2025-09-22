                    ### $$$ sample separation from Xenium slide before batch correction $$$ ###

library(future)
plan("multicore", workers = 12)   

library(Seurat)
# loading seurat objs

MASH_obj_1 <- readRDS('D:/Bioinformatics Files/XeniumR/experiment_1/experiment_1_SH_MASH_obj.RDS')

MASH_obj_4 <- readRDS('D:/Bioinformatics Files/XeniumR/experiment_4/experiment_4_SH_MASH_obj.RDS')

# assigning each slide with a new column
MASH_obj_1$slide <- 1
MASH_obj_4$slide <- 4

MASH_obj_1 <- subset(MASH_obj_1, 
                  subset = nFeature_Xenium > 5 & nFeature_Xenium < 200 & 
                    nCount_Xenium > 10 & nCount_Xenium < 1500)

MASH_obj_4 <- subset(MASH_obj_4, 
                     subset = nFeature_Xenium > 5 & nFeature_Xenium < 200 & 
                       nCount_Xenium > 10 & nCount_Xenium < 1500)

MASH_obj_1 <- SCTransform(MASH_obj_1, assay = "Xenium", verbose = FALSE)
MASH_obj_4 <- SCTransform(MASH_obj_4, assay = "Xenium", verbose = FALSE)

# Merge two Seurat objects
mash_obj <- merge(
  x = MASH_obj_1,
  y = MASH_obj_4,
  add.cell.ids = c("slide1", "slide4")  # prefixes to keep track of origins
)
# mash_ob <- SCTransform(mash_ob, assay = "Xenium", verbose = FALSE)

mash_obj <- PrepSCTFindMarkers(mash_obj)

mash_obj <- RunPCA(mash_obj, assay = "SCT", features = row.names(mash_obj))

################# ploting PCA dim 
DimPlot(mash_obj, group.by = 'slide')
ElbowPlot(mash_obj)

# compute correlations of top PCs with common covariates
pc_emb <- Embeddings(mash_obj, "pca")[, 1:10]
meta <- mash_obj@meta.data

# Example correlations 
cor(pc_emb[,1], meta$nCount_Xenium, use = "complete.obs")
cor(pc_emb[,1], meta$nFeature_Xenium, use = "complete.obs")
# checked for other PCs and it shows no correlation with technical artifacts

# running UMAP
mash_obj <- RunUMAP(mash_obj, dims = 1:30)

# ploting umap dim 
DimPlot(mash_obj, group.by = 'slide', reduction = "umap")


mash_obj <- FindNeighbors(mash_obj, reduction = "pca", dims = 1:30)
mash_obj <- FindClusters(mash_obj, resolution = 0.3)
DimPlot(mash_obj)


library(ggplot2)

# Count cells per sample × cluster
df <- table(mash_obj$sample_ID, Idents(mash_obj))
df <- as.data.frame(df)
colnames(df) <- c("sample_ID", "cluster", "n_cells")

# Plot
ggplot(df, aes(x = sample_ID, y = n_cells, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  labs(title = "Number of cells per cluster per sample",
       y = "Cell count", x = "Sample ID")
