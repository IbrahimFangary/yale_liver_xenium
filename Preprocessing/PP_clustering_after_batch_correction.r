### $$$ sample separation from Xenium slide before batch correction $$$ ###

library(future)
plan("multicore", workers = 12)   
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
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

# mash_obj <- PrepSCTFindMarkers(mash_obj)


mash_obj <- RunPCA(mash_obj, assay = "SCT", features = rownames(mash_obj))


### running harmony
mash_obj <- mash_obj %>% RunHarmony("sample_ID")

harmony.embeddings <- Embeddings(mash_obj, reduction = "harmony")
p1 <- DimPlot(object = mash_obj, reduction = "harmony", pt.size = .1, group.by = "sample_ID")
p2 <- VlnPlot(object = mash_obj, features = "harmony_1", group.by = "sample_ID",  pt.size = .1)
plot_grid(p1,p2)


mash_obj <- RunUMAP(mash_obj, dims = 1:30, reduction = "harmony")

mash_obj <- mash_obj %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.3) 

# ploting umap after harmonizing
p1 <- DimPlot(mash_obj, reduction = "umap", group.by = "sample_ID")
p2 <- DimPlot(mash_obj, reduction = "umap", group.by = "seurat_clusters")
plot_grid(p1, p2)

DimPlot(mash_obj, reduction = "umap", group.by = "seurat_clusters")


mash_obj <- PrepSCTFindMarkers(mash_obj)
all.markers <- FindAllMarkers(mash_obj)

saveRDS(all.markers, 'MASH_allmarkers.rds')
saveRDS(mash_obj, 'MASH_after_harmony.rds')

mash_obj<- readRDS('MASH_after_harmony.rds')

library(dplyr)

top10_markers <- all.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

print(top10_markers)

genes_per_cluster <- top10_markers %>%
  group_by(cluster) %>%
  summarize(genes = paste(gene, collapse = ", "))

print(genes_per_cluster)


mash_obj$celltype <- plyr::mapvalues(
  mash_obj$seurat_clusters,
  from = 0:16,  
  to = c("Endothelial Cell / Macrophage",        
         "Fibroblast-like Cell",
         "Macrophages",
         "plasma cells",
         "Proliferating cells",
         "Likely tumor/stromal cells",
         "Cytotoxic T cells",
         "liver sinusoidal endothelial cells (LSECs)",
         "Activated hepatic stellate cell",
         "stromal cells/NMI",
         "Mast cells",
         "Plasma cells 2",
         "Cancer-associated fibroblasts",
         "Hepatocytes",
         "Gamma delta T cells",
         "Fibroblasts with immune modulation",
         "Tumor Cells"
  )
)


meta1 <- subset(mash_obj, subset = slide ==1)@meta.data


MASH_obj_1$celltype <- meta1$celltype


meta4 <- subset(mash_obj, subset = slide ==4)@meta.data

MASH_obj_4$celltype <- meta4$celltype


meta4$cell_id <- gsub("^slide1_", "", rownames(meta4))
head(meta4)

meta4 <- meta4[, c("cell_id","celltype")]
meta <- rbind(meta1, meta4)



library(cowplot)

plots <- list()

for (i in unique(mash_obj$sample_ID)) {
  
  # Initialize empty plots
  p1 <- NULL
  p2 <- NULL
  
  # Only plot if sample_ID exists in MASH_obj_1
  if (i %in% MASH_obj_1$sample_ID) {
    p1 <- ImageDimPlot(subset(MASH_obj_1, subset = sample_ID == i), group.by = 'celltype') 
  }
  
  # Only plot if sample_ID exists in MASH_obj_4
  if (i %in% HCV_obj_4$sample_ID) {
    p2 <- ImageDimPlot(subset(MASH_obj_4, subset = sample_ID == i), group.by = 'celltype') 
  }
  
  # Store both plots for this sample
  plots[[i]] <- plot_grid(p1, p2, ncol = 2)
}

# Combine all samples into one big plot grid
plots[[7]]