### $$$ sample separation from Xenium slide $$$ ###

library(Seurat)

setwd("xeniumR")

# Load Xenium object
experiment_4_obj <- LoadXenium("experiment_4")

# Path to CSVs exported from Xenium Explorer
path <- "experiment_4/sample_separation"
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# Initialize sample_ID column
experiment_4_obj$sample_ID <- NA

# Loop over sample CSVs
for (file in files) {
  df <- read.csv(file, skip = 2)
  cell_ids <- df[[1]]
  
  # Parse sample name from filename
  filename <- basename(file)
  name_parts <- strsplit(filename, "_")[[1]]
  key <- paste(name_parts[3:length(name_parts)], collapse = "_") 
  key <- sub("\\.csv$", "", key)
  
  # Assign sample_ID to the cells
  experiment_4_obj$sample_ID[colnames(experiment_4_obj) %in% cell_ids] <- key
}

# Check assignment
table(experiment_4_obj$sample_ID, useNA = "ifany")


# Subset by disease type
mash_obj <- subset(experiment_4_obj, subset = grepl("MASH", sample_ID))
hcv_obj  <- subset(experiment_4_obj, subset = grepl("HCV", sample_ID))

# adding another column for etiology
mash_obj$etiology <- 'MASH'
hcv_obj$etiology <- 'HCV'

# saving .RDS
saveRDS(mash_obj, 'experiment_4/experiment_4_SH_MASH_obj.RDS')
saveRDS(hcv_obj, 'experiment_4/experiment_4_SH_HCV_obj.RDS')