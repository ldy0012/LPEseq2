# Supplementary Note S1 / Figure 2B
# GSE142647 SUM149PT: prepare QC-filtered Seurat objects.

library(Seurat)

# Input and output paths
data_dir <- file.path("S1_SUM149PT", "data")
raw_dir <- file.path(data_dir, "GSE142647_RAW")
tar_file <- file.path(data_dir, "GSE142647_RAW.tar")
output_file <- file.path(data_dir, "sum149pt_qc.rds")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# Find 10x count matrix files
matrix_files <- list.files(
  raw_dir,
  pattern = "_matrix\\.mtx\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

# Select biological replicates
sample_ids <- list(
  DMSO = c(Rep1 = "GSM4235302", Rep2 = "GSM4235310"),
  TSA  = c(Rep1 = "GSM4235306", Rep2 = "GSM4235314")
)

read_sample <- function(gsm, condition, replicate) {
  
  sample_files <- matrix_files[
    startsWith(basename(matrix_files), paste0(gsm, "_"))
  ]
  
  prefix <- sub("_matrix\\.mtx\\.gz$", "", sample_files)
  
  barcode_file <- paste0(prefix, "_barcodes.tsv.gz")
  feature_file <- paste0(prefix, "_genes.tsv.gz")
  
  counts <- ReadMtx(
    mtx = sample_files,
    cells = barcode_file,
    features = feature_file,
    feature.column = 2,
    unique.features = TRUE
  )
  
  obj <- CreateSeuratObject(
    counts,
    min.cells = 3,
    min.features = 200
  )
  
  obj$replicate <- replicate
  obj$sample <- paste0(condition, "_", replicate)
  
  return(obj)
}


# Merge replicates and apply cell QC
prepared <- list()

for (condition in names(sample_ids)) {
  
  ids <- sample_ids[[condition]]
  
  rep1 <- read_sample(ids["Rep1"], condition, "Rep1")
  rep2 <- read_sample(ids["Rep2"], condition, "Rep2")
  
  obj <- merge(
    rep1, rep2,
    add.cell.ids = c("Rep1", "Rep2")
  )
  
  obj <- JoinLayers(obj)
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj, pattern = "^MT-"
  )
  
  obj <- subset(
    obj,
    subset = nFeature_RNA >= 200 &
      nFeature_RNA <= 7500 &
      percent.mt <= 15
  )
  
  obj$celltype <- "all"
  
  n_rep1 <- sum(obj$replicate == "Rep1")
  n_rep2 <- sum(obj$replicate == "Rep2")
  
  prepared[[condition]] <- obj
  
  message(
    condition,
    ": QC-passed Rep1 = ", n_rep1,
    ", Rep2 = ", n_rep2
  )
}


# Save QC filtered raw count objects
saveRDS(prepared, output_file)
message("Saved: ", output_file)
