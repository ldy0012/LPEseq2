# Supplementary Note S2: prepare U2OS DRUG-seq input data (GSE176150)
# Input: GSE176150_RAW.tar and GSE176150_metadata.csv.gz in S2_U2OS/data/
# Output: S2_U2OS/data/u2os_counts_meta.rds

input_dir <- file.path("S2_U2OS", "data")
tar_file <- file.path(input_dir, "GSE176150_RAW.tar")
metadata_file <- file.path(input_dir, "GSE176150_metadata.csv.gz")
raw_dir <- file.path(input_dir, "GSE176150_RAW")
output_file <- file.path(input_dir, "u2os_counts_meta.rds")

# 1. Extract the raw data archive.
if (!dir.exists(raw_dir)) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  untar(tar_file, exdir = raw_dir)
}

# 2. Read all nine plate-level UMI count matrices.
files <- sort(list.files(raw_dir, pattern = "UMI_results_.*\\.txt\\.gz$",
                         recursive = TRUE, full.names = TRUE))
if (length(files) != 9) stop("Expected nine UMI files; found ", length(files))

matrices <- vector("list", length(files))
for (i in seq_along(files)) {
  plate_table <- read.delim(gzfile(files[i]), check.names = FALSE)
  if (!("GeneId" %in% names(plate_table))) stop("GeneId column missing in ", files[i])

  plate_matrix <- as.matrix(plate_table[, names(plate_table) != "GeneId", drop = FALSE])
  rownames(plate_matrix) <- plate_table$GeneId
  matrices[[i]] <- plate_matrix
}

# 3. Keep genes shared by all plates, then combine the wells.
shared_genes <- rownames(matrices[[1]])
for (i in 2:length(matrices)) {
  shared_genes <- intersect(shared_genes, rownames(matrices[[i]]))
}

aligned_matrices <- vector("list", length(matrices))
for (i in seq_along(matrices)) {
  aligned_matrices[[i]] <- matrices[[i]][shared_genes, , drop = FALSE]
}
counts <- do.call(cbind, aligned_matrices)

# 4. Match each count-matrix well to its metadata row.
meta <- read.csv(gzfile(metadata_file), stringsAsFactors = FALSE,
                 check.names = FALSE)
well_ids <- paste(meta$plate_barcode, meta$well_index, sep = "_")
positions <- match(colnames(counts), well_ids)
if (anyNA(positions)) stop("Some count-matrix wells have no metadata")

meta <- meta[positions, , drop = FALSE]
rownames(meta) <- colnames(counts)
meta$dose <- factor(signif(as.numeric(meta$treatment_dose), 3))
meta$batch <- factor(meta$batch)

# 5. Save raw counts and aligned metadata (before normalization or DE testing).
saveRDS(list(counts = counts, meta = meta), output_file)
message("Saved: ", output_file)
