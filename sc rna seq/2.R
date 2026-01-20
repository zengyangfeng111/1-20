# --------------------------
# 1. Define all paths
# --------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

# 2. Get all sample files from healthy group and DFU group
health_files <- list.files(health_path, pattern = "*.csv.gz", full.names = TRUE)
dfu_files <- list.files(dfu_path, pattern = "*.csv.gz", full.names = TRUE)

cat("Number of healthy group samples: ", length(health_files), "\n")
cat("Number of DFU group samples: ", length(dfu_files), "\n")

# 3. Merge all sample files
all_files <- c(health_files, dfu_files)

# 4. Create grouping information
sample_groups <- c(rep("Normal", length(health_files)), rep("DFU", length(dfu_files)))
names(sample_groups) <- basename(all_files)

cat("\nTotal number of samples: ", length(all_files), "\n")
cat("Grouping information:\n")
print(table(sample_groups))

# --------------------------
# Step 2: Identify problematic sample (7th sample)
# --------------------------
problem_index <- 7  # According to your error message, the 7th sample has issues
problem_file <- all_files[problem_index]
problem_sample_name <- gsub(".csv.gz", "", basename(problem_file))

# Remove problematic sample from file list and create normal sample list
normal_files <- all_files[-problem_index]
normal_sample_names <- gsub(".csv.gz", "", basename(normal_files))
normal_sample_groups <- sample_groups[basename(normal_files)]

# --------------------------
# Step 3: Processing function (for normal samples)
# --------------------------
process_sample_normal <- function(file, sample_name) {
  cat("Processing sample:", sample_name, "\n")
  
  # 1. Read file
  df <- read_csv(file, show_col_types = FALSE, na = c("", "NA"))
  
  # 2. Extract gene names
  gene_names <- df[[1]]
  
  # 3. Check and handle duplicate genes (simplified version, keep first occurrence)
  if (any(duplicated(gene_names))) {
    dup_count <- sum(duplicated(gene_names))
    cat("Sample", sample_name, "has", dup_count, "duplicate genes, keeping first occurrence\n")
    
    # Keep first occurrence of genes
    keep_rows <- !duplicated(gene_names)
    df <- df[keep_rows, ]
    gene_names <- gene_names[keep_rows]
  }
  
  # 4. Extract expression matrix
  count_data <- df[, -1, drop = FALSE]
  
  # Convert to numeric matrix
  count_mat <- as.matrix(count_data)
  mode(count_mat) <- "numeric"  # Force conversion to numeric type
  
  # Set row and column names
  rownames(count_mat) <- gene_names
  colnames(count_mat) <- paste(sample_name, colnames(count_mat), sep = "_")
  
  # 5. Handle NA values
  count_mat[is.na(count_mat)] <- 0
  
  cat("Processing completed, dimensions:", dim(count_mat), "\n")
  
  return(count_mat)
}

# --------------------------
# Step 4: Process normal samples
# --------------------------
# Create a list to store matrices of all normal samples
normal_matrices <- list()

# Process each normal sample
for (i in 1:length(normal_files)) {
  sample_base <- normal_sample_names[i]
  group <- normal_sample_groups[basename(normal_files[i])]
  
  cat(sprintf("\n[%d/%d] Processing normal sample: %s (Group %s)\n", 
              i, length(normal_files), sample_base, group))
  
  tryCatch({
    mat <- process_sample_normal(normal_files[i], sample_base)
    
    # Check if matrix is valid
    if (nrow(mat) > 0 && ncol(mat) > 0) {
      normal_matrices[[sample_base]] <- mat
      cat(sprintf("Successfully processed: %d genes, %d cells\n", nrow(mat), ncol(mat)))
    } else {
      cat(sprintf("Warning: Processing result for sample %s is empty matrix\n", sample_base))
    }
    
    # Release memory
    rm(mat)
    gc()
    
  }, error = function(e) {
    cat(sprintf("Error processing sample %s: %s\n", sample_base, e$message))
  })
}

# --------------------------
# Step 5: Process problematic sample separately (using your provided code)
# --------------------------
cat("\n" + "="*50 + "\n")
cat("Start processing problematic sample separately:", problem_sample_name, "\n")
cat("="*50 + "\n")

# Process problematic sample separately
problem_sample <- all_files[problem_index]
problem_sample_name <- gsub(".csv.gz", "", basename(problem_sample))

cat("Special processing for problematic sample:", problem_sample_name, "\n")

# Try different reading methods
# Method 1: Check file structure
df_test <- read_csv(problem_sample, n_max = 10, show_col_types = FALSE)
cat("Preview of first 10 rows:\n")
print(df_test)
cat("Number of columns:", ncol(df_test), "\n")
cat("Column names:", colnames(df_test), "\n")

# Check column types
col_types <- sapply(df_test, class)
cat("Column types:\n")
print(col_types)

# Method 2: Specify all columns as numeric (except first column)
cat("Non-numeric columns found, attempting forced conversion...\n")

# Read entire file
df_full <- read_csv(problem_sample, show_col_types = FALSE)

# Convert all columns except first to numeric
for (i in 2:ncol(df_full)) {
  df_full[[i]] <- as.numeric(as.character(df_full[[i]]))
}

# Handle duplicate genes (simplified version)
gene_names <- df_full[[1]]
if (any(duplicated(gene_names))) {
  cat("Sample", problem_sample_name, "has", sum(duplicated(gene_names)), "duplicate genes, processing...\n")
  
  # Calculate row sums
  row_sums <- rowSums(df_full[, -1], na.rm = TRUE)
  df_full$row_sum <- row_sums
  
  # Group by gene, keep row with maximum sum
  df_unique <- df_full %>%
    group_by(across(1)) %>%  # Group by first column (gene names)
    slice_max(order_by = row_sum, n = 1) %>%
    ungroup() %>%
    select(-row_sum)
  
  df_full <- df_unique
}

# Create matrix
count_mat <- as.matrix(df_full[, -1])
rownames(count_mat) <- df_full[[1]]
colnames(count_mat) <- paste(problem_sample_name, colnames(count_mat), sep = "_")
count_mat[is.na(count_mat)] <- 0

# Save this specially processed matrix
problem_matrix <- count_mat

cat("Problematic sample processing completed, dimensions:", dim(problem_matrix), "\n")
cat("Data type:", class(problem_matrix[1,1]), "\n")

# --------------------------
# --------------------------
# Step 6: Save each sample matrix to disk (avoid memory overflow)
# --------------------------
cat("\n" + "="*50 + "\n")
cat("Saving each sample matrix to disk...\n")
cat("="*50 + "\n")

# Create directory to store intermediate results
temp_dir <- "temp_sample_matrices"
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

# Save normal sample matrices (without modifying original list)
normal_sample_names <- names(normal_matrices)
cat("Number of normal samples to save：", length(normal_sample_names), "\n")

for (i in 1:length(normal_sample_names)) {
  sample_name <- normal_sample_names[i]
  saveRDS(normal_matrices[[sample_name]], file = file.path(temp_dir, paste0(sample_name, ".rds")))
  cat(sprintf("Saved normal sample %d/%d: %s\n", i, length(normal_sample_names), sample_name))
  
  if (i %% 5 == 0) gc()
}

# Save problematic sample matrix
if (exists("problem_matrix")) {
  saveRDS(problem_matrix, file = file.path(temp_dir, paste0(problem_sample_name, ".rds")))
  cat(sprintf("Saved problematic sample: %s\n", problem_sample_name))
} else {
  cat("Warning: problem_matrix does not exist\n")
}

# Check which matrix objects exist
cat("\nObject list before cleanup：\n")
print(ls())

# Clear all matrices from memory
if (exists("normal_matrices")) {
  rm(normal_matrices)
  cat("Removed normal_matrices\n")
}
if (exists("problem_matrix")) {
  rm(problem_matrix)
  cat("Removed problem_matrix\n")
}
if (exists("all_matrices")) {
  rm(all_matrices)
  cat("Removed all_matrices\n")
}

gc()

cat("\nAll sample matrices have been saved to", temp_dir, "directory\n")
cat("Current memory usage after releasing memory：\n")
print(gc())
# --------------------------
# Step 7: Merge sample matrices in batches (avoid memory overflow)
# --------------------------
cat("\n" + "="*50 + "\n")
cat("Start merging sample matrices in batches...\n")
cat("="*50 + "\n")

# Get all saved sample files
sample_files <- list.files(temp_dir, pattern = "\\.rds$", full.names = TRUE)
sample_names <- gsub("\\.rds$", "", basename(sample_files))

cat("Found", length(sample_files), "sample matrix files\n")

# Batch processing parameters
batch_size <- 5  # Process 5 samples each time, adjust according to memory
n_batches <- ceiling(length(sample_files) / batch_size)

cat("Will merge in", n_batches, "batches,", batch_size, "samples per batch\n")

# Initialize empty list to store batch merge results
batch_results <- list()

# Batch processing loop
for (batch in 1:n_batches) {
  cat(sprintf("\n--- Processing batch %d/%d ---\n", batch, n_batches))
  
  # Calculate sample indices for current batch
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(sample_files))
  current_batch_files <- sample_files[start_idx:end_idx]
  current_batch_names <- sample_names[start_idx:end_idx]
  
  cat("Samples processed in this batch：", paste(current_batch_names, collapse = ", "), "\n")
  
  # Load all sample matrices for current batch
  batch_matrices <- list()
  for (i in 1:length(current_batch_files)) {
    batch_matrices[[current_batch_names[i]]] <- readRDS(current_batch_files[i])
    cat(sprintf("  Loaded sample %d/%d: %s (dimensions: %d x %d)\n", 
                i, length(current_batch_files), 
                current_batch_names[i],
                nrow(batch_matrices[[current_batch_names[i]]]),
                ncol(batch_matrices[[current_batch_names[i]]])))
  }
  
  # Merge all samples within current batch
  cat("Merging samples within current batch...\n")
  if (length(batch_matrices) == 1) {
    # If only one sample, use directly
    batch_combined <- batch_matrices[[1]]
  } else {
    # Method 1: Use fixed index without modifying list
    # Merge multiple samples
    batch_combined <- batch_matrices[[1]]
    
    # Get all matrix names
    matrix_names <- names(batch_matrices)
    
    for (i in 2:length(matrix_names)) {
      sample_name <- matrix_names[i]
      cat(sprintf("  Merging: %s\n", sample_name))
      
      # Merge using merge function
      temp_merge <- merge(
        x = batch_combined,
        y = batch_matrices[[sample_name]],
        by = "row.names",
        all = TRUE
      )
      
      rownames(temp_merge) <- temp_merge$Row.names
      batch_combined <- temp_merge[, -1]
      batch_combined[is.na(batch_combined)] <- 0
      
      # Note: Do not delete batch_matrices[[sample_name]] here
      # Release temp_merge memory
      rm(temp_merge)
      gc()
    }
    
    
  }
  
  # Save batch merge results
  batch_file <- file.path(temp_dir, sprintf("batch_%d_combined.rds", batch))
  saveRDS(batch_combined, file = batch_file)
  batch_results[[batch]] <- batch_file
  
  cat(sprintf("Batch %d merging completed, dimensions: %d x %d, saved to: %s\n", 
              batch, nrow(batch_combined), ncol(batch_combined), batch_file))
  
  # Clean up current batch matrices to release memory
  rm(batch_matrices, batch_combined)
  gc()
}

# --------------------------
# Step 8: Merge results from all batches
# --------------------------
cat("\n" + "="*50 + "\n")
cat("Merging results from all batches...\n")
cat("="*50 + "\n")

cat("Number of batch files：", length(batch_results), "\n")

# If only one batch, use directly
if (length(batch_results) == 1) {
  cat("Only one batch, loading directly...\n")
  final_combined <- readRDS(batch_results[[1]])
} else {
  # Load first batch as starting point
  cat("Loading batch 1 as starting point...\n")
  final_combined <- readRDS(batch_results[[1]])
  
  # Merge other batches step by step
  for (i in 2:length(batch_results)) {
    cat(sprintf("Merging batch %d/%d...\n", i, length(batch_results)))
    
    # Load current batch
    current_batch <- readRDS(batch_results[[i]])
    
    # Merge
    temp_merge <- merge(
      x = final_combined,
      y = current_batch,
      by = "row.names",
      all = TRUE
    )
    
    rownames(temp_merge) <- temp_merge$Row.names
    final_combined <- temp_merge[, -1]
    final_combined[is.na(final_combined)] <- 0
    
    # Clean up
    rm(current_batch, temp_merge)
    gc()
    
    # Save intermediate results periodically
    if (i %% 2 == 0) {
      saveRDS(final_combined, file = file.path(temp_dir, "intermediate_combined.rds"))
      cat(sprintf("  Saved intermediate merge results (batch %d)\n", i))
    }
  }
}


# Save final merged matrix
saveRDS(final_combined, file = "final_combined_matrix.rds")
cat("Final merged matrix saved to：final_combined_matrix.rds\n")

# --------------------------
# Step 9: Convert to numeric matrix
# --------------------------
cat("\nConverting to numeric matrix...\n")

# Check data type
if (!is.numeric(final_combined[1,1])) {
  cat("Non-numeric type detected, converting...\n")
  
  # Convert in chunks to avoid memory overflow
  n_cols <- ncol(final_combined)
  chunk_size <- 1000  # Process 1000 columns each time
  n_chunks <- ceiling(n_cols / chunk_size)
  
  cat(sprintf("Will convert in %d chunks, %d columns per chunk\n", n_chunks, chunk_size))
  
  # Create empty matrix to store results
  final_numeric <- matrix(0, nrow = nrow(final_combined), ncol = n_cols)
  rownames(final_numeric) <- rownames(final_combined)
  colnames(final_numeric) <- colnames(final_combined)
  
  for (chunk in 1:n_chunks) {
    start_col <- (chunk - 1) * chunk_size + 1
    end_col <- min(chunk * chunk_size, n_cols)
    
    cat(sprintf("  Converting chunk %d/%d (columns %d-%d)\n", chunk, n_chunks, start_col, end_col))
    
    # Convert current chunk
    chunk_data <- final_combined[, start_col:end_col, drop = FALSE]
    chunk_numeric <- apply(chunk_data, 2, function(col) as.numeric(col))
    chunk_numeric[is.na(chunk_numeric)] <- 0
    
    # Store in result matrix
    final_numeric[, start_col:end_col] <- chunk_numeric
    
    # Clean up
    rm(chunk_data, chunk_numeric)
    gc()
  }
  
  # Replace original matrix
  combined_mat_num <- final_numeric
  
} else {
  # If already numeric matrix, use directly
  cat("Matrix is already numeric type, no conversion needed\n")
  combined_mat_num <- final_combined
}

# Check and handle possible NA values
na_count <- sum(is.na(combined_mat_num))
if (na_count > 0) {
  cat("Found", na_count, "NA values, replaced with 0\n")
  combined_mat_num[is.na(combined_mat_num)] <- 0
}

cat("Final numeric matrix dimensions：", dim(combined_mat_num), "\n")
cat("Data type：", class(combined_mat_num[1,1]), "\n")

# Save numeric matrix
saveRDS(combined_mat_num, file = "final_numeric_matrix.rds")
cat("Numeric matrix saved to：final_numeric_matrix.rds\n")

# Clean up memory
rm(final_combined, final_numeric)
gc()

# --------------------------
# Step 10: Create metadata (including grouping information)
# --------------------------
cat("\nCreating metadata...\n")

# Create cell metadata
cell_ids <- colnames(combined_mat_num)

# Extract sample information from cell IDs
get_sample_from_cellname <- function(cell_id) {
  # Simple extraction logic: take first two parts as sample name
  parts <- unlist(strsplit(cell_id, "_"))
  # For most samples, first two parts are sample identifiers
  if (length(parts) >= 2) {
    return(paste(parts[1:2], collapse = "_"))
  } else {
    return(parts[1])
  }
}

sample_ids <- sapply(cell_ids, get_sample_from_cellname)

# Get grouping information based on sample name prefix
get_group_from_sample <- function(sample_id) {
  # Check if it's problematic sample
  if (grepl(problem_sample_name, sample_id)) {
    return("Normal")  # Problematic sample belongs to Normal group
  }
  
  # For other samples, judge based on filename
  for (i in 1:length(all_files)) {
    base_name <- gsub(".csv.gz", "", basename(all_files[i]))
    if (grepl(base_name, sample_id)) {
      return(sample_groups[basename(all_files[i])])
    }
  }
  
  # If no match, try to judge from sample name prefix
  if (grepl("GSM50505[0-9]{2}", sample_id)) {
    # Judge based on GSM number range
    gsm_num <- as.numeric(gsub("GSM50505", "", substr(sample_id, 1, 10)))
    if (gsm_num >= 34 && gsm_num <= 56) {
      # According to your data, this may be Normal group
      return("Normal")
    } else if (gsm_num >= 23 && gsm_num <= 33) {
      # This may be DFU group
      return("DFU")
    }
  }
  
  return("Unknown")
}

cell_groups <- sapply(sample_ids, get_group_from_sample)

# Check grouping results
cat("Grouping result statistics：\n")
print(table(cell_groups))

# Create metadata frame
metadata <- data.frame(
  row.names = cell_ids,
  cell_id = cell_ids,
  sample_id = sample_ids,
  group = cell_groups,
  stringsAsFactors = FALSE
)

# Save metadata
write.csv(metadata, file = "cell_metadata_full.csv", row.names = TRUE)
cat("Metadata saved to：cell_metadata_full.csv\n")

# --------------------------
# Step 11: Create Seurat object (with metadata)
# --------------------------
cat("\nCreating Seurat object...\n")

# Use sparse matrix to save memory
library(Matrix)

# Convert numeric matrix to sparse matrix
cat("Converting to sparse matrix...\n")
combined_mat_num_sparse <- as(combined_mat_num, "sparseMatrix")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = combined_mat_num_sparse,
  project = "GSE165816_full_20samples",
  meta.data = metadata,
  min.cells = 3,
  min.features = 200
)

# View Seurat object information
cat("\nSeurat object information：\n")
print(seurat_obj)

# View metadata
cat("\nMetadata column names：\n")
print(colnames(seurat_obj@meta.data))

# View grouping statistics
cat("\nFinal grouping statistics：\n")
group_table <- table(seurat_obj$group)
print(group_table)

# --------------------------
# Step 12: Save final results
# --------------------------
cat("\nSaving final results...\n")

# Save Seurat object
saveRDS(seurat_obj, file = "seurat_obj_full_20samples_preQC.rds")
cat("Seurat object saved to：seurat_obj_full_20samples_preQC.rds\n")

# Save statistical information
stats <- list(
  total_samples = length(all_files),
  normal_samples = length(health_files),
  dfu_samples = length(dfu_files),
  total_genes = nrow(seurat_obj),
  total_cells = ncol(seurat_obj),
  normal_cells = sum(seurat_obj$group == "Normal"),
  dfu_cells = sum(seurat_obj$group == "DFU"),
  unknown_cells = sum(seurat_obj$group == "Unknown"),
  problem_sample = problem_sample_name,
  problem_sample_cells = sum(grepl(problem_sample_name, seurat_obj$sample_id)),
  memory_usage = format(object.size(seurat_obj), units = "MB")
)

capture.output(print(stats), file = "final_processing_stats.txt")
cat("Statistical information saved to：final_processing_stats.txt\n")

# Clean up temporary files (optional)
cat("\nCleaning up temporary files...\n")
if (dir.exists(temp_dir)) {
  # Can choose to keep or delete temporary files
  # If disk space is sufficient, it's recommended to keep for debugging
  cat("Temporary files are saved in：", temp_dir, "\n")
  cat("To clean up, please manually delete this directory\n")
}

# --------------------------
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

cat("\nNormalization completed\n")

# ============================================================================
# Single-cell RNA-seq Analysis - Complete Restart Pipeline
# Starting from saved intermediate files
# ============================================================================

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(Matrix)

# Set graphic parameters
options(stringsAsFactors = FALSE)
theme_set(theme_classic(base_size = 12))

# --------------------------
# 1. Check and load data
# --------------------------
cat(strrep("=", 70), "\n")
cat("Step 1: Check and load data\n")
#

# Check if files exist
files_to_check <- c(
  "final_numeric_matrix.rds",
  "cell_metadata_full.csv",
  "seurat_obj_full_20samples_preQC.rds"
)

for (file in files_to_check) {
  if (file.exists(file)) {
    cat(sprintf("✓ File exists: %s\n", file))
  } else {
    cat(sprintf("✗ File missing: %s\n", file))
  }
}

# Load Seurat object (most convenient way)
cat("\nLoading Seurat object...\n")
seurat_obj <- readRDS("seurat_obj_full_20samples_preQC.rds")

cat("\nBasic information of Seurat object:\n")
print(seurat_obj)

# View grouping statistics
if ("group" %in% colnames(seurat_obj@meta.data)) {
  cat("\nGroup statistics:\n")
  group_stats <- table(seurat_obj$group)
  print(group_stats)
  
  # Plot grouped bar chart
  p_group <- ggplot(data.frame(group = names(group_stats), count = as.numeric(group_stats)),
                    aes(x = group, y = count, fill = group)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.5) +
    scale_fill_manual(values = c("Normal" = "#2E86AB", "DFU" = "#A23B72", "Unknown" = "#8A8D91")) +
    labs(title = "Cell group distribution", x = "Group", y = "Number of cells") +
    theme(legend.position = "none")
  
  ggsave("group_distribution.png", p_group, width = 8, height = 6, dpi = 300)
}

# Check if quality control metrics already exist
cat("\nChecking existing quality control metrics:\n")
qc_columns <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
existing_qc <- qc_columns[qc_columns %in% colnames(seurat_obj@meta.data)]
cat("Existing QC metrics:", paste(existing_qc, collapse = ", "), "\n")

# --------------------------
# 2. Calculate mitochondrial gene percentage (if not already calculated)
# --------------------------
cat("Step 2: Calculate mitochondrial gene percentage\n")
#

if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  cat("Calculating mitochondrial gene percentage...\n")
  
  # Try multiple mitochondrial gene patterns
  mt_patterns <- c("^MT-", "^mt-", "^Mt-")
  mt_genes_found <- FALSE
  
  for (pattern in mt_patterns) {
    mt_genes <- rownames(seurat_obj)[grep(pattern, rownames(seurat_obj))]
    if (length(mt_genes) > 0) {
      cat(sprintf("Found %d mitochondrial genes using pattern '%s'\n", pattern, length(mt_genes)))
      cat("First 10 mitochondrial genes:", head(mt_genes, 10), "...\n")
      
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = pattern)
      mt_genes_found <- TRUE
      break
    }
  }
  
  if (!mt_genes_found) {
    cat("Warning: Failed to automatically detect mitochondrial genes\n")
    cat("Possible reasons: mitochondrial genes use other naming conventions, or data does not contain mitochondrial genes\n")
    
    # Manually check common mitochondrial genes
    common_mt_genes <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", 
                         "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6")
    
    existing_mt <- common_mt_genes[common_mt_genes %in% rownames(seurat_obj)]
    
    if (length(existing_mt) > 0) {
      cat("Found the following mitochondrial genes:", paste(existing_mt, collapse = ", "), "\n")
      # Calculate percentage using these genes
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features = existing_mt)
    } else {
      cat("No common mitochondrial genes found, setting percent.mt to 0\n")
      seurat_obj[["percent.mt"]] <- 0
    }
  }
  
  # View mitochondrial gene percentage statistics
  cat("\nMitochondrial gene percentage statistics:\n")
  print(summary(seurat_obj$percent.mt))
} else {
  cat("Mitochondrial gene percentage already exists, skipping calculation\n")
  cat("Statistical information:\n")
  print(summary(seurat_obj$percent.mt))
}

# --------------------------
# 3. Quality control and filtering
# --------------------------
cat("Step 3: Data quality control and filtering\n")
#

cat("Number of cells before filtering:", ncol(seurat_obj), "\n")
cat("Number of genes before filtering:", nrow(seurat_obj), "\n")

# View current QC metric distribution
cat("\nCurrent QC metric statistics:\n")
qc_summary <- data.frame(
  Metric = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  Mean = c(mean(seurat_obj$nFeature_RNA), 
           mean(seurat_obj$nCount_RNA), 
           mean(seurat_obj$percent.mt)),
  Median = c(median(seurat_obj$nFeature_RNA), 
             median(seurat_obj$nCount_RNA), 
             median(seurat_obj$percent.mt)),
  Min = c(min(seurat_obj$nFeature_RNA), 
          min(seurat_obj$nCount_RNA), 
          min(seurat_obj$percent.mt)),
  Max = c(max(seurat_obj$nFeature_RNA), 
          max(seurat_obj$nCount_RNA), 
          max(seurat_obj$percent.mt))
)
print(qc_summary)

# Plot QC plots before filtering
cat("\nGenerating QC plots before filtering...\n")

# Function: plot QC violin plot
plot_qc_violin <- function(seurat_obj, title) {
  p <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1,
    group.by = "group"
  ) + 
    ggtitle(title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  return(p)
}

# Function: plot QC scatter plot
plot_qc_scatter <- function(seurat_obj, title) {
  # Create data frame
  qc_data <- data.frame(
    nFeature_RNA = seurat_obj$nFeature_RNA,
    nCount_RNA = seurat_obj$nCount_RNA,
    percent.mt = seurat_obj$percent.mt,
    group = seurat_obj$group
  )
  
  p1 <- ggplot(qc_data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_gradient(low = "blue", high = "red", name = "MT%") +
    labs(title = "Number of cells vs Number of genes", x = "UMI count", y = "Number of genes") +
    theme_minimal()
  
  p2 <- ggplot(qc_data, aes(x = nCount_RNA, y = percent.mt, color = group)) +
    geom_point(alpha = 0.5, size = 0.5) +
    scale_color_manual(values = c("Normal" = "#2E86AB", "DFU" = "#A23B72", "Unknown" = "#8A8D91")) +
    labs(title = "Number of cells vs Mitochondrial%", x = "UMI count", y = "Mitochondrial gene percentage") +
    theme_minimal()
  
  combined <- p1 + p2 + plot_layout(ncol = 2)
  return(combined)
}

# Generate and save QC plots before filtering
p_qc_before_violin <- plot_qc_violin(seurat_obj, "QC metric distribution before filtering")
ggsave("qc_before_filtering_violin.png", p_qc_before_violin, width = 14, height = 6, dpi = 300)

p_qc_before_scatter <- plot_qc_scatter(seurat_obj, "QC scatter plot before filtering")
ggsave("qc_before_filtering_scatter.png", p_qc_before_scatter, width = 14, height = 6, dpi = 300)

# Set filtering thresholds
cat("\nSetting filtering thresholds:\n")
cat("  nFeature_RNA: 200-6000\n")
cat("  nCount_RNA: 500-8000\n")
cat("  percent.mt: <10%\n")

# Perform filtering
cat("\nPerforming filtering...\n")
seurat_obj_filtered <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    nCount_RNA >= 500 &
    nCount_RNA <= 8000 &
    percent.mt < 10
)

cat("Number of cells after filtering:", ncol(seurat_obj_filtered), "\n")
cat("Number of genes after filtering:", nrow(seurat_obj_filtered), "\n")
cat("Percentage of cells filtered out:", round((1 - ncol(seurat_obj_filtered)/ncol(seurat_obj)) * 100, 1), "%\n")

# Generate QC plots after filtering
cat("\nGenerating QC plots after filtering...\n")
p_qc_after_violin <- plot_qc_violin(seurat_obj_filtered, "QC metric distribution after filtering")
ggsave("qc_after_filtering_violin.png", p_qc_after_violin, width = 14, height = 6, dpi = 300)

p_qc_after_scatter <- plot_qc_scatter(seurat_obj_filtered, "QC scatter plot after filtering")
ggsave("qc_after_filtering_scatter.png", p_qc_after_scatter, width = 14, height = 6, dpi = 300)

# Save filtered object
saveRDS(seurat_obj_filtered, file = "seurat_obj_postQC.rds")
cat("\n✓ Seurat object after QC saved to: seurat_obj_postQC.rds\n")

# Update object
seurat_obj <- seurat_obj_filtered
rm(seurat_obj_filtered)
gc()

# --------------------------
# 4. Identify and remove mitochondrial genes
# --------------------------
# Supplementary Step 4: Identify and remove mitochondrial genes (corrected version)

# --------------------------
# 4. Identify and remove mitochondrial genes
# --------------------------
cat(strrep("-", 70), "\n", sep="")
cat("Step 4: Identify and remove mitochondrial genes\n")
cat(strrep("-", 70), "\n", sep="")

# Identify mitochondrial genes
cat("Identifying mitochondrial genes...\n")

# Try multiple mitochondrial gene naming patterns
mt_patterns <- c("^MT-", "^mt-", "^Mt-", "^Mito-", "^mito-", "^MT\\.", "^mt\\.")
mt_genes <- character(0)

for (pattern in mt_patterns) {
  genes <- rownames(seurat_obj)[grep(pattern, rownames(seurat_obj))]
  if (length(genes) > 0) {
    cat(sprintf("Found %d mitochondrial genes using pattern '%s'\n", pattern, length(genes)))
    mt_genes <- union(mt_genes, genes)
  }
}

# If not found by pattern matching, try common mitochondrial gene list
if (length(mt_genes) == 0) {
  cat("Failed to find mitochondrial genes by pattern matching, trying common mitochondrial gene list...\n")
  
  common_mt_genes <- c(
    # Human mitochondrial genes
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8",
    "MT-CYB", "MT-RNR1", "MT-RNR2", "MT-TA", "MT-TC", "MT-TD", "MT-TE",
    "MT-TF", "MT-TG", "MT-TH", "MT-TI", "MT-TK", "MT-TL1", "MT-TL2",
    "MT-TM", "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1", "MT-TS2",
    "MT-TT", "MT-TV", "MT-TW", "MT-TY",
    # Mouse mitochondrial genes
    "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6",
    "mt-Co1", "mt-Co2", "mt-Co3", "mt-Atp6", "mt-Atp8",
    "mt-Cytb", "mt-Rnr1", "mt-Rnr2"
  )
  
  mt_genes <- common_mt_genes[common_mt_genes %in% rownames(seurat_obj)]
  
  if (length(mt_genes) > 0) {
    cat(sprintf("Found %d mitochondrial genes from common gene list\n", length(mt_genes)))
  }
}

# Statistical results
if (length(mt_genes) > 0) {
  cat(sprintf("\nTotal of %d mitochondrial genes found\n", length(mt_genes)))
  
  # Display first 20 mitochondrial genes
  if (length(mt_genes) <= 20) {
    cat("Mitochondrial gene list:", paste(mt_genes, collapse = ", "), "\n")
  } else {
    cat("First 20 mitochondrial genes:", paste(mt_genes[1:20], collapse = ", "), "...\n")
  }
  
  # Calculate expression ratio of mitochondrial genes
  counts_matrix <- GetAssayData(seurat_obj, slot = "counts")
  mt_expression <- Matrix::rowSums(counts_matrix[mt_genes, , drop = FALSE])
  total_expression <- Matrix::colSums(counts_matrix)
  
  # Calculate ratio of mitochondrial gene expression to total expression
  mt_percent_expression <- sum(mt_expression) / sum(total_expression) * 100
  
  cat(sprintf("Ratio of mitochondrial gene expression to total expression: %.2f%%\n", mt_percent_expression))
  
  # Plot mitochondrial gene expression heatmap (if not too many)
  if (length(mt_genes) <= 50 && length(mt_genes) > 0) {
    cat("\nGenerating mitochondrial gene expression heatmap...\n")
    
    # Randomly select some cells to display to avoid oversized heatmap
    set.seed(123)
    if (ncol(seurat_obj) > 200) {
      cells_to_plot <- sample(colnames(seurat_obj), 200)
    } else {
      cells_to_plot <- colnames(seurat_obj)
    }
    
    # Check if there are enough cells
    if (length(cells_to_plot) > 1 && length(mt_genes) > 1) {
      # Create temporary object for heatmap plotting
      temp_obj <- subset(seurat_obj, cells = cells_to_plot)
      
      # Check if normalized data already exists
      if (!"RNA" %in% Assays(temp_obj) || length(GetAssayData(temp_obj, slot = "data")) == 0) {
        cat("Normalized data not found, normalizing subset...\n")
        temp_obj <- NormalizeData(temp_obj, normalization.method = "LogNormalize", scale.factor = 10000)
      }
      
      # Create heatmap
      tryCatch({
        p_heatmap <- DoHeatmap(
          temp_obj,
          features = mt_genes,
          slot = "data",
          group.by = "group",
          angle = 45,
          size = 3
        ) +
          scale_fill_gradientn(colors = c("blue", "white", "red")) +
          ggtitle("Mitochondrial gene expression heatmap") +
          theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
          )
        
        ggsave("mitochondrial_genes_heatmap.png", p_heatmap, width = 12, height = 8, dpi = 300)
        cat("✓ Mitochondrial gene expression heatmap saved to: mitochondrial_genes_heatmap.png\n")
      }, error = function(e) {
        cat("Error generating heatmap:", e$message, "\n")
        cat("Trying to use count data instead of normalized data...\n")
        
        # Try using count data
        p_heatmap <- DoHeatmap(
          temp_obj,
          features = mt_genes,
          slot = "counts",
          group.by = "group",
          angle = 45,
          size = 3
        ) +
          scale_fill_gradientn(colors = c("blue", "white", "red")) +
          ggtitle("Mitochondrial gene expression heatmap (raw counts)") +
          theme(
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
          )
        
        ggsave("mitochondrial_genes_heatmap_counts.png", p_heatmap, width = 12, height = 8, dpi = 300)
        cat("✓ Mitochondrial gene expression heatmap (raw counts) saved to: mitochondrial_genes_heatmap_counts.png\n")
      })
      
      # Also create simple heatmap
      cat("\nGenerating simplified mitochondrial gene expression heatmap...\n")
      
      # Extract expression matrix
      expr_matrix <- as.matrix(GetAssayData(temp_obj, slot = "counts")[mt_genes, ])
      
      # Create heatmap
      pheatmap::pheatmap(
        expr_matrix,
        scale = "row",  # Normalize by row
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        main = "Mitochondrial gene expression heatmap",
        filename = "mitochondrial_genes_heatmap_simple.png",
        width = 10,
        height = 8
      )
      cat("✓ Simplified mitochondrial gene expression heatmap saved to: mitochondrial_genes_heatmap_simple.png\n")
      
      rm(temp_obj)
    } else {
      cat("Too few cells or mitochondrial genes, skipping heatmap generation\n")
    }
  } else if (length(mt_genes) > 50) {
    cat("Too many mitochondrial genes (>50), skipping heatmap generation\n")
    cat("Generating heatmap for top 50 mitochondrial genes...\n")
    
    # Select top 50 mitochondrial genes by expression level
    mt_means <- rowMeans(GetAssayData(seurat_obj, slot = "counts")[mt_genes, ])
    top_mt_genes <- names(sort(mt_means, decreasing = TRUE))[1:min(50, length(mt_means))]
    
    # Create heatmap
    tryCatch({
      # Create heatmap using pheatmap
      expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "counts")[top_mt_genes, sample(colnames(seurat_obj), min(100, ncol(seurat_obj)))])
      
      pheatmap::pheatmap(
        expr_matrix,
        scale = "row",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        main = "Top 50 Mitochondrial gene expression heatmap",
        filename = "mitochondrial_genes_top50_heatmap.png",
        width = 12,
        height = 10
      )
      cat("✓ Heatmap of top 50 mitochondrial genes saved to: mitochondrial_genes_top50_heatmap.png\n")
    }, error = function(e) {
      cat("Error generating heatmap:", e$message, "\n")
    })
  }
  
  # Remove mitochondrial genes
  cat("\nRemoving mitochondrial genes...\n")
  
  # Create backup (including mitochondrial genes)
  seurat_obj_with_mt <- seurat_obj
  
  # Remove mitochondrial genes
  non_mt_genes <- setdiff(rownames(seurat_obj), mt_genes)
  seurat_obj <- subset(seurat_obj, features = non_mt_genes)
  
  cat(sprintf("Removed %d mitochondrial genes\n", length(mt_genes)))
  cat(sprintf("Number of remaining genes after removing mitochondrial genes: %d\n", nrow(seurat_obj)))
  cat(sprintf("Number of cells after removing mitochondrial genes: %d\n", ncol(seurat_obj)))
  
  # Check mitochondrial percentage after removing mitochondrial genes
  # Recalculate mitochondrial percentage (should be 0 or very close to 0)
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$percent.mt <- 0  # Because mitochondrial genes have been removed
    cat("Note: Mitochondrial genes have been removed, percent.mt set to 0\n")
  }
  
  # Save mitochondrial gene list
  writeLines(mt_genes, "mitochondrial_genes_list.txt")
  cat("✓ Mitochondrial gene list saved to: mitochondrial_genes_list.txt\n")
  
  # Save object after removing mitochondrial genes
  saveRDS(seurat_obj, file = "seurat_obj_noMT.rds")
  cat("✓ Seurat object after removing mitochondrial genes saved to: seurat_obj_noMT.rds\n")
  
  # Optional: Save object containing mitochondrial genes
  saveRDS(seurat_obj_with_mt, file = "seurat_obj_with_MT.rds")
  cat("✓ Seurat object containing mitochondrial genes saved to: seurat_obj_with_MT.rds\n")
  
  # Clean up memory
  rm(seurat_obj_with_mt)
  
} else {
  cat("No mitochondrial genes found, no need to remove\n")
  
  # If percent.mt exists but is non-zero, check possible reasons
  if ("percent.mt" %in% colnames(seurat_obj@meta.data) && 
      any(seurat_obj$percent.mt > 0)) {
    cat("\nWarning: percent.mt > 0 but no mitochondrial genes detected\n")
    cat("Possible reasons:\n")
    cat("1. Mitochondrial genes use non-standard naming conventions\n")
    cat("2. percent.mt was calculated by other methods\n")
    cat("3. Data may contain mitochondrial pseudogenes in nuclear genome\n")
    
    # Display percent.mt statistics
    cat("\nCurrent percent.mt statistics:\n")
    print(summary(seurat_obj$percent.mt))
  }
}

# Run garbage collection to free memory
gc()

cat("\n", strrep("=", 70), "\n", sep="")
cat("Step 4 completed: Mitochondrial genes identified and removed\n")
cat(strrep("=", 70), "\n", sep="")

# Display current object information
cat("\nCurrent Seurat object information:\n")
print(seurat_obj)
cat(sprintf("Number of genes: %d\n", nrow(seurat_obj)))
cat(sprintf("Number of cells: %d\n", ncol(seurat_obj)))
# 5. Data normalization
# --------------------------
cat("Step 5: Data normalization\n")
#

cat("Performing LogNormalize normalization...\n")
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("Normalization completed\n")

# --------------------------
# 6. Identify highly variable genes
# --------------------------
cat("Step 6: Identify highly variable genes\n")
#

cat("Identifying highly variable genes using vst method...\n")
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

# Get highly variable genes
var_features <- VariableFeatures(seurat_obj)
cat(sprintf("Identified %d highly variable genes\n", length(var_features)))

# Display first 10 highly variable genes
cat("\nFirst 10 highly variable genes:\n")
print(head(var_features, 10))

# Remove possible residual mitochondrial genes from highly variable genes
if (exists("all_mt_genes") && length(all_mt_genes) > 0) {
  var_features <- var_features[!var_features %in% all_mt_genes]
  VariableFeatures(seurat_obj) <- var_features
  cat(sprintf("Removed mitochondrial genes from highly variable genes, %d highly variable genes remaining\n", length(var_features)))
}

# Plot highly variable genes
cat("\nGenerating highly variable genes plot...\n")
p_var_features <- VariableFeaturePlot(seurat_obj)
top_genes <- head(var_features, 10)
p_labeled <- LabelPoints(plot = p_var_features, points = top_genes, repel = TRUE)
ggsave("variable_features_labeled.png", p_labeled, width = 10, height = 8, dpi = 300)

# Save current state
saveRDS(seurat_obj, file = "seurat_obj_normalized.rds")
cat("✓ Normalized Seurat object saved to: seurat_obj_normalized.rds\n")

# --------------------------
# 7. Data scaling (batch processing to avoid memory issues)
# --------------------------
#
cat("Step 7: Data scaling (batch processing)\n")
#

# Check memory usage
obj_size <- format(object.size(seurat_obj), units = "MB")
cat("Current object size:", obj_size, "\n")

# Get highly variable genes
var_features <- VariableFeatures(seurat_obj)
cat("Number of highly variable genes to scale:", length(var_features), "\n")

# Batch scaling parameters
batch_size <- 300  # Adjust based on memory, smaller is more memory-efficient but slower
n_batches <- ceiling(length(var_features) / batch_size)

cat(sprintf("Will scale in %d batches, %d genes per batch\n", n_batches, batch_size))

# Batch scaling
for (batch in 1:n_batches) {
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, length(var_features))
  batch_genes <- var_features[start_idx:end_idx]
  
  cat(sprintf("Scaling batch %d/%d (genes %d-%d)\n", 
              batch, n_batches, start_idx, end_idx))
  
  seurat_obj <- ScaleData(
    seurat_obj,
    features = batch_genes,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
}

cat("Data scaling completed\n")

# View scaled data
cat("\nPreview of scaled data (first 5 genes and first 5 cells):\n")
scaled_data <- GetAssayData(seurat_obj, slot = "scale.data")
if (!is.null(scaled_data)) {
  print(scaled_data[1:5, 1:5])
} else {
  cat("Scaled data is empty, may need to check scaling step\n")
}

# Save scaled object
saveRDS(seurat_obj, file = "seurat_obj_scaled.rds")
cat("✓ Scaled Seurat object saved to: seurat_obj_scaled.rds\n")

# --------------------------
# 8. PCA analysis
# --------------------------
#
cat("Step 8: PCA analysis\n")
#

cat("Running PCA analysis...\n")
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = 50,
  verbose = FALSE
)

cat("PCA completed\n")
cat("Available dimensionality reduction methods:", names(seurat_obj@reductions), "\n")

# Plot PCA results
cat("\nGenerating PCA visualization plots...\n")

# PCA colored by group
p_pca_group <- DimPlot(seurat_obj, reduction = "pca", group.by = "group") +
  ggtitle("PCA - Colored by group") +
  theme(plot.title = element_text(hjust = 0.5))

# PCA colored by features
p_pca_features <- FeaturePlot(seurat_obj, 
                              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                              reduction = "pca",
                              ncol = 3) &
  theme(plot.title = element_text(size = 10))

# Combined plot
p_pca_combined <- p_pca_group / p_pca_features + plot_layout(heights = c(1, 0.8))
ggsave("pca_visualization.png", p_pca_combined, width = 12, height = 10, dpi = 300)

# Elbow plot
p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
ggsave("elbow_plot.png", p_elbow, width = 8, height = 6, dpi = 300)

# Check variance explained by PCs
cat("\nVariance explained by first 10 principal components:\n")
pca_variance <- seurat_obj@reductions$pca@stdev^2
pca_variance_percent <- pca_variance / sum(pca_variance) * 100
print(data.frame(
  PC = 1:10,
  Variance = pca_variance[1:10],
  Percent_Variance = round(pca_variance_percent[1:10], 2)
))

# Save object after PCA
saveRDS(seurat_obj, file = "seurat_obj_pca.rds")
cat("✓ PCA analysis completed, object saved to: seurat_obj_pca.rds\n")

# --------------------------
# 9. Cell clustering
# --------------------------
#
cat("Step 9: Cell clustering\n")
#

# Select number of PCs based on elbow plot
n_pcs <- 14  # Can adjust based on elbow plot
cat(sprintf("Using first %d principal components for clustering\n", n_pcs))

# Find neighbors
cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:n_pcs,
  verbose = FALSE
)

# Cluster
cat("Performing clustering...\n")
# Try multiple resolutions
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = c(0.2, 0.5, 0.8),
  verbose = FALSE
)

# View number of clusters at different resolutions
cat("\nNumber of clusters at different resolutions:\n")
for (res in c(0.2, 0.5, 0.8)) {
  cluster_col <- paste0("RNA_snn_res.", res)
  if (cluster_col %in% colnames(seurat_obj@meta.data)) {
    n_clusters <- length(unique(seurat_obj@meta.data[[cluster_col]]))
    cat(sprintf("Resolution %.1f: %d clusters\n", res, n_clusters))
  }
}

# Select a resolution (0.2 selected here, can adjust as needed)
selected_res <- 0.5
Idents(seurat_obj) <- paste0("RNA_snn_res.", selected_res)

cat(sprintf("\nSelected resolution %.1f, resulting in %d clusters\n", 
            selected_res, length(unique(Idents(seurat_obj)))))

# View cluster sizes
cluster_sizes <- table(Idents(seurat_obj))
cat("\nNumber of cells per cluster:\n")
print(cluster_sizes)

# --------------------------
# 10. t-SNE dimensionality reduction
# --------------------------
#
cat("Step 10: t-SNE dimensionality reduction\n")
#

cat("Running t-SNE...\n")
set.seed(123)  # Set random seed for reproducibility
cat("\nChecking for duplicate points in PCA embedding...\n")

# Get PCA embedding
pca_embedding <- Embeddings(seurat_obj, "pca")

# Check for duplicate rows
duplicate_rows <- duplicated(pca_embedding)
num_duplicates <- sum(duplicate_rows)

if (num_duplicates > 0) {
  cat(sprintf("Found %d duplicate PCA points, adding small noise...\n", num_duplicates))
  
  # Add small random noise to all points (avoid identical points)
  set.seed(123)  # Set random seed for reproducibility
  noise <- matrix(rnorm(nrow(pca_embedding) * ncol(pca_embedding), 
                        mean = 0, sd = 1e-10), 
                  nrow = nrow(pca_embedding), ncol = ncol(pca_embedding))
  
  # Only add noise to duplicate points (or add to all points to ensure no duplicates)
  # Here we simply add small noise to all points
  pca_embedding <- pca_embedding + noise
  
  # Put modified PCA embedding back into object
  seurat_obj[["pca"]] <- CreateDimReducObject(
    embeddings = pca_embedding,
    loadings = Loadings(seurat_obj[["pca"]]),
    key = "PC_",
    assay = "RNA"
  )
}

# Rerun t-SNE
seurat_obj <- RunTSNE(
  seurat_obj,
  dims = 1:n_pcs,
  perplexity = 30
)

cat("t-SNE completed\n")

cat("t-SNE completed\n")

# Plot t-SNE plots
cat("\nGenerating t-SNE visualization plots...\n")

# t-SNE colored by cluster
p_tsne_cluster <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) +
  ggtitle("t-SNE - Colored by cluster") +
  theme(plot.title = element_text(hjust = 0.5))

# t-SNE colored by group
p_tsne_group <- DimPlot(seurat_obj, reduction = "tsne", group.by = "group") +
  ggtitle("t-SNE - Colored by group") +
  theme(plot.title = element_text(hjust = 0.5))

# Combined plot
p_tsne_combined <- p_tsne_cluster + p_tsne_group
ggsave("tsne_visualization.png", p_tsne_combined, width = 14, height = 6, dpi = 300)

# --------------------------
# 11. Save final results
# --------------------------
#
cat("Step 11: Save final results\n")
#

# Save complete Seurat object
saveRDS(seurat_obj, file = "seurat_analysis_complete.rds")
cat("✓ Complete Seurat object saved to: seurat_analysis_complete.rds\n")

# Save metadata
write.csv(seurat_obj@meta.data, file = "final_cell_metadata.csv", row.names = TRUE)
cat("✓ Cell metadata saved to: final_cell_metadata.csv\n")

# Save clustering information
cluster_info <- data.frame(
  Cell = colnames(seurat_obj),
  Cluster = Idents(seurat_obj),
  Group = seurat_obj$group,
  Sample = seurat_obj$sample_id
)

write.csv(cluster_info, file = "cluster_assignments.csv", row.names = FALSE)
cat("✓ Cluster assignment information saved to: cluster_assignments.csv\n")

# Generate analysis report
cat("\nGenerating analysis report...\n")
analysis_report <- data.frame(
  Metric = c(
    "Original number of cells",
    "Number of cells after QC",
    "Original number of genes",
    "Number of genes after removing mitochondrial genes",
    "Number of highly variable genes",
    "Number of PCA principal components",
    "Number of clusters",
    "Number of Normal group cells",
    "Number of DFU group cells",
    "Number of Unknown group cells"
  ),
  Value = c(
    "Obtained from saved object",  # Original number of cells
    ncol(seurat_obj),    # Number of cells after QC
    "Obtained from saved object",  # Original number of genes
    nrow(seurat_obj),    # Number of genes after removing mitochondrial genes
    length(VariableFeatures(seurat_obj)),  # Number of highly variable genes
    n_pcs,               # Number of PCA principal components
    length(unique(Idents(seurat_obj))),  # Number of clusters
    sum(seurat_obj$group == "Normal"),   # Number of Normal group cells
    sum(seurat_obj$group == "DFU"),      # Number of DFU group cells
    sum(seurat_obj$group == "Unknown")   # Number of Unknown group cells
  )
)

write.csv(analysis_report, file = "analysis_summary.csv", row.names = FALSE)
cat("✓ Analysis summary saved to: analysis_summary.csv\n")

# Print final summary
cat("Analysis completed!\n")
#

# Clean up memory
rm(list = setdiff(ls(), "seurat_obj"))
gc()

cat("\n✓ All steps completed! You can now proceed with downstream analysis.\n")
# --------------------------
# 11. UMAP dimensionality reduction (for cell annotation)
# --------------------------
cat("Step 11: UMAP dimensionality reduction\n")
#
n_pcs <- 14  # Can adjust based on elbow plot

cat("Running UMAP...\n")
set.seed(123)
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:n_pcs,
  verbose = FALSE
)

cat("UMAP completed\n")

# Plot UMAP plots
cat("\nGenerating UMAP visualization plots...\n")

# UMAP colored by cluster
p_umap_cluster <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP - Colored by cluster") +
  theme(plot.title = element_text(hjust = 0.5))

# UMAP colored by group
p_umap_group <- DimPlot(seurat_obj, reduction = "umap", group.by = "group", pt.size = 0.5) +
  ggtitle("UMAP - Colored by group") +
  theme(plot.title = element_text(hjust = 0.5))

# Combined plot
p_umap_combined <- p_umap_cluster + p_umap_group
ggsave("umap_visualization.png", p_umap_combined, width = 14, height = 6, dpi = 300)

# Save object containing UMAP
saveRDS(seurat_obj, file = "seurat_obj_with_umap.rds")
cat("✓ UMAP analysis completed, object saved to: seurat_obj_with_umap.rds\n")

# Save seurat_obj to .RData file (can specify path, e.g., "~/data/seurat_obj.RData")
save(seurat_obj, file = "seurat_obj.RData")

# --------------------------
# 12. Cell type annotation
# --------------------------
cat("Step 12: Cell type annotation\n")
#

# Define marker genes for common cell types
marker_genes <- list(
  "T_cells" = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B"),
  "B_cells" = c("CD19", "CD79A", "MS4A1", "CD79B"),
  "NK_cells" = c("NKG7", "GNLY", "NCR1", "KLRD1"),
  "Monocytes_Macrophages" = c("CD14", "CD68", "CD163", "FCGR3A", "ITGAM"),
  "Dendritic_cells" = c("CD1C", "FCER1A", "CLEC9A", "CD209"),
  "Neutrophils" = c("FCGR3B", "CSF3R", "S100A8", "S100A9"),
  "Erythrocytes" = c("HBB", "HBA1", "HBA2"),
  "Platelets" = c("PPBP", "PF4", "GP9"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  "Fibroblasts" = c("COL1A1", "COL3A1", "DCN", "LUM", "PDGFRA"),
  "Epithelial" = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  "Keratinocytes" = c("KRT5", "KRT14", "KRT10", "KRT1"),
  "Melanocytes" = c("MLANA", "PMEL", "TYR", "DCT"),
  "Mast_cells" = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
  "Epithelial_cells" = c("CDH1", "TPSB2", "CPA3", "KIT"),
  
  "Plasma_cells" = c("SDC1", "CD38", "MZB1", "IGHA1", "IGHG1")
)
marker_genes <- list(
  # ========== T cells and subsets ==========
  "T_NK_cell" = c("CD2", "CD3D", "CD4", "CD8A", "CD3E", "CD7", "CD52", "CD3G"),
  "Tfh_cell" = c("B3GAT1", "CD57", "BTLA", "CTLA4", "CXCL13", "FAS", "CD95", "IL2RA", "CD25", "IL21", "SH2D1A", "SAP", "TNFRSF4", "OX40"),
  "Naive_CD4_T" = c("CD3D", "CD4", "CCR7", "SELL"),
  "Active_CD4_T" = c("CD3D", "CD4", "LTB", "GATA3"),
  "cTh1" = c("CD3D", "CD4", "CXCR3", "TBX21", "STAT4", "IFNG"),
  "ex_cTh17" = c("CD3D", "CD4", "CXCR3", "TBX21", "STAT4", "IFNG"),
  "Treg" = c("IL2RA", "FOXP3", "IKZF2"),
  "Naive_CD8_T" = c("CD3D", "CD8A", "CD8B", "IL7R", "CCR7"),
  "CD8_gdT_Tscm" = c("CD3D", "CD8A", "CD8B", "TCF7", "LEF1", "CXCR3", "XCL1"),
  "CX3CR1_CD8Teff" = c("CD3D", "CD8A", "CD8B", "CX3CR1", "GZMA", "GZMB"),
  "CD8Trm" = c("CD3D", "CD8A", "CD8B", "CD69", "CD103", "IFNG", "GZMA", "GZMM"),
  "Mitotic_CD8_T" = c("CD3D", "CD8A", "CD8B", "CD69", "MKI67"),
  "gdT_1" = c("CD3D", "CD8A", "FCER1G", "CXCR3"),
  "gdT_2" = c("CD3D", "PDCD1", "GZMM"),
  "CX3CR1_gdT" = c("CD3D", "TRDV2", "TRGV9", "CD27", "CX3CR1"),
  "ex_gdT17" = c("CD3D", "CD8A", "TRDV2", "TRGV9", "IL23R", "CCR6", "IL7R", "TBX21", "STAT4", "EOMES"),
  "MAIT_1" = c("CD3D", "CD8A"),
  "MAIT_2" = c("CD3D", "CD8A", "TRAV1-2", "RORC", "RORA", "IL23R", "CCR6"),
  "ex_MAIT17" = c("CD3D", "TRAV1-2", "IL23R", "CCR6", "IL7R", "TBX21", "STAT4", "EOMES"),
  "MAIT" = c("PLZF", "MR1", "ZBTB16", "CD103", "CD69"),
  
  # ========== NK cells and subsets ==========
  "NK_cell" = c("CD16", "FCGR3A", "CD56", "NCAM1", "NKG7", "KLRD1", "KLRF1"),
  "CD56dim_CD16hi_NK" = c("CD56", "NCAM1", "CD16"),
  "cir_NK" = c("KLRD1", "NCAM1", "CX3CR1"),
  "ILC1_like_NK" = c("IL7R", "TCF7", "SELL", "KLRD1", "NCAM1"),
  "lr_NK" = c("KLRD1", "NCAM1", "CD69"),
  "lr_NK_IFNG_pos" = c("KLRD1", "NCAM1", "CD69", "IFNG"),
  "LTi" = c("IL7R", "TCF7", "SELL", "KIT", "AHR", "LTB"),
  "Cir_NKT" = c("KLRD1", "CX3CR1", "CD3D"),
  "Lr_NKT" = c("KLRD1", "CD3D", "CD69"),
  
  # ========== Neutrophils ==========
  "Neutrophils" = c("CD16", "FCGR3B", "FCGR3A", "CD66b", "CEACAM8", "CD15", "FUT9", "S100A8", "S100A9", "CXCL8", "MMP9"),
  
  # ========== B cells and subsets ==========
  "B_cell" = c("CD79A", "CD79B", "MS4A1", "CD20", "VPREB3", "CD19"),
  "Activated_B" = c("CD27", "CD86", "HLA-DRA"),
  "Immature_B" = c("CD19", "CD10"),
  "T1_B" = c("CD19", "CD24", "CD38"),
  "T2_B" = c("CD19", "CD24", "CD38", "IGHD"),
  "Naive_B" = c("CD19", "IGHD", "IGHM"),
  "Transitional_B_TBCs" = c("CD27", "CD24", "IGHM", "IGHD"),
  "USMBs_Unswitched_Memory_B" = c("CD27", "CD38", "IGHM", "IGHD"),
  "SMBs_Switched_Memory_B" = c("CD27", "IGHA1", "IGHG1", "IGHE"),
  "CD27_pos_Memory_B" = c("CD27", "IGD"),
  "CD27_neg_Memory_B" = c("IGD"),
  "Plasmablast" = c("CD27", "CD38", "MKI67", "IGHG1", "SDC1"),
  "Plasma_cells" = c("CD27", "CD38", "TNFRSF17", "MZB1"),
  "Pre_pro_B" = c("CD34", "CD10"),
  "Pro_B" = c("CD19", "CD34", "CD10"),
  "Pre_B" = c("CD19", "CD10"),
  
  # ========== Myeloid cells (monocytes/macrophages/mast cells) ==========
  "Myeloid" = c("ITGAX", "CSF1R", "CD68", "FCER1G", "CD14", "AIF1", "TYROBP", "LYZ", "MS4A6A", "CD1E", "IL3RA"),
  "Activated_HSC" = c("ACTA2", "MFAP4", "CCN2", "MMP10"),
  "Monocytes" = c("VCAN", "FCN1", "S100A12"),
  "Classical_Mon" = c("CD14"),
  "Non_classical_Mon" = c("CD16"),
  "Intermediate_Mon" = c("CD14", "CD16"),
  "Pro_inflammatory_Mon" = c("IL1B"),
  "Lr_CD14pp_CD16p_Mon" = c("CD14", "CD16"),
  "Macrophages" = c("CD14", "CD68", "CD163", "APOE", "C1QA", "C1QB", "CD11B", "CD64", "MERTK"),
  "M1_Macrophages" = c("IL1B", "IL6", "TNFA"),
  "M2a_Macrophages" = c("IL10", "TGFB"),
  "M2b_Macrophages" = c("IL1", "IL6"),
  "M2c_Macrophages" = c("CD163"),
  "M2d_TAM" = c(), # No exclusive markers for tumor-associated macrophages, leave blank for manual addition
  "Mast_cells" = c("CPA3"),
  
  # ========== DC cells and subsets ==========
  "pDC" = c("LILRA4", "IL3RA", "TCF4", "TCL1A", "CLEC4C", "GZMB", "CLIC3"),
  "cDC1" = c("CLEC9A", "IRF8", "IDO1", "XCR1", "BATF3", "C1orf54"),
  "cDC2" = c("CD1E", "FCER1A", "CLEC10A", "CD1A", "CD1C", "TIMP1", "CEBPD", "LST1"),
  "DC3" = c("CCR7", "LAMP3", "FSCN1", "CCL19", "CCL22", "BIRC3"),
  "Liver_resident_pDC" = c("KLRB1", "CCL5", "NKG7"),
  
  # ========== Epithelial and fibrosis-related cells ==========
  "Epithelial" = c("EPCAM", "KLR5", "KRT6", "SPRR2", "SPRR3", "KRT19", "KRT18", "KRT15"),
  "BEC_Biliary_Epithelial" = c("AMBP", "KRT18", "KRT8", "KRT9"),
  "Fibrosis_related" = c("PDGFRA", "RGS5"),
  
  # ========== Fibroblasts/stromal cells ==========
  "CAF_Fibroblasts" = c("DCN", "TAGLN", "COL3A1", "COL1A1", "COL6A1", "ACTA2"),
  "HSCs" = c("PDGFRB", "PDGFRA", "CD140a", "DCN", "ACTA2", "COL3A1", "CD34"),
  "Stromal_cells" = c("PRPX1", "TWIST1"),
  
  # ========== Endothelial cell subsets ==========
  "Lymphatic_Endothelial" = c("MEOX1", "TBX1", "DTX1", "PROX1", "PDPN", "LYVE1", "FLT4", "HOXD3", "NR2F1", "NR2F2", "GPR182", "TEK"),
  "Arterial_ECs" = c("EFNB2", "SOX17", "BMX", "SEMA3G", "HEY1", "LTBP4", "FBLN5", "GJA5", "GJA4"),
  "Capillary_ECs" = c("CA4", "PRX", "RGCC", "SPARC", "SGK1"),
  "Venous_ECs" = c("NR2F2", "VCAM1", "ACKR1", "FCN3", "SELP"),
  
  # ========== Other cell types ==========
  "Neural_tube_cells" = c("SOX2", "SOX3", "HES5"),
  "Motor_neurons" = c("OLIG2", "LHX1", "LHX3", "ONECUT1", "ONECUT2", "MNX1"),
  "Sensory_neurons" = c("NEUROD4", "NEUROD1", "POU4F1", "ISL1", "EYA2", "SIX1"),
  "Cardiomyocytes" = c("MYL4", "TNNI2", "MYL7", "TNNI1"),
  "Hepatocytes" = c("AFR", "HNF4A", "TTR"),
  "Capsule_cells" = c("TCF21", "RSPO3", "WT1", "GATA4"),
  
  # ========== Function-related gene sets ==========
  "Proliferation_related" = c("TUBB", "STMN1", "HIST1H4C", "PCNA", "TK1"),
  "IFN_gamma_response" = c("CD40", "CD74", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "IRF8", "BST2"),
  "T_exhaustion_markers" = c("TIGIT", "BTLA", "LAG3", "PDCD1", "CTLA4", "HAVCR2"),
  "Naive_T_markers" = c("SELL", "TCF7", "LEF1", "CCR7"),
  "Effector_function_genes" = c("GZMA", "PRF1", "IL2", "GNLY", "GZMB", "GZMK", "IFNG", "NKG7"),
  "Inhibitory_markers" = c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2"),
  "Cytotoxic_genes" = c("IL2", "GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "IFNG", "NKG7"),
  "Stromal_cell_markers" = c("ENG", "CD105", "CD73", "NT5E", "CD90", "THY1", "CD51", "ITGAV", "CD146", "MCAM"),
  "Endothelial_markers" = c("CD34", "PECAM1", "CD31")
)

# Check marker genes present in dataset
cat("Checking marker genes present in dataset...\n")
available_markers <- list()
for (cell_type in names(marker_genes)) {
  available_genes <- marker_genes[[cell_type]][marker_genes[[cell_type]] %in% rownames(seurat_obj)]
  if (length(available_genes) > 0) {
    available_markers[[cell_type]] <- available_genes
    cat(sprintf("%s: Found %d/%d marker genes\n", 
                cell_type, length(available_genes), length(marker_genes[[cell_type]])))
  }
}

# Plot expression dot plots for all marker genes
cat("\nGenerating marker gene expression dot plots...\n")
for (cell_type in names(available_markers)) {
  genes <- available_markers[[cell_type]]
  if (length(genes) <= 6) {  # Limit number of genes per plot
    tryCatch({
      p <- FeaturePlot(seurat_obj, 
                       features = genes,
                       reduction = "umap",
                       ncol = min(3, length(genes)),
                       pt.size = 0.3,
                       order = TRUE,
                       combine = FALSE)
      
      # Add title to each subplot
      for (i in seq_along(p)) {
        p[[i]] <- p[[i]] + ggtitle(genes[i])
      }
      
      # Combine plots
      p_combined <- patchwork::wrap_plots(p, ncol = min(3, length(genes)))
      p_combined <- p_combined + 
        patchwork::plot_annotation(title = paste(cell_type, "Marker Genes"),
                                   theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      
      filename <- paste0("AAAAAAAAumap_markers_", gsub("_", "", cell_type), ".png")
      ggsave(filename, p_combined, width = 4 * min(3, length(genes)), height = 4, dpi = 300)
      
      cat(sprintf("✓ %s marker gene plot saved: %s\n", cell_type, filename))
    }, error = function(e) {
      cat(sprintf("Error generating %s marker gene plot: %s\n", cell_type, e$message))
    })
  }
}

# Plot marker gene DoHeatmap
cat("\nGenerating marker gene heatmap...\n")
# Get all available marker genes
all_available_markers <- unlist(available_markers)
all_available_markers <- unique(all_available_markers)

if (length(all_available_markers) > 0) {
  # Select top 30 marker genes (avoid oversized heatmap)
  markers_to_plot <- head(all_available_markers, 30)
  
  # Generate heatmap
  p_heatmap <- DoHeatmap(
    seurat_obj,
    features = markers_to_plot,
    slot = "data",  # Use normalized data
    group.by = "ident",
    group.colors = scales::hue_pal()(length(unique(Idents(seurat_obj)))),
    angle = 45,
    size = 3,
    label = TRUE
  ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    ggtitle("Marker Genes Expression Heatmap") +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  ggsave("marker_genes_heatmap.png", p_heatmap, width = 14, height = 10, dpi = 300)
  cat("✓ Marker gene heatmap saved to: marker_genes_heatmap.png\n")
}

# Create interactive annotation tool function
create_annotation_guide <- function(seurat_obj) {
  # Calculate average marker gene expression for each cluster
  cluster_ids <- levels(Idents(seurat_obj))
  marker_matrix <- matrix(0, nrow = length(cluster_ids), ncol = length(all_available_markers))
  rownames(marker_matrix) <- paste0("Cluster_", cluster_ids)
  colnames(marker_matrix) <- all_available_markers
  
  # Get normalized data
  norm_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Calculate average expression of each marker gene in each cluster
  for (i in seq_along(cluster_ids)) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster_ids[i])
    if (length(cluster_cells) > 0) {
      for (j in seq_along(all_available_markers)) {
        gene <- all_available_markers[j]
        if (gene %in% rownames(norm_data)) {
          marker_matrix[i, j] <- mean(norm_data[gene, cluster_cells])
        }
      }
    }
  }
  
  # Save marker gene expression matrix
  write.csv(marker_matrix, file = "cluster_marker_expression.csv")
  cat("✓ Cluster marker gene expression matrix saved to: cluster_marker_expression.csv\n")
  
  # Recommend possible cell types for each cluster
  cat("\nMarker gene expression profile for each cluster (top 5 most highly expressed genes):\n")
  for (i in seq_along(cluster_ids)) {
    cluster_name <- paste0("Cluster_", cluster_ids[i])
    top_genes <- names(sort(marker_matrix[i, ], decreasing = TRUE))[1:5]
    top_scores <- sort(marker_matrix[i, ], decreasing = TRUE)[1:5]
    
    cat(sprintf("\n%s (%d cells):\n", cluster_name, sum(Idents(seurat_obj) == cluster_ids[i])))
    for (j in 1:5) {
      if (top_scores[j] > 0) {
        cat(sprintf("  %s: %.3f\n", top_genes[j], top_scores[j]))
      }
    }
    
    # Suggest cell types based on marker genes
    cat("  Possible cell types: ")
    suggestions <- c()
    
    # Check T cell markers
    t_markers <- c("CD3D", "CD3E", "CD4", "CD8A")
    t_score <- mean(marker_matrix[i, t_markers[t_markers %in% colnames(marker_matrix)]])
    if (t_score > 0.5) suggestions <- c(suggestions, "T cells")
    
    # Check B cell markers
    b_markers <- c("CD19", "MS4A1", "CD79A")
    b_score <- mean(marker_matrix[i, b_markers[b_markers %in% colnames(marker_matrix)]])
    if (b_score > 0.5) suggestions <- c(suggestions, "B cells")
    
    # Check macrophage markers
    macro_markers <- c("CD14", "CD68", "CD163")
    macro_score <- mean(marker_matrix[i, macro_markers[macro_markers %in% colnames(marker_matrix)]])
    if (macro_score > 0.5) suggestions <- c(suggestions, "Macrophages")
    
    # Check fibroblast markers
    fibro_markers <- c("COL1A1", "COL3A1", "DCN")
    fibro_score <- mean(marker_matrix[i, fibro_markers[fibro_markers %in% colnames(marker_matrix)]])
    if (fibro_score > 0.5) suggestions <- c(suggestions, "Fibroblasts")
    
    # Check endothelial cell markers
    endo_markers <- c("PECAM1", "VWF", "CDH5")
    endo_score <- mean(marker_matrix[i, endo_markers[endo_markers %in% colnames(marker_matrix)]])
    if (endo_score > 0.5) suggestions <- c(suggestions, "Endothelial cells")
    
    if (length(suggestions) > 0) {
      cat(paste(suggestions, collapse = ", "))
    } else {
      cat("Unknown")
    }
    cat("\n")
  }
  
  return(marker_matrix)
}

# Run annotation guide
marker_matrix <- create_annotation_guide(seurat_obj)

# Save current object
saveRDS(seurat_obj, file = "seurat_obj_pre_annotation.rds")
cat("✓ Pre-annotation Seurat object saved to: seurat_obj_pre_annotation.rds\n")

# Define target gene to plot (CD125 corresponds to gene IL5RA)
target_gene <- "PDGFRB"  # Modify here if gene name is in other format in your data

# Plot UMAP annotation for this gene separately
tryCatch({
  # Plot single gene FeaturePlot
  p <- FeaturePlot(
    object = seurat_obj,
    features = target_gene,  # Pass single gene directly
    reduction = "umap",      # Same dimensionality reduction as in your original code
    pt.size = 0.3,           # Same point size as in your original code
    order = TRUE,            # Order points by expression level, high expression points on top
    label = FALSE            # Set to TRUE to display cell type labels if Idents(seurat_obj) has been defined
  ) +
    # Add gene title
    ggtitle(paste0("CD125 (", target_gene, ") Expression")) +
    # Center and bold title
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  # Save plot
  filename <- "AAAAAAA.png"
  ggsave(
    filename = filename,
    plot = p,
    width = 6,        # Width doesn't need to be too large for single gene plot
    height = 5,
    dpi = 300,
    bg = "white"      # White background to avoid transparent background after saving
  )
  
  cat(sprintf("✓ CD125(%s) annotation plot saved: %s\n", target_gene, filename))
}, error = function(e) {
  cat(sprintf("Error generating CD125(%s) annotation plot: %s\n", target_gene, e$message))
})

### Manual annotation
# --------------------------
# 14. Manual cell type annotation (corrected version)
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 14: Manual cell type annotation (corrected version)\n")
cat("=" * 70, "\n")

# View current clustering results
cat("Current clustering results:\n")
cluster_counts <- table(Idents(seurat_obj))
print(cluster_counts)

# Ensure cluster numbers match your annotations
all_clusters <- as.character(sort(unique(as.numeric(Idents(seurat_obj)))))
cat(sprintf("Total of %d clusters, numbered from %s to %s\n", 
            length(all_clusters), min(all_clusters), max(all_clusters)))

# Perform manual annotation according to your provided annotations
manual_annotation <- list(
  "0" = "Smooth_muscle_cells",      # Smooth muscle cells
  "1" = "Fibroblasts",              # Fibroblasts
  "2" = "B_cells",                  # B cells
  "3" = "Smooth_muscle_cells",      # Smooth muscle cells
  "4" = "Endothelial_cells",        # Endothelial cells
  "5" = "Endothelial_cells",        # Endothelial cells
  "6" = "Keratinocytes",            # Keratinocytes
  "7" = "Monocytes_Macrophages",    # Monocytes/Macrophages
  "8" = "Keratinocytes",            # Keratinocytes
  "9" = "Fibroblasts",              # Fibroblasts
  "10" = "Endothelial_cells",       # Endothelial cells
  "11" = "Melanocytes",             # Melanocytes
  "12" = "Mast_cells",              # Mast cells
  "13" = "Endothelial_cells"        # Endothelial cells
)

# Check for unannotated clusters
missing_clusters <- setdiff(all_clusters, names(manual_annotation))

if (length(missing_clusters) > 0) {
  cat("Warning: The following clusters have no annotations: ", paste(missing_clusters, collapse = ", "), "\n")
  cat("Assigning 'Unknown' label to these clusters\n")
  
  for (clust in missing_clusters) {
    manual_annotation[[clust]] <- "Unknown"
  }
}

# Apply annotations - safer method
cat("\nApplying cell type annotations...\n")

# Method: Create a new cell type vector
cell_types <- character(length(colnames(seurat_obj)))
names(cell_types) <- colnames(seurat_obj)

# Get clustering information for each cell
cluster_ids <- as.character(Idents(seurat_obj))

# Map clusters to cell types
for (i in seq_along(cell_types)) {
  cluster_id <- cluster_ids[i]
  if (cluster_id %in% names(manual_annotation)) {
    cell_types[i] <- manual_annotation[[cluster_id]]
  } else {
    cell_types[i] <- "Unknown"
  }
}

# Add to metadata
seurat_obj$cell_type <- cell_types

# Ensure no NA values
if (any(is.na(seurat_obj$cell_type))) {
  cat("NA values found, converting to 'Unknown'\n")
  seurat_obj$cell_type[is.na(seurat_obj$cell_type)] <- "Unknown"
}

# Convert cell type to factor, ordered by frequency
celltype_counts <- table(seurat_obj$cell_type)
celltype_levels <- names(sort(celltype_counts, decreasing = TRUE))
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = celltype_levels)

# View annotation results
cat("\nCell type annotation statistics:\n")
celltype_stats <- table(seurat_obj$cell_type)
print(celltype_stats)

# Create annotation summary
annotation_summary <- data.frame(
  Cluster = character(),
  Cell_Type = character(),
  n_Cells = numeric(),
  stringsAsFactors = FALSE
)

for (cluster_id in all_clusters) {
  if (cluster_id %in% names(manual_annotation)) {
    cell_type <- manual_annotation[[cluster_id]]
    n_cells <- sum(cluster_ids == cluster_id)
    annotation_summary <- rbind(annotation_summary, 
                                data.frame(Cluster = cluster_id,
                                           Cell_Type = cell_type,
                                           n_Cells = n_cells))
  }
}

write.csv(annotation_summary, file = "celltype_annotation_summary.csv", row.names = FALSE)
cat("✓ Cell type annotation summary saved to: celltype_annotation_summary.csv\n")

# 15. Plot annotated UMAP (corrected version)
cat("\n" + "=" * 70, "\n")
cat("Step 15: Generate annotated visualization plots\n")
cat("=" * 70, "\n")

# Check and ensure all cells have UMAP coordinates
cat("\nChecking UMAP coordinate integrity...\n")
if ("umap" %in% names(seurat_obj@reductions)) {
  # Get cells with UMAP coordinates
  cells_with_umap <- rownames(seurat_obj@reductions$umap@cell.embeddings)
  all_cells <- colnames(seurat_obj)
  missing_umap <- setdiff(all_cells, cells_with_umap)
  
  if (length(missing_umap) > 0) {
    cat(sprintf("Warning: %d cells missing UMAP coordinates\n", length(missing_umap)))
    cat("This may be due to previous calculation issues, attempting to recalculate UMAP...\n")
    
    # Rerun UMAP
    set.seed(123)
    seurat_obj <- RunUMAP(
      seurat_obj,
      dims = 1:n_pcs,
      verbose = FALSE
    )
    cat("UMAP recalculation completed\n")
  } else {
    cat("All cells have UMAP coordinates\n")
  }
} else {
  cat("Warning: UMAP dimensionality reduction results not found, recalculating UMAP...\n")
  set.seed(123)
  seurat_obj <- RunUMAP(
    seurat_obj,
    dims = 1:n_pcs,
    verbose = FALSE
  )
  cat("UMAP calculation completed\n")
}

# Plot annotated UMAP
cat("\nGenerating annotated UMAP plot...\n")

# Method 1: Using DimPlot
p_umap_celltype <- DimPlot(seurat_obj, 
                           reduction = "umap", 
                           group.by = "cell_type",
                           label = TRUE,
                           pt.size = 0.5,
                           repel = TRUE,
                           label.size = 4) +
  ggtitle("UMAP - Cell Type Annotation") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave("umap_celltype_annotated.png", p_umap_celltype, 
       width = 14, height = 10, dpi = 300, bg = "white")
cat("✓ Annotated UMAP plot saved to: umap_celltype_annotated.png\n")

# Method 2: Direct plotting with ggplot2 for more control
cat("\nGenerating more detailed UMAP plot using ggplot2...\n")

# Extract UMAP coordinates and cell types
umap_coords <- as.data.frame(seurat_obj@reductions$umap@cell.embeddings)
colnames(umap_coords) <- c("UMAP1", "UMAP2")
umap_coords$cell_type <- seurat_obj$cell_type
umap_coords$cell_id <- rownames(umap_coords)

# Calculate centroids for each cell type for labeling
label_positions <- umap_coords %>%
  group_by(cell_type) %>%
  summarize(
    UMAP1 = median(UMAP1),
    UMAP2 = median(UMAP2),
    n_cells = n()
  )

# Plot UMAP
p_umap_ggplot <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.5, alpha = 0.6) +
  geom_text(data = label_positions,
            aes(label = cell_type, x = UMAP1, y = UMAP2),
            size = 4, fontface = "bold", color = "black",
            check_overlap = FALSE) +
  scale_color_viridis_d(option = "plasma", name = "Cell Type") +
  labs(
    title = "UMAP - Cell Type Annotation",
    subtitle = sprintf("Total %d cells, %d cell types", nrow(umap_coords), length(unique(umap_coords$cell_type)))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey50"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave("umap_celltype_ggplot.png", p_umap_ggplot, 
       width = 14, height = 10, dpi = 300, bg = "white")
cat("✓ ggplot2 version UMAP plot saved to: umap_celltype_ggplot.png\n")

# Faceted display by sample group
if ("group" %in% colnames(seurat_obj@meta.data)) {
  cat("\nGenerating UMAP plot faceted by sample group...\n")
  
  # Check for NA values
  if (any(is.na(seurat_obj$group))) {
    cat("Warning: NA values found in group column, converting to 'Unknown'\n")
    seurat_obj$group[is.na(seurat_obj$group)] <- "Unknown"
  }
  
  p_umap_split <- DimPlot(seurat_obj,
                          reduction = "umap",
                          group.by = "cell_type",
                          split.by = "group",
                          ncol = min(3, length(unique(seurat_obj$group))),
                          pt.size = 0.3,
                          label = FALSE) +
    ggtitle("UMAP - Cell Type Distribution by Sample Group") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  ggsave("umap_celltype_by_group.png", p_umap_split, 
         width = 5 * min(3, length(unique(seurat_obj$group))), 
         height = 5, dpi = 300, bg = "white")
  cat("✓ UMAP plot faceted by sample group saved to: umap_celltype_by_group.png\n")
}

# 16. Save annotated results
cat("\n" + "=" * 70, "\n")
cat("Step 16: Save final annotation results\n")
cat("=" * 70, "\n")

# Save complete annotated object
saveRDS(seurat_obj, file = "seurat_obj_fully_annotated.rds")
cat("✓ Fully annotated Seurat object saved to: seurat_obj_fully_annotated.rds\n")

# Save metadata containing annotations
write.csv(seurat_obj@meta.data, file = "final_metadata_with_celltypes.csv", row.names = TRUE)
cat("✓ Metadata with cell type annotations saved to: final_metadata_with_celltypes.csv\n")

cat("\n" + "=" * 70, "\n")
cat("Manual annotation completed!\n")
cat("=" * 70, "\n")

cat("\nAnnotation summary:\n")
for (i in 1:nrow(annotation_summary)) {
  cat(sprintf("Cluster %s: %s (%d cells, %.1f%%)\n",
              annotation_summary$Cluster[i],
              annotation_summary$Cell_Type[i],
              annotation_summary$n_Cells[i],
              annotation_summary$n_Cells[i] / sum(annotation_summary$n_Cells) * 100))
}

# --------------------------
# 13. Plot expression dot plots for target genes
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 13: Generate expression dot plots for core genes VEGFA, THBS1, MTHFR\n")
cat("=" * 70, "\n")

# Target gene list
target_genes <- c("VEGFA", "THBS1", "MTHFR")

# Check if target genes are present in data
cat("Checking if target genes are present in data...\n")
available_genes <- target_genes[target_genes %in% rownames(seurat_obj)]
if (length(available_genes) == 0) {
  stop("Error: No target genes found! Please check if gene names are correct.")
} else {
  cat(sprintf("Found %d/%d target genes: %s\n", 
              length(available_genes), length(target_genes), 
              paste(available_genes, collapse = ", ")))
}

# Check if cell type annotations exist
if (!"cell_type" %in% colnames(seurat_obj@meta.data)) {
  cat("Warning: Cell type annotations not found, using clustering results as cell types\n")
  seurat_obj$cell_type <- Idents(seurat_obj)
}

# Check if sample group information exists
if (!"group" %in% colnames(seurat_obj@meta.data)) {
  cat("Warning: Sample group information not found, setting to 'Unknown'\n")
  if (!"group" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$group <- "Unknown"
  }
}

# 1. Calculate summary data by cell type
cat("\nCalculating summary data by cell type...\n")

# Get all cell types
cell_types <- unique(seurat_obj$cell_type)

# Create summary data frame
cell_type_summary <- data.frame()

# Get normalized expression data
expr_data <- GetAssayData(seurat_obj, slot = "data")

# Calculate statistics for each gene in each cell type
for (cell_type in cell_types) {
  # Get cells of this cell type
  cell_type_cells <- colnames(seurat_obj)[seurat_obj$cell_type == cell_type]
  
  if (length(cell_type_cells) > 0) {
    for (gene in available_genes) {
      if (gene %in% rownames(expr_data)) {
        # Get expression values of this gene in this cell type
        gene_expr <- expr_data[gene, cell_type_cells]
        
        # Calculate average expression (mean of non-zero values)
        avg_expr <- mean(gene_expr[gene_expr > 0])
        if (is.na(avg_expr)) avg_expr <- 0
        
        # Calculate percentage of cells expressing this gene
        expr_percent <- sum(gene_expr > 0) / length(gene_expr) * 100
        
        # Add a row of data
        cell_type_summary <- rbind(cell_type_summary, data.frame(
          CellType = cell_type,
          Gene = gene,
          AvgExpression = avg_expr,
          PercentExpressed = expr_percent,
          nCells = length(cell_type_cells)
        ))
      }
    }
  }
}

# 2. Calculate summary data by sample group
cat("\nCalculating summary data by sample group...\n")

# Get all sample groups
sample_groups <- unique(seurat_obj$group)

# Create summary data frame
group_summary <- data.frame()

# Calculate statistics for each gene in each sample group
for (group in sample_groups) {
  # Get cells of this group
  group_cells <- colnames(seurat_obj)[seurat_obj$group == group]
  
  if (length(group_cells) > 0) {
    for (gene in available_genes) {
      if (gene %in% rownames(expr_data)) {
        # Get expression values of this gene in this group
        gene_expr <- expr_data[gene, group_cells]
        
        # Calculate average expression
        avg_expr <- mean(gene_expr[gene_expr > 0])
        if (is.na(avg_expr)) avg_expr <- 0
        
        # Calculate percentage of cells expressing this gene
        expr_percent <- sum(gene_expr > 0) / length(gene_expr) * 100
        
        # Add a row of data
        group_summary <- rbind(group_summary, data.frame(
          Group = group,
          Gene = gene,
          AvgExpression = avg_expr,
          PercentExpressed = expr_percent,
          nCells = length(group_cells)
        ))
      }
    }
  }
}

# Save summary data
write.csv(cell_type_summary, file = "target_genes_celltype_summary.csv", row.names = FALSE)
write.csv(group_summary, file = "target_genes_group_summary.csv", row.names = FALSE)
cat("✓ Summary data saved to: target_genes_celltype_summary.csv and target_genes_group_summary.csv\n")

# 3. Plot dot plot by cell type
cat("\nGenerating dot plot by cell type...\n")

# Set color and size ranges
# Color: darker color indicates higher expression (using viridis color palette)
library(viridis)

# Ensure Gene order
cell_type_summary$Gene <- factor(cell_type_summary$Gene, levels = available_genes)

# Create dot plot
p_celltype <- ggplot(cell_type_summary, 
                     aes(x = Gene, y = reorder(CellType, desc(CellType)))) +
  geom_point(aes(size = PercentExpressed, color = AvgExpression), alpha = 0.8) +
  scale_color_viridis(
    name = "Average Expression",
    option = "plasma",
    direction = -1,
    limits = c(0, max(cell_type_summary$AvgExpression) * 1.1)
  ) +
  scale_size_continuous(
    name = "Percent of Expressing Cells (%)",
    range = c(2, 10),
    breaks = c(20, 40, 60, 80, 100),
    labels = c("20%", "40%", "60%", "80%", "100%")
  ) +
  labs(
    title = "Core Gene Expression Across Different Cell Types",
    x = "Target Genes",
    y = "Cell Types"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Add cell count labels
  geom_text(
    aes(label = paste0("n=", nCells)),
    size = 3,
    color = "grey40",
    nudge_x = 0.3
  )

# Save dot plot by cell type
ggsave("dotplot_target_genes_by_celltype.png", p_celltype, 
       width = 10, height = 8, dpi = 300, bg = "white")
cat("✓ Dot plot by cell type saved to: dotplot_target_genes_by_celltype.png\n")

# 18. Generate target gene expression dot plots for disease and normal groups separately
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 18: Generate target gene expression dot plots for disease and normal groups separately\n")
cat("=" * 70, "\n")

# Target gene list
target_genes <- c("VEGFA", "THBS1", "MTHFR")

# Check if target genes are present in data
cat("Checking if target genes are present in data...\n")
available_genes <- target_genes[target_genes %in% rownames(seurat_obj)]
if (length(available_genes) == 0) {
  stop("Error: No target genes found! Please check if gene names are correct.")
} else {
  cat(sprintf("Found %d/%d target genes: %s\n", 
              length(available_genes), length(target_genes), 
              paste(available_genes, collapse = ", ")))
}

# Check if grouping information exists
if (!"group" %in% colnames(seurat_obj@meta.data)) {
  cat("Warning: Grouping information not found, using 'Unknown' as group\n")
  seurat_obj$group <- "Unknown"
}

# Check if Normal and DFU groups exist
if (!all(c("Normal", "DFU") %in% seurat_obj$group)) {
  cat("Warning: Normal and DFU groups not found. Current groups are:", paste(unique(seurat_obj$group), collapse = ", "), "\n")
  cat("Will plot using existing groups.\n")
}

# Get normalized expression data
expr_data <- GetAssayData(seurat_obj, slot = "data")

# 1. Calculate cell type summary data for Normal and DFU groups separately
cat("\nCalculating cell type summary data for Normal and DFU groups separately...\n")

# Function: Calculate cell type summary data for specified group
calculate_group_celltype_summary <- function(group_name) {
  # Get cells of this group
  group_cells <- colnames(seurat_obj)[seurat_obj$group == group_name]
  
  if (length(group_cells) == 0) {
    cat(sprintf("Warning: No cells in %s group\n", group_name))
    return(NULL)
  }
  
  # Create subset for this group
  seurat_group <- subset(seurat_obj, cells = group_cells)
  
  # Get cell types of this group
  cell_types <- unique(seurat_group$cell_type)
  
  # Create summary data frame
  summary_df <- data.frame()
  
  # Calculate statistics for each gene in each cell type
  for (cell_type in cell_types) {
    # Get cells of this cell type
    cell_type_cells <- colnames(seurat_group)[seurat_group$cell_type == cell_type]
    
    if (length(cell_type_cells) > 0) {
      for (gene in available_genes) {
        if (gene %in% rownames(expr_data)) {
          # Get expression values of this gene in this cell type
          gene_expr <- expr_data[gene, cell_type_cells]
          
          # Calculate average expression (mean of non-zero values)
          avg_expr <- mean(gene_expr[gene_expr > 0])
          if (is.na(avg_expr)) avg_expr <- 0
          
          # Calculate percentage of cells expressing this gene
          expr_percent <- sum(gene_expr > 0) / length(gene_expr) * 100
          
          # Add a row of data
          summary_df <- rbind(summary_df, data.frame(
            Group = group_name,
            CellType = cell_type,
            Gene = gene,
            AvgExpression = avg_expr,
            PercentExpressed = expr_percent,
            nCells = length(cell_type_cells)
          ))
        }
      }
    }
  }
  
  return(summary_df)
}

# Calculate summary data for Normal group
cat("Calculating summary data for Normal group...\n")
normal_summary <- calculate_group_celltype_summary("Normal")

# Calculate summary data for DFU group
cat("Calculating summary data for DFU group...\n")
dfu_summary <- calculate_group_celltype_summary("DFU")

# Merge summary data of two groups
group_celltype_summary <- rbind(normal_summary, dfu_summary)

# Save summary data
write.csv(group_celltype_summary, file = "target_genes_group_celltype_summary.csv", row.names = FALSE)
cat("✓ Grouped cell type summary data saved to: target_genes_group_celltype_summary.csv\n")

# 2. Generate cell type dot plot for Normal group
if (!is.null(normal_summary) && nrow(normal_summary) > 0) {
  cat("\nGenerating cell type dot plot for Normal group...\n")
  
  # Ensure Gene order
  normal_summary$Gene <- factor(normal_summary$Gene, levels = available_genes)
  
  # Create dot plot
  p_normal_celltype <- ggplot(normal_summary, 
                              aes(x = Gene, y = reorder(CellType, desc(CellType)))) +
    geom_point(aes(size = PercentExpressed, color = AvgExpression), alpha = 0.8) +
    scale_color_viridis(
      name = "Average Expression",
      option = "plasma",
      direction = -1,
      limits = c(0, max(normal_summary$AvgExpression) * 1.1)
    ) +
    scale_size_continuous(
      name = "Percent of Expressing Cells (%)",
      range = c(2, 10),
      breaks = c(20, 40, 60, 80, 100),
      labels = c("20%", "40%", "60%", "80%", "100%")
    ) +
    labs(
      title = "Normal Group - Core Gene Expression Across Different Cell Types",
      x = "Target Genes",
      y = "Cell Types",
      caption = "Normal Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.caption = element_text(hjust = 1, size = 12, face = "italic", color = "#2E86AB"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption.position = "plot"
    ) +
    # Add cell count labels
    geom_text(
      aes(label = paste0("n=", nCells)),
      size = 3,
      color = "grey40",
      nudge_x = 0.3
    )
  
  # Save Normal group dot plot
  ggsave("dotplot_target_genes_normal_celltype.png", p_normal_celltype, 
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("✓ Normal group cell type dot plot saved to: dotplot_target_genes_normal_celltype.png\n")
} else {
  cat("Warning: No data or empty data for Normal group, skipping Normal group dot plot generation\n")
}

# 3. Generate cell type dot plot for DFU group
if (!is.null(dfu_summary) && nrow(dfu_summary) > 0) {
  cat("\nGenerating cell type dot plot for DFU group...\n")
  
  # Ensure Gene order
  dfu_summary$Gene <- factor(dfu_summary$Gene, levels = available_genes)
  
  # Create dot plot
  p_dfu_celltype <- ggplot(dfu_summary, 
                           aes(x = Gene, y = reorder(CellType, desc(CellType)))) +
    geom_point(aes(size = PercentExpressed, color = AvgExpression), alpha = 0.8) +
    scale_color_viridis(
      name = "Average Expression",
      option = "plasma",
      direction = -1,
      limits = c(0, max(dfu_summary$AvgExpression) * 1.1)
    ) +
    scale_size_continuous(
      name = "Percent of Expressing Cells (%)",
      range = c(2, 10),
      breaks = c(20, 40, 60, 80, 100),
      labels = c("20%", "40%", "60%", "80%", "100%")
    ) +
    labs(
      title = "DFU Group - Core Gene Expression Across Different Cell Types",
      x = "Target Genes",
      y = "Cell Types",
      caption = "DFU Group"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.caption = element_text(hjust = 1, size = 12, face = "italic", color = "#A23B72"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(face = "bold", size = 14),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "grey90", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.caption.position = "plot"
    ) +
    # Add cell count labels
    geom_text(
      aes(label = paste0("n=", nCells)),
      size = 3,
      color = "grey40",
      nudge_x = 0.3
    )
  
  # Save DFU group dot plot
  ggsave("dotplot_target_genes_dfu_celltype.png", p_dfu_celltype, 
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("✓ DFU group cell type dot plot saved to: dotplot_target_genes_dfu_celltype.png\n")
} else {
  cat("Warning: No data or empty data for DFU group, skipping DFU group dot plot generation\n")
}

# 4. Generate dot plot by sample group
cat("\nGenerating dot plot by sample group...\n")

# Ensure Gene order
group_summary$Gene <- factor(group_summary$Gene, levels = available_genes)

# Set group colors (if only Normal and DFU groups)
if (all(c("Normal", "DFU") %in% sample_groups)) {
  group_colors <- c("Normal" = "#2E86AB", "DFU" = "#A23B72", "Unknown" = "#8A8D91")
} else {
  group_colors <- scales::hue_pal()(length(sample_groups))
  names(group_colors) <- sample_groups
}

# Create dot plot - modified parameter version
p_group <- ggplot(group_summary, 
                  aes(x = Gene, y = reorder(Group, desc(Group)))) +
  geom_point(aes(size = PercentExpressed, color = AvgExpression), alpha = 0.8) +
  scale_color_viridis(
    name = "Average Expression",
    option = "plasma",
    direction = -1,
    limits = c(1.2, 1.7),  # Limit color range to 1.2-1.7 to enhance discrimination in this interval
    na.value = "lightgrey"  # Display values outside range in grey
  ) +
  scale_size_continuous(
    name = "Percent of Expressing Cells (%)",
    range = c(3, 10),  # Appropriately adjust point size range
    limits = c(5, 25)   # Limit size range to 5%-25%
  ) +
  labs(
    title = "Core Gene Expression Across Different Sample Groups",
    x = "Target Genes",
    y = "Sample Groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Add cell count labels
  geom_text(
    aes(label = paste0("n=", nCells)),
    size = 4,
    color = "grey40",
    nudge_x = 0.3
  )

# Save dot plot by sample group
ggsave("dotplot_target_genes_by_group.png", p_group, 
       width = 10, height = 6, dpi = 300, bg = "white")
cat("✓ Dot plot by sample group saved to: dotplot_target_genes_by_group.png\n")

# 5. Create combined plot (including two dot plots)
cat("\nCreating combined plot...\n")

# Combine two dot plots
library(patchwork)

p_combined <- (p_celltype + theme(legend.position = "none")) / 
  (p_group + theme(legend.position = "right")) + 
  plot_annotation(
    title = "Core Gene Expression Analysis of VEGFA, THBS1, MTHFR",
    subtitle = "Point size indicates percentage of expressing cells, color indicates average expression level",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey50")
    )
  ) +
  plot_layout(heights = c(2, 1))

ggsave("dotplot_target_genes_combined.png", p_combined, 
       width = 16, height = 12, dpi = 300, bg = "white")
cat("✓ Combined dot plot saved to: dotplot_target_genes_combined.png\n")

# 6. Generate detailed heatmap-style dot plot
cat("\nGenerating heatmap-style dot plot...\n")

# Cell type heatmap
p_heatmap_celltype <- ggplot(cell_type_summary, 
                             aes(x = Gene, y = reorder(CellType, desc(CellType)))) +
  geom_tile(aes(fill = AvgExpression, alpha = PercentExpressed/100), color = "white") +
  scale_fill_viridis(
    name = "Average Expression",
    option = "plasma",
    direction = -1
  ) +
  scale_alpha_continuous(
    name = "Expression Percentage",
    range = c(0.3, 1),
    breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)
  ) +
  geom_text(aes(label = sprintf("%.2f\n(%.1f%%)", AvgExpression, PercentExpressed)), 
            size = 3, color = "white", fontface = "bold") +
  labs(
    title = "Core Gene Expression Heatmap (By Cell Type)",
    x = "Target Genes",
    y = "Cell Types",
    caption = "Value format: Average Expression\n(Percent of Expressing Cells%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    plot.caption = element_text(size = 10, color = "grey50", hjust = 0)
  )

ggsave("heatmap_target_genes_by_celltype.png", p_heatmap_celltype, 
       width = 10, height = 8, dpi = 300, bg = "white")
cat("✓ Cell type heatmap saved to: heatmap_target_genes_by_celltype.png\n")

# 7. Generate UMAP plots for target genes
cat("\nGenerating UMAP plots for target genes...\n")

# Generate UMAP plot for each target gene individually
for (gene in available_genes) {
  tryCatch({
    p_umap <- FeaturePlot(
      seurat_obj,
      features = gene,
      reduction = "umap",
      pt.size = 0.5,
      order = TRUE
    ) +
      scale_color_gradientn(
        colors = c("lightgrey", "yellow", "orange", "red"),
        name = "Expression"
      ) +
      ggtitle(sprintf("UMAP: %s Expression", gene)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right"
      )
    
    filename <- sprintf("umap_%s_expression.png", gene)
    ggsave(filename, p_umap, width = 10, height = 8, dpi = 300, bg = "white")
    cat(sprintf("✓ UMAP plot for %s saved to: %s\n", gene, filename))
    
  }, error = function(e) {
    cat(sprintf("Error generating UMAP plot for %s: %s\n", gene, e$message))
  })
}

# 19. Generate UMAP expression plots for core genes in disease group
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 19: Generate UMAP expression plots for core genes in disease group\n")
cat("=" * 70, "\n")

# Check if DFU group exists
if (!"DFU" %in% unique(seurat_obj$group)) {
  stop("Error: DFU group not found! Please check grouping information.")
}

# Extract DFU group cells
dfu_cells <- colnames(seurat_obj)[seurat_obj$group == "DFU"]
cat(sprintf("Total %d cells in DFU group\n", length(dfu_cells)))

if (length(dfu_cells) == 0) {
  stop("Error: No cells in DFU group!")
}

# Create Seurat subset for DFU group
cat("Creating Seurat subset for DFU group...\n")
seurat_dfu <- subset(seurat_obj, cells = dfu_cells)

# Check if target genes are present in data
target_genes <- c("VEGFA", "THBS1", "MTHFR")
available_genes <- target_genes[target_genes %in% rownames(seurat_dfu)]
cat(sprintf("Found %d/%d target genes in DFU group: %s\n", 
            length(available_genes), length(target_genes), 
            paste(available_genes, collapse = ", ")))

if (length(available_genes) == 0) {
  stop("Error: No target genes found in DFU group!")
}

# 1. Generate cell type UMAP plot for DFU group
cat("\n1. Generating cell type UMAP plot for DFU group...\n")

# Define color scheme for DFU group
dfu_celltype_colors <- c(
  "Smooth_muscle_cells" = "#1f77b4",
  "Fibroblasts" = "#ff7f0e",
  "B_cells" = "#2ca02c",
  "Endothelial_cells" = "#d62728",
  "Keratinocytes" = "#9467bd",
  "Monocytes_Macrophages" = "#8c564b",
  "Melanocytes" = "#e377c2",
  "Mast_cells" = "#7f7f7f"
)

# Generate cell type UMAP plot for DFU group
p_dfu_celltype <- DimPlot(seurat_dfu,
                          reduction = "umap",
                          group.by = "cell_type",
                          pt.size = 0.8,
                          label = TRUE,
                          repel = TRUE,
                          cols = dfu_celltype_colors) +
  ggtitle("DFU Group - Cell Type Distribution") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave("umap_dfu_celltype.png", p_dfu_celltype, 
       width = 12, height = 9, dpi = 300, bg = "white")
cat("✓ DFU group cell type UMAP plot saved to: umap_dfu_celltype.png\n")

# 2. Generate expression UMAP plot for each core gene in DFU group
cat("\n2. Generating expression UMAP plots for core genes in DFU group...\n")

# Generate individual UMAP plot for each gene
for (gene in available_genes) {
  cat(sprintf("  Generating expression UMAP plot for %s...\n", gene))
  
  # Check if gene is expressed
  gene_expr <- GetAssayData(seurat_dfu, slot = "data")[gene, ]
  if (sum(gene_expr > 0) == 0) {
    cat(sprintf("  Warning: %s not expressed in DFU group, skipping plot generation\n", gene))
    next
  }
  
  # Create expression UMAP plot
  p_gene <- FeaturePlot(seurat_dfu,
                        features = gene,
                        reduction = "umap",
                        pt.size = 0.8,
                        order = TRUE,
                        min.cutoff = "q10",
                        max.cutoff = "q90") +
    scale_color_gradientn(
      colors = c("lightgrey", "#FFD700", "#FF8C00", "#FF4500", "#8B0000"),
      name = "Expression"
    ) +
    ggtitle(sprintf("DFU Group - %s Expression Distribution", gene)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Save individual gene plot
  filename <- sprintf("umap_dfu_%s_expression.png", gene)
  ggsave(filename, p_gene, width = 10, height = 8, dpi = 300, bg = "white")
  cat(sprintf("  ✓ UMAP plot for %s saved to: %s\n", gene, filename))
}

# 8. Generate grouped UMAP plots
cat("\nGenerating grouped UMAP plots for target genes...\n")

if ("umap" %in% names(seurat_obj@reductions)) {
  # Generate grouped UMAP plot for each target gene
  for (gene in available_genes) {
    tryCatch({
      p_split <- FeaturePlot(
        seurat_obj,
        features = gene,
        reduction = "umap",
        split.by = "group",
        pt.size = 0.5,
        order = TRUE,
        combine = FALSE
      )
      
      # Add titles
      for (i in seq_along(p_split)) {
        p_split[[i]] <- p_split[[i]] + 
          ggtitle(sprintf("%s - %s", gene, names(p_split)[i])) +
          theme(plot.title = element_text(size = 10, hjust = 0.5))
      }
      
      # Combine plots
      n_groups <- length(unique(seurat_obj$group))
      n_cols <- min(3, n_groups)
      p_combined_split <- patchwork::wrap_plots(p_split, ncol = n_cols)
      p_combined_split <- p_combined_split + 
        plot_annotation(
          title = sprintf("%s Expression by Group", gene),
          theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        )
      
      filename <- sprintf("umap_%s_by_group.png", gene)
      ggsave(filename, p_combined_split, 
             width = 5 * n_cols, height = 5, dpi = 300, bg = "white")
      cat(sprintf("✓ Grouped UMAP plot for %s saved to: %s\n", gene, filename))
      
    }, error = function(e) {
      cat(sprintf("Error generating grouped UMAP plot for %s: %s\n", gene, e$message))
    })
  }
}

# 10. Create statistical analysis report
cat("\nCreating statistical analysis report...\n")

# Calculate some basic statistics
stats_report <- data.frame(
  Metric = character(),
  VEGFA = character(),
  THBS1 = character(),
  MTHFR = character(),
  stringsAsFactors = FALSE
)

# Add statistical data
for (gene in available_genes) {
  # Get expression of this gene in all cells
  gene_expr <- expr_data[gene, ]
  
  stats_report <- rbind(stats_report, data.frame(
    Metric = "Overall Average Expression",
    VEGFA = ifelse(gene == "VEGFA", sprintf("%.3f", mean(gene_expr[gene_expr > 0])), ""),
    THBS1 = ifelse(gene == "THBS1", sprintf("%.3f", mean(gene_expr[gene_expr > 0])), ""),
    MTHFR = ifelse(gene == "MTHFR", sprintf("%.3f", mean(gene_expr[gene_expr > 0])), "")
  ))
  
  stats_report <- rbind(stats_report, data.frame(
    Metric = "Overall Percentage of Expressing Cells",
    VEGFA = ifelse(gene == "VEGFA", sprintf("%.1f%%", sum(gene_expr > 0)/length(gene_expr)*100), ""),
    THBS1 = ifelse(gene == "THBS1", sprintf("%.1f%%", sum(gene_expr > 0)/length(gene_expr)*100), ""),
    MTHFR = ifelse(gene == "MTHFR", sprintf("%.1f%%", sum(gene_expr > 0)/length(gene_expr)*100), "")
  ))
  
  stats_report <- rbind(stats_report, data.frame(
    Metric = "Expression Range",
    VEGFA = ifelse(gene == "VEGFA", sprintf("[%.3f, %.3f]", min(gene_expr[gene_expr > 0]), max(gene_expr)), ""),
    THBS1 = ifelse(gene == "THBS1", sprintf("[%.3f, %.3f]", min(gene_expr[gene_expr > 0]), max(gene_expr)), ""),
    MTHFR = ifelse(gene == "MTHFR", sprintf("[%.3f, %.3f]", min(gene_expr[gene_expr > 0]), max(gene_expr)), "")
  ))
  
  # Add group-specific statistics
  for (group in c("Normal", "DFU")) {
    if (group %in% seurat_obj$group) {
      group_cells <- colnames(seurat_obj)[seurat_obj$group == group]
      group_expr <- expr_data[gene, group_cells]
      
      stats_report <- rbind(stats_report, data.frame(
        Metric = paste0(group, " Average Expression"),
        VEGFA = ifelse(gene == "VEGFA", sprintf("%.3f", mean(group_expr[group_expr > 0])), ""),
        THBS1 = ifelse(gene == "THBS1", sprintf("%.3f", mean(group_expr[group_expr > 0])), ""),
        MTHFR = ifelse(gene == "MTHFR", sprintf("%.3f", mean(group_expr[group_expr > 0])), "")
      ))
      
      stats_report <- rbind(stats_report, data.frame(
        Metric = paste0(group, " Percentage of Expressing Cells"),
        VEGFA = ifelse(gene == "VEGFA", sprintf("%.1f%%", sum(group_expr > 0)/length(group_expr)*100), ""),
        THBS1 = ifelse(gene == "THBS1", sprintf("%.1f%%", sum(group_expr > 0)/length(group_expr)*100), ""),
        MTHFR = ifelse(gene == "MTHFR", sprintf("%.1f%%", sum(group_expr > 0)/length(group_expr)*100), "")
      ))
    }
  }
}

# Remove empty rows and reorganize
stats_report <- stats_report[rowSums(stats_report != "") > 1, ]
rownames(stats_report) <- NULL

# Save statistical report
write.csv(stats_report, file = "target_genes_statistical_report.csv", row.names = FALSE)
cat("✓ Statistical report saved to: target_genes_statistical_report.csv\n")

# Print statistical summary to console
cat("\n" + "=" * 70, "\n")
cat("Target Gene Statistical Summary\n")
cat("=" * 70, "\n")
print(stats_report, row.names = FALSE)

# --------------------------
# 20. Group-wise differential expression analysis of target genes
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 20: Differential expression analysis between Normal and DFU groups\n")
cat("=" * 70, "\n")

# Check if both groups exist
if (!all(c("Normal", "DFU") %in% unique(seurat_obj$group))) {
  cat("Warning: Missing Normal or DFU group, skipping differential expression analysis\n")
} else {
  # Prepare data for differential analysis
  de_results_list <- list()
  
  for (gene in available_genes) {
    cat(sprintf("Analyzing differential expression of %s...\n", gene))
    
    # Extract expression data by group
    normal_expr <- expr_data[gene, seurat_obj$group == "Normal"]
    dfu_expr <- expr_data[gene, seurat_obj$group == "DFU"]
    
    # Perform Wilcoxon rank-sum test (non-parametric test suitable for single-cell data)
    wilcox_test <- wilcox.test(normal_expr, dfu_expr, alternative = "two.sided")
    
    # Calculate fold change (log2)
    normal_mean <- mean(normal_expr[normal_expr > 0])
    dfu_mean <- mean(dfu_expr[dfu_expr > 0])
    log2_fc <- ifelse(normal_mean > 0 & dfu_mean > 0, log2(dfu_mean / normal_mean), NA)
    
    # Calculate expression ratio
    expr_ratio <- ifelse(normal_mean > 0, dfu_mean / normal_mean, NA)
    
    # Store results
    de_results_list[[gene]] <- data.frame(
      Gene = gene,
      Normal_Average = normal_mean,
      DFU_Average = dfu_mean,
      Log2_FoldChange = log2_fc,
      FoldChange_Ratio = expr_ratio,
      P_value = wilcox_test$p.value,
      P_adjusted = p.adjust(wilcox_test$p.value, method = "fdr", n = length(available_genes)),
      Test_Method = "Wilcoxon rank-sum test"
    )
  }
  
  # Combine all differential analysis results
  de_results <- do.call(rbind, de_results_list)
  rownames(de_results) <- NULL
  
  # Add significance classification
  de_results$Significance <- case_when(
    de_results$P_adjusted < 0.001 ~ "***",
    de_results$P_adjusted < 0.01 ~ "**",
    de_results$P_adjusted < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  # Save differential analysis results
  write.csv(de_results, file = "target_genes_differential_expression.csv", row.names = FALSE)
  cat("✓ Differential expression analysis results saved to: target_genes_differential_expression.csv\n")
  
  # Print differential analysis summary
  cat("\nDifferential Expression Analysis Summary:\n")
  print(de_results[, c("Gene", "Log2_FoldChange", "P_adjusted", "Significance")], row.names = FALSE)
}

# --------------------------
# 21. Violin plot of target gene expression (by cell type and group)
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 21: Generate violin plots of target gene expression\n")
cat("=" * 70, "\n")

# Prepare expression data for plotting
plot_data_list <- list()

for (gene in available_genes) {
  # Extract expression and metadata
  expr_values <- as.vector(expr_data[gene, ])
  plot_df <- data.frame(
    CellType = seurat_obj$cell_type,
    Group = seurat_obj$group,
    Expression = expr_values
  )
  
  # Filter out zero-expression cells (optional, can be adjusted)
  plot_df <- plot_df[plot_df$Expression > 0, ]
  
  if (nrow(plot_df) > 0) {
    plot_df$Gene <- gene
    plot_data_list[[gene]] <- plot_df
  } else {
    cat(sprintf("Warning: No positive expression data for %s, skipping violin plot\n", gene))
  }
}

# Combine plot data
if (length(plot_data_list) > 0) {
  plot_data <- do.call(rbind, plot_data_list)
  
  # 1. Violin plot by cell type
  cat("\nGenerating violin plot by cell type...\n")
  p_violin_celltype <- ggplot(plot_data, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_jitter(size = 0.5, alpha = 0.3, color = "black") +
    facet_wrap(~Gene, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Target Gene Expression Distribution Across Cell Types",
      x = "Cell Type",
      y = "Normalized Expression"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_fill_viridis_d(option = "plasma")
  
  ggsave("violinplot_target_genes_by_celltype.png", p_violin_celltype, 
         width = 14, height = 8, dpi = 300, bg = "white")
  cat("✓ Violin plot by cell type saved to: violinplot_target_genes_by_celltype.png\n")
  
  # 2. Violin plot by group (split by cell type)
  if (all(c("Normal", "DFU") %in% plot_data$Group)) {
    cat("\nGenerating violin plot by group (split by cell type)...\n")
    
    # Select top 5 cell types with the most cells for clarity
    top_celltypes <- names(sort(table(plot_data$CellType), decreasing = TRUE))[1:5]
    plot_data_sub <- plot_data[plot_data$CellType %in% top_celltypes, ]
    
    p_violin_group <- ggplot(plot_data_sub, aes(x = Group, y = Expression, fill = Group)) +
      geom_violin(alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
      facet_grid(Gene ~ CellType, scales = "free_y") +
      theme_minimal(base_size = 12) +
      labs(
        title = "Target Gene Expression in Normal vs DFU Groups (Top 5 Cell Types)",
        x = "Group",
        y = "Normalized Expression"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "none",
        strip.text = element_text(face = "bold", size = 10)
      ) +
      scale_fill_manual(values = c("Normal" = "#2E86AB", "DFU" = "#A23B72"))
    
    ggsave("violinplot_target_genes_by_group_celltype.png", p_violin_group, 
           width = 16, height = 10, dpi = 300, bg = "white")
    cat("✓ Violin plot by group and cell type saved to: violinplot_target_genes_by_group_celltype.png\n")
  }
}

# --------------------------
# 22. Final result summary and save all objects
# --------------------------
cat("\n" + "=" * 70, "\n")
cat("Step 22: Final result summary and save all analysis objects\n")
cat("=" * 70, "\n")

# Create analysis summary text
summary_text <- paste0(
  "Single-Cell Transcriptome Analysis Summary for DFU Research\n",
  "=========================================================\n",
  "Analysis Date: ", Sys.Date(), "\n",
  "Total Cells: ", ncol(seurat_obj), "\n",
  "Total Cell Types: ", length(unique(seurat_obj$cell_type)), "\n",
  "Sample Groups: ", paste(unique(seurat_obj$group), collapse = ", "), "\n",
  "Target Genes Analyzed: ", paste(available_genes, collapse = ", "), "\n\n",
  "Cell Type Distribution:\n"
)

# Add cell type distribution
celltype_dist <- as.data.frame(table(seurat_obj$cell_type))
colnames(celltype_dist) <- c("CellType", "Count")
celltype_dist$Percentage <- round(celltype_dist$Count / sum(celltype_dist$Count) * 100, 2)

for (i in 1:nrow(celltype_dist)) {
  summary_text <- paste0(summary_text,
                         "  - ", celltype_dist$CellType[i], ": ", 
                         celltype_dist$Count[i], " cells (", 
                         celltype_dist$Percentage[i], "%)\n")
}

# Add differential expression summary if available
if (exists("de_results") && nrow(de_results) > 0) {
  summary_text <- paste0(summary_text, "\nDifferential Expression Results (Normal vs DFU):\n")
  for (i in 1:nrow(de_results)) {
    summary_text <- paste0(summary_text,
                           "  - ", de_results$Gene[i], ": Log2FC = ", 
                           round(de_results$Log2_FoldChange[i], 3), 
                           ", P-adjusted = ", format(de_results$P_adjusted[i], scientific = TRUE),
                           " (", de_results$Significance[i], ")\n")
  }
}

# Save summary text
writeLines(summary_text, con = "analysis_summary.txt")
cat("✓ Analysis summary saved to: analysis_summary.txt\n")

# Save all key objects into a single RData file
save(
  seurat_obj, seurat_dfu, marker_matrix, cell_type_summary,
  group_summary, de_results, stats_report,
  file = "dfu_singlecell_analysis_results.RData"
)
cat("✓ All analysis objects saved to: dfu_singlecell_analysis_results.RData\n")

# --------------------------
# Analysis completion message
