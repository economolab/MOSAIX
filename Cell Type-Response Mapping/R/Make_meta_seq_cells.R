# make_meta_seq_cells.R
# Full Stream 2 pipeline: create meta cells from scRNA-seq, fit GMMs, binarize.
# Outputs: intermediate/All_genes_bin.csv, intermediate/GMM_Thresholds.csv,
#          intermediate/meta_log_counts.csv, intermediate/meta_norm_counts.csv

library(dplyr)
library(Seurat)
library(ggplot2)

source('R/config_loader.R')
source('R/fit_gmm_thresholds.R')
source('R/classifier_functions.R')

# --- Configuration ---
paths <- load_paths()

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Slco2a1',
               'Fam84b', 'Npnt', 'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')

# Per-gene component counts matching the original Python GMM fits.
# Derived from the existing GMM_Thresholds.csv: count non-NA means per gene.
n_components_list <- list(
  Penk = 3, Calb1 = 3, Lamp5 = 3, Rorb = 3, Ccdc80 = 3, Pamr1 = 4,
  Dkkl1 = 2, Slco2a1 = 2, Fam84b = 3, Npnt = 3, Myl4 = 3,
  Syt6 = 2, Ctgf = 2, Slc32a1 = 2, Pvalb = 3, Sst = 2, Vip = 3
)

# Optional: provide initial means for specific genes that need guided fitting
# (set to NULL for automatic selection)
means_init_list <- list(
  Penk = c(0, 1.8, 3.6),
  Calb1 = c(0.05, 2.9, 4.4),
  Lamp5 = c(0.2, 2.2, 3.6),
  Npnt = c(0.1, 1.6, 3.0),
  Myl4 = c(0.25, 3, 4.5))
# Example:
# means_init_list <- list(
#   Ctgf = c(0.05, 2.2, 3.2)
# )

# --- Load reference data ---
cat('Loading scRNA-seq reference object...\n')
mop_object <- readRDS(file.path(paths$reference_data_dir, 'sub_mop_neurons.rds'))
merfish_data <- read.table(file.path(paths$reference_data_dir, 'merfish_depth_metadata.csv'),
                           sep = ',', header = TRUE)

# Create the my_supertype labels in the RNA-Seq data
mop_object[['my_supertype']] <- mop_object$allen_supertype
mop_object$my_supertype[!grepl('L\\d', mop_object$my_supertype)] <-
  mop_object$allen_subclass[!grepl('L\\d', mop_object$my_supertype)]

# --- Step 1: Create meta cells ---
cat('Finding neighbors...\n')
mop_object <- FindNeighbors(mop_object, reduction = 'pca', dims = 1:30, return.neighbor = TRUE, k = 500)

neighbor_dists <- mop_object@neighbors$RNA.nn@nn.dist[, 2:500]
rownames(neighbor_dists) <- colnames(mop_object)
neighbor_idx <- mop_object@neighbors$RNA.nn@nn.idx[, 2:500]
rownames(neighbor_idx) <- colnames(mop_object)

# Only include clusters with sufficient cells and present in MERFISH data
type_counts <- table(mop_object$my_supertype)
clusters_to_include <- names(type_counts[(type_counts > 80) &
                                          (names(type_counts) %in% merfish_data$my_supertype)])
n_clusters <- length(clusters_to_include)

cat(sprintf('Creating meta cells for %d clusters...\n', n_clusters))
meta_norm_counts <- matrix(nrow = length(hcr_genes), ncol = sum(type_counts[clusters_to_include]))

count <- 0
types_kept <- rep(NA_character_, sum(type_counts[clusters_to_include]))

for (cluster_type in clusters_to_include) {

  sub_counts <- mop_object@assays$RNA@data[, mop_object$my_supertype == cluster_type]
  sub_counts <- expm1(sub_counts)
  sub_neighbors <- neighbor_idx[mop_object$my_supertype == cluster_type, ]

  set.seed(10)
  sample_meta_cell <- sample(rownames(sub_neighbors), type_counts[cluster_type], replace = FALSE)

  for (j in 1:type_counts[cluster_type]) {
    count <- count + 1
    neighbors <- neighbor_idx[sample_meta_cell[j], ]
    neighbor_names <- colnames(mop_object)[neighbors]
    neighbors_to_keep <- neighbor_names[neighbor_names %in% colnames(sub_counts)]

    if (length(neighbors_to_keep) < 9) {
      warning(sprintf("Cluster '%s': only %d same-type neighbors found (minimum 9 required)",
                      cluster_type, length(neighbors_to_keep)))
    } else {
      mean_counts <- apply(as.matrix(sub_counts[hcr_genes, neighbors_to_keep[1:9]]), 1, sum)
      mean_counts <- mean_counts + sub_counts[hcr_genes, sample_meta_cell[j]]
      types_kept[count] <- cluster_type
      meta_norm_counts[, count] <- mean_counts
    }
  }
}

meta_norm_counts <- t(meta_norm_counts)
meta_norm_counts <- as.data.frame(meta_norm_counts)
colnames(meta_norm_counts)[1:17] <- hcr_genes
meta_norm_counts[, 'my_supertype'] <- types_kept
meta_norm_counts <- meta_norm_counts[!is.na(meta_norm_counts$my_supertype), ]

# --- Step 2: Add taxonomy metadata ---
all_types <- unique(meta_norm_counts$my_supertype)
cluster_to_subclass <- character(length(all_types))
cluster_to_class <- character(length(all_types))
cluster_to_nt <- character(length(all_types))

for (i in 1:length(all_types)) {
  cluster_type <- all_types[i]
  idxs <- mop_object$my_supertype == cluster_type
  subclasses <- unique(mop_object$allen_subclass[idxs])
  classes <- unique(mop_object$allen_class[idxs])
  nt <- unique(mop_object$allen_NT[idxs])

  if (length(nt) > 1) {
    warning(sprintf("Cluster '%s' has multiple neurotransmitter types: %s", cluster_type, paste(nt, collapse = ', ')))
  }

  cluster_to_subclass[i] <- subclasses[1]
  cluster_to_class[i] <- classes[1]
  cluster_to_nt[i] <- nt[1]
}
names(cluster_to_subclass) <- all_types
names(cluster_to_class) <- all_types
names(cluster_to_nt) <- all_types

meta_norm_counts[, 'allen_subclass'] <- cluster_to_subclass[meta_norm_counts$my_supertype]
meta_norm_counts[, 'allen_class'] <- cluster_to_class[meta_norm_counts$my_supertype]
meta_norm_counts[, 'allen_NT'] <- cluster_to_nt[meta_norm_counts$my_supertype]

# Filter out ambiguous types
meta_norm_counts <- meta_norm_counts[!meta_norm_counts$allen_NT == 'Glut-GABA', ]
meta_norm_counts <- meta_norm_counts[!meta_norm_counts$allen_class == '08 CNU-MGE GABA', ]

# --- Step 3: Log-transform ---
meta_log_counts <- meta_norm_counts
meta_log_counts[, 1:17] <- log(meta_log_counts[, 1:17])

# --- Step 4: Fit GMMs and binarize all genes ---
cat('Fitting GMMs for all genes...\n')
gmm_results <- fit_all_genes_gmm(meta_log_counts, hcr_genes,
                                  means_init_list = means_init_list,
                                  n_components_list = n_components_list)

all_genes_bin <- gmm_results$all_genes_bin
gmm_thresholds <- gmm_results$gmm_thresholds

# --- Step 5: Save outputs ---
dir.create(paths$intermediate_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(meta_norm_counts, file.path(paths$intermediate_dir, 'meta_norm_counts.csv'), row.names = FALSE)
cat('Saved:', file.path(paths$intermediate_dir, 'meta_norm_counts.csv'), '\n')

write.csv(meta_log_counts, file.path(paths$intermediate_dir, 'meta_log_counts.csv'), row.names = FALSE)
cat('Saved:', file.path(paths$intermediate_dir, 'meta_log_counts.csv'), '\n')

write.csv(all_genes_bin, file.path(paths$intermediate_dir, 'All_genes_bin.csv'), row.names = FALSE)
cat('Saved:', file.path(paths$intermediate_dir, 'All_genes_bin.csv'), '\n')

write.csv(gmm_thresholds, file.path(paths$intermediate_dir, 'GMM_Thresholds.csv'))
cat('Saved:', file.path(paths$intermediate_dir, 'GMM_Thresholds.csv'), '\n')

cat(sprintf('\nDone. Processed %d meta cells across %d clusters for %d genes.\n',
            nrow(meta_norm_counts), n_clusters, length(hcr_genes)))
cat(sprintf('GMM components per gene:\n'))
for (gene in hcr_genes) {
  cat(sprintf('  %s: %d components, thresholds = %s\n',
              gene, gmm_results$per_gene_results[[gene]]$n_components,
              paste(round(gmm_results$per_gene_results[[gene]]$thresholds, 3), collapse = ', ')))
}
