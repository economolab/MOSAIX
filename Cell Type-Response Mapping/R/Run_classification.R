# run_classification.R
# Main entry point for running HCR cell type classification on an experiment.
#
# Usage: Set experiment_config_file below and source this script in RStudio,
# or run from command line: Rscript R/run_classification.R config/experiments/2024-07-30_thalamus.json

library(dplyr)
library(Seurat)
library(ggplot2)

# --- Configuration ---
# Set this to the experiment you want to classify:
experiment_config_file <- 'config/experiments/2024-07-30_thalamus.json'

# Override with command line argument if provided
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  experiment_config_file <- args[1]
}

# --- Source pipeline modules ---
source('R/config_loader.R')
source('R/classifier_functions.R')
source('R/hcr_data_loader.R')

# --- Load configuration ---
config <- load_experiment_config(experiment_config_file)
paths <- config$paths

cat('Running classification for experiment:', config$experiment_label, '\n')

# --- Load reference data ---
cat('Loading reference data...\n')
#mop_neurons <- readRDS(file.path(paths$reference_data_dir, 'sub_mop_neurons.rds'))
merfish_metadata <- read.table(file.path(paths$reference_data_dir, 'merfish_depth_metadata.csv'),
                               sep = ',', header = TRUE)

# Load GMM-binarized reference gene expression
meta_gene_matrix <- read.table(file.path(paths$intermediate_dir, 'meta_log_counts.csv'),
                                sep = ',', header = TRUE)
meta_gene_matrix <- CombineGenes(meta_gene_matrix, config$genes_combined)

mop_object_hcr <- meta_gene_matrix
mop_object_hcr[, 'all_cells'] <- rep('all', nrow(mop_object_hcr))

meta_bin <- read.table(file.path(paths$intermediate_dir, 'All_genes_bin.csv'),
                        sep = ',', header = TRUE)
meta_bin <- as.data.frame(meta_bin)
meta_bin <- meta_bin[, config$hcr_genes[config$hcr_genes != 'JF585']]

# Attach taxonomy metadata from the log-counts matrix
meta_bin <- cbind(meta_bin, mop_object_hcr[, (ncol(mop_object_hcr)-4):ncol(mop_object_hcr)])

# Convert multi-level GMM bins to binary (on/off)
# Genes with 3+ Gaussian components get special handling:
# some use the highest bin, others use the second-highest
high_bin_genes <- c('Rorb', 'Fam84b', 'Ccdc80', 'Pvalb', 'Vip', 'Npnt')
for (i in 1:(ncol(meta_bin) - 5)) {
  if (2 %in% meta_bin[, i]) {
    if (colnames(meta_bin)[i] %in% high_bin_genes) {
      meta_bin[, i] <- meta_bin[, i] == max(meta_bin[, i])
    } else {
      meta_bin[, i] <- meta_bin[, i] == max(meta_bin[, i]) - 1
    }
  } else {
    meta_bin[, i] <- meta_bin[, i] > 0
  }
}
meta_bin[, 1:(ncol(meta_bin) - 5)] <- meta_bin[, 1:(ncol(meta_bin) - 5)] * 1
meta_bin <- CombineGenes(meta_bin, config$genes_combined)

mop_object_hcr <- meta_bin
mop_object_hcr[, 'all_cells'] <- rep('all', nrow(mop_object_hcr))

# Build subclass lookup dictionary
subclass_rows <- unique(meta_bin[, c('my_supertype', 'allen_subclass')])
subclass_dict <- subclass_rows$allen_subclass
names(subclass_dict) <- subclass_rows$my_supertype

# --- Load HCR imaging data ---
cat('Loading HCR imaging data...\n')
hcr_data <- load_hcr_data(config)
mean.vals <- hcr_data$mean.vals
ilastik.anno <- hcr_data$ilastik.anno
voltron.anno <- hcr_data$voltron.anno

# Align ilastik annotation columns with meta_bin
ilastik.anno <- ilastik.anno[, colnames(meta_bin)[1:(ncol(meta_bin) - 5)]]
ilastik.anno$newy <- mean.vals$newy

# --- Run classification ---
cat('Running HierarchicalGMMBayes classifier...\n')
bayes.classif <- HierarchicalGMMBayes(
  ilastik.anno, meta_bin,
  merfish_data = merfish_metadata,
  mop_object_hcr = mop_object_hcr,
  use.lamina = TRUE,
  depth.round = 3
)

bayes.classif$classif.label[is.na(bayes.classif$classif.label)] <- 'Unknown'
names(bayes.classif$classif.label) <- mean.vals$cell.id
names(bayes.classif$classif.prob) <- mean.vals$cell.id

# Sort classification labels
order_labels <- sub('\\d{3,4} (.*)', '\\1', unique(bayes.classif$classif.label))
order_nums <- as.numeric(sub('(\\d{3,4}) .*', '\\1', unique(bayes.classif$classif.label)))
label_order <- order(order_labels, order_nums)
bayes.classif$classif.label <- factor(bayes.classif$classif.label,
                                       levels = unique(bayes.classif$classif.label)[label_order])

## Colors definition:

excit.types <- unique(mop_object_hcr$my_supertype)[grepl('L\\d', unique(mop_object_hcr$my_supertype)) & unique(mop_object_hcr$my_supertype) %in% merfish_metadata$my_supertype] # 'L\\d( |/)' Exclude L6b cells by adding the space after the number
inhib.types <- unique(mop_object_hcr$allen_subclass)[!(grepl('L\\d', unique(mop_object_hcr$allen_subclass)))]
all.types <- c(excit.types, inhib.types)
col_order <- sub('\\d{3,4} (.*)', '\\1', all.types)
col_order2 <- as.numeric(sub('(\\d{3,4}) .*', '\\1', all.types))
col_order <- order(col_order, col_order2)
col_names <- factor(all.types, levels =  all.types[col_order])

cluster_colors <- c(colorRampPalette(c('#d4775e', '#a33e59', '#765d8c', '#79a4b9', '#afd435', '#fdec28'), bias = 1.3)(length(col_names[grepl('L\\d', col_names)])),
                    colorRampPalette(c('#ffce74', '#ff9fa1', '#ffa9e1', '#ca94ff', '#94acff'), bias = 1)(length(col_names[!grepl('L\\d', col_names)])))
names(cluster_colors) <- levels(col_names)
other_colors <- c('#b4b4b4', '#faebcf', '#edd2ca', '#d4bcc2', '#d5cade', '#d5caee', '#dbe3c1', '#ede9be', '#e8e8e8')
names(other_colors) <- c('9990 VVgat', '9991 ZNo markers L1', '9992 ZNo markers L2/3', 
                         '9993 ZNo markers L4/5a', '9994 ZNo markers UL5b', '9995 ZNo markers LL5b', 
                         '9996 ZNo markers L6a', '9997 ZNo markers L6b', '9999 ZNo Good Match')
cluster_colors <- c(cluster_colors, other_colors)
colors <- list(cell_type = cluster_colors)

# --- Write output ---
subclasses <- sapply(as.character(bayes.classif$classif.label),
                     function(x) ifelse(x %in% names(subclass_dict), subclass_dict[x], x))

table.to.write <- data.frame(
  cell_ids = as.numeric(names(bayes.classif$classif.label)),
  slice = mean.vals$slice,
  depth = mean.vals$newy,
  supertype = as.character(bayes.classif$classif.label),
  subclass = subclasses,
  classif_cols = cluster_colors[as.character(bayes.classif$classif.label)],
  slice_depth = mean.vals$Z
)

output_path <- file.path(config$output_dir, 'HCR_classification.csv')
write.table(table.to.write, output_path, sep = ',', col.names = TRUE, row.names = FALSE)
cat('Classification output written to:', output_path, '\n')

# Also write to the experiment's annotation directory for MATLAB compatibility
legacy_path <- file.path(config$anno_dir, 'HCR_classification.csv')
write.table(table.to.write, legacy_path, sep = ',', col.names = TRUE, row.names = FALSE)
cat('Legacy output written to:', legacy_path, '\n')

cat('Classification complete.\n')
cat('Classified', nrow(table.to.write), 'cells into', length(unique(bayes.classif$classif.label)), 'types.\n')
