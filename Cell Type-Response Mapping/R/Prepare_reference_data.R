# prepare_reference_data.R
# One-time script to prepare the scRNA-seq reference dataset from the Allen Brain Cell Atlas.
# Run this once when the reference data needs to be recreated.
# Outputs: merfish_depth_metadata.csv, mop_neurons.rds, and sub_mop_neurons.rds in the reference data directory.

library(Seurat)
library(dplyr)
library(hdf5r)
library(SeuratDisk)

ref_dir <- 'input'  # relative to the repo root; run this script from there

# Load raw Allen Brain Atlas isocortex data
ctx_object <- LoadH5Seurat(file.path(ref_dir, 'raw_abc/WMB-10Xv2-Isocortex-1-raw.h5seurat'),
                           meta.data = FALSE, commands = FALSE, misc = FALSE, tools = FALSE)

metadata <- read.csv(file.path(ref_dir, 'raw_abc/cell_metadata.csv'))
rownames(metadata) <- metadata$cell_label

# Subset to MOp (primary motor cortex) cells
metadata_mop <- metadata[metadata$region_of_interest_acronym == 'MOp',]
mop_object <- ctx_object[, colnames(ctx_object) %in% metadata$cell_label]

# Map cluster annotations from the Allen taxonomy
library_data <- read.csv(file.path(ref_dir, 'raw_abc/cluster_to_cluster_annotation_membership.csv'))
anno_library <- matrix(
  c(library_data[library_data$cluster_annotation_term_set_name == 'cluster', 'cluster_annotation_term_name'],
    library_data[library_data$cluster_annotation_term_set_name == 'supertype', 'cluster_annotation_term_name'],
    library_data[library_data$cluster_annotation_term_set_name == 'subclass', 'cluster_annotation_term_name'],
    library_data[library_data$cluster_annotation_term_set_name == 'class', 'cluster_annotation_term_name'],
    library_data[library_data$cluster_annotation_term_set_name == 'neurotransmitter', 'cluster_annotation_term_name']),
  nrow = length(unique(library_data$cluster_alias)),
  ncol = 5,
  dimnames = list(unique(library_data$cluster_alias),
                  c('cluster', 'supertype', 'subclass', 'class', 'neurotransmitter')))

mop_object[['allen_cluster']] <- as.character(anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'cluster'])
mop_object[['allen_supertype']] <- as.character(anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'supertype'])
mop_object[['allen_subclass']] <- as.character(anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'subclass'])
mop_object[['allen_class']] <- as.character(anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'class'])
mop_object[['allen_NT']] <- as.character(anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'neurotransmitter'])

# Normalize, cluster, and embed
mop_object <- NormalizeData(mop_object) %>%
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims=1:30) %>%
  FindClusters(., resolution=1) %>%
  RunUMAP(., dims=1:30)

# Filter to cortical neuron types only
types_wanted <- unique(mop_object$allen_subclass)
types_wanted <- c(types_wanted[grepl('CTX', types_wanted)],
                  types_wanted[grepl('0(53)|(56)|(47)|(49)|(52)|(51)|(46)|(50) .* Gaba', types_wanted)])
types_wanted <- types_wanted[!grepl('(CLA)|(STR)', types_wanted)]
mop_neurons <- mop_object[, mop_object$allen_subclass %in% types_wanted]

# Re-process the neuron subset
mop_neurons <- NormalizeData(mop_neurons) %>%
  FindVariableFeatures(.) %>%
  ScaleData(.) %>%
  RunPCA(.) %>%
  FindNeighbors(., dims=1:30) %>%
  FindClusters(., resolution=1) %>%
  RunUMAP(., dims=1:30)

# Optional: Create custom type labels
#mop_neurons[['my_supertype']] <- mop_neurons$allen_supertype
#mop_neurons$my_supertype[!grepl('L\\d', mop_neurons$my_supertype)] <- mop_neurons$allen_subclass[!grepl('L\\d', mop_neurons$my_supertype)]

# Generate color palette
# supertypes <- sort(unique(mop_neurons$my_supertype))
# cluster_colors <- c(
#   colorRampPalette(c('#d4775e', '#a33e59', '#765d8c', '#79a4b9', '#afd435', '#fdec28'))(
#     length(supertypes[grepl('L\\d', supertypes)])),
#   colorRampPalette(c('#ffce74', '#ff9fa1', '#ffa9e1', '#ca94ff', '#94acff'))(
#     length(supertypes[!grepl('L\\d', supertypes)])))
# names(cluster_colors) <- c(supertypes[grepl('L\\d', supertypes)],
#                            supertypes[!grepl('L\\d', supertypes)])

# Process MERFISH metadata
merfish_metadata <- read.table(file.path(ref_dir, 'raw_abc/merfish_metadata.csv'), sep = ',', header = T)

merfish_metadata[, 'allen_class'] <- anno_library[as.character(merfish_metadata$cluster_alias), 'class']
merfish_metadata[, 'allen_subclass'] <- anno_library[as.character(merfish_metadata$cluster_alias), 'subclass']
merfish_metadata[, 'allen_supertype'] <- anno_library[as.character(merfish_metadata$cluster_alias), 'supertype']
merfish_metadata[, 'allen_NT'] <- anno_library[as.character(merfish_metadata$cluster_alias), 'neurotransmitter']
merfish_metadata[, 'my_supertype'] <- merfish_metadata$allen_supertype
merfish_metadata$my_supertype[!grepl('L\\d', merfish_metadata$my_supertype)] <- merfish_metadata$allen_subclass[!grepl('L\\d', merfish_metadata$my_supertype)]
  
ctx_merdata <- merfish_metadata[(merfish_metadata$z == '10.2' | merfish_metadata$z == '10.4'),] # These slices contain MOp data
excit_merdata <- ctx_merdata[grepl('Glut', ctx_merdata$allen_class),]
mop_merdata_right <- ctx_merdata[(ctx_merdata$x < 7.75 & ctx_merdata$x > 6.625 & ctx_merdata$y > 2.5 & ctx_merdata$y < 4.75),]
mop_merdata_left <- ctx_merdata[(ctx_merdata$x < 4.3 & ctx_merdata$x > 3.4 & ctx_merdata$y > 2.75 & ctx_merdata$y < 5.0),]

# Rotate MERFISH cells to make axial cortical depths
origin.xr <- mean(mop_merdata_right$x)
origin.yr <- mean(mop_merdata_right$y)
angler <- -18*pi/180

origin.xl <- mean(mop_merdata_left$x)
origin.yl <- mean(mop_merdata_left$y)
anglel <- 28*pi/180

new.xr <- (mop_merdata_right$x-origin.xr)*cos(angler) - (mop_merdata_right$y-origin.yr)*sin(angler) + origin.xr
new.yr <- (mop_merdata_right$x-origin.xr)*sin(angler) + (mop_merdata_right$y-origin.yr)*cos(angler) + origin.yr
new.xl <- (mop_merdata_left$x-origin.xl)*cos(anglel) - (mop_merdata_left$y-origin.yl)*sin(anglel) + origin.xl
new.yl <- (mop_merdata_left$x-origin.xl)*sin(anglel) + (mop_merdata_left$y-origin.yl)*cos(anglel) + origin.yl

new.xr <- (new.xr - min(new.xr))/(max(new.xr) - min(new.xr))
new.yr <- (new.yr - min(new.yr))/(max(new.yr) - min(new.yr))
new.xl <- (-new.xl - min(-new.xl))/(max(-new.xl) - min(-new.xl))
new.yl <- (new.yl - min(new.yl))/(max(new.yl) - min(new.yl))

mop_merdata_right$newx <- new.xr
mop_merdata_right$newy <- new.yr
mop_merdata_left$newx <- new.xl
mop_merdata_left$newy <- new.yl

mop_merdata <- rbind(mop_merdata_right, mop_merdata_left)

# Save MERFISH data

write.table(mop_merdata, file.path(ref_dir, 'processed_abc/merfish_depth_metadata.csv'), sep = ',', col.names = T, row.names = F)

# Save a smaller RNA-Seq subsample for quick visualization
sub_mop_neurons <- mop_neurons[, sample(1:ncol(mop_neurons), ncol(mop_neurons)/4, replace = FALSE)]

saveRDS(mop_neurons, file.path(ref_dir, 'processed_abc/mop_neurons.rds'))
saveRDS(sub_mop_neurons, file.path(ref_dir, 'processed_abc/sub_mop_neurons.rds'))

cat("Reference data saved to:", ref_dir, "\n")
