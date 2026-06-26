library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(reshape2)

# --- Configuration ---
# Set this to the experiment you want to generate figures for:
experiment_config_file <- 'config/experiments/2024-07-30_thalamus.json'

source('R/config_loader.R')
source('R/classifier_functions.R')
source('R/hcr_data_loader.R')

config <- load_experiment_config(experiment_config_file)
paths <- config$paths

# --- Load HCR classification table ---
classif_path <- file.path(config$output_dir, 'HCR_classification.csv')
hcr_classif <- read.csv(classif_path)

# --- Load HCR imaging data (for spatial coordinates and gene expression) ---
hcr_data <- load_hcr_data(config)
mean.vals <- hcr_data$mean.vals
ilastik.anno <- hcr_data$ilastik.anno

# Align ilastik annotation columns with reference genes
meta_bin <- read.table(file.path(paths$intermediate_dir, 'All_genes_bin.csv'),
                       sep = ',', header = TRUE)
meta_bin <- as.data.frame(meta_bin)
meta_bin <- meta_bin[, config$hcr_genes[config$hcr_genes != 'JF585']]

meta_gene_matrix <- read.table(file.path(paths$intermediate_dir, 'meta_log_counts.csv'),
                               sep = ',', header = TRUE)
meta_gene_matrix <- CombineGenes(meta_gene_matrix, config$genes_combined)
mop_object_hcr <- meta_gene_matrix
mop_object_hcr[, 'all_cells'] <- rep('all', nrow(mop_object_hcr))

# Convert multi-level GMM bins to binary
high_bin_genes <- c('Rorb', 'Fam84b', 'Ccdc80', 'Pvalb', 'Vip', 'Npnt')
for (i in 1:(ncol(meta_bin))) {
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
meta_bin[, 1:(ncol(meta_bin))] <- meta_bin[, 1:(ncol(meta_bin))] * 1
meta_bin <- CombineGenes(meta_bin, config$genes_combined)

ilastik.anno <- ilastik.anno[, colnames(meta_bin)]
ilastik.anno$newy <- mean.vals$newy

# --- Load MERFISH reference data ---
merfish_metadata <- read.table(file.path(paths$reference_data_dir, 'merfish_depth_metadata.csv'),
                               sep = ',', header = TRUE)

# --- Build classification label vector (replaces bayes.classif$classif.label) ---
classif_label <- hcr_classif$supertype
names(classif_label) <- hcr_classif$cell_ids

# Sort classification labels
order_labels <- sub('\\d{3,4} (.*)', '\\1', unique(classif_label))
order_nums <- as.numeric(sub('(\\d{3,4}) .*', '\\1', unique(classif_label)))
label_order <- order(order_labels, order_nums)
classif_label <- factor(classif_label,
                        levels = unique(classif_label)[label_order])

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

#### Scatterplot of types: #####

excit.types <- grepl('L\\d', classif_label) # 'L\\d( |/)' Exclude L6b cells by adding the space after the number
inhib.types <- !grepl('L\\d', classif_label) | !grepl('ZNo', classif_label)
custom.types <- grepl('ET', classif_label)
type <- c(excit.types|inhib.types)
big_groups <- names(table(classif_label))[table(classif_label) > 10]
scatter_data <- data.frame(classif=classif_label, X=mean.vals$X, Y=-mean.vals$Y, Z=mean.vals$Z)
scatter_data <- scatter_data[scatter_data$classif %in% big_groups,]

ggplot(scatter_data[!grepl('ZNo ', scatter_data$classif),], aes(x=X, y=Y, col=classif)) +
  geom_point() +
  ylim(c(min(-config$img_lims$ymax), max(-config$img_lims$ymin))) +
  theme_minimal() +
  scale_color_manual(values = cluster_colors) +
  guides(color = 'none')

## Plotting laminar position of cells, violin. Only works if you did the rotation at the start

excit.types <- grepl('L\\d', classif_label)
inhib.types <- !grepl('L\\d', classif_label) + grepl('ZNo|Err|all', classif_label)
custom.types <- grepl('L6b', classif_label)
type <- c(excit.types | inhib.types)
big_groups <- names(table(classif_label))[table(classif_label) > 10]

rel.depths <- data.frame(Y = ilastik.anno$newy[type])
rel.depths$use_type <- classif_label[type]
rel.depths <- rel.depths[rel.depths$use_type %in% big_groups,]

ylims <- rel.depths %>%
  group_by(use_type) %>%
  summarise(Q1 = quantile(Y, 0.05, na.rm = T), Q9 = quantile(Y, 0.95, na.rm = T)) # Only plot the middle 80% of cells

rel.depths.mod <- data.frame()
for (i in 1:nrow(ylims)) {
  max.Y <- as.numeric(ylims[i,'Q9'])
  min.Y <- as.numeric(ylims[i,'Q1'])
  rel.depths.rows <- rel.depths$use_type == ylims[i,][['use_type']]
  rel.depths.inquant <- rel.depths[(rel.depths.rows & (rel.depths$Y >= min.Y) & (rel.depths$Y < max.Y)),]
  rel.depths.mod <- rbind(rel.depths.mod, rel.depths.inquant)
}

ggplot(rel.depths.mod, aes(x=use_type, y=Y, fill=use_type)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill ='white'), axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
  geom_hline(yintercept = c(-0.075, -0.25, -0.4, -0.55, -0.7, -0.92, -1), color = 'grey80') +
  geom_violin(trim=T, scale='width') +
  ylim(c(-1.0, 0)) +
  scale_fill_manual(values=cluster_colors) +
  guides(fill = 'none')

#### Comparison of HCR population to MERFISH: ####

merfish.cells <- ExtractMerfishLamina(merfish_metadata, 'my_supertype', mop_object_hcr, bin = F, bin_to = 2)
type_densities <- merfish.cells[[1]]
rel_depths <- merfish.cells[[2]]
prop_depths <- merfish.cells[[3]]

merfish_type_props <- table(rel_depths$use_type)/nrow(rel_depths)
all_types <- classif_label[!grepl('(ZNo)|(VV)', classif_label)]
all_type_props <- table(all_types)/length(all_types)

voltron.cells <- as.character(classif_label)
voltron.cells <- voltron.cells[!grepl('(ZNo)|(VV)', voltron.cells)]
unique.voltron <- unique(voltron.cells)
voltron.counts <- c()
for (i in 1:length(unique.voltron)) {
  voltron.counts <- c(voltron.counts, sum(voltron.cells == unique.voltron[i]))
}
names(voltron.counts) <- unique.voltron
voltron.props <- voltron.counts/length(voltron.cells)

both.props <- data.frame(voltron = c(rep('Voltron', length(voltron.props)), rep('Merfish', length(merfish_type_props))),
                         class = factor(c(names(voltron.props), names(merfish_type_props)), levels = names(cluster_colors)[1:39]),
                         prop = c(voltron.props, merfish_type_props))

ggplot(both.props, aes(x = voltron, y = prop, fill = class)) +
  theme_minimal() +
  geom_bar(position = 'stack', stat = 'identity', color = 'white') +
  scale_fill_manual(values=cluster_colors)

## MERFISH depth prior bar plot (or density plot):

prop_cast <- dcast(data = prop_depths, formula = Y ~ use_type, fun.aggregate = sum, value.var = 'prop')
prop_cast_adj <- prop_cast + 0.002
prop_melt <- melt(prop_cast_adj, id.vars = 'Y', variable.name = 'use_type', value.name = 'prop')

smooth_prop <- prop_melt %>% group_by(use_type) %>% arrange(Y) %>%
  mutate(smooth_prop = loess(prop ~ Y, span = 0.16)$fitted) %>% ungroup() %>%
  group_by(Y) %>% mutate(smooth_prop = pmax(smooth_prop, 0.002),
                         smooth_prop = smooth_prop/sum(smooth_prop))

ggplot(smooth_prop, aes(x = Y, y = smooth_prop, group = use_type, fill = use_type)) +
  geom_area(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  xlim(-1, 0) +
  coord_cartesian(ylim = c(0, 1), expand = F) +
  guides(fill = 'none')

#### Dot plot showing proportion of HCR cells in each type that express each gene. To compare to above

bayes.dot.df <- cbind(ilastik.anno[1:(ncol(ilastik.anno)-1)]*1, classif_label)
names(bayes.dot.df)[ncol(bayes.dot.df)] <- 'classif.label'
bayes.dot.tidy <- bayes.dot.df %>% group_by(classif.label) %>% summarise(across(colnames(bayes.dot.df)[1:16], sum))
bayes.dot.n <- bayes.dot.df %>% group_by(classif.label) %>% summarise(n = n())
bayes.dot.prop <- mapply(`/`, data.frame(bayes.dot.tidy[,2:ncol(bayes.dot.tidy)]), bayes.dot.n[,2])
colnames(bayes.dot.prop) <- sub('\\.', '\\+', colnames(bayes.dot.prop))

bayes.dot.prop <- data.frame(bayes.dot.prop[,colnames(ilastik.anno[,1:(ncol(ilastik.anno)-1)])],
                             id = bayes.dot.n[,1])
bayes.dot.melt <- melt(bayes.dot.prop)

ggplot(bayes.dot.melt, aes(x = classif.label, y = variable, size = value, color = classif.label)) +
  geom_point(stroke = 0) +
  theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, hjust =1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = cluster_colors) +
  scale_size(range = c(-0.05, 8)) +
  guides(color = 'none')

## Plotting the HCR heatmap:

order <- order(classif_label)
anno.df <- as.data.frame(cbind(ilastik.anno[,1:16], classif_label, rownames(ilastik.anno)))
colnames(anno.df)[17:18] <- c('type', 'cell.id')
anno.melt <- melt(anno.df, value.name = 'value', variable.name = 'gene', id.vars = c('cell.id', 'type'))
anno.melt$type <- factor(anno.melt$type, levels = levels(classif_label))
anno.melt <- anno.melt %>% arrange(factor(anno.melt$type, levels = levels(anno.melt$type)))
anno.melt$cell.id <- factor(anno.melt$cell.id, levels = unique(anno.melt$cell.id))
anno.melt$gene <- factor(anno.melt$gene, levels = rev(colnames(ilastik.anno)[1:16]))

ggplot(anno.melt, aes(x=cell.id, y=gene, group=type)) +
  geom_tile(mapping=aes(fill=value)) + scale_fill_gradient(low='white', high = '#00368b') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# The scRNA-Seq heatmap, with number of cells matched to the HCR data

sc_table <- cbind(meta_bin, mop_object_hcr$my_supertype)
colnames(sc_table)[ncol(sc_table)] <- 'my_supertype'
sc_matched <- sc_table[sc_table$my_supertype %in% unique(classif_label),]
to_compare <- table(classif_label[classif_label %in% unique(sc_matched$my_supertype)])

sc_compile <- matrix(nrow = sum(to_compare), ncol = 17)
count <- 1
for (i in names(to_compare)[names(to_compare) %in% unique(sc_matched$my_supertype)]) {
  which_cells_match <- which(sc_matched$my_supertype == i)
  cells_match <- as.matrix(sc_matched[sample(which_cells_match, to_compare[i], replace = T),])
  sc_compile[count:(count+to_compare[i]-1),1:17] <- cells_match[,colnames(cells_match)[1:17]]
  count <- count+as.numeric(to_compare[i])
}

colnames(sc_compile) <- colnames(sc_matched)[1:17]
sc_compile <- as.data.frame(sc_compile)
sc_compile$my_supertype <- factor(sc_compile$my_supertype, levels = levels(classif_label)[levels(classif_label) %in% unique(sc_compile$my_supertype)])
sc_compile[,'cell.id'] <- rownames(sc_compile)
sc_melt <- melt(sc_compile, variable.name = 'gene', value.name = 'value', id.vars = c('my_supertype', 'cell.id'), na.rm = T)
sc_melt <- sc_melt %>% arrange(factor(sc_melt$my_supertype, levels = levels(sc_melt$my_supertype)))
sc_melt$cell.id <- factor(sc_melt$cell.id, levels = unique(sc_melt$cell.id))
sc_melt$gene <- factor(sc_melt$gene, levels = rev(colnames(ilastik.anno)[1:16]))
sc_melt$value <- as.numeric(sc_melt$value)

ggplot(sc_melt, aes(x=cell.id, y=gene, group=my_supertype)) +
  geom_tile(mapping=aes(fill=value)) + scale_fill_gradient(low='white', high = '#00368b') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

