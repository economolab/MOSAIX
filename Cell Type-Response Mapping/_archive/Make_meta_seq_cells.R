library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)

#setwd('~/Desktop/Transcriptomic data/')
#mop_object <- readRDS('~/Desktop/Transcriptomic data/CTX/mop_object.rds')
mop_object <- readRDS('~/Desktop/Local code/HCR Annotation/sub_mop_neurons.rds')
merfish_data <- read.table('~/Desktop/Local code/HCR Annotation/merfish_depth_metadata.csv', sep = ',', header = T)
mop_object[['my_supertype']] <- mop_object$allen_supertype
mop_object$my_supertype[!grepl('L\\d', mop_object$my_supertype)] <- mop_object$allen_subclass[!grepl('L\\d', mop_object$my_supertype)]

mop_object <- FindNeighbors(mop_object, reduction = 'pca', dims = 1:30, return.neighbor = T, k = 500)

neighbor_dists <- mop_object@neighbors$RNA.nn@nn.dist[,2:500]
rownames(neighbor_dists) <- colnames(mop_object)# rows are cells, cols are neighbor dists for the 19 nearest - first column is empty?
neighbor_idx <- mop_object@neighbors$RNA.nn@nn.idx[,2:500] # rows are cells, cols are column indices of neighbors in the object structure
rownames(neighbor_idx) <- colnames(mop_object)

# only include clusters with at least 80 cells in it

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Slco2a1', 'Fam84b', 'Npnt', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')

type_counts <- table(mop_object$my_supertype)
clusters_to_include <- names(type_counts[(type_counts > 80) & (names(type_counts) %in% merfish_data$my_supertype)])
n_clusters <- length(clusters_to_include)

#mop_object@assays$RNA@meta.features$gene_symbol[mop_object@assays$RNA@meta.features$gene_symbol == 'Lratd2'] <- 'Fam84b'
#mop_object@assays$RNA@meta.features$gene_symbol[mop_object@assays$RNA@meta.features$gene_symbol == 'Ccn2'] <- 'Ctgf'
#mop_object@assays$RNA@counts@Dimnames[[1]] <- mop_object@assays$RNA@meta.features[rownames(mop_object),1]
#mop_object@assays$RNA@data@Dimnames[[1]] <- mop_object@assays$RNA@meta.features[rownames(mop_object),1]

meta_norm_counts <- matrix(nrow = length(hcr_genes), ncol = sum(type_counts[clusters_to_include])) #50

count <- 0
types_kept <- c(rep(NULL, sum(type_counts[clusters_to_include])))
for (i in clusters_to_include) {

  sub_counts <- mop_object@assays$RNA@data[,mop_object$my_supertype == i]
  sub_counts <- expm1(sub_counts)
  sub_neighbors <- neighbor_idx[mop_object$my_supertype == i,]
  
  set.seed(10)
  sample_meta_cell <- sample(rownames(sub_neighbors), type_counts[i], replace = F)
  for (j in 1:type_counts[i]) {
    count <- count + 1
    neighbors <- neighbor_idx[sample_meta_cell[j],]
    neighbor_names <- colnames(mop_object)[neighbors]
    neighbor_types <- mop_object$my_supertype[neighbors]
    neighbors_to_keep <- neighbor_names[neighbor_names %in% colnames(sub_counts)]
    
    if (length(neighbors_to_keep) < 9) {
      print(i)
      #print(neighbor_types)
      #print(length(neighbors_to_keep))
      print('Shit')
    }
    else{
      mean_counts <- apply(as.matrix(sub_counts[hcr_genes, neighbors_to_keep[1:9]]), 1, sum)
      mean_counts <- mean_counts + sub_counts[hcr_genes, sample_meta_cell[j]]
      types_kept[count] <- i

      meta_norm_counts[,count] <- mean_counts
    }
    
  }
  
}

meta_norm_counts <- t(meta_norm_counts)
meta_norm_counts <- as.data.frame(meta_norm_counts)
colnames(meta_norm_counts)[1:17] <- hcr_genes
meta_norm_counts[,'my_supertype'] <- types_kept
meta_norm_counts <- meta_norm_counts[!is.na(meta_norm_counts$my_supertype),]

##### If file already exists: #####
meta_norm_counts <- read.csv("meta_norm_counts.csv", sep = ',', header = T)
all_types <- unique(meta_norm_counts$my_supertype)
#cluster_to_supertype <- vector(mode = 'character', length = length(all_types))
cluster_to_subclass <- vector(mode = 'character', length = length(all_types))
cluster_to_class <- vector(mode = 'character', length = length(all_types))
cluster_to_nt <- vector(mode = 'character', length = length(all_types))
for (i in 1:length(all_types)){
  
  cluster_type <- all_types[i]
  idxs <- mop_object$my_supertype == cluster_type
  #supertypes <- unique(mop_object$allen_supertype[idxs])
  subclasses <- unique(mop_object$allen_subclass[idxs])
  classes <- unique(mop_object$allen_class[idxs])
  nt <- unique(mop_object$allen_NT[idxs])
  if (length(nt) > 1) {
    print(nt)
    print("Wtf")
  }
  cluster_to_subclass[i] <- subclasses
  cluster_to_class[i] <- classes
  #cluster_to_supertype[i] <- supertypes
  if (length(nt) > 1) {
    cluster_to_nt[i] <- nt[1]
  }
  else{
    cluster_to_nt[i] <- nt
  }
    
}
names(cluster_to_subclass) <- all_types
names(cluster_to_class) <- all_types
names(cluster_to_nt) <- all_types
#names(cluster_to_supertype) <- all_types

#meta_norm_counts[,'allen_supertype'] <- cluster_to_supertype[meta_norm_counts$allen_cluster]
meta_norm_counts[,'allen_subclass'] <- cluster_to_subclass[meta_norm_counts$my_supertype]
meta_norm_counts[,'allen_class'] <- cluster_to_class[meta_norm_counts$my_supertype]
meta_norm_counts[,'allen_NT'] <- cluster_to_nt[meta_norm_counts$my_supertype]
meta_norm_counts <- meta_norm_counts[!meta_norm_counts$allen_NT == 'Glut-GABA',]
meta_norm_counts <- meta_norm_counts[!meta_norm_counts$allen_class == '08 CNU-MGE GABA',]

meta_log_counts <- meta_norm_counts
meta_log_counts[,1:17] <- log(meta_log_counts[,1:17]) # This is the natural log without adding 1 (to keep things)

hist(meta_log_counts[,'Pvalb'], 200, ylim = c(0,700), xlim = c(-2, 10))
hist(mop_object@assays$RNA@data['Fam84b',], 200, ylim = c(0,500), xlim = c(0.1, 5))

sub_log_counts <- meta_log_counts[sample(1:nrow(meta_log_counts), 5000),]
hist(sub_log_counts[,'Slco2a1'], 100, ylim = c(0,100), xlim = c(-2, 8))
hist(meta_norm_counts[,'Fam84b'], 200, ylim = c(0,2000), xlim = c(0, 200))

to_save <- meta_log_counts$Vip
to_save[to_save == -Inf] <- -5
# Save each gene to load into python as a np array
write.table(to_save, 
            '~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/Vip.csv', 
            sep = ',', col.names = F, row.names = F)
#write.table(meta_norm_counts, '~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/meta_norm_counts.csv', sep = ',', col.names = T, row.names = F)

##### After running the GMM binarization in python: #####
meta_bin <- read.table('~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/All_genes_bin.csv', sep = ',', header = T)
meta_bin <- as.data.frame(meta_bin)
meta_bin[,18:ncol(meta_log_counts)] <- meta_log_counts[,18:ncol(meta_log_counts)]


### Plot patterns: ###
meta.dot.df <- meta_bin[,1:18]
for (i in 1:17){
  if (2 %in% meta.dot.df[,i]){
    if (colnames(meta.dot.df)[i] %in% c('Rorb', 'Fam84b', 'Ccdc80', 'Pvalb', 'Vip', 'Npnt')){
      meta.dot.df[,i] <- meta.dot.df[,i] == max(meta.dot.df[,i])
    }
    else{
      meta.dot.df[,i] <- meta.dot.df[,i] == max(meta.dot.df[,i])-1
      print(colnames(meta.dot.df)[i])
    }
  }
  else{
    meta.dot.df[,i] <- meta.dot.df[,i] > 0
  }
}
meta.dot.df[,1:17] <- meta.dot.df[,1:17] * 1

for (i in 1:length(genes_combined)) {
  
  gene.1 <- genes_combined[[i]][1]
  gene.2 <- genes_combined[[i]][2]
  meta.dot.df[,gene.1] <- pmax(meta.dot.df[,gene.1], meta.dot.df[,gene.2])
  meta.dot.df <- meta.dot.df[,colnames(meta.dot.df) != gene.2]
  colnames(meta.dot.df)[colnames(meta.dot.df) == gene.1] <- paste(gene.1, gene.2, sep='+')
  
}

meta.dot.tidy <- meta.dot.df %>% group_by(my_supertype) %>% summarise(across(colnames(meta.dot.df)[1:16], sum))
meta.dot.n <- meta.dot.df %>% group_by(my_supertype) %>% count()
meta.dot.prop <- mapply(`/`, data.frame(meta.dot.tidy[,2:ncol(meta.dot.tidy)]), meta.dot.n[,2])

meta.dot.prop <- data.frame(meta.dot.prop[,colnames(meta.dot.df[,1:16])],
                             id = meta.dot.n[,1])
meta.dot.prop <- meta.dot.prop[meta.dot.prop$my_supertype %in% merfish.rel.depths$use_type,]
meta.dot.melt <- melt(meta.dot.prop)
meta.dot.melt$my_supertype <- factor(meta.dot.melt$my_supertype, levels = names(cluster_colors))

ggplot(meta.dot.melt, aes(x = my_supertype, y = variable, size = value, color = my_supertype)) + 
  geom_point(stroke = 0) + #color = '#0000ff'
  theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, hjust =1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = cluster_colors) +
  scale_size_continuous(range = c(0, 8)) +
  guides(color = 'none')


#### Plot gaussians: #####
gene_thresholds_table <- read.table('~/Desktop/Local code/HCR Annotation/Classification with Gaussian mixture/GMM_Thresholds.csv', sep = ',', header = T)
rownames(gene_thresholds_table) <- gene_thresholds_table$X
gene <- 'Slc32a1'

p3 <- ggplot(log(expm1(data.frame(Slc32a1=mop_object@assays$RNA@data[gene,]))), aes(Slc32a1)) +
  theme_minimal() +
  geom_histogram(binwidth = 0.05) +
  ylim(c(0,500)) +
  xlim(c(-1.5,5)) #5.5
p1 <- ggplot(meta_log_counts, aes(Slc32a1)) +
  theme_minimal() +
  geom_histogram(binwidth = 0.05) +
  ylim(c(0,1100)) +
  xlim(c(-1.5,5))
p2 <- ggplot(meta_log_counts, aes(Slc32a1)) +
  theme_minimal() +
  #geom_histogram(binwidth = 0.1) +
  #ylim(c(0,20)) +
  stat_function(fun = function(x){dnorm(x, mean = gene_thresholds_table[gene, 'Mean.1'], 
                                        sd = gene_thresholds_table[gene, 'Stddev.1'])*gene_thresholds_table[gene, 'Weight.1']}) +
  stat_function(fun = function(x){dnorm(x, mean = gene_thresholds_table[gene, 'Mean.2'], 
                                        sd = gene_thresholds_table[gene, 'Stddev.2'])*gene_thresholds_table[gene, 'Weight.2']}) +
  stat_function(fun = function(x){dnorm(x, mean = gene_thresholds_table[gene, 'Mean.3'], 
                                        sd = gene_thresholds_table[gene, 'Stddev.3'])*gene_thresholds_table[gene, 'Weight.3']}) +
  stat_function(fun = function(x){dnorm(x, mean = gene_thresholds_table[gene, 'Mean.4'], 
                                        sd = gene_thresholds_table[gene, 'Stddev.4'])*gene_thresholds_table[gene, 'Weight.4']}) +
  #scale_y_continuous(breaks = NULL) +
  ylim(c(0, 0.75)) +
  xlim(c(-1.5,5)) +
  geom_vline(xintercept = c(gene_thresholds_table[gene, 'Threshold.1'], 
                            gene_thresholds_table[gene, 'Threshold.2'],
                            gene_thresholds_table[gene, 'Threshold.3']))
p3/p1/p2




