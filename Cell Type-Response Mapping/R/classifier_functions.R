# classifier_functions.R
# Core classification functions for the MOSAIX pipeline.
# All functions are self-contained — no global variable dependencies.

library(dplyr)
library(reshape2)

#' Primary classifier: Hierarchical GMM-based Bayesian classifier.
#' Walks the Allen taxonomy using GMM-binarized expression tables with optional
#' laminar depth smoothing from MERFISH data.
#'
#' @param hcr.table Data frame of binarized HCR expression with newy column
#' @param bin_table Binarized meta-cell gene expression table with taxonomy columns
#' @param merfish_data MERFISH depth metadata
#' @param mop_object_hcr Reference expression matrix (defaults to bin_table with all_cells column)
#' @param use.lamina Whether to use laminar depth priors (default TRUE)
#' @param depth.round Decimal places for depth binning (default 2)
#' @return List with classif.label and classif.prob vectors
HierarchicalGMMBayes <- function(hcr.table, bin_table, merfish_data = NULL, mop_object_hcr = NULL, use.lamina = TRUE, depth.round = 2) {

  if (is.null(mop_object_hcr)) {
    mop_object_hcr <- bin_table
    mop_object_hcr[,'all_cells'] <- rep('all', nrow(mop_object_hcr))
  }

  all.labels <- character(nrow(hcr.table))
  label.prob <- numeric(nrow(hcr.table))
  
  print('Performing cell-by-cell classification...')
  pb <- txtProgressBar(min = 0, max = nrow(hcr.table), initial = 0, style = 3)
  
  hcr_depths <- round(hcr.table[,'newy'], depth.round)
  classes <- c('allen_NT', 'allen_class', 'allen_subclass', 'my_supertype')
  used_mer_depths <- list()
  used_percents <- list()
  
  for (i in 1:nrow(hcr.table)) {
    
    cell <- hcr.table[i,1:ncol(hcr.table)-1]
    cell.depth <- hcr_depths[i]
    cell.pos <- subset(cell, select=(cell == 1 | cell == 0.5))
    cell.neg <- subset(cell, select=(cell == 0))
    
    if (ncol(cell.pos) == 0) {
      depth_of_neg <- hcr.table[i,'newy']
      
      if (depth_of_neg > -0.075) {label <- '9991 ZNo markers L1'}
      else if (depth_of_neg <= -0.075 & depth_of_neg > -0.22) {label <- '9992 ZNo markers L2/3'}
      else if (depth_of_neg <= -0.22 & depth_of_neg > -0.4) {label <- '9993 ZNo markers L4/5a'}
      else if (depth_of_neg <= -0.4 & depth_of_neg > -0.55) {label <- '9994 ZNo markers UL5b'}
      else if (depth_of_neg <= -0.55 & depth_of_neg > -0.7) {label <- '9995 ZNo markers LL5b'}
      else if (depth_of_neg <= -0.7 & depth_of_neg > -0.9) {label <- '9996 ZNo markers L6a'}
      else if (depth_of_neg <= -0.9) {label <- '9997 ZNo markers L6b'}
      else {label <- '9998 ZNo markers'}
      
      prob.closest <- 1
    }
    else if (ncol(cell.pos) == 1 & 'Slc32a1' %in% colnames(cell.pos)) {
      label <- '9990 VVgat'
      prob.closest <- 1
    }
    else {
      
      label <- 'all'
      
      for (j in 1:4) {
        
        if (j == 1) {
          next_class <- classes[j]
          
          if (label %in% names(used_percents)){
            each.percent <- used_percents[[label]]
          }
          else{
            
            each.percent <- data.frame(matrix(nrow = length(unique(bin_table[,next_class])), ncol = ncol(bin_table)-5))
            rownames(each.percent) <- unique(bin_table[,next_class])
            colnames(each.percent) <- colnames(bin_table)[1:(ncol(bin_table)-5)]
            
            input_matrix <- as.matrix(bin_table[bin_table[,next_class] %in% rownames(each.percent),])
            for (k in 1:nrow(each.percent)) {
              
              input <- as.matrix(bin_table[bin_table[,next_class] == rownames(each.percent)[k], 1:(ncol(bin_table)-5)])
              percents <- apply(input, 2, function(col){(sum(col)+1)/length(col)})
              each.percent[k,] <- percents
              
            }
            
            each.percent <- cbind(rownames(each.percent), each.percent)
            names(each.percent)[1] <- 'id'
            used_percents[[label]] <- each.percent
          }
        }
        else {
          prev_class <- next_class
          next_class <- classes[j]
          
          if (label %in% names(used_percents)){
            each.percent <- used_percents[[label]]
          }
          else{
            
            labels_of_prev <- unique(bin_table[bin_table[,prev_class] == label, next_class])
            each.percent <- data.frame(matrix(nrow = length(labels_of_prev), ncol = ncol(bin_table)-5))
            rownames(each.percent) <- labels_of_prev
            colnames(each.percent) <- colnames(bin_table)[1:(ncol(bin_table)-5)]
            
            input_matrix <- as.matrix(bin_table[bin_table[,next_class] %in% rownames(each.percent),])
            gene_prior <- character(length = ncol(each.percent))
            for (k in 1:nrow(each.percent)) {
              
              input <- as.matrix(bin_table[bin_table[,next_class] == rownames(each.percent)[k], 1:(ncol(bin_table)-5)])
              percents <- apply(input, 2, function(col){(sum(col)+1)/length(col)})
              each.percent[k,] <- percents
              
            }
            
            each.percent <- cbind(rownames(each.percent), each.percent)
            names(each.percent)[1] <- 'id'
            used_percents[[label]] <- each.percent
          }
        }
        
        # Merfish depth extraction:
        if (label %in% names(used_mer_depths)){
          depth_data <- used_mer_depths[[label]]
        }
        else{
          depth_data <- ExtractMerfishLamina(merfish_data, next_class, mop_object_hcr, bin = T, bin_to = depth.round, smooth = T)
          used_mer_depths[[label]] <- depth_data
        }
        
        type_densities <- depth_data[[1]]
        rel_depths <- depth_data[[2]]
        prop_depths <- depth_data[[3]]
        merfish_type_props <- table(rel_depths$use_type)/nrow(rel_depths)
        
        prob.closest <- 0
        label <- '9999 ZNo Good Match'
        each.percent <- each.percent[(rownames(each.percent) %in% names(merfish_type_props) & rownames(each.percent) %in% colnames(type_densities)),]
        
        for (k in 1:nrow(each.percent)) {
          
          # Bernoulli NB:
          test.type <- each.percent[k,'id']
          prop.depth <- prop_depths$prop[prop_depths$use_type == as.character(test.type) & prop_depths$Y == cell.depth]
          
          if (use.lamina & j != 1){
            if(length(prop.depth)) {
              type.depth.prob <- prop.depth
            }
            else{
              type.depth.prob <- min(prop_depths$prop)
            }
          }
          else {
            type.depth.prob <- 1
          }
          # print(type.depth.prob)
          # if (j == 1){
          #   print((each.percent[k,colnames(cell.pos)]))
          #   print(1-(as.numeric(each.percent[k,colnames(cell.neg)])))
          # }
          #fnr <- 0.05
          pos.genes <- prod(as.numeric(each.percent[k,colnames(cell.pos)]))
          neg.genes <- prod(1-(as.numeric(each.percent[k,colnames(cell.neg)])))
          #neg.genes <- neg.genes*(1-fnr) + pos.genes*fnr
          
          prob <- type.depth.prob*pos.genes*neg.genes
          
          if (prob > prob.closest) {
            prob.closest <- prob
            label <- as.character(test.type)
          }
          
        }
        
        if (prob.closest == 0){
          break
        }
      }
      
    }
    if (prob.closest == 0) {
      label <- '9999 ZNo Good Match'
    }
    all.labels[i] <- label
    label.prob[i] <- prob.closest
    setTxtProgressBar(pb, i)
    
  }
  close(pb)
  classifs <- list('classif.label' = all.labels, 'classif.prob' = label.prob)
  if(length(unique(classifs$classif.label)) == 1){
    print('ONLY ONE CELL TYPE DETECTED! You may have forgotten to turn the Ilastik labels into booleans')
  }
  return(classifs)
}


#' Calculate percent of cells above threshold for each gene.
#'
#' @param scale.object Gene x cell matrix
#' @param threshold Named vector of thresholds per gene, or a single quantile value
#' @param quant If TRUE, use quantile-based thresholding instead of absolute
#' @return Named numeric vector of percentages per gene
PercentAboveVM <- function(scale.object, threshold, quant = FALSE) {

  if (!quant) {

    percents <- vector('numeric', length = nrow(scale.object))

    for (i in 1:nrow(scale.object)) {

      thresh <- threshold[,rownames(scale.object)[i]]
      percent <- sum(scale.object[i,] > thresh)/ncol(scale.object)
      percents[i] <- percent

    }

    names(percents) <- rownames(scale.object)
  }

  else {
    percents <- vector('numeric', length = nrow(scale.object))

    for (i in 1:nrow(scale.object)) {
      qVal <- quantile(scale.object[i,scale.object[i,] > 0], threshold)
      percent <- length(scale.object[i,scale.object[i,] > qVal])/length(scale.object[i,])
      percents[i] <- percent
    }
    names(percents) <- rownames(scale.object)

  }
  return(percents)
}


#' Extract MERFISH laminar depth distributions per cell type.
#' Bins and optionally smooths depth data for use as laminar priors.
#'
#' @param mer_metadata MERFISH metadata data frame with newy and type columns
#' @param selected_class Column name in metadata to group by (e.g., 'my_supertype')
#' @param mop_object_hcr Reference expression matrix to filter types against
#' @param quant_low Lower quantile for depth normalization (default 0.02)
#' @param quant_high Upper quantile for depth normalization (default 0.91)
#' @param bin Whether to bin depths (default FALSE)
#' @param bin_to Decimal places for binning (default 2)
#' @param smooth Whether to apply loess smoothing to depth distributions (default TRUE)
#' @return List of (type_densities matrix, rel.depths data frame, prop_depths data frame)
ExtractMerfishLamina <- function(mer_metadata, selected_class, mop_object_hcr, quant_low = 0.02, quant_high = 0.91, bin = FALSE, bin_to = 2, smooth = TRUE) {

  mop_merdata <- mer_metadata

  # Redo min-max based on specified quants:
  new.y <- (mop_merdata$newy - quantile(mop_merdata$newy, quant_low))/(quantile(mop_merdata$newy, quant_high) - quantile(mop_merdata$newy, quant_low))
  mop_merdata$new.quant.y <- new.y

  rel.depths <- data.frame(X = mop_merdata$newx, Y = -mop_merdata$new.quant.y)
  rel.depths[,'use_type'] <- mop_merdata[[selected_class]]
  rel.depths <- rel.depths[order(rel.depths$use_type),]
  rel.depths <- rel.depths[rel.depths$use_type %in% mop_object_hcr[,selected_class],]
  rel.depths <- rel.depths[rel.depths$Y >= -1.0,]
  rel.depths <- rel.depths[!(grepl('L6 IT', rel.depths$use_type) & rel.depths$Y > -0.5),]

  ylims <- rel.depths %>%
    group_by(use_type) %>%
    summarise(Q1 = quantile(Y, 0.01), Q9 = quantile(Y, 0.99))

  rel.depths.mod <- data.frame()
  for (i in 1:nrow(ylims)) {
    max.Y <- as.numeric(ylims[i,'Q9'])
    min.Y <- as.numeric(ylims[i,'Q1'])
    rel.depths.rows <- rel.depths$use_type == as.character(ylims[i,'use_type'])
    rel.depths.inquant <- rel.depths[(rel.depths.rows & (rel.depths$Y >= min.Y) & (rel.depths$Y < max.Y)),]
    rel.depths.mod <- rbind(rel.depths.mod, rel.depths.inquant)
  }
  rel.depths.mod[,'simple_type'] <- sub('(.*)_\\d+', '\\1', rel.depths.mod$use_type)

  if (bin) {
    rel.depths.mod$Y <- round(rel.depths.mod$Y, bin_to) # Round the depths to "bin" them to 100 values
  }
  
  # Compute proportional depth densities per type
  densities <- matrix(nrow = 101, ncol = length(unique(rel.depths.mod$use_type)), dimnames = list(c(0:100*-0.01), unique(rel.depths.mod$use_type)))
  for (i in 1:ncol(densities)) {
    
    type <- unique(rel.depths.mod$use_type)[i]
    type_Ys <- rel.depths.mod[rel.depths.mod$use_type == type, 'Y']
    for (j in rownames(densities)) {
      
      depth_prop <- sum(type_Ys == as.numeric(j))/length(type_Ys)
      densities[j,i] <- depth_prop
      
    }
  }
  n_depth <- rel.depths.mod %>% 
    group_by(use_type, Y) %>% 
    mutate(n_type = n()) %>% 
    ungroup() %>%
    mutate(n_depth = n(), .by = Y) %>%
    mutate(prop = n_type/n_depth)
  
  prop_depth <- as.data.frame(unique(n_depth[,c('Y', 'use_type', 'prop')]))
  
  if (smooth) {
    prop_cast <- dcast(data = prop_depth, formula = Y ~ use_type, fun.aggregate = sum, value.var = 'prop')
    prop_cast_adj <- prop_cast + 0.002
    prop_melt <- melt(prop_cast_adj, id.vars = 'Y', variable.name = 'use_type', value.name = 'prop')
    
    smooth_prop <- prop_melt %>% group_by(use_type) %>% arrange(Y) %>%
      mutate(smooth_prop = loess(prop ~ Y, span = 0.16)$fitted) %>% ungroup() %>%
      group_by(Y) %>% mutate(smooth_prop = pmax(smooth_prop, 0.002),
                             smooth_prop = smooth_prop/sum(smooth_prop))
    
    prop_depth <- data.frame(Y = smooth_prop$Y,
                             use_type = smooth_prop$use_type,
                             prop = smooth_prop$smooth_prop)
  }
  
  return(list(densities, rel.depths.mod, prop_depth))
}


#' Combine co-expressed gene pairs by taking the max expression value.
#'
#' @param table Data frame containing gene columns
#' @param genes_combined List of character vectors, each of length 2
#' @return Modified data frame with combined gene columns
CombineGenes <- function(table, genes_combined) {
  if (length(genes_combined) == 0) return(table)
  for (i in 1:length(genes_combined)) {

    gene.1 <- genes_combined[[i]][1]
    gene.2 <- genes_combined[[i]][2]
    table[,gene.1] <- pmax(table[,gene.1], table[,gene.2])
    table <- table[,colnames(table) != gene.2]
    colnames(table)[colnames(table) == gene.1] <- paste(gene.1, gene.2, sep='+')

  }
  return(table)
}


#' Dot product-based cell type classification.
#' Classifies cells by computing dot product similarity against average
#' expression profiles per cell type.
#'
#' @param mop_object Seurat object with gene expression data
#' @param hcr.table Data frame of binarized HCR expression with newy column
#' @param merfish_data MERFISH depth metadata
#' @param meta_gene_matrix Reference log-normalized gene expression matrix
#' @param genes_combined List of gene pairs to combine
#' @return List with classif.label and classif.prob vectors
DotProductABC <- function(mop_object, hcr.table, merfish_data, meta_gene_matrix, genes_combined) {

  mop_object[['my_supertype']] <- mop_object$allen_supertype
  mop_object$my_supertype[!grepl('L\\d', mop_object$my_supertype)] <- mop_object$allen_subclass[!grepl('L\\d', mop_object$my_supertype)]
  gene_average <- AverageExpression(mop_object, features = colnames(meta_gene_matrix)[1:17],
                                    group.by = 'my_supertype', layer = 'counts')
  mm_gene_average <- apply(as.data.frame(gene_average), 1, function(x){(x-min(x))/(max(x)-min(x))})
  comb_gene_average <- CombineGenes(mm_gene_average, genes_combined)
  type_names <- sub('g(\\d{3,4} .*(Glut-\\d|Gaba+))', '\\1', colnames(gene_average$RNA))
  type_names <- sub('-', '_', type_names)
  rownames(comb_gene_average) <- type_names

  comb_gene_average <- comb_gene_average[rownames(comb_gene_average) %in% merfish_data$my_supertype,]

  all.labels <- c()
  all.prods <- c()
  for(i in 1:nrow(hcr.table)){
    dot.prods <- apply(comb_gene_average, 1, function(x){sum((hcr.table[i,1:(ncol(hcr.table)-1)]+0.0001)*x)})
    winning.type <- names(dot.prods)[dot.prods == max(dot.prods)]
    all.labels <- c(all.labels, winning.type)
    all.prods <- c(all.prods, max(dot.prods))
  }

  classifs <- list('classif.label' = all.labels, 'classif.prob' = all.prods)
  return(classifs)
}
