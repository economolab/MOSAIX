library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(cowplot)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(Polychrome)
library(reshape2)
library(clipr)
library(scatterplot3d)
library(ggraph)
library(igraph)
library(scales)

BayesClassifier <- function(hcr.table, bin_table, all.probs, lamina = NULL, mer_probs = NULL) {
  
  all.labels <- character(nrow(hcr.table))
  label.prob <- numeric(nrow(hcr.table))
  
  print('Performing cell-by-cell classification...')
  pb <- txtProgressBar(min = 0, max = nrow(hcr.table), initial = 0, style = 3)
  
  if (is.null(lamina)) {
    for (i in 1:nrow(hcr.table)) {
      
      cell <- hcr.table[i,1:ncol(hcr.table)-1]
      cell.pos <- subset(cell, select=(cell == 1))
      cell.neg <- subset(cell, select=(cell == 0 | cell == 0.5))
      if (ncol(cell.pos) == 0) {
        
        label <- 'ZNo markers'
        prob.closest <- 1
        
      }
      else {
        
        prob.closest <- 0
        label <- 0
        prop.table <- meta_bin %>% group_by(my_supertype) %>% summarise(across(colnames(meta_bin)[1:16], sum))
        hcr.table[,'type.prob'] <- rep(1, nrow(hcr.table))

        for (j in 1:nrow(hcr.table)) {
          
          #if (sc.table[j, 'id'] %in% excit.types) {
          #  type.prob <- 0.78
          #}
          #else {
          #  type.prob <- 0.14
          #}
          type.prob <- hcr.table[j,'type.prob']
          pos.numerator <- prod((as.numeric(hcr.table[j,colnames(cell.pos)])+0.00001))
          neg.numerator <- prod(1-((as.numeric(hcr.table[j,colnames(cell.neg)])+0.00001)))#(cell.type.probs[as.character(sc.table[j,'id'])])
          pos.denominator <- prod((as.numeric(all.probs[colnames(cell.pos)])+0.00001))
          neg.denominator <- prod(1-((as.numeric(all.probs[colnames(cell.neg)])+0.00001)))
          prob <- (pos.numerator*neg.numerator*type.prob)/(pos.denominator*neg.denominator)
          
          if (prob > prob.closest) {
            prob.closest <- prob
            label <- as.character(sc.table[j,'id'])
            
          }
          
          else {
            
            next
            
          }
        }
      }
      
      all.labels[i] <- label
      label.prob[i] <- prob.closest
      setTxtProgressBar(pb, i)  
    }
  }
  else {
    
    hcr_depths <- round(hcr.table[,'newy'], 2)
    
    for (i in 1:nrow(hcr.table)) {
      
      cell <- hcr.table[i,1:ncol(hcr.table)-1]
      cell.depth <- hcr_depths[i]
      cell.pos <- subset(cell, select=(cell == 1 | cell == 0.5))
      cell.neg <- subset(cell, select=(cell == 0))
      if (ncol(cell.pos) == 0) {
        
        label <- 'ZNo markers'
        prob.closest <- 0
        
      }
      else {
        
        prob.closest <- 0
        label <- 0
        sc.table[,'type.prob'] <- mer_probs[rownames(sc.table)]
        
        for (j in 1:nrow(sc.table)) {
          
          # Bernoulli NB:
          type.prob <- sc.table[j,'type.probs']
          type.depth.prob <- lamina[as.character(cell.depth),sc.table[j,'id']]
          pos.genes <- prod((as.numeric(sc.table[j,colnames(cell.pos)])))
          neg.genes <- prod(1-((as.numeric(sc.table[j,colnames(cell.neg)]))))
          prob <- type.prob*type.depth.prob*pos.genes*neg.genes
          #print(prob)
          
          if (prob > prob.closest) {
            prob.closest <- prob
            label <- as.character(sc.table[j,'id'])
            
          }
          
          else {
            
            next
            
          }
        }
      }
      
      all.labels[i] <- label
      label.prob[i] <- prob.closest
      setTxtProgressBar(pb, i)  
    }
    
  }
  close(pb)
  classifs <- list('classif.label' = all.labels, 'classif.prob' = label.prob)
  if(length(unique(classifs$classif.label)) == 1){
    print('ONLY ONE CELL TYPE DETECTED! You may have forgotten to turn the Ilastik labels into booleans')
  }
  return(classifs)
}

HierarchicalBayes <- function(hcr.table, thresholds, merfish_data = NULL) {
  
  all.labels <- character(nrow(hcr.table))
  label.prob <- numeric(nrow(hcr.table))
  
  print('Performing cell-by-cell classification...')
  pb <- txtProgressBar(min = 0, max = nrow(hcr.table), initial = 0, style = 3)
  
  hcr_depths <- round(hcr.table[,'newy'], 2)
  classes <- c('allen_NT', 'allen_class', 'allen_subclass', 'allen_supertype')
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
      else if (depth_of_neg <= -0.075 & depth_of_neg > -0.25) {label <- '9992 ZNo markers L2/3'}
      else if (depth_of_neg <= -0.25 & depth_of_neg > -0.4) {label <- '9993 ZNo markers L4/5a'}
      else if (depth_of_neg <= -0.4 & depth_of_neg > -0.55) {label <- '9994 ZNo markers UL5b'}
      else if (depth_of_neg <= -0.55 & depth_of_neg > -0.7) {label <- '9995 ZNo markers LL5b'}
      else if (depth_of_neg <= -0.7 & depth_of_neg > -0.9) {label <- '9996 ZNo markers L6a'}
      else if (depth_of_neg <= -0.9) {label <- '9997 ZNo markers L6b'}
      else {label <- '9998 ZNo markers'}
      
      prob.closest <- 1
    }
    else {
      
      label <- 'all'
      
      for (j in 1:4) {
        
        if (grepl('Pvalb|Sst|Lamp5|Vip|Sncg|Error', label)){
          break
        }
        if (j == 1) {
          next_class <- classes[j]
          sub_thresholds <- thresholds[label, 2:ncol(thresholds)]
          if (label %in% names(used_percents)){
            each.percent <- used_percents[[label]]
          }
          else{
            each.percent <- data.frame(matrix(nrow = length(unique(mop_object_hcr[,next_class])), ncol = length(hcr_genes)-length(genes_combined)))
            rownames(each.percent) <- unique(mop_object_hcr[,next_class])
            colnames(each.percent) <- colnames(sub_thresholds)
            
            # gene.1 <- genes_combined[[1]][1]
            # gene.2 <- genes_combined[[1]][2]
            for (k in 1:nrow(each.percent)) {
              
              input <- as.matrix(mop_object_hcr[mop_object_hcr[,next_class] == rownames(each.percent)[k], 1:16])
              input <- t(input)
              # input[gene.1,] <- pmax(input[gene.1,],input[gene.2,])
              # rownames(input)[rownames(input) == gene.1] <- paste(gene.1, gene.2, sep='+')
              # input <- input[!(rownames(input) == gene.2),]
              percents <- PercentAboveVM(input, threshold=sub_thresholds, quant = F) # 0.09 when quantiles are not used
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
          sub_thresholds <- thresholds[label, 2:ncol(thresholds)]
          if (label %in% names(used_percents)){
            each.percent <- used_percents[[label]]
          }
          else{
            labels_of_prev <- unique(mop_object_hcr[mop_object_hcr[,prev_class] == label, next_class])
            each.percent <- data.frame(matrix(nrow = length(labels_of_prev), ncol = length(hcr_genes)-length(genes_combined)))
            rownames(each.percent) <- labels_of_prev
            colnames(each.percent) <- colnames(sub_thresholds)
            
            #gene.1 <- genes_combined[[1]][1]
            #gene.2 <- genes_combined[[1]][2]
            for (k in 1:nrow(each.percent)) {
              
              input <- as.matrix(mop_object_hcr[mop_object_hcr[,next_class] == rownames(each.percent)[k], 1:16])
              input <- t(input)
              # input[gene.1,] <- pmax(input[gene.1,],input[gene.2,])
              # rownames(input)[rownames(input) == gene.1] <- paste(gene.1, gene.2, sep='+')
              # input <- input[!(rownames(input) == gene.2),]
              percents <- PercentAboveVM(input, threshold=sub_thresholds, quant = F) # 0.09 when quantiles are not used
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
          depth_data <- ExtractMerfishLamina(merfish_data, next_class)
          used_mer_depths[[label]] <- depth_data
        }
        
        type_densities <- depth_data[[1]]
        rel_depths <- depth_data[[2]]
        prop_depths <- depth_data[[3]]
        merfish_type_props <- table(rel_depths$use_type)/nrow(rel_depths)
        
        prob.closest <- 0
        label <- '9999 ZNo Good Match'
        each.percent <- each.percent[(rownames(each.percent) %in% names(merfish_type_props) & rownames(each.percent) %in% colnames(type_densities)),]
        each.percent[,'type.probs'] <- merfish_type_props[rownames(each.percent)]
        each.percent <- each.percent[sample(1:nrow(each.percent), nrow(each.percent), replace = F),]
        
        for (k in 1:nrow(each.percent)) {
          
          # Bernoulli NB:
          test.type <- each.percent[k,'id']
          #type.prob <- each.percent[k,'type.probs']
          #type.depth.prob <- type_densities[as.character(cell.depth),test.type] + 1e-4
          prop.depth <- prop_depths$prop[prop_depths$use_type == as.character(test.type) & prop_depths$Y == cell.depth]
          if(length(prop.depth) == 1) {
            type.depth.prob <- prop.depth
          }
          else{
            type.depth.prob <- 0
          }
          pos.genes <- prod((as.numeric(each.percent[k,colnames(cell.pos)])) + 1e-15)
          neg.genes <- prod((1-(as.numeric(each.percent[k,colnames(cell.neg)]))) + 1e-15) 
          #prob <- type.prob*type.depth.prob*pos.genes*neg.genes
          prob <- type.depth.prob*pos.genes*neg.genes
          
          if (prob > prob.closest) {
            prob.closest <- prob
            label <- as.character(test.type)
             # if (grepl('NP |CT ', label)){
             #   print(cell.pos)
             #   print(cell.neg)
             #   print(pos.genes)
             #   print(neg.genes)
             #   print(label)
             # }
            
          }
          
          else {
            
            next
            
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

HierarchicalGMMBayes <- function(hcr.table, bin_table, merfish_data = NULL, use.lamina = T, depth.round = 2) {
  
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
          depth_data <- ExtractMerfishLamina(merfish_data, next_class, bin = T, bin_to = depth.round, smooth = T)
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

ExtractMerfishLamina <- function(mer_metadata, selected_class, quant_low = 0.02, quant_high = 0.91, bin = F, bin_to = 2, smooth = T) {
  
  # ggplot(mop_merdata[(mop_merdata$z %in% c('10.2', '10.4')),], aes(x = newx, y = -new.quant.y, color = my_supertype)) +
  #  geom_point(size = 2.0) +
  #  scale_color_manual(values = cluster_colors) +
  #  #ylim(quantile(-mop_merdata$newy,0.10), quantile(-mop_merdata$newy,0.99)) +
  # #geom_point(data = excit_merdata[excit_merdata$allen_subclass == '029 L6b CTX Glut',], aes(x = x, y = -y), color = 'red', size = 1.0) +
  #  theme_minimal() +
  #  theme(legend.position = "none")
  
  # Extract depths
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
    summarise(Q1 = quantile(Y, 0.01), Q9 = quantile(Y, 0.99)) # Only plot the middle 80% of cells
  
  rel.depths.mod <- data.frame()
  for (i in 1:nrow(ylims)) {
    max.Y <- as.numeric(ylims[i,'Q9'])
    min.Y <- as.numeric(ylims[i,'Q1'])
    rel.depths.rows <- rel.depths$use_type == as.character(ylims[i,'use_type'])
    rel.depths.inquant <- rel.depths[(rel.depths.rows & (rel.depths$Y >= min.Y) & (rel.depths$Y < max.Y)),]
    rel.depths.mod <- rbind(rel.depths.mod, rel.depths.inquant)
  }
  rel.depths.mod[,'simple_type'] <- sub('(.*)_\\d+', '\\1', rel.depths.mod$use_type)
  
  #big_groups <- names(table(rel.depths.mod$use_type)[table(rel.depths.mod$use_type) > 10])
  #rel.depths.mod <- rel.depths.mod[rel.depths.mod$use_type %in% big_groups,]
  
   # ggplot(rel.depths.mod, aes(x=use_type, y=Y, fill=use_type)) +
   #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill ='white'), axis.text.x = element_text(angle = 45, vjust =1, hjust=1)) +
   #   ggtitle('MERFISH depths') +
   #   geom_violin(trim=T, scale='width') +
   #   scale_x_discrete(name = names(cluster_colors)) +
   #   #geom_boxplot() +
   #   #geom_point()
   #   scale_fill_manual(values=cluster_colors) +
   #   geom_hline(yintercept = c(-0.075, -0.22, -0.4, -0.55, -0.7, -0.9, -0.95)) +
   #   guides(fill = 'none')
   #   #ylim(1,0)
  
  if (bin) {
    rel.depths.mod$Y <- round(rel.depths.mod$Y, bin_to) # Round the depths to "bin" them to 100 values
  }
  
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

CombineGenes <- function(table, genes_combined) {
  for (i in 1:length(genes_combined)) {
    
    gene.1 <- genes_combined[[i]][1]
    gene.2 <- genes_combined[[i]][2]
    table[,gene.1] <- pmax(table[,gene.1], table[,gene.2])
    table <- table[,colnames(table) != gene.2]
    colnames(table)[colnames(table) == gene.1] <- paste(gene.1, gene.2, sep='+')
    
  }
  return(table)
}

DotProductABC <- function(mop_object, hcr.table, merfish_data) {
  
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

#### Load in sequencing related files: ####

setwd('~/Desktop/Local code/HCR Annotation/')

#alm_object <- readRDS('alm_object.rds')
mop_neurons <- readRDS('mop_my_neurons.rds') # or mop_object if you want the full thing

merfish_metadata <- read.table('merfish_depth_metadata.csv', sep = ',', header = T)

gene_thresholds <- read.table('Classification for meta Seq data/gene_thresholds/gene_thresholds_hier.csv', sep = ',', header = T)
colnames(gene_thresholds)[colnames(gene_thresholds) == 'Ctgf.Penk'] <- 'Ctgf+Penk'
rownames(gene_thresholds) <- gene_thresholds$class
colnames(gene_thresholds)[4:ncol(gene_thresholds)] <- sub('X0', '0', colnames(gene_thresholds)[4:ncol(gene_thresholds)])

sub_gene_matrix <- read.table('Classification for meta Seq data/binary_sub_gene_matrix_my_type.csv', sep = ',', header = T)
num_sub_gene_matrix <- sub_gene_matrix[,1:(ncol(sub_gene_matrix)-1)]
colnames(num_sub_gene_matrix)[colnames(num_sub_gene_matrix) == 'depth'] <- "newy"
colnames(num_sub_gene_matrix)[colnames(num_sub_gene_matrix) == 'Ctgf.Penk'] <- "Ctgf+Penk"
num_sub_gene_matrix <- num_sub_gene_matrix[num_sub_gene_matrix$newy < 0,]

# accuracies <- read.table('Classification for meta Seq data/Classification_accuracy_by_edge.csv', sep = ',', header = T)
# 
# graph_df <- data.frame(all = mop_object_hcr$all_cells,
#                        nt = mop_object_hcr$allen_NT,
#                        class = mop_object_hcr$allen_class,
#                        subclass = mop_object_hcr$allen_subclass,
#                        supertype = mop_object_hcr$allen_supertype)
# 
# graph_df <- unique(graph_df)
# graph_df <- graph_df[order(graph_df$supertype, decreasing = F),]
# rownames(graph_df) <- 1:nrow(graph_df)
# 
# edges_all_nt <- graph_df %>% select(all, nt) %>% unique %>% rename(from=all, to=nt)
# edges_nt_class <- graph_df %>% select(nt, class) %>% unique %>% rename(from=nt, to=class)
# edges_class_subclass <- graph_df %>% select(class, subclass) %>% unique %>% rename(from=class, to=subclass)
# edges_subclass_supertype <- graph_df %>% select(subclass, supertype) %>% unique %>% rename(from=subclass, to=supertype)
# edge_list <- rbind(edges_all_nt, edges_nt_class, edges_class_subclass, edges_subclass_supertype)
# 

#### Preparing the new single cell seq dataset to extract MOp specfic cells: ####

# ctx_object <- LoadH5Seurat('WMB-10Xv2-Isocortex-1-raw.h5seurat', meta.data = F, commands = F, misc = F, tools = F)
# 
# metadata <- read.csv('singlecell_metadata.csv')
# rownames(metadata) <- metadata$cell_label
# metadata_mop <- metadata[metadata$region_of_interest_acronym == 'MOp',]
# 
# #mop2 <- ctx_object[,colnames(ctx_object) %in% metadata$cell_label]
# mop_object <- ctx_object[,colnames(ctx_object) %in% metadata$cell_label]
# #saveRDS(trimmed_ctx, 'mop_object.rds')
# 
# mop_object <- readRDS('mop_object.rds')
# 
# library <- read.csv('cluster_to_cluster_annotation_membership.csv')
# anno_library <- matrix(c(library[library$cluster_annotation_term_set_name == 'cluster', 'cluster_annotation_term_name'],
#                             library[library$cluster_annotation_term_set_name == 'supertype', 'cluster_annotation_term_name'],
#                             library[library$cluster_annotation_term_set_name == 'subclass', 'cluster_annotation_term_name'],
#                             library[library$cluster_annotation_term_set_name == 'class', 'cluster_annotation_term_name'],
#                             library[library$cluster_annotation_term_set_name == 'neurotransmitter', 'cluster_annotation_term_name']),
#                        nrow = length(unique(library$cluster_alias)),
#                        ncol = 5,
#                        dimnames = list(unique(library$cluster_alias), c('cluster', 'supertype', 'subclass', 'class', 'neurotransmitter')))
# 
# 
# mop_object[['allen_cluster']] <- anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'cluster']
# mop_object[['allen_supertype']] <- anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'supertype']
# mop_object[['allen_subclass']] <- anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'subclass']
# mop_object[['allen_class']] <- anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'class']
# mop_object[['allen_NT']] <- anno_library[as.character(metadata[colnames(mop_object), 'cluster_alias']), 'neurotransmitter']
# 
# # Perform clustering:
# 
# mop_object <- NormalizeData(mop_object) %>%
#   FindVariableFeatures(.) %>%
#   ScaleData(.) %>%
#   RunPCA(.) %>%
#   FindNeighbors(., dims=1:30) %>%
#   FindClusters(., resolution=1) %>%
#   RunUMAP(., dims=1:30)
# 
# DimPlot(mop_object, reduction = 'umap', cols = cluster_colors, label = T, repel = T, group.by = 'allen_supertype') + NoLegend()
# 
# types_wanted <- unique(mop_object$allen_subclass)
# types_wanted <- c(types_wanted[grepl('CTX', types_wanted)], types_wanted[grepl('0(53)|(56)|(47)|(49)|(52)|(51)|(46)|(50) .* Gaba', types_wanted)])
# types_wanted <- types_wanted[!grepl('(CLA)|(STR)', types_wanted)]
# mop_neurons <- mop_object[,mop_object$allen_subclass %in% types_wanted]
# 
# mop_neurons <- NormalizeData(mop_neurons) %>%
#   FindVariableFeatures(.) %>%
#   ScaleData(.) %>%
#   RunPCA(.) %>%
#   FindNeighbors(., dims=1:30) %>%
#   FindClusters(., resolution=1) %>%
#   RunUMAP(., dims=1:30)
# 
# 
# mop_neurons[['my_subclass']] <- sub('\\d{3} (.*)', '\\1', sub_mop_neurons$allen_subclass)
# mop_neurons[['my_supertype']] <- sub('\\d{4} (.*)', '\\1', sub_mop_neurons$allen_supertype)
# 
# supertypes <- unique(mop_neurons$my_supertype)
# supertypes <- supertypes[order(supertypes)]
# cluster_colors <- c(colorRampPalette(c('#d4775e', '#a33e59', '#765d8c', '#79a4b9', '#afd435', '#fdec28'))(length(supertypes[grepl('L\\d', supertypes)])),
#                     colorRampPalette(c('#ffce74', '#ff9fa1', '#ffa9e1', '#ca94ff', '#94acff'))(length(supertypes[!grepl('L\\d', supertypes)])))
# names(cluster_colors) <- c(supertypes[grepl('L\\d', supertypes)], supertypes[!grepl('L\\d', supertypes)])
# 
# 
# sub_mop_neurons <- mop_neurons[,sample(1:ncol(mop_neurons), ncol(mop_neurons)/4, replace = F)]
# DimPlot(sub_mop_neurons, reduction = 'umap', group.by = 'my_supertype', cols = cluster_colors, label = T, repel = T) + NoLegend()
# 
# saveRDS(sub_mop_neurons, 'sub_mop_neurons.rds')
# saveRDS(mop_neurons, 'mop_neurons.rds')


# Change gene names from ensembl to commong gene symbols

###

################## 04/16/2022 SynMap Cortical dataset: ######

anno.dir <- '/Volumes/SynMap Data/2022-04-16 Full Cortical Probe set/Will/'
reg.file.list <- list.files(anno.dir, full.names = T)
n.slices <- 1
gene.file <- sub('.*/(.*)_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(1900), 
                 ymin = c(70),
                 zmin = c(0),
                 zmax = c(24),
                 angle = c(25*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Slco2a1', 'Fam84b', 'Npnt', 
               'Myl4', 'Syt6', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Dkkl1', 'Slco2a1'))

file.colnames <- c(cell_id = 'labelimage_oid', X = 'Center.of.the.object_0', Y = 'Center.of.the.object_1', Z = 'Center.of.the.object_2',
                   Roi_mean = 'Mean.Intensity', Neighbor_mean = 'Mean.Intensity.in.neighborhood', Label = 'Predicted.Class')

################## 10/22/2022 SynMap Cortical dataset with Voltron cells: ######

anno.dir <- '/Volumes/SynMap Data/2022-10-22 Puromycin test Ctrl03/'
reg.file.list <- list.files(anno.dir, full.names = T)
n.slices <- 1
gene.file <- sub('.*/(.+)_\\d{3}_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(2007), 
                 ymin = c(0),
                 zmin = c(7),
                 zmax = c(23),
                 angle = c(10*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Fam84b', 'Npnt', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'JF585')
genes_combined <- list(c('Rorb', 'Vip'))

file.colnames <- c(cell_id = 'labelimage_oid', X = 'Center.of.the.object_0', Y = 'Center.of.the.object_1', Z = 'Center.of.the.object_2',
                   Roi_mean = 'Mean.Intensity', Neighbor_mean = 'Mean.Intensity.in.neighborhood', Label = 'Predicted.Class')

################## 12/05/2023 SynMap Cortical dataset with Voltron cells: ######

anno.dir <- '/Volumes/SynMap Data/2023-12-05 SynMap/s02/'
reg.file.list <- list.files(anno.dir, full.names = T)
n.slices <- 1
gene.file <- sub('.*/(.+)_\\d{3}_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(2550), 
                 ymin = c(300),
                 zmin = c(1),
                 zmax = c(34),
                 angle = c(-15*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 
               'Slco2a1', 'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip', 'JF585')
genes_combined <- list()

file.colnames <- c(cell_id = 'labelimage_oid', X = 'Center.of.the.object_0', Y = 'Center.of.the.object_1', Z = 'Center.of.the.object_2',
                   Roi_mean = 'Mean.Intensity', Neighbor_mean = 'Mean.Intensity.in.neighborhood', Label = 'Predicted.Class')

################## 07/30/2024 SynMap Thalamus dataset: ######

anno.dir <- '/Volumes/SynMap Data/2024-07-30 SynMap Thalamus/Registered/'
reg.file.list <- c()
n.slices <- sum(grepl('Slice ', list.files(anno.dir)))
for (i in 1:n.slices){
  
  cur.dir <- paste(anno.dir, 'Slice ', as.character(i+2), '/Gene expression tables/', sep = '')
  reg.file.list <- c(reg.file.list, list.files(cur.dir, pattern = 'table', full.names = T))
  
}

gene.file <- sub('.*/(s\\d{2}_.*)\\d{3}_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(2050, 2050, 1890, 1860, 1865), #These values come from the rotated image in Fiji, without image size expansion upon rotation
                 ymin = c(384, 480, 370, 430, 405),
                 zmin = c(0, 1, 1, 1, 1),
                 zmax = c(24, 35, 35, 35, 35),
                 angle = c(0*pi/180, 2*pi/180, 12*pi/180, 6*pi/180, 4*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 'Slco2a1', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Ctgf', 'Penk'))

file.colnames <- c(cell_id = 'CellID', X = 'Center_X', Y = 'Center_Y', Z = 'Center_Z',
                   Roi_mean = 'MeanROIIntensity', Neighbor_mean = 'MeanNeighborhoodIntensity', Label = 'UserLabel')

#file.colnames <- c(cell_id = 'labelimage_oid', X = 'Center.of.the.object_0', Y = 'Center.of.the.object_1', Z = 'Center.of.the.object_2',
#                   Roi_mean = 'Mean.Intensity', Neighbor_mean = 'Mean.Intensity.in.neighborhood', Label = 'User.Label')

################## 08/02/2024 SynMap Thalamus dataset: ######

anno.dir <- '/Volumes/SynMap Data/2024-08-02 SynMap Thalamus/Registered/'
reg.file.list <- c()
n.slices <- sum(grepl('Slice ', list.files(anno.dir)))
for (i in 1:n.slices){
  
  cur.dir <- paste(anno.dir, 'Slice ', as.character(i+1), '/Gene expression tables/', sep = '')
  reg.file.list <- c(reg.file.list, list.files(cur.dir, pattern = 'table', full.names = T))
  
}

gene.file <- sub('.*/(s\\d{2}_.*)\\d{3}_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(1860, 1770, 1600, 1800, 1800),
                 ymin = c(155, 208, 180, 440, 50),
                 zmin = c(0, 1, 1, 0, 1),
                 zmax = c(35, 40, 40, 40, 35),
                 angle = c(9*pi/180, 6*pi/180, 10*pi/180, 22*pi/180, -7*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 'Slco2a1', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Ctgf', 'Penk'))

file.colnames <- c(cell_id = 'CellID', X = 'Center_X', Y = 'Center_Y', Z = 'Center_Z',
                   Roi_mean = 'MeanROIIntensity', Neighbor_mean = 'MeanNeighborhoodIntensity', Label = 'UserLabel')

#file.colnames <- c(cell_id = 'labelimage_oid', X = 'Center.of.the.object_0', Y = 'Center.of.the.object_1', Z = 'Center.of.the.object_2',
#                   Roi_mean = 'Mean.Intensity', Neighbor_mean = 'Mean.Intensity.in.neighborhood', Label = 'User.Label')

################## 08/11/2024 SynMap Contra MCtx dataset: ######

anno.dir <- '/Volumes/SynMap Data/2024-08-11 SynMap Contra MCtx/Registered/'
reg.file.list <- c()
n.slices <- sum(grepl('Slice ', list.files(anno.dir)))
for (i in 1:n.slices){
  
  cur.dir <- paste(anno.dir, 'Slice ', as.character(i), '/Gene expression tables/', sep = '')
  reg.file.list <- c(reg.file.list, list.files(cur.dir, pattern = 'table', full.names = T))
  
}

gene.file <- sub('.*/(s\\d{2}_.*)\\d{3}_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(1674, 1662, 1791, 1629, 1749), #ymax = c(1611, 1866, 1752, 2144),
                 ymin = c(258, 300, 192, 135, 228),
                 zmin = c(0, 0, 0, 0, 0),
                 zmax = c(35, 36, 40, 35, 35),
                 angle = c(7*pi/180, 17*pi/180, 12*pi/180, 9*pi/180, 7*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 'Slco2a1', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Ctgf', 'Penk'))

file.colnames <- c(cell_id = 'CellID', X = 'Center_X', Y = 'Center_Y', Z = 'Center_Z',
                   Roi_mean = 'MeanROIIntensity', Neighbor_mean = 'MeanNeighborhoodIntensity', Label = 'UserLabel')

################## 08/14/2024 SynMap Contra MCtx dataset: ######
anno.dir <- '/Volumes/SynMap Data/2024-08-14 SynMap Contra MCtx/Registered/'
reg.file.list <- c()
n.slices <- sum(grepl('Slice ', list.files(anno.dir)))
for (i in 1:n.slices){
  
  cur.dir <- paste(anno.dir, 'Slice ', as.character(i+1), '/Gene expression tables/', sep = '')
  reg.file.list <- c(reg.file.list, list.files(cur.dir, pattern = 'table', full.names = T))
  
}

gene.file <- sub('.*/(s\\d{2}_.*)\\d{3}_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(1646, 1823, 1931, 1655), #ymax = c(1611, 1866, 1752, 2144),
                 ymin = c(191, 122, 118, 169),
                 zmin = c(0, 0, 0, 0),
                 zmax = c(46, 46, 46, 46),
                 angle = c(18*pi/180, 2*pi/180, 12*pi/180, 12*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 'Slco2a1', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Ctgf', 'Penk'))

file.colnames <- c(cell_id = 'CellID', X = 'Center_X', Y = 'Center_Y', Z = 'Center_Z',
                   Roi_mean = 'MeanROIIntensity', Neighbor_mean = 'MeanNeighborhoodIntensity', Label = 'UserLabel')


################## 10/08/2024 SynMap Contra MCtx dataset: ######
anno.dir <- '/Volumes/SynMap Data/2024-10-08 SynMap Contra MCtx/Registered/'
reg.file.list <- c()
n.slices <- sum(grepl('Slice ', list.files(anno.dir)))
for (i in 1:n.slices){
  
  cur.dir <- paste(anno.dir, 'Slice ', as.character(i), '/Gene expression tables/', sep = '')
  reg.file.list <- c(reg.file.list, list.files(cur.dir, pattern = 'table', full.names = T))
  
}

gene.file <- sub('.*/(s\\d{2}_.*)\\d{3}_.*_table\\.csv', '\\1', reg.file.list)

img.lims <- list(ymax = c(1711, 1966, 1902, 2244), #ymax = c(1611, 1866, 1752, 2144),
                 ymin = c(213, 123, 78, 250),
                 zmin = c(0, 0, 0, 0),
                 zmax = c(40, 40, 36, 39),
                 angle = c(-5*pi/180, -12*pi/180, -3*pi/180, -5*pi/180))

hcr_genes <- c('Penk', 'Calb1', 'Lamp5', 'Rorb', 'Ccdc80', 'Pamr1', 'Dkkl1', 'Npnt', 'Fam84b', 'Slco2a1', 
               'Myl4', 'Syt6', 'Ctgf', 'Slc32a1', 'Pvalb', 'Sst', 'Vip')
genes_combined <- list(c('Ctgf', 'Penk'))

file.colnames <- c(cell_id = 'CellID', X = 'Center_X', Y = 'Center_Y', Z = 'Center_Z',
                   Roi_mean = 'MeanROIIntensity', Neighbor_mean = 'MeanNeighborhoodIntensity', Label = 'UserLabel')

#################### Read in the files ###########################

#first.file <- read.csv(reg.file.list[1], sep = ',', header = TRUE)
#max.size <- quantile(first.file$Size.in.pixels, 0.85)
#min.size <- quantile(first.file$Size.in.pixels, 0.15)
#first.file <- first.file[first.file$Size.in.pixels < max.size & first.file$Size.in.pixels > min.size,]
mean.vals <- data.frame(row.names = c(), slice = c(), X = c(), Y = c(), Z = c(), volt = c())
gene.file.len <- length(hcr_genes) - length(genes_combined)

for (j in 1:n.slices){
  if (j == 1){
    print('Compiling Ilastik annotations...')
    pb <- txtProgressBar(min = 0, max = length(reg.file.list), initial = 0, style = 3)
  }
  name = reg.file.list[j*gene.file.len-(gene.file.len-1)]
  current.file <- read.csv(name, sep = ',', header = TRUE)
  cur.mean.vals <- data.frame(cell.id = as.character(current.file[,file.colnames['cell_id']]),
                              slice = as.numeric(sub('s(\\d{2})_.*', '\\1', gene.file[j*gene.file.len])),
                              X = current.file[,file.colnames['X']],
                              Y = current.file[,file.colnames['Y']],
                              Z = current.file[,file.colnames['Z']],
                              volt = rep(FALSE, nrow(current.file)))
  
  for (i in 1:gene.file.len) {
    idx <- i + (j*gene.file.len-gene.file.len)
    name = reg.file.list[idx]
    current.file <- read.csv(name, sep = ',', header = TRUE)
    if ('Label 1' %in% current.file[,file.colnames['Label']]) {
      user_labels <- (current.file[,file.colnames['Label']] == 'Label 1')*1
    }
    else {
      
      #if (grepl('JF585', name)){user_labels <- current.file[,'User.Label']}
      #else{user_labels <- current.file[,file.colnames['Label']]}
      user_labels <- current.file[,file.colnames['Label']]
    }
    gene.name <- sub('s\\d{2}_(.*)', '\\1', gene.file[idx])
    cur.mean.vals[[paste(gene.name, '.mean.intensity', sep='')]] <- current.file[,file.colnames['Roi_mean']]
    cur.mean.vals[[paste(gene.name, '.mean.neighborhood.intensity', sep='')]] <- current.file[,file.colnames['Neighbor_mean']]
    cur.mean.vals[[paste(gene.name, '.object.classification', sep='')]] <- user_labels
    setTxtProgressBar(pb, idx)
  }
  mean.vals <- rbind(mean.vals, cur.mean.vals)
  
}
close(pb)


#mean.vals <- mean.vals[((mean.vals$X > xmin) & (mean.vals$X < xmax) & (mean.vals$Y < ymax) & (mean.vals$Y > ymin) & (mean.vals$Z >= zmin)),]

### Rotating cell positions for laminar structure: #####

new.x <- c()
new.y <- c()

for (i in 1:n.slices){
  
  slice.oi <- unique(mean.vals$slice)[i]
  rows.oi <- mean.vals[mean.vals$slice == slice.oi,]
  origin.x <- mean(rows.oi$X)
  origin.y <- mean(rows.oi$Y)
  
  cur.new.x <- (rows.oi$X-origin.x)*cos(img.lims$angle[i]) - (rows.oi$Y-origin.y)*sin(img.lims$angle[i]) + origin.x
  cur.new.y <- (rows.oi$X-origin.x)*sin(img.lims$angle[i]) + (rows.oi$Y-origin.y)*cos(img.lims$angle[i]) + origin.y
  
  #new.x <- (new.x - min(new.x))/(max(new.x) - min(new.x))
  #new.y <- (new.y - min(new.y))/(max(new.y) - min(new.y))
  norm.new.x <- (cur.new.x - min(cur.new.x))/(max(cur.new.x) - min(cur.new.x))
  norm.new.y <- (cur.new.y - img.lims$ymin[i])/(img.lims$ymax[i] - img.lims$ymin[i])
  #plot(new.x, -new.y, pch = 16, cex = 1.0)
  new.x <- c(new.x, norm.new.x)
  new.y <- c(new.y, norm.new.y)
  
}

mean.vals$newx <- new.x
mean.vals$newy <- -new.y
plot(mean.vals$newx, mean.vals$newy, col = mean.vals$slice, pch = 16, cex = 1.0)

#mop_neurons <- SetIdent(mop_neurons, value = 'my_supertype')
#mop_neurons@assays$RNA@meta.features$gene_symbol[mop_neurons@assays$RNA@meta.features$gene_symbol == 'Lratd2'] <- 'Fam84b'
#mop_neurons@assays$RNA@meta.features$gene_symbol[mop_neurons@assays$RNA@meta.features$gene_symbol == 'Ccn2'] <- 'Ctgf'
# alm_object <- SetIdent(alm_object, value = 'my_type')
# alm_object <- BuildClusterTree(alm_object, dims = 1:30)
# PlotClusterTree(alm_object, edge.color = 'blue')
# data.tree <- Tool(object = alm_object, slot = "BuildClusterTree")
# ape::plot.phylo(x = data.tree, direction = "downwards")
# tree.types <- read.table('alm_tree_node.txt', sep='\t', header=TRUE)
# my_type_1 <- tree.types$New.node.name.1
# names(my_type_1) <- tree.types$Cluster.Label..my_type.
# translated_my_type_1 <- my_type_1[alm_object$my_type]
# names(translated_my_type_1) <- colnames(alm_object)
# alm_object[['my_type_1']] <- translated_my_type_1
# #my_type_2 <- tree.types$New.node.name.2
# #names(my_type_2) <- tree.types$Cluster.Label..my_type.
# #alm_object[['my_type_2']] <- my_type_2[alm_object$my_type]
# #my_type_3 <- tree.types$New.node.name.3
# #names(my_type_3) <- tree.types$Cluster.Label..my_type.
# #alm_object[['my_type_3']] <- my_type_3[alm_object$my_type]
# 
# type_level <- my_type_1 # my_type_1 my_type_2
# select_type <- 'my_type_1'
# 
# orig_types <- unique(alm_object$cluster_label)
# cell.counts <- table(alm_object[[select_type]])
# final.cell.types <- cell.counts[cell.counts >= 10] # Exclude cell types that have fewer than 10 cells in the dataset
#names(final.cell.types) <- type.translation[names(final.cell.types), 'my_type']

#mop_neurons@assays$RNA@counts@Dimnames[[1]] <- mop_neurons@assays$RNA@meta.features[rownames(mop_neurons),1]
#mop_neurons@assays$RNA@data@Dimnames[[1]] <- mop_neurons@assays$RNA@meta.features[rownames(mop_neurons),1]

#mop_object_hcr <- mop_neurons[hcr_genes, mop_neurons$allen_supertype %in% merfish_metadata$allen_supertype]
#mop_object_hcr[['all_cells']] <- rep('cells', ncol(mop_object_hcr))

meta_gene_matrix <- read.table('Classification with Gaussian mixture/meta_log_counts.csv', sep = ',', header = T)
#meta_gene_matrix[,'my_supertype'] <- meta_gene_matrix$allen_supertype
#meta_gene_matrix$my_supertype[!grepl('L\\d', meta_gene_matrix$my_supertype)] <- meta_gene_matrix$allen_subclass[!grepl('L\\d', meta_gene_matrix$my_supertype)]
meta_gene_matrix <- CombineGenes(meta_gene_matrix, genes_combined)

mop_object_hcr <- meta_gene_matrix
mop_object_hcr[,'all_cells'] <- rep('all', nrow(mop_object_hcr))

meta_bin <- read.table('Classification with Gaussian mixture/All_genes_bin.csv', sep = ',', header = T)
meta_bin <- as.data.frame(meta_bin)
meta_bin <- meta_bin[,hcr_genes[hcr_genes != 'JF585']]
meta_bin <- cbind(meta_bin, mop_object_hcr[,(ncol(mop_object_hcr)-4):ncol(mop_object_hcr)])
for (i in 1:(ncol(meta_bin)-5)){
  if (2 %in% meta_bin[,i]){
    if (colnames(meta_bin)[i] %in% c('Rorb', 'Fam84b', 'Ccdc80', 'Pvalb', 'Vip', 'Npnt')){
      meta_bin[,i] <- meta_bin[,i] == max(meta_bin[,i])
    }
    else{
      meta_bin[,i] <- meta_bin[,i] == max(meta_bin[,i])-1
      print(colnames(meta_bin)[i])
    }
  }
  else{
    meta_bin[,i] <- meta_bin[,i] > 0
  }
}
meta_bin[,1:(ncol(meta_bin)-5)] <- meta_bin[,1:(ncol(meta_bin)-5)] * 1
meta_bin <- CombineGenes(meta_bin, genes_combined)
#meta_bin <- meta_bin[!grepl('0013|0015|0016|0017|0027|0092|0021|0111|0114|0116|0119|0113', meta_bin$my_supertype),]
mop_object_hcr <- meta_bin
mop_object_hcr[,'all_cells'] <- rep('all', nrow(mop_object_hcr))
subclass_rows <- unique(meta_bin[,c('my_supertype', 'allen_subclass')])
subclass_dict <- subclass_rows$allen_subclass
names(subclass_dict) <- subclass_rows$my_supertype

##### Attempting Bayesian classification: ####

ilastik.anno <- mean.vals[,grepl('.*\\.object\\.classification', colnames(mean.vals))]
colnames(ilastik.anno) <- unique(sub('s\\d{2}_(.*)', '\\1', gene.file))
if ('JF585' %in% colnames(ilastik.anno)) {
  voltron.anno <- ilastik.anno[,'JF585']
  names(voltron.anno) <- rownames(ilastik.anno)
  ilastik.anno <- ilastik.anno[,colnames(ilastik.anno)[colnames(ilastik.anno) != 'JF585']]
} else if (TRUE %in% mean.vals$volt) {
  voltron.anno <- mean.vals$volt
  names(voltron.anno) <- rownames(ilastik.anno)
}

ilastik.anno <- ilastik.anno[,colnames(meta_bin)[1:(ncol(meta_bin)-5)]]
ilastik.anno$newy <- mean.vals$newy

#### HCR data classification ####

#bayes.classif <- HierarchicalBayes(ilastik.anno, gene_thresholds, merfish_data = merfish_metadata)
bayes.classif <- HierarchicalGMMBayes(ilastik.anno, meta_bin, merfish_data = merfish_metadata, use.lamina = T, depth.round = 3)
bayes.classif <- DotProductABC(mop_object, ilastik.anno, merfish_metadata)

bayes.classif$classif.label[is.na(bayes.classif$classif.label)] <- 'Unknown'
names(bayes.classif$classif.label) <- mean.vals$cell.id
names(bayes.classif$classif.prob) <- mean.vals$cell.id

order <- sub('\\d{3,4} (.*)', '\\1', unique(bayes.classif$classif.label))
order2 <- as.numeric(sub('(\\d{3,4}) .*', '\\1', unique(bayes.classif$classif.label)))
order <- order(order, order2)
bayes.classif$classif.label <- factor(bayes.classif$classif.label, levels = unique(bayes.classif$classif.label)[order])

#### Test data classification ####

merfish.cells <- ExtractMerfishLamina(merfish_metadata, 'my_supertype')
merfish.rel.depths <- merfish.cells[[2]]
newy <- c()
mop_for_test <- mop_object_hcr[1,1:(ncol(ilastik.anno))]
for (i in bayes.classif$classif.label) {
  if (i %in% mop_object_hcr$my_supertype) {
    new_cell <- mop_object_hcr[sample(which(mop_object_hcr$my_supertype == i),1),1:(ncol(ilastik.anno))]
    mop_for_test <- rbind(mop_for_test, new_cell)
  }
  else {
    next
  }
}
#mop_for_test <- mop_object_hcr[mop_object_hcr$my_supertype %in% merfish.rel.depths$use_type,]
for (i in mop_for_test$my_supertype){
  s_depth <- sample(merfish.rel.depths$Y[merfish.rel.depths$use_type == i], 1)
  newy <- c(newy, s_depth)
}
mop_for_test_bin <- mop_for_test[,1:(ncol(mop_for_test)-1)] #mop_for_test[,1:(ncol(mop_for_test)-5)]
mop_for_test_bin <- cbind(mop_for_test_bin, newy)

# For shuffling binary values (to test false negatives):
# for (i in 1:nrow(mop_for_test_bin)){
#   pos_ones <- which(mop_for_test_bin[i,] == 1)
#   if (length(pos_ones) == 0){
#     next
#   }
#   gene_to_shuffle <- sample(pos_ones, 1)
#   #gene_to_shuffle <- sample(1:(ncol(mop_for_test_bin)-1), 2)
#   #mop_for_test_bin[i,gene_to_shuffle] <- (mop_for_test_bin[i,gene_to_shuffle] == 0)*1
#   mop_for_test_bin[i,gene_to_shuffle] <- 0
# }

#test.bayes.classif <- HierarchicalBayes(as.data.frame(num_sub_gene_matrix), gene_thresholds, merfish_data = merfish_metadata)
test.bayes.classif <- HierarchicalGMMBayes(mop_for_test_bin, meta_bin, merfish_data = merfish_metadata, use.lamina = T, depth.round = 3)
# orig.test.bayes.classif <- test.bayes.classif

names(test.bayes.classif$classif.label) <- rownames(mop_for_test)
names(test.bayes.classif$classif.prob) <- rownames(mop_for_test)
test.bayes.classif$classif.label[test.bayes.classif$classif.label == 'Errored'] <- 'Unknown'
test.bayes.classif$real.label <- mop_for_test$my_supertype
names(test.bayes.classif$real.label) <- rownames(mop_for_test)

order <- sub('\\d{3,4} (.*)', '\\1', unique(test.bayes.classif$real.label))
order2 <- sub('(\\d{3,4}) .*', '\\1', unique(test.bayes.classif$real.label))
order <- order(order, order2)
test.bayes.classif$classif.label <- factor(test.bayes.classif$classif.label, levels = names(cluster_colors))
test.bayes.classif$real.label <- factor(test.bayes.classif$real.label, levels = names(cluster_colors))

test.bayes.classif$depth <- mop_for_test_bin$newy

##### Spatial plots for Bayesian classifs #######

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

##### Write simple classif output: #####
subclasses <- sapply(as.character(bayes.classif$classif.label), function(x) ifelse(x %in% names(subclass_dict), subclass_dict[x], x))
table.to.write <- data.frame(cell_ids = as.numeric(names(bayes.classif$classif.label)),
                             slice = mean.vals$slice,
                             depth = mean.vals$newy,
                             supertype = as.character(bayes.classif$classif.label),
                             subclass = subclasses,
                             classif_cols = cluster_colors[as.character(bayes.classif$classif.label)],
                             slice_depth = mean.vals$Z)
write.table(table.to.write, paste(anno.dir, 'HCR_classification.csv', sep = ''), sep = ',', col.names = T, row.names = F)

#### Scatterplot of types: #####

#par(mar=c(5, 4, 4, 15), xpd=TRUE)
excit.types <- grepl('L\\d', bayes.classif$classif.label) # 'L\\d( |/)' Exclude L6b cells by adding the space after the number
inhib.types <- !grepl('L\\d', bayes.classif$classif.label) | !grepl('ZNo', bayes.classif$classif.label)
custom.types <- grepl('ET', bayes.classif$classif.label)
type <- c(excit.types|inhib.types)
big_groups <- names(table(bayes.classif$classif.label))[table(bayes.classif$classif.label) > 10]
scatter_data <- data.frame(classif=bayes.classif$classif.label, X=mean.vals$X, Y=-mean.vals$Y, Z=mean.vals$Z)
scatter_data <- scatter_data[scatter_data$classif %in% big_groups,]

ggplot(scatter_data[!grepl('ZNo ', scatter_data$classif),], aes(x=X, y=Y, col=classif)) +
  geom_point() +
  ylim(c(min(-img.lims$ymax), max(-img.lims$ymin))) +
  theme_minimal() +
  scale_color_manual(values = cluster_colors) +
  guides(color = 'none')#+
  #geom_point(data = scatter_data[voltron.anno == 1 & !grepl('ZNo ', scatter_data$classif),], 
  #           color = 'red', shape = 1, stroke = 2)

## Plotting laminar position of cells, violin. Only works if you did the rotation at the start

excit.types <- grepl('L\\d', bayes.classif$classif.label)
inhib.types <- !grepl('L\\d', bayes.classif$classif.label) + grepl('ZNo|Err|all', bayes.classif$classif.label)
custom.types <- grepl('L6b', bayes.classif$classif.label)
type <- c(excit.types | inhib.types)
big_groups <- names(table(bayes.classif$classif.label))[table(bayes.classif$classif.label) > 10]
#rel.depths <- mean.vals[type, c('X','Y')]
rel.depths <- data.frame(Y = ilastik.anno$newy[type])
rel.depths$use_type <- bayes.classif$classif.label[type]
rel.depths <- rel.depths[rel.depths$use_type %in% big_groups,]

rel.depths <- data.frame(Y = mop_for_test_bin$newy[type])
rel.depths$use_type <- test.bayes.classif$classif.label[type]
#rel.depths <- rel.depths[,]

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
  geom_violin(trim=T, scale='width') + #area, count
  ylim(c(-1.0, 0)) +
  #geom_vline(xintercept = c(-0.075, -0.25, -0.4, -0.55, -0.7, -0.92, -1), color = 'grey80') +
  #ggridges::geom_density_ridges(mapping = aes(y = use_type, x = Y)) +
  #xlim(c(-1,0)) +
  scale_fill_manual(values=cluster_colors) +
  guides(fill = 'none')

ggplot(rel.depths.mod, aes(y=Y, group=use_type, fill=use_type, color = use_type)) + 
  geom_density(aes(x=..count.., y=Y), alpha = 0.3, adjust = 0.5) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal()
#xlim(-1.05, 0.05)

#### Comparison of HCR population to MERFISH: ####

merfish.cells <- ExtractMerfishLamina(merfish_metadata, 'my_supertype', bin = F, bin_to = 2)
type_densities <- merfish.cells[[1]]
rel_depths <- merfish.cells[[2]]
prop_depths <- merfish.cells[[3]]
#rel_depths <- rel_depths[rel_depths$use_type %in% bayes.classif$classif.label,]
merfish_type_props <- table(rel_depths$use_type)/nrow(rel_depths)
all_types <- bayes.classif$classif.label[!grepl('(ZNo)|(VV)', bayes.classif$classif.label)]
all_type_props <- table(all_types)/length(all_types)

voltron.cells <- as.character(bayes.classif$classif.label)#[voltron.anno == 1]
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
both.props <- data.frame(voltron = c(rep('Voltron', length(voltron.props)), rep('All', length(all_type_props)), rep('Merfish', length(merfish_type_props))),
                         class = factor(c(names(voltron.props), names(all_type_props), names(merfish_type_props)), levels = names(cluster_colors)[1:39]),
                         prop = c(voltron.props, all_type_props, merfish_type_props))

ggplot(both.props, aes(x = voltron, y = prop, fill = class)) +
  theme_minimal() +
  geom_bar(position = 'stack', stat = 'identity', color = 'white') +
  scale_fill_manual(values=cluster_colors)

#### Dot plot showing proportion of HCR cells in each type that express each gene. To compare to above

bayes.dot.df <- cbind(mop_for_test_bin*1, test.bayes.classif$real.label)
bayes.dot.df <- cbind(ilastik.anno*1, bayes.classif$classif.label)
names(bayes.dot.df)[ncol(bayes.dot.df)] <- 'classif.label'
bayes.dot.tidy <- bayes.dot.df %>% group_by(classif.label) %>% summarise(across(colnames(bayes.dot.df)[1:16], sum))
bayes.dot.n <- bayes.dot.df %>% group_by(classif.label) %>% count()
bayes.dot.prop <- mapply(`/`, data.frame(bayes.dot.tidy[,2:ncol(bayes.dot.tidy)]), bayes.dot.n[,2])
colnames(bayes.dot.prop) <- sub('\\.', '\\+', colnames(bayes.dot.prop))

bayes.dot.prop <- data.frame(bayes.dot.prop[,colnames(mop_for_test_bin[,1:(ncol(mop_for_test_bin)-1)])],
                             id = bayes.dot.n[,1])
bayes.dot.prop <- data.frame(bayes.dot.prop[,colnames(ilastik.anno[,1:(ncol(ilastik.anno)-1)])],
                             id = bayes.dot.n[,1])
bayes.dot.melt <- melt(bayes.dot.prop)

ggplot(bayes.dot.melt, aes(x = classif.label, y = variable, size = value, color = classif.label)) + 
  geom_point(stroke = 0) + #color = '#0000ff'
  theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, hjust =1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = cluster_colors) +
  scale_size(range = c(-0.05, 8)) + 
  guides(color = 'none')

#### Difference in expression pattern between predicted HCR and test seq samples: ####
#### Do the above dotplot first! ###

test.bayes.dot.df <- cbind(mop_for_test_bin, test.bayes.classif$real.label)
names(test.bayes.dot.df)[ncol(test.bayes.dot.df)] <- 'classif.label'
test.bayes.dot.tidy <- test.bayes.dot.df %>% group_by(classif.label) %>% summarise(across(colnames(test.bayes.dot.df)[1:16], sum))
test.bayes.dot.n <- test.bayes.dot.df %>% group_by(classif.label) %>% count()
test.bayes.dot.prop <- mapply(`/`, data.frame(test.bayes.dot.tidy[,2:ncol(test.bayes.dot.tidy)]), test.bayes.dot.n[,2])
colnames(test.bayes.dot.prop) <- sub('\\.', '\\+', colnames(test.bayes.dot.prop))

test.bayes.dot.prop <- data.frame(test.bayes.dot.prop[,colnames(mop_for_test_bin[,1:(ncol(mop_for_test_bin)-1)])],
                                  id = test.bayes.dot.n[,1])
test.bayes.dot.prop <- test.bayes.dot.prop[test.bayes.dot.prop$classif.label %in% bayes.dot.prop$classif.label,]

test.bayes.dot.melt <- melt(test.bayes.dot.prop)

ggplot(bayes.dot.melt, aes(x = classif.label, y = variable, size = value, color = classif.label)) + 
  geom_point(stroke = 0) + #color = '#0000ff'
  theme(panel.background = element_rect(fill = 'white'), axis.text.x = element_text(angle = 45, hjust =1)) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values = cluster_colors) +
  geom_point(data = test.bayes.dot.melt, shape = 21, color = 'black', stroke = 1,
             aes(x = classif.label, y = variable, size = value)) +
  scale_size(range = c(-0.05, 8)) + 
  guides(color = 'none')

#### Mismatched test classification matrix: ####

all.cells <- as.character(test.bayes.classif$classif.label)
all.types <- as.character(test.bayes.classif$real.label)

all.cells <- sapply(as.character(test.bayes.classif$classif.label), function(x) ifelse(x %in% names(subclass_dict), subclass_dict[x], x))
all.types <- sapply(as.character(test.bayes.classif$real.label), function(x) ifelse(x %in% names(subclass_dict), subclass_dict[x], x))

assignments <- matrix(nrow = length(unique(all.cells)), ncol = length(unique(all.types)))
rownames(assignments) <- unique(all.cells)
colnames(assignments) <- unique(all.types)

for (i in all.types) {
  correct <- which(all.types == i)
  preds <- all.cells[correct]
  preds_table <- table(preds)
  summary <- as.numeric(preds_table/length(correct))
  assignments[names(preds_table),i] <- summary
}

assignments.df <- melt(assignments)
colnames(assignments.df) <- c('pred.label', 'real.label', 'occur')
order1 <- sub('\\d{3,4} (.*)', '\\1', unique(c(all.cells, all.types)))
order2 <- sub('(\\d{3,4}) .*', '\\1', unique(c(all.cells, all.types)))
order_custom <- order(order1, order2)
level_types <- unique(c(all.cells, all.types))[order_custom]
assignments.df$real.label <- factor(assignments.df$real.label, levels = level_types)
assignments.df$pred.label <- factor(assignments.df$pred.label, levels = level_types)
assignments.df$occur[is.na(assignments.df$occur)] <- 0

ggplot(assignments.df, aes(y = real.label, x = pred.label, fill = occur)) +
  geom_tile() +
  scale_fill_viridis(option = 'D', na.value = 'black', limits = c(0,1), oob = squish) +
  theme(axis.text.x = element_text(angle = 45, vjust =1, hjust=1))
  #geom_abline(slope = 0.8, color = 'white') +
  #geom_vline(lwd = 0.3, color = 'white', xintercept = c(2.5, 6.5, 9.5, 13.5, 16.5, 21.5, 25.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5, 36.5, 37.5)) +
  #geom_hline(lwd = 0.3, color = 'white', yintercept = c(2.5, 6.5, 9.5, 13.5, 13.5, 16.5, 21.5, 25.5, 29.5, 30.5, 31.5, 32.5, 33.5, 34.5, 35.5))

## Dendritic tree of classification accuracies for the optimized thresholds:

accuracies <- read.table('Classification for meta Seq data/Classification_accuracy_by_edge.csv', sep = ',', header = T)
accuracies$accuracy[is.na(accuracies$accuracy)] <- 0

graph_df <- data.frame(all = mop_object_hcr$all_cells,
                       nt = mop_object_hcr$allen_NT,
                       class = mop_object_hcr$allen_class,
                       subclass = mop_object_hcr$allen_subclass,
                       supertype = mop_object_hcr$my_supertype)

graph_df <- unique(graph_df)
graph_df <- graph_df[graph_df$supertype %in% accuracies$to,]
graph_df <- graph_df[order(graph_df$supertype, decreasing = F),]
rownames(graph_df) <- 1:nrow(graph_df)

edges_all_nt <- graph_df %>% select(all, nt) %>% unique %>% rename(from=all, to=nt)
edges_nt_class <- graph_df %>% select(nt, class) %>% unique %>% rename(from=nt, to=class)
edges_class_subclass <- graph_df %>% select(class, subclass) %>% unique %>% rename(from=class, to=subclass)
edges_subclass_supertype <- graph_df %>% select(subclass, supertype) %>% unique %>% rename(from=subclass, to=supertype)
edge_list <- rbind(edges_all_nt, edges_nt_class, edges_class_subclass, edges_subclass_supertype)

edge_graph <- graph_from_data_frame(edge_list)
ggraph(edge_graph, layout = 'dendrogram', circular = FALSE) + 
  geom_edge_diagonal() +
  geom_node_point(aes(size = c(1, accuracies$accuracy), color = c(1, accuracies$accuracy))) +
  geom_node_text(aes(label = name, angle = 45, hjust = -0.1, vjust = -0.05, size = 0.4)) +
  scale_color_viridis(option = 'H', na.value = 'black', limits = c(0,1), oob = squish) +
  #scale_color_gradient2(low = 'black', high = 'red', midpoint = 0.65) +
  theme_void()

## Unique expression patterns in the HCR data: ###

seq_patterns <- matrix(ncol = 16)
for (i in unique(meta_bin$my_supertype)[grepl('0090', unique(meta_bin$my_supertype))]) {
  rows_to_sample <- meta_bin[meta_bin$my_supertype == i, 1:16]
  if (!type %in% names(merfish_type_props)){
    next
  }
  else{
    cells_sampled <- sample(1:nrow(rows_to_sample), ceiling(nrow(meta_bin)*merfish_type_props[i]), replace = T)
    cells_sampled <- rows_to_sample[cells_sampled,]
    seq_patterns <- rbind(seq_patterns, as.matrix(cells_sampled)) 
  }
}

colnames(seq_patterns) <- colnames(meta_bin[,1:16])
n_rows <- nrow(seq_patterns)
shuffle_seq_patterns <- seq_patterns
seq_patterns <- as.data.frame(seq_patterns) %>% group_by_all() %>% summarise(COUNT = n())
seq_patterns <- as.data.frame(seq_patterns[order(seq_patterns$COUNT, decreasing = T),])
seq_patterns[,'cumsum'] <- cumsum(seq_patterns$COUNT)/(n_rows)
seq_patterns[,'pattern_id'] <- 1:nrow(seq_patterns)

shuffle_seq_patterns <- apply(shuffle_seq_patterns, 2, function(col){sample(col, length(col), replace = F)})
shuffle_seq_patterns <- as.data.frame(shuffle_seq_patterns)
shuffle_seq_patterns <- shuffle_seq_patterns %>% group_by_all() %>% summarise(COUNT = n())
shuffle_seq_patterns <- as.data.frame(shuffle_seq_patterns[order(shuffle_seq_patterns$COUNT, decreasing = T),])
shuffle_seq_patterns[,'cumsum'] <- cumsum(shuffle_seq_patterns$COUNT)/(n_rows)
shuffle_seq_patterns[,'pattern_id'] <- 1:nrow(shuffle_seq_patterns)

bin_patterns <- ilastik.anno[grepl('0090',bayes.classif$classif.label),1:15]
shuffle_patterns <- bin_patterns
bin_patterns <- bin_patterns %>% group_by_all() %>% summarise(COUNT = n())
bin_patterns <- as.data.frame(bin_patterns[order(bin_patterns$COUNT, decreasing = T),])
n_zeros <- bin_patterns$COUNT[1]
bin_patterns <- bin_patterns[2:nrow(bin_patterns),]
bin_patterns[,'cumsum'] <- cumsum(bin_patterns$COUNT)/(nrow(ilastik.anno) - n_zeros)
bin_patterns[,'pattern_id'] <- 1:nrow(bin_patterns)

n_zeros <- sum(rowSums(shuffle_patterns) == 0)
shuffle_patterns <- shuffle_patterns[!rowSums(shuffle_patterns) == 0,]
shuffle_patterns <- apply(shuffle_patterns, 2, function(col){sample(col, length(col), replace = F)})
shuffle_patterns <- as.data.frame(shuffle_patterns)
shuffle_patterns <- shuffle_patterns %>% group_by_all() %>% summarise(COUNT = n())
shuffle_patterns <- as.data.frame(shuffle_patterns[order(shuffle_patterns$COUNT, decreasing = T),])
shuffle_patterns[,'cumsum'] <- cumsum(shuffle_patterns$COUNT)/(sum(shuffle_patterns$COUNT))
shuffle_patterns[,'pattern_id'] <- 1:nrow(shuffle_patterns)

ggplot(bin_patterns, aes(x = pattern_id, y = cumsum)) +
  theme_minimal() +
  geom_step(color = 'red') +
  geom_step(data = seq_patterns, aes(x = pattern_id, y = cumsum), color = '#44a548') +
  geom_step(data = shuffle_patterns, aes(x = pattern_id, y = cumsum), color = 'blue') +
  geom_step(data = shuffle_seq_patterns, aes(x = pattern_id, y = cumsum), color = 'black') +
  xlim(c(0,200))

## MERFISH depth prior bar plot (or density plot):

merfish.cells <- ExtractMerfishLamina(merfish_metadata, 'my_supertype', bin = T, bin_to = 5)
type_densities <- merfish.cells[[1]]
rel_depths <- merfish.cells[[2]]
prop_depths <- merfish.cells[[3]]

prop_cast <- dcast(data = prop_depths, formula = Y ~ use_type, fun.aggregate = sum, value.var = 'prop')
prop_cast_adj <- prop_cast + 0.002
prop_melt <- melt(prop_cast_adj, id.vars = 'Y', variable.name = 'use_type', value.name = 'prop')

smooth_prop <- prop_melt %>% group_by(use_type) %>% arrange(Y) %>%
  mutate(smooth_prop = loess(prop ~ Y, span = 0.16)$fitted) %>% ungroup() %>%
  group_by(Y) %>% mutate(smooth_prop = pmax(smooth_prop, 0.002),
                         smooth_prop = smooth_prop/sum(smooth_prop))

ggplot(smooth_prop, aes(x = Y, y = smooth_prop, group = use_type, fill = use_type)) +
  geom_area(position = 'stack', stat = 'identity') +
  #geom_density(position = 'fill') +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  xlim(-1, 0) +
  coord_cartesian(ylim = c(0, 1), expand = F) +
  guides(fill = 'none')

## Plotting expression:

gene <- 'Fam84b'
max.val <- max(ilastik.rel.mean[,gene])
min.val <- min(ilastik.rel.mean[,gene])
min.max <- ceiling(abs(ilastik.rel.mean[,gene]-min.val)/(max.val-min.val)*10000)
exp.colors <- colorRampPalette(c("#ffffff", "#ff0000"))(max(min.max))
pos.w.exp <- cbind(mean.vals[,c('X', 'Y')], min.max)
pos.w.exp <- cbind(pos.w.exp, exp.colors[pos.w.exp[,'min.max']+1])
names(pos.w.exp) <- c('X', 'Y', 'val', 'col')
anno.colors <- colorRampPalette(c("#ffffff", "#0000ff"))(2)
pos.w.anno <- cbind(mean.vals[mean.vals$Z < zmax,c('X', 'Y')], ilastik.anno[mean.vals$Z < zmax,gene])
pos.w.anno <- cbind(pos.w.anno, anno.colors[pos.w.anno[,3]+1])
names(pos.w.anno) <- c('X', 'Y', 'val', 'col')

par(mar=c(5, 4, 4, 15), xpd=TRUE)
plot(pos.w.exp[,'X'], -pos.w.exp[,'Y'], 
     col = pos.w.exp[,'col'], pch = 16, cex = 1.0, xlim = c(xmin, xmax), ylim=c(-ymax,-ymin))
points(pos.w.anno[,'X'], -pos.w.anno[,'Y'], 
       col = pos.w.anno[,'col'], pch = 1, cex = 1.0, xlim = c(xmin, xmax), ylim=c(-ymax,-ymin))


## Plot cells spatially colored by certainty of classification:
round.probs <- bayes.classif$classif.prob
round.probs <- ((round.probs-min(round.probs))/(quantile(round.probs, 0.75)-min(round.probs)))
round.probs[round.probs > 1] <- 1
round.probs <- round.probs*100
prob.colors.pallete <- colorRampPalette(c('#e0e0e0', '#ff0000'))(max(round.probs))
prob.colors <- prob.colors.pallete[round(round.probs, digits = 0)]

par(mar=c(5, 4, 4, 15), xpd=TRUE)
excit.types <- grepl('L\\d', bayes.classif$classif.label)
inhib.types <- !grepl('L\\d', bayes.classif$classif.label) + grepl('ZNo', bayes.classif$classif.label)
custom.types <- grepl('L4/5 IT ', bayes.classif$classif.label)
type <- custom.types
Z.pos <- mean.vals$Z > zmin & mean.vals$Z < zmax
plot(mean.vals[type & Z.pos,'X'], -mean.vals[type & Z.pos,'Y'], 
     col = prob.colors[type & Z.pos], pch = 16, cex = 1.0, xlim = c(xmin+10, xmax-10), ylim=c(-ymax+50,-ymin+50))
legend(x = xmax+100, y = 0, cex = 0.55, title = 'Classification \nconfidence', xjust = 0,
       fill = c(prob.colors.pallete[1], 
                prob.colors.pallete[length(prob.colors.pallete)/4], 
                prob.colors.pallete[length(prob.colors.pallete)/2], 
                prob.colors.pallete[length(prob.colors.pallete)/4*3], prob.colors.pallete[length(prob.colors.pallete)]),
       legend = c(paste('0%', length(round.probs[type & round.probs >= 0 & round.probs < 12.5]), sep = ', n = '),
                  paste('25%', length(round.probs[type & round.probs >= 12.5 & round.probs < 37.5]), sep = ', n = '),
                  paste('50%', length(round.probs[type & round.probs >= 37.5 & round.probs < 62.5]), sep = ', n = '), 
                  paste('75%', length(round.probs[type & round.probs >= 62.5 & round.probs < 87.5]), sep = ', n = '), 
                  paste('100%', length(round.probs[type & round.probs >= 87.5]), sep = ', n = ')), 
       col = c(prob.colors.pallete[1], 
               prob.colors.pallete[length(prob.colors.pallete)/4], 
               prob.colors.pallete[length(prob.colors.pallete)/2], 
               prob.colors.pallete[length(prob.colors.pallete)/4*3], 
               prob.colors.pallete[length(prob.colors.pallete)]),)


## Plot confidences by cell type:
type.confidence <- data.frame(type = bayes.classif$classif.label,
                              conf = round.probs)

ggplot(type.confidence, aes(x = type, y = conf, fill = type)) +
  geom_boxplot() + scale_fill_manual(values = cluster_colors)

## Plotting the HCR heatmap:

order <- order(bayes.classif$classif.label)
#colors_list <- list(cell_type = cluster_colors[unique(bayes.classif$classif.label)])
#names(colors_list$cell_type)[is.na(names(colors_list$cell_type))] <- 'ZNo markers'
#colors_list$cell_type <- colors_list$cell_type[order(names(colors_list$cell_type))]
anno.df <- as.data.frame(cbind(ilastik.anno[,1:16], bayes.classif$classif.label, rownames(ilastik.anno)))
colnames(anno.df)[17:18] <- c('type', 'cell.id')
anno.melt <- melt(anno.df, value.name = 'value', variable.name = 'gene', id.vars = c('cell.id', 'type'))
anno.melt$type <- factor(anno.melt$type, levels = levels(bayes.classif$classif.label))
anno.melt <- anno.melt %>% arrange(factor(anno.melt$type, levels = levels(anno.melt$type)))
anno.melt$cell.id <- factor(anno.melt$cell.id, levels = unique(anno.melt$cell.id))
anno.melt$gene <- factor(anno.melt$gene, levels = rev(colnames(ilastik.anno)[1:16]))
#anno.df <- as.data.frame(cbind(ilastik.anno[order,], voltron.anno[order])) * 1
#names(anno.df)[names(anno.df) == 'voltron.anno[order]'] <- 'Voltron'
ggplot(anno.melt, aes(x=cell.id, y=gene, group=type)) +
  geom_tile(mapping=aes(fill=value)) + scale_fill_gradient(low='white', high = '#00368b') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  #geom_bar(position = 'stack', stat = 'identity', mapping = aes(x = type, y = 17, fill = type)) + 
  #scale_fill_manual(values = cluster_colors)

pheatmap(t(anno.df), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c('white', '#00368b'))(2), treeheight_col = 0,
         border_color = NA, annotation_col = data.frame(cell_type = factor(bayes.classif$classif.label[order]), row.names = colnames(t(anno.df))), 
         show_colnames = FALSE, annotation_colors = colors_list, fontsize = 7)

ordered.classif <- bayes.classif$classif.label[order]
prop.df <- anno.df
for (i in 1:length(unique(ordered.classif))) {
  now.type <- unique(ordered.classif)[i]
  type.index <- ordered.classif == now.type
  total.type <- sum(type.index)
  sum.type <- apply(prop.df[type.index,], 2, sum)
  prop.df[type.index,] <- as.list(sum.type/total.type)
}

pheatmap(t(prop.df), cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c('white', '#00368b'))(101), treeheight_col = 0,
         border_color = NA, annotation_col = data.frame(cell_type = factor(bayes.classif$classif.label[order]), row.names = colnames(t(anno.df))), 
         show_colnames = FALSE, annotation_colors = colors_list, fontsize = 7)

# The scRNA-Seq heatmap, with number of cells matched to the HCR data

sc_matched <- meta_bin[meta_bin$my_supertype %in% unique(bayes_all),]
to_compare <- table(bayes_all)

sc_compile <- matrix(nrow = length(bayes_all), ncol = 17)
count <- 1
for (i in names(to_compare)[names(to_compare) %in% unique(sc_matched$my_supertype)]) {
  which_cells_match <- which(sc_matched$my_supertype == i)
  cells_match <- as.matrix(sc_matched[sample(which_cells_match, to_compare[i], replace = T),])
  sc_compile[count:(count+to_compare[i]-1),1:17] <- cells_match[,colnames(cells_match)[1:17]]
  count <- count+as.numeric(to_compare[i])
}

colnames(sc_compile) <- colnames(sc_matched)[1:17]
sc_compile <- as.data.frame(sc_compile)
sc_compile$my_supertype <- factor(sc_compile$my_supertype, levels = levels(bayes_all)[levels(bayes_all) %in% unique(sc_compile$my_supertype)])
sc_compile[,'cell.id'] <- rownames(sc_compile)
sc_melt <- melt(sc_compile, variable.name = 'gene', value.name = 'value', id.vars = c('my_supertype', 'cell.id'), na.rm = T)
sc_melt <- sc_melt %>% arrange(factor(sc_melt$my_supertype, levels = levels(sc_melt$my_supertype)))
sc_melt$cell.id <- factor(sc_melt$cell.id, levels = unique(sc_melt$cell.id))
sc_melt$gene <- factor(sc_melt$gene, levels = rev(colnames(ilastik_all)[1:16]))
sc_melt$value <- as.numeric(sc_melt$value)

ggplot(sc_melt, aes(x=cell.id, y=gene, group=my_supertype)) +
  geom_tile(mapping=aes(fill=value)) + scale_fill_gradient(low='white', high = '#00368b') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pheatmap(heatmap_values, 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c('white', '#00368b'))(101), border_color = NA,
         annotation_col = data.frame(cell_type = sc_matched$id, row.names = rownames(sc_matched)),
         annotation_colors = list(cell_type = colors_list$cell_type[unique(sc_matched$id)]), show_colnames = FALSE)

# The two heatmaps, but now single block per type

anno.df.label <- ilastik.anno[,colnames(each.percent[,2:ncol(each.percent)])]*1
anno.df.label[,'label'] <- bayes.classif$classif.label
anno.df.group <- anno.df.label %>% group_by(label) %>% summarise(across(1:ncol(anno.df.label)-1, sum))
total.group <- table(bayes.classif$classif.label)
anno.df.prop <- anno.df.group[,2:ncol(anno.df.group)]/total.group * 100
rownames(anno.df.prop) <- names(total.group)

pheatmap(t(anno.df.prop), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c('white', '#00368b'))(101), border_color = NA,
         annotation_col = data.frame(cell_type = rownames(anno.df.prop), row.names = rownames(anno.df.prop)),
         annotation_colors = list(cell_type = colors_list$cell_type[rownames(anno.df.prop)]), show_colnames = FALSE,
         breaks = c(seq(0,100,length.out = 101)))

matched.percent <- each.percent[rownames(anno.df.prop),colnames(anno.df.prop)]*100
matched.percent <- matched.percent[1:nrow(matched.percent)-1,]

pheatmap(t(matched.percent), 
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c('white', '#00368b'))(101), border_color = NA,
         annotation_col = data.frame(cell_type = rownames(matched.percent), row.names = rownames(matched.percent)),
         annotation_colors = list(cell_type = colors_list$cell_type[rownames(matched.percent)]), show_colnames = FALSE,
         breaks = c(seq(0,100,length.out = 101)))


# Plotting a heatmap of average expression per cluster

AvgAll <- AverageExpression(alm_object_hcr, return.seurat = TRUE, group.by = 'my_type')
AvgAll@meta.data <- cbind(AvgAll@meta.data, rownames(AvgAll@meta.data))
colnames(AvgAll@meta.data)[4] <- 'my_type'

DoHeatmap(AvgAll,
          features = rownames(AvgAll),
          draw.lines = FALSE,
          group.by = 'my_type',
          slot = 'counts', 
          disp.min = 0, 
          disp.max = 0.5) + NoLegend() + scale_fill_gradientn(colors = c("black", "orange", "red")) #+ NoLegend()


### Plus responders!! ###
responders <- read.table('~/Desktop/Slice.05 Responders IoU.txt', sep = '\t', header = TRUE)
responders <- responders[responders$Ilastik.cell.ID > 0,]
voltron.classif <- bayes.classif$classif.label[voltron.anno == 1]
#names(voltron.classif) <- as.character(as.numeric(names(voltron.classif))-1e6)
voltron.sizes <- first.file$Size.in.pixels[voltron.anno == 1]
names(voltron.sizes) <- names(voltron.classif)
responders[['type']] <- voltron.classif[as.character(responders$Ilastik.cell.ID)]
responders[['size']] <- voltron.sizes[as.character(responders$Ilastik.cell.ID)]
#responders <- responders[responders$size > 100,]
responder.yes.indices <- responders$Response.SNR > 2.0 & responders$Baseline.sd < 0.065
responder.no.indices <- responders$Response.SNR <= 2.0 & responders$Response.SNR > 0 & responders$Baseline.sd < 0.065
responders.yes.ids <- responders$Ilastik.cell.ID[responder.yes.indices] #2.4
responders.no.ids <- responders$Ilastik.cell.ID[responder.no.indices] #2.4

responder.yes.types <- responders$type[responder.yes.indices]
responder.no.types <- responders$type[responder.no.indices]

responder.yes <- c()
responder.no <- c()
for (i in 1:length(unique.all)) {
  responder.yes <- c(responder.yes, sum(responder.yes.types == unique.all[i]))
  responder.no <- c(responder.no, sum(responder.no.types == unique.all[i]))
}
names(responder.yes) <- unique.all
names(responder.no) <- unique.all
responder.yes.props <- responder.yes/length(responder.yes.types)
responder.no.props <- responder.no/length(responder.no.types)

whole.props <- data.frame(voltron = c(rep('Voltron', length(voltron.props)), rep('All', length(all.props)), 
                                      rep('Responding', length(responder.yes.props)), rep('Not responding', length(responder.no.props))),
                          class = c(names(voltron.props), names(all.props), names(responder.yes.props), names(responder.no.props)),
                          prop = c(voltron.props, all.props, responder.yes.props, responder.no.props))

whole.prop.colors <- cluster_colors[whole.props$class]
names(whole.prop.colors)[is.na(whole.prop.colors)] <- 'ZNo markers'
whole.prop.colors[is.na(whole.prop.colors)] <- '#FFFFFF'
whole.props$voltron <- factor(whole.props$voltron, levels = c('All', 'Voltron', 'Responding', 'Not responding'))

ggplot(whole.props, aes(x = voltron, y = prop, fill = class)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values=whole.prop.colors)

## Test whether a particular group is over represented over chance

n_voltron <- sum(voltron.anno)
thousand_randos <- matrix(nrow = n_voltron, ncol = 1000)

for (i in 1:1000) {
  
  randos <- sample(bayes.classif$classif.label, size=n_voltron)
  thousand_randos[,i] <- randos
}
colnames(thousand_randos) <- c(1:1000)
randos_df <- melt(thousand_randos)
colnames(randos_df) <- c('cell.num', 'sample.group', 'cell.type')

summ_randos <- randos_df %>% group_by(sample.group) %>% count(cell.type)
mean_randos <- summ_randos %>% group_by(cell.type) %>% mutate(means = mean(n))
std_randos <- mean_randos %>% group_by(cell.type) %>% mutate(std = sd(n))

std_randos <- std_randos[std_randos$sample.group == '1',]
std_randos <- std_randos[std_randos$cell.type %in% unique(bayes.classif$classif.label[voltron.anno == 1]),]
std_randos[,'voltron.z'] <- apply(std_randos, 1, function(x){(sum(bayes.classif$classif.label[voltron.anno == 1] == x[['cell.type']])-as.numeric(x[['means']]))/as.numeric(x[['std']])})

sig_randos <- std_randos[abs(std_randos$voltron.z) >= 1.96,]
sig_randos

# Now for responders:
n_responders <- length(responder.yes.types)
hundred_randos <- matrix(nrow = n_responders, ncol = 1000)

for (i in 1:1000) {
  #randos <- sample(responders$type[responder.no.types], size=n_responders)
  randos <- sample(responders$type[responders$Baseline.sd < 0.065 & responders$Response.SNR > 0], size=n_responders)
  #randos <- sample(bayes.classif$classif.label[voltron.anno == 1], size=n_responders)
  hundred_randos[,i] <- randos
}
colnames(hundred_randos) <- c(1:1000)
responder_randos_df <- melt(hundred_randos)
colnames(responder_randos_df) <- c('cell.num', 'sample.group', 'cell.type')

responder_summ_randos <- responder_randos_df %>% group_by(sample.group) %>% count(cell.type)
responder_mean_randos <- responder_summ_randos %>% group_by(cell.type) %>% mutate(means = mean(n))
responder_std_randos <- responder_mean_randos %>% group_by(cell.type) %>% mutate(std = sd(n))

responder_std_randos <- responder_std_randos[responder_std_randos$sample.group == '1',]
responder_std_randos <- responder_std_randos[responder_std_randos$cell.type %in% unique.voltron,]
responder_std_randos[,'voltron.z'] <- apply(responder_std_randos, 1, function(x){(sum(responder.yes.types == x[['cell.type']])-as.numeric(x[['means']]))/as.numeric(x[['std']])})

responder_sig_randos <- responder_std_randos[abs(responder_std_randos$voltron.z) >= 1.96,]
responder_sig_randos
# Do something with sig_randos

#### Response amplitudes by cell type: #####
responders.data <- responders
responders.data <- responders[responders$Baseline.sd < 0.065 & responders$Response.SNR > 0,] #0.065
ggplot(responders.data, aes(y = Response.Amp, x = type, color = type)) + 
  geom_point(size = 3) + 
  scale_color_manual(values = cluster_colors) +
  stat_summary(geom = 'point', fun.y = 'mean', color = '#000000') +
  geom_hline(yintercept = c(0, 0.13)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

responders.mean <- responders.data %>% group_by(type) %>% mutate(means = mean(Response.Amp), std = sd(Response.Amp))
responders.mean <- responders.mean[,c('type', 'means', 'std')]
responders.mean <- unique(responders.mean)
responders.mean <- responders.mean[order(responders.mean$type),]

ggplot(responders.mean, aes(x = type, y = 1, fill = means)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradient2(low = 'white', mid = 'white', high = '#c82c2d', midpoint = 0, limits = c(min(responders$Response.Amp), max(responders$Response.Amp))) # or c(NA, NA)

#### Tau by cell type: #####
ggplot(responders.data, aes(y = Response.Tau, x = type, fill = type)) + 
  geom_bar(stat = 'summary', fun.y = 'mean') +
  #geom_errorbar(aes(ymin = mean(Response.Tau)-sd(Response.Tau), ymax = mean(Response.Tau)+sd(Response.Tau)))+
  geom_point(color = '#000000') +
  scale_fill_manual(values = cluster_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##### Response amplitude plotted spatially #####

par(mar=c(5, 4, 4, 15), xpd=TRUE)
plt.responders <- responders$Ilastik.cell.ID[responders$Baseline.sd < 0.065 & responders$FOV.number == 37]#+1e6 # Add 1e6 to the version with a separate file for JF585 annotations 
rounded.amps <- round(responders$Response.Amp*10)
positive.amps <- rounded.amps+abs(min(rounded.amps))
amp.palette <- c(rep('#ffffff', abs(min(rounded.amps))), colorRampPalette(c('#ffffff', '#c82c2d'))(max(positive.amps)-min(positive.amps)))
lgd_ <- rep(NA, length(amp.palette))
lgd_[c(1,length(amp.palette)%/%2, length(amp.palette))] <- c(min(responders$Response.Amp), median(responders$Response.Amp), max(responders$Response.Amp))
plot(mean.vals[as.character(plt.responders), 'X'], -mean.vals[as.character(plt.responders), 'Y'],
     pch = 21, cex = 1.5, lwd = 3, col = cluster_colors[bayes.classif$classif.label[as.character(plt.responders)]], 
     bg = amp.palette, 
     xlim = c(xmin+10, xmax-10), ylim=c(-ymax+50,-ymin+50))
#points(mean.vals[as.character(plt.responders), 'X'], -mean.vals[as.character(plt.responders), 'Y'],
#       pch = 21, cex = 1.5, lwd = 3, col = cluster_colors[bayes.classif$classif.label[as.character(plt.responders)]], 
#       bg = amp.palette, 
#       xlim = c(xmin+10, xmax-10), ylim=c(-ymax+50,-ymin+50))
legend(x = xmax+100, y = 0,  cex = 0.3, 
       fill = unique(cluster_colors[bayes.classif$classif.label[as.character(plt.responders)]])[order(unique(bayes.classif$classif.label[as.character(plt.responders)]))], 
       legend = unique(bayes.classif$classif.label[as.character(plt.responders)])[order(unique(bayes.classif$classif.label[as.character(plt.responders)]))], 
       col = unique(cluster_colors[bayes.classif$classif.label[as.character(plt.responders)]])[order(unique(bayes.classif$classif.label[as.character(plt.responders)]))])
legend(x = xmax+100, y = -1000, cex = 0.5,
       fill = rev(amp.palette),
       legend = rev(lgd_), y.intersp = 0.25, x.intersp = 0.4,
       border = NA)

#### Response amplitude as a function of depth? #####

y.depths <- -mean.vals[as.character(responders$Ilastik.cell.ID),'Y']
responders.y <- responders
responders.y[['Y']] <- y.depths

ggplot(responders.y, aes(x = -Y, y = Response.Amp, color = type)) + 
  geom_point() + 
  scale_color_manual(values = cluster_colors)

