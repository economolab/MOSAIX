# hcr_data_loader.R
# Functions for reading HCR imaging data and performing spatial corrections.

#' Load HCR expression data from Ilastik-generated CSV files.
#' Reads per-gene expression tables, extracts cell positions, gene intensities,
#' and Ilastik binary labels. Performs spatial rotation for laminar alignment.
#'
#' @param config Experiment configuration list from load_experiment_config()
#' @return List with mean.vals (full data frame) and ilastik.anno (binarized expression matrix with newy)
load_hcr_data <- function(config) {

  hcr_genes <- config$hcr_genes
  genes_combined <- config$genes_combined
  file_colnames <- config$file_colnames
  reg_file_list <- config$reg_file_list
  gene_file <- config$gene_file
  n_slices <- config$n_slices
  img_lims <- config$img_lims

  gene_file_len <- length(hcr_genes) - length(genes_combined)

  mean.vals <- data.frame(row.names = c(), slice = c(), X = c(), Y = c(), Z = c(), volt = c())

  print('Compiling Ilastik annotations...')
  pb <- txtProgressBar(min = 0, max = length(reg_file_list), initial = 0, style = 3)

  for (j in 1:n_slices) {

    name <- reg_file_list[j * gene_file_len - (gene_file_len - 1)]
    current.file <- read.csv(name, sep = ',', header = TRUE)
    cur.mean.vals <- data.frame(
      cell.id = as.character(current.file[, file_colnames[['cell_id']]]),
      slice = as.numeric(sub('s(\\d{2})_.*', '\\1', gene_file[j * gene_file_len])),
      X = current.file[, file_colnames[['X']]],
      Y = current.file[, file_colnames[['Y']]],
      Z = current.file[, file_colnames[['Z']]],
      volt = rep(FALSE, nrow(current.file))
    )

    for (i in 1:gene_file_len) {
      idx <- i + (j * gene_file_len - gene_file_len)
      name <- reg_file_list[idx]
      current.file <- read.csv(name, sep = ',', header = TRUE)

      if ('Label 1' %in% current.file[, file_colnames[['Label']]]) {
        user_labels <- (current.file[, file_colnames[['Label']]] == 'Label 1') * 1
      } else {
        user_labels <- current.file[, file_colnames[['Label']]]
      }

      gene.name <- sub('s\\d{2}_(.*)', '\\1', gene_file[idx])
      cur.mean.vals[[paste0(gene.name, '.mean.intensity')]] <- current.file[, file_colnames[['Roi_mean']]]
      cur.mean.vals[[paste0(gene.name, '.mean.neighborhood.intensity')]] <- current.file[, file_colnames[['Neighbor_mean']]]
      cur.mean.vals[[paste0(gene.name, '.object.classification')]] <- user_labels
      setTxtProgressBar(pb, idx)
    }
    mean.vals <- rbind(mean.vals, cur.mean.vals)
  }
  close(pb)

  # Rotate cell positions for laminar structure alignment
  new.x <- c()
  new.y <- c()

  for (i in 1:n_slices) {
    slice.oi <- unique(mean.vals$slice)[i]
    rows.oi <- mean.vals[mean.vals$slice == slice.oi, ]
    origin.x <- mean(rows.oi$X)
    origin.y <- mean(rows.oi$Y)

    cur.new.x <- (rows.oi$X - origin.x) * cos(img_lims$angle[i]) - (rows.oi$Y - origin.y) * sin(img_lims$angle[i]) + origin.x
    cur.new.y <- (rows.oi$X - origin.x) * sin(img_lims$angle[i]) + (rows.oi$Y - origin.y) * cos(img_lims$angle[i]) + origin.y

    norm.new.x <- (cur.new.x - min(cur.new.x)) / (max(cur.new.x) - min(cur.new.x))
    norm.new.y <- (cur.new.y - img_lims$ymin[i]) / (img_lims$ymax[i] - img_lims$ymin[i])

    new.x <- c(new.x, norm.new.x)
    new.y <- c(new.y, norm.new.y)
  }

  mean.vals$newx <- new.x
  mean.vals$newy <- -new.y

  # Extract Ilastik binary annotations into classification-ready matrix
  ilastik.anno <- mean.vals[, grepl('.*\\.object\\.classification', colnames(mean.vals))]
  colnames(ilastik.anno) <- unique(sub('s\\d{2}_(.*)', '\\1', gene_file))

  # Handle JF585/Voltron channel separately if present
  voltron.anno <- NULL
  if ('JF585' %in% colnames(ilastik.anno)) {
    voltron.anno <- ilastik.anno[, 'JF585']
    names(voltron.anno) <- rownames(ilastik.anno)
    ilastik.anno <- ilastik.anno[, colnames(ilastik.anno) != 'JF585']
  } else if (TRUE %in% mean.vals$volt) {
    voltron.anno <- mean.vals$volt
    names(voltron.anno) <- rownames(ilastik.anno)
  }

  ilastik.anno$newy <- mean.vals$newy

  return(list(
    mean.vals = mean.vals,
    ilastik.anno = ilastik.anno,
    voltron.anno = voltron.anno
  ))
}
