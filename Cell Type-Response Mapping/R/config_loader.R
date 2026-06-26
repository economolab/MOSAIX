# config_loader.R
# Loads experiment configurations from JSON files.

library(jsonlite)

#' Get the project root directory (where config/ lives).
#' @return Absolute path to the MOSAIX pipeline root
get_project_root <- function() {
  # If sourced from within R/, go up one level
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  if (basename(script_dir) == 'R') {
    return(dirname(script_dir))
  }
  return(script_dir)
}

#' Load shared pipeline paths from config/paths.json.
#' Resolves relative paths against the project root.
#'
#' @param project_root Optional override for project root directory
#' @return Named list of resolved directory paths
load_paths <- function(project_root = NULL) {
  if (is.null(project_root)) project_root <- get_project_root()

  paths <- fromJSON(file.path(project_root, 'config', 'paths.json'))

  # Expand ~ in paths
  paths <- lapply(paths, path.expand)

  # Resolve relative paths against project root
  for (name in names(paths)) {
    if (!grepl('^[/~]', paths[[name]])) {
      paths[[name]] <- file.path(project_root, paths[[name]])
    }
  }

  return(paths)
}

#' Load an experiment configuration from a JSON file.
#' Converts angle_deg to radians and builds the file listing.
#'
#' @param config_path Path to experiment JSON file (absolute or relative to config/experiments/)
#' @param project_root Optional override for project root directory
#' @return Named list with all experiment parameters
load_experiment_config <- function(config_path, project_root = NULL) {
  if (is.null(project_root)) project_root <- get_project_root()
  paths <- load_paths(project_root)

  # Allow short names like "2024-07-30_thalamus.json"
  if (!file.exists(config_path)) {
    config_path <- file.path(project_root, 'config', 'experiments', config_path)
  }

  config <- fromJSON(config_path)

  # Convert angles from degrees to radians
  config$img_lims$angle <- config$img_lims$angle_deg * pi / 180

  # Build the annotation directory path
  config$anno_dir <- file.path(paths$mosaix_data_dir, config$data_subdir)

  # Build the output directory
  config$output_dir <- file.path(paths$output_dir, config$experiment_id)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build the file listing based on enumeration type
  if (config$slice_enumeration == 'flat') {
    config$reg_file_list <- list.files(config$anno_dir, full.names = TRUE)
    config$gene_file <- sub(config$gene_file_regex, '\\1', config$reg_file_list)
  } else {
    # Folder-based: enumerate Slice directories
    if (is.null(config$n_slices)) {
      config$n_slices <- sum(grepl('Slice ', list.files(config$anno_dir)))
    }
    config$reg_file_list <- c()
    for (i in 1:config$n_slices) {
      cur_dir <- paste0(config$anno_dir, 'Slice ', as.character(i + config$slice_offset), '/Gene expression tables/')
      config$reg_file_list <- c(config$reg_file_list, list.files(cur_dir, pattern = 'table', full.names = TRUE))
    }
    config$gene_file <- sub(config$gene_file_regex, '\\1', config$reg_file_list)
  }

  # Convert genes_combined from list-of-vectors to proper R list
  if (is.matrix(config$genes_combined)) {
    gc_list <- list()
    for (i in 1:nrow(config$genes_combined)) {
      gc_list[[i]] <- config$genes_combined[i,]
    }
    config$genes_combined <- gc_list
  } else if (is.character(config$genes_combined) && length(config$genes_combined) == 2) {
    config$genes_combined <- list(config$genes_combined)
  } else if (length(config$genes_combined) == 0) {
    config$genes_combined <- list()
  }

  # Attach shared paths
  config$paths <- paths

  return(config)
}
