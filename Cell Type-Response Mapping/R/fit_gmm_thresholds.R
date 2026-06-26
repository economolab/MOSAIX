# fit_gmm_thresholds.R
# GMM-based gene expression thresholding, ported from Find_binary_thresholds_GMM.py.
# Uses sklearn.mixture.GaussianMixture via reticulate to match the original Python results exactly.

library(reticulate)

# Import sklearn once at load time
sklearn_mixture <- import("sklearn.mixture")


#' Find the weighted intersection point between two Gaussians.
#' Solves the quadratic equation for the crossing point of two weighted
#' normal density functions.
#'
#' @param mu1 Mean of first Gaussian
#' @param sigma1 Standard deviation of first Gaussian
#' @param w1 Weight of first Gaussian
#' @param mu2 Mean of second Gaussian
#' @param sigma2 Standard deviation of second Gaussian
#' @param w2 Weight of second Gaussian
#' @return Intersection point between mu1 and mu2 (numeric scalar)
find_weighted_intersection <- function(mu1, sigma1, w1, mu2, sigma2, w2) {
  a <- 1 / (2 * sigma1^2) - 1 / (2 * sigma2^2)
  b <- mu2 / (sigma2^2) - mu1 / (sigma1^2)
  c_coef <- mu1^2 / (2 * sigma1^2) - mu2^2 / (2 * sigma2^2) - log((sigma2 / sigma1) * (w1 / w2))

  roots <- polyroot(c(c_coef, b, a))
  roots <- Re(roots[abs(Im(roots)) < 1e-6])

  # Return the root that lies between the two means
  valid <- roots[roots > mu1 & roots < mu2]
  if (length(valid) == 0) {
    # Fallback: return midpoint
    return((mu1 + mu2) / 2)
  }
  return(valid[1])
}


#' Classify expression values into discrete levels based on thresholds.
#'
#' @param values Numeric vector of expression values
#' @param thresholds Sorted numeric vector of thresholds
#' @return Integer vector of expression levels (0, 1, 2, ...)
classify_expression <- function(values, thresholds) {
  levels <- rep(0L, length(values))
  for (t in thresholds) {
    levels <- levels + as.integer(values > t)
  }
  return(levels)
}


#' Fit a sklearn GaussianMixture to a single gene's expression distribution.
#' Selects optimal number of components via BIC elbow (matching the original
#' Python workflow), finds thresholds at weighted Gaussian intersection points.
#'
#' @param expression_vector Numeric vector of (log-transformed) expression values
#' @param gene_name Name of the gene (for output labeling)
#' @param max_components Maximum number of GMM components to test (default 6)
#' @param bic_gain_threshold Minimum relative BIC improvement to add a component (default 0.01)
#' @param n_components Optional fixed number of components (overrides BIC selection)
#' @param means_init Optional initial means for the GMM (numeric vector).
#'   Passed to sklearn as means_init. If provided, n_components is inferred from its length.
#' @return List with: thresholds, means, stddevs, weights, n_components, binarized
fit_gene_gmm <- function(expression_vector, gene_name, max_components = 6,
                          bic_gain_threshold = 0.01, n_components = NULL,
                          means_init = NULL) {

  # Remove zero-count cells: the original pipeline replaces log(0)=-Inf with -5,
  # then excludes values <= -5 before fitting.
  expr_clean <- expression_vector[expression_vector > -5 & is.finite(expression_vector)]
  # sklearn expects a 2D column vector
  expr_2d <- array_reshape(as.numeric(expr_clean), c(length(expr_clean), 1L))

  # If means_init is provided, infer n_components from it
  if (!is.null(means_init) && is.null(n_components)) {
    n_components <- length(means_init)
  }

  # Determine optimal_n: use fixed n_components if provided, otherwise BIC selection
  if (!is.null(n_components)) {
    optimal_n <- n_components
  } else {
    # Fit GMMs for 1:max_components and compute BIC (sklearn: lower is better)
    bic_scores <- numeric(max_components)
    for (n in 1:max_components) {
      gmm <- sklearn_mixture$GaussianMixture(
        n_components = as.integer(n), random_state = 42L,
        max_iter = 500L, init_params = 'k-means++', n_init = 5L
      )
      gmm$fit(expr_2d)
      bic_scores[n] <- gmm$bic(expr_2d)
    }

    # Select optimal n via BIC elbow (same logic as original Python script)
    if (max_components > 1) {
      percent_bic_gain <- abs(diff(bic_scores)) / abs(bic_scores[1:(max_components - 1)])
      percent_bic_gain[is.nan(percent_bic_gain)] <- 0

      optimal_n <- 1
      for (n in 1:(max_components - 1)) {
        # sklearn BIC is lower-is-better, so improvement means bic decreases
        if (percent_bic_gain[n] > bic_gain_threshold && bic_scores[n + 1] < bic_scores[n]) {
          optimal_n <- n + 1
        } else {
          break
        }
      }
    } else {
      optimal_n <- 1
    }
  }

  # Fit the final model with the selected number of components
  gmm_args <- list(
    n_components = as.integer(optimal_n), random_state = 42L,
    max_iter = 500L, init_params = 'k-means++', n_init = 5L
  )
  if (!is.null(means_init)) {
    gmm_args$means_init <- array_reshape(as.numeric(means_init), c(length(means_init), 1L))
  }
  final_gmm <- do.call(sklearn_mixture$GaussianMixture, gmm_args)
  final_gmm$fit(expr_2d)

  # Extract and sort parameters
  means <- as.numeric(final_gmm$means_)
  stddevs <- as.numeric(sqrt(final_gmm$covariances_))
  weights <- as.numeric(final_gmm$weights_)

  sorted_idx <- order(means)
  means <- means[sorted_idx]
  stddevs <- stddevs[sorted_idx]
  weights <- weights[sorted_idx]

  # Find thresholds at intersection points between adjacent components
  thresholds <- numeric(0)
  if (optimal_n >= 2) {
    for (k in 1:(optimal_n - 1)) {
      thresh <- find_weighted_intersection(
        means[k], stddevs[k], weights[k],
        means[k + 1], stddevs[k + 1], weights[k + 1]
      )
      thresholds <- c(thresholds, thresh)
    }
  }

  # Classify the full expression vector
  binarized <- classify_expression(expression_vector, thresholds)

  return(list(
    thresholds = thresholds,
    means = means,
    stddevs = stddevs,
    weights = weights,
    n_components = optimal_n,
    binarized = binarized
  ))
}


#' Fit GMMs to all genes and produce combined output tables.
#' Replaces the manual per-gene Python workflow.
#'
#' @param meta_log_counts Data frame with log-transformed gene expression (genes as columns)
#' @param hcr_genes Character vector of gene names to process
#' @param means_init_list Optional named list of initial means per gene
#' @param n_components_list Optional named list of fixed n_components per gene
#'   (overrides BIC selection for those genes)
#' @param max_components Maximum components to test per gene (default 6)
#' @return List with: all_genes_bin (data frame), gmm_thresholds (data frame), per_gene_results (list)
fit_all_genes_gmm <- function(meta_log_counts, hcr_genes, means_init_list = NULL,
                               n_components_list = NULL, max_components = 6) {

  all_genes_bin <- data.frame(matrix(nrow = nrow(meta_log_counts), ncol = 0))
  threshold_rows <- list()
  per_gene_results <- list()

  for (gene in hcr_genes) {
    cat(sprintf("Fitting GMM for %s...\n", gene))

    expr <- meta_log_counts[[gene]]
    # Replace -Inf with -5 (same as original pipeline)
    expr[expr == -Inf] <- -5

    gene_means_init <- if (!is.null(means_init_list) && gene %in% names(means_init_list)) {
      means_init_list[[gene]]
    } else {
      NULL
    }

    gene_n_components <- if (!is.null(n_components_list) && gene %in% names(n_components_list)) {
      n_components_list[[gene]]
    } else {
      NULL
    }

    result <- fit_gene_gmm(expr, gene, max_components = max_components,
                            n_components = gene_n_components, means_init = gene_means_init)
    per_gene_results[[gene]] <- result

    all_genes_bin[[gene]] <- result$binarized

    # Build threshold row (pad with NA for consistent column count)
    max_thresh <- 3
    max_gauss <- 4
    thresh_padded <- c(result$thresholds, rep(NA, max_thresh - length(result$thresholds)))
    means_padded <- c(result$means, rep(NA, max_gauss - length(result$means)))
    stddevs_padded <- c(result$stddevs, rep(NA, max_gauss - length(result$stddevs)))
    weights_padded <- c(result$weights, rep(NA, max_gauss - length(result$weights)))

    threshold_rows[[gene]] <- c(thresh_padded, means_padded, stddevs_padded, weights_padded)
  }

  # Assemble GMM_Thresholds table
  gmm_thresholds <- do.call(rbind, threshold_rows)
  colnames(gmm_thresholds) <- c(
    paste0('Threshold.', 1:3),
    paste0('Mean.', 1:4),
    paste0('Stddev.', 1:4),
    paste0('Weight.', 1:4)
  )
  gmm_thresholds <- as.data.frame(gmm_thresholds)
  rownames(gmm_thresholds) <- hcr_genes

  return(list(
    all_genes_bin = all_genes_bin,
    gmm_thresholds = gmm_thresholds,
    per_gene_results = per_gene_results
  ))
}
