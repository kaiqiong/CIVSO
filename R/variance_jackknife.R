#' Internal Engine: Generic Block Jackknife
#'
#' @description
#' Estimates the variance using the delete-d-jackknife method.
#' * If 'blocks' are provided (GLS), it performs a "Leave-One-Block-Out" jackknife,
#'   ensuring LD structure is preserved.
#' * If 'blocks' are NULL (Diagonal), it creates contiguous chunks automatically.
#'
#' @param estimator_func The engine function (.civso_engine).
#' @param n_blocks Number of jackknife blocks (default 200). Ignored if 'blocks' is provided.
#' @param total_snps Total number of SNPs (M).
#' @param ... Additional arguments (betaX, seX, blocks, etc.) passed to the engine.
#'
#' @return A list containing variance, se, and mean_est.
#' @keywords internal
.compute_jackknife_variance <- function(estimator_func, n_blocks = 200, total_snps, ...) {

  # Capture the data arguments
  args <- list(...)
  provided_blocks <- args$blocks

  # --- 1. Define Block IDs for every SNP ---

  if (!is.null(provided_blocks)) {
    # CASE A: User provided LD Blocks (GLS)
    # Use actual block sizes
    blk_sizes <- vapply(provided_blocks, function(b) length(b$betaX), numeric(1))

    # Create map: SNP -> Block ID
    block_ids <- rep(seq_along(provided_blocks), blk_sizes)
    n_blocks <- length(provided_blocks)

  } else {
    # CASE B: Diagonal Mode (No LD blocks)
    block_size <- floor(total_snps / n_blocks)
    block_ids <- rep(1:n_blocks, each = block_size)

    # Handle remainder
    if (length(block_ids) < total_snps) {
      remainder <- total_snps - length(block_ids)
      block_ids <- c(block_ids, rep(n_blocks, remainder))
    }
  }

  # --- 2. Run Jackknife Loop ---
  pseudovalues <- numeric(n_blocks)

  for (i in 1:n_blocks) {

    # A. Identify SNPs to KEEP
    idx_keep <- which(block_ids != i)

    # B. Prepare Arguments
    iter_args <- args
    iter_args$idx_subset <- idx_keep

    # C. Handle Blocks List Slicing (GLS)
    if (!is.null(provided_blocks)) {
      iter_args$blocks <- provided_blocks[-i]
    }

    # D. Call the Engine
    # Result is a LIST (beta, params, etc.)
    raw_result <- do.call(estimator_func, iter_args)

    # --- FIX: Extract only the beta value ---
    # We check if it's a list (standard engine output) or raw number
    if (is.list(raw_result)) {
      pseudovalues[i] <- raw_result$beta
    } else {
      pseudovalues[i] <- raw_result
    }
  }

  # --- 3. Calculate Variance ---
  mean_jack <- mean(pseudovalues)

  # Variance formula
  var_jack <- ((n_blocks - 1) / n_blocks) * sum((pseudovalues - mean_jack)^2)

  return(list(
    variance = var_jack,
    se = sqrt(var_jack),
    mean_est = mean_jack
  ))
}
