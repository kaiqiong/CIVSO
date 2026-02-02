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
    # CASE A: User provided LD Blocks (GLS / HAPNEST)
    # We must respect these boundaries exactly.

    # Calculate the size of each provided block
    # We assume the blocks cover the data in order.
    blk_sizes <- vapply(provided_blocks, function(b) length(b$betaX), numeric(1))

    # Create a map: SNP 1..N -> Block ID
    # e.g. [1, 1, 1, 2, 2, 3, 3, 3...]
    block_ids <- rep(seq_along(provided_blocks), blk_sizes)

    # Update n_blocks to match the actual number of LD blocks
    n_blocks <- length(provided_blocks)

  } else {
    # CASE B: No LD Blocks (Diagonal Mode)
    # Create arbitrary contiguous chunks
    block_size <- floor(total_snps / n_blocks)
    block_ids <- rep(1:n_blocks, each = block_size)

    # Handle remainder
    if (length(block_ids) < total_snps) {
      remainder <- total_snps - length(block_ids)
      block_ids <- c(block_ids, rep(n_blocks, remainder))
    }
  }

  # Double check length (Safety)
  if (length(block_ids) != total_snps) {
    # This happens if 'blocks' don't sum up to 'total_snps' (e.g., missing SNPs dropped)
    # In that case, we trust the 'block_ids' derived from the blocks and ignore total_snps mismatch warning
    # provided total_snps was just passed as metadata.
  }

  # --- 2. Run Jackknife Loop ---
  pseudovalues <- numeric(n_blocks)

  for (i in 1:n_blocks) {

    # A. Identify SNPs to KEEP (Everything except block i)
    idx_keep <- which(block_ids != i)

    # B. Prepare Arguments for this iteration
    # We copy the original arguments
    iter_args <- args

    # Update 'idx_subset' in the arguments
    iter_args$idx_subset <- idx_keep

    # C. Handle Blocks List Slicing (CRITICAL for GLS)
    # If we are doing GLS, we must also remove the block from the list
    # so the engine doesn't try to loop over a block whose data we just removed.
    if (!is.null(provided_blocks)) {
      iter_args$blocks <- provided_blocks[-i]
    }

    # D. Call the Engine
    # do.call is robust way to pass a list of arguments
    pseudovalues[i] <- do.call(estimator_func, iter_args)
  }

  # --- 3. Calculate Variance ---
  mean_jack <- mean(pseudovalues)

  # Standard Jackknife Variance Formula
  var_jack <- ((n_blocks - 1) / n_blocks) * sum((pseudovalues - mean_jack)^2)

  return(list(
    variance = var_jack,
    se = sqrt(var_jack),
    mean_est = mean_jack
  ))
}
