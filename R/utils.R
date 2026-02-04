#' Internal: Pre-calculate Trace Constants for LD Blocks
#'
#' @description
#' Enriches the 'blocks' list with scalar trace values (tr_R2, tr_R3, tr_R4)
#' required for the Analytic Variance calculation.
#' Uses the fast 'sum(A*B)' trick for symmetric matrices.
#'
#' @param blocks List of LD blocks (containing $R matrices).
#' @return The same list, but each block has $tr_R2, $tr_R3, $tr_R4 added.
#' @keywords internal
.enrich_blocks_with_traces <- function(blocks) {
  if (is.null(blocks)) return(NULL)

  # Check if first block already has traces (Assume all do if first does)
  if (!is.null(blocks[[1]]$tr_R4)) {
    return(blocks) # Already enriched, skip work.
  }

  for (i in seq_along(blocks)) {
    # Extract R
    R <- blocks[[i]]$R

    # 1. R^2 Matrix (Needed once)
    LD2 <- R %*% R

    # 2. Compute Scalar Traces (The "Trace Trick")
    # Tr(A B) = sum(A * B) for symmetric matrices
    tr_R2 <- sum(R * R)       # Tr(R^2)
    tr_R3 <- sum(LD2 * R)     # Tr(R^3) --- double sumation of l_jm and rho_jm
    tr_R4 <- sum(LD2 * LD2)   # Tr(R^4)

    # 3. Store
    blocks[[i]]$tr_R2 <- tr_R2
    blocks[[i]]$tr_R3 <- tr_R3
    blocks[[i]]$tr_R4 <- tr_R4
    blocks[[i]]$size  <- nrow(R)
  }

  return(blocks)
}
