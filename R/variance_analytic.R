#' Internal Engine: Analytic Variance (Optimized)
#'
#' @description
#' Computes the asymptotic variance using pre-calculated trace scalars.
#' Complexity: O(1) relative to matrix size (instant).
#'
#' @keywords internal
.compute_analytic_variance <- function(blocks,
                                       beta_hat,
                                       h_xx, v_x,
                                       h_yy, v_y,
                                       h_xy, v_xy,
                                       n_snp, n_x, n_y, overlap_prop) {

  # 1. Strict Requirement
  if (is.null(blocks)) return(NA)

  # 2. Aggregate Traces (Sum across all blocks)
  Sum_R2 <- 0; Sum_R3 <- 0; Sum_R4 <- 0; Sum_M  <- 0

  for(b in seq_along(blocks)) {
    # Note: Blocks must be enriched via .enrich_blocks_with_traces() before calling this
    Sum_R2 <- Sum_R2 + blocks[[b]]$tr_R2
    Sum_R3 <- Sum_R3 + blocks[[b]]$tr_R3
    Sum_R4 <- Sum_R4 + blocks[[b]]$tr_R4
    Sum_M  <- Sum_M  + blocks[[b]]$size
  }

  # 3. Define Structural Coefficients (Scalars)
  # Constants
  c1_x <- 1 + 1/n_x;  c2_x <- n_snp/n_x
  c1_y <- 1 + 1/n_y;  c2_y <- n_snp/n_y

  d1 <- 1 + overlap_prop/n_y
  d2 <- (n_snp * overlap_prop) / n_y

  # Structural Params (A * R^2 + B * R)
  A_x <- h_xx * c1_x;  B_x <- h_xx * c2_x + v_x
  A_y <- h_yy * c1_y;  B_y <- h_yy * c2_y + v_y
  A_xy <- h_xy * d1;   B_xy <- h_xy * d2 + v_xy

  # 4. Global Variance Calculation

  # Denominator Trace: Tr(Xi)
  sum_xi_trace <- d1 * Sum_R2 + d2 * Sum_M

  # Numerator Components
  # Helper to compute Tr((A1*R^2 + B1*R)(A2*R^2 + B2*R))
  calc_trace_prod <- function(A1, B1, A2, B2) {
    term4 <- A1 * A2 * Sum_R4
    term3 <- (A1 * B2 + A2 * B1) * Sum_R3
    term2 <- B1 * B2 * Sum_R2
    return(term4 + term3 + term2)
  }

  T1 <- calc_trace_prod(A_x, B_x, A_y, B_y)     # Tr(S_XX * S_YY)
  T2 <- calc_trace_prod(A_xy, B_xy, A_xy, B_xy) # Tr(S_XY * S_XY)
  T3 <- calc_trace_prod(A_x, B_x, A_x, B_x)     # Tr(S_XX * S_XX)
  T4 <- calc_trace_prod(A_x, B_x, A_xy, B_xy)   # Tr(S_XX * S_XY)

  b2 <- beta_hat^2
  sum_variance_terms <- T1 + T2 + 2 * b2 * T3 - 4 * beta_hat * T4

  # 5. Final Scaling
  if (abs(h_xx) < 1e-8 || abs(sum_xi_trace) < 1e-8) return(NA)

  mu_D <- h_xx * sum_xi_trace
  final_var <- sum_variance_terms / (mu_D^2)

  if (final_var < 0) return(NA)
  return(final_var)
}
