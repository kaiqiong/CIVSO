#' Internal Engine: Analytic Variance (Delta Method)
#'
#' @description
#' Computes the asymptotic variance of the CIVSO estimator using the
#' closed-form expression (Equation \ref{varexp}) derived via the Delta Method.
#'
#' V_CIVSO = (1 / mu_D^2) * Sum [ sigma_X*sigma_Y + sigma_XY^2 + 2*beta*sigma_X*(beta*sigma_X - 2*sigma_XY) ]
#'
#' @keywords internal
.compute_analytic_variance <- function(blocks,
                                       beta_hat,
                                       h_xx, v_x,
                                       h_yy, v_y,
                                       h_xy, v_xy,
                                       n_snp, n_x, n_y, overlap_prop) {

  # Accumulators
  sum_variance_terms <- 0
  sum_xi_trace <- 0      # Needed for theoretical denominator (mu_D)

  b2 <- beta_hat^2

  for(b in 1:length(blocks)) {
    blk <- blocks[[b]]
    R <- blk$R

    # --- 1. Construct Sigmas (Covariance Structures) ---
    LD2_mat <- R %*% R

    # Sigma XX (Covariance of Gamma_hat)
    k_scale_X <- 1 + 1/n_x
    k_noise_X <- n_snp / n_x
    Kappa_X   <- k_scale_X * LD2_mat + k_noise_X * R
    Sigma_XX  <- h_xx * Kappa_X + v_x * R

    # Sigma YY (Covariance of Gamma_hat_outcome)
    k_scale_Y <- 1 + 1/n_y
    k_noise_Y <- n_snp / n_y
    Kappa_Y   <- k_scale_Y * LD2_mat + k_noise_Y * R
    Sigma_YY  <- h_yy * Kappa_Y + v_y * R

    # Sigma XY (Cross-Covariance)
    # Xi Matrix construction matches geometric constant xi_j
    term_N    <- (n_snp * overlap_prop) / n_y
    xi_scale  <- 1 + overlap_prop/n_y
    Xi_mat    <- xi_scale * LD2_mat + term_N * R

    Sigma_XY  <- h_xy * Xi_mat + v_xy * R

    # --- 2. Calculate Theoretical Expectation of Denominator ---
    # Text: mu_D = h_xx * Sum(xi_j)
    # The diagonal of Xi_mat contains the vector xi_j for this block
    sum_xi_trace <- sum_xi_trace + sum(diag(Xi_mat))

    # --- 3. Calculate Variance Components (Traces) ---
    # We use sum(A * B) which is equivalent to Tr(AB) for symmetric matrices
    # Corresponds to the double summation Sum_{j,m}

    # (i) Var(A) components
    # T1 = Sum(sigma_X,jm * sigma_Y,jm)
    T1 <- sum(Sigma_XX * Sigma_YY)
    # T2 = Sum(sigma_XY,jm^2)
    T2 <- sum(Sigma_XY * Sigma_XY)

    # (ii) Var(D) component
    # T3 = Sum(sigma_X,jm^2) -> Var(D) = 2 * T3
    T3 <- sum(Sigma_XX * Sigma_XX)

    # (iii) Cov(A,D) component
    # T4 = Sum(sigma_X,jm * sigma_XY,jm) -> Cov(A,D) = 2 * T4
    T4 <- sum(Sigma_XX * Sigma_XY)

    # --- 4. Combine Block Contribution ---
    # Formula: Var(A) + beta^2*Var(D) - 2*beta*Cov(A,D)
    # Var(A) = T1 + T2
    # Var(D) = 2 * T3
    # Cov(A,D) = 2 * T4
    #
    # Expansion: (T1 + T2) + beta^2*(2*T3) - 2*beta*(2*T4)
    #          = T1 + T2 + 2*beta^2*T3 - 4*beta*T4

    block_val <- T1 + T2 + 2 * b2 * T3 - 4 * beta_hat * T4

    sum_variance_terms <- sum_variance_terms + block_val
  }

  # --- 5. Final Scaling (Equation 13) ---
  # Denominator is mu_D^2 = (h_xx * Sum_xi)^2
  mu_D <- h_xx * sum_xi_trace

  final_var <- sum_variance_terms / (mu_D^2)

  return(final_var)
}
