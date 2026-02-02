#' Internal Engine: Full-Moment GLS (Block-wise)
#'
#' @description
#' Implements Formulation 2 (GLS) for a SINGLE LD block.
#' Calculates the components required for the global summation:
#' Numerator = A' * Omega_inv * M
#' Denominator = A' * Omega_inv * A
#'
#' @param beta Vector of effect sizes for the block (length M_b).
#' @param se Vector of standard errors for the block.
#' @param R The LD correlation matrix for the block (M_b x M_b).
#' @param n_total Total number of SNPs in the FULL genome (M_global).
#' @param n_samp Sample size (N).
#' @param h_init Initial heritability estimate (from Diagonal WLS Step 1).
#' @param v_init Initial intercept estimate (from Diagonal WLS Step 1).
#'
#' @return A list containing the numerator and denominator scalars for this block.
#' @keywords internal
.estimate_block_gls_h2 <- function(beta, se, R, n_total, n_samp, h_init, v_init) {

  M_b <- length(beta)


  # --- 1. Construct Structural Matrix (Kappa) ---
  # Raw Bivariate LD Score: l_jm^(2) = sum(r_jk * r_mk)
  # This represents the expected correlation of effects through LD
  LD2_mat <- R %*% R

  # Full Kappa Definition (matching paper derivation):
  # kappa_jm = (1 + 1/N) * l_jm^(2) + (M/N) * r_jm
  term_scaling <- 1 + 1/n_samp
  term_noise   <- n_total / n_samp

  Kappa_mat <- term_scaling * LD2_mat + term_noise * R

  # --- 2. Construct Regressors (A) ---
  # A_jm = Kappa_jm - (M-1)/N * r_jm
  # Note: This simplifies algebraically!
  # If we plug Kappa into A:
  # A = [(1+1/N)L^2 + (M/N)R] - [(M-1)/N]R
  # A = (1+1/N)L^2 + (1/N)R

  # --- 2. Construct Regressors (A) ---
  # Formula: A_jk = Kappa_jk - (M_global - 1)/N * rho_jk
  # Note: The structural constant A must be vectorized to match the moments
  correction_factor <- (n_total - 1) / n_samp
  A_mat <- Kappa_mat - correction_factor * R

  # Vectorize lower triangle (including diagonal)
  # We use lower.tri(diag=TRUE) because the matrix is symmetric
  # but GLS usually uses the full vector or a half-vectorization.
  # Using full vectorization (M*M) is easier for Kronecker math.
  a_vec <- as.vector(A_mat)

  # --- 3. Construct Observed Moments (M) ---
  # Raw cross-products: gamma_j * gamma_k
  Gamma_cross <- tcrossprod(beta) # beta %*% t(beta)

  # Bias correction matrix
  # Diagonal: sigma_j^2
  # Off-diagonal: rho_jk * (sigma_j^2 + sigma_k^2)/2
  Sigma_sq <- se^2

  # OPTIMIZATION: Use outer() instead of a loop for speed
  # This creates the matrix of (sigma_j^2 + sigma_k^2) / 2
  Avg_Sigma_Mat <- outer(Sigma_sq, Sigma_sq, "+") / 2

  # Multiply by R to get the rho-scaled bias
  # Note: On diagonal, R_jj=1 and (s^2+s^2)/2 = s^2, so this works perfectly for diag too.
  Bias_mat <- R * Avg_Sigma_Mat


  M_mat <- Gamma_cross - Bias_mat
  m_vec <- as.vector(M_mat)

  # --- 4. Construct Weighting Matrix (Omega) ---
  # Sigma_jk = h_init * Kappa_jk + v_init * R_jk
  Sigma_marginal <- h_init * Kappa_mat + v_init * R

  # Kronecker Product: Sigma (x) Sigma
  # Size becomes M^2 x M^2. WARNING: Grows fast!
  # M=100 -> 10,000 x 10,000 matrix (800MB RAM)
  Kron_SS <- kronecker(Sigma_marginal, Sigma_marginal)

  # Commutation Matrix part: (Sigma (x) Sigma) * K
  # For symmetric Sigma, Omega = (I + K) * (Sigma (x) Sigma)
  # But constructing K is expensive.
  # Shortcut: For symmetric matrices A and B, K(A (x) B) = (B (x) A)K?
  # The formula is Omega = Sigma (x) Sigma + (Sigma (x) Sigma) K

  # We construct K explicitly (Commutation Matrix)
  # K_mn maps vec(X) -> vec(X')
  # For M=100, this is feasible.
  K_perm <- .commutation_matrix(M_b)

  # Omega = Kron + Kron %*% K
  # Note: This is equivalent to summing the Parallel and Crossed terms
  Omega <- Kron_SS + Kron_SS[, K_perm]

  # --- 5. Calculate Components ---
  # We need: inv(Omega)
  # Use Cholesky or SVD for stability if possible, but solve() is standard
  Omega_inv <- tryCatch(solve(Omega), error = function(e) MASS::ginv(Omega))

  # Numerator: a' W m
  num_part <- as.numeric(t(a_vec) %*% Omega_inv %*% m_vec)

  # Denominator: a' W a
  den_part <- as.numeric(t(a_vec) %*% Omega_inv %*% a_vec)

  return(list(
    num = num_part,
    den = den_part,
    M_dim = M_b # For checking
  ))
}


#' Internal Engine: Full-Moment GLS Cross-Trait (Block-wise)
#'
#' @keywords internal
.estimate_block_gls_hxy <- function(betaX, betaY, R,
                                    n_total, n_x, n_y, overlap_prop,
                                    covXY_theory,
                                    res_XX_init,# <--- Output from Diagonal WLS (Trait 1) -- list(slope=, incpt= )
                                    res_YY_init,# <--- Output from Diagonal WLS (Trait 2)
                                    h_xy_init, # <--- Output from Diagonal WLS (Cross)
                                    v_xy_init # <--- Output from Diagonal WLS (Cross)
                                    ) {

  M_b <- length(betaX)


  # --- 1. Construct Constants ---
  # Raw Bivariate LD Score
  LD2_mat <- R %*% R


  # A. Kappa for X (Uses n_x)
  # kappa_X = (1 + 1/Nx)*L^2 + (M/Nx)*R
  k_scale_X <- 1 + 1/n_x
  k_noise_X <- n_total / n_x
  Kappa_X   <- k_scale_X * LD2_mat + k_noise_X * R

  # B. Kappa for Y (Uses n_y)
  # kappa_Y = (1 + 1/Ny)*L^2 + (M/Ny)*R
  k_scale_Y <- 1 + 1/n_y
  k_noise_Y <- n_total / n_y
  Kappa_Y   <- k_scale_Y * LD2_mat + k_noise_Y * R

  # C. Xi for Cross-Trait (Uses Overlap & Ny)
  # xi = (1 + No/Ny)*L^2 + (M*No/Ny)*R
  term_N    <- (n_total * overlap_prop) / n_y
  xi_scale  <- 1 + overlap_prop/n_y
  Xi_mat    <- xi_scale * LD2_mat + term_N * R


  # --- 3. Construct Regressors (A) and Moments (M) ---
  # A_xy = Xi - term_N * rho_jk
  # Note: The term_N * R cancels out exactly here!
  # A_xy becomes just: term_scaling * LD2_mat
  # However, keeping the full subtraction is safer for code readability/debugging.
  A_mat <- Xi_mat - term_N * R
  a_vec <- as.vector(A_mat)



  # --- 2. Construct Moments ---
  # Gamma_j * Gamma_k
  Cross_prod <- tcrossprod(betaX, betaY) # betaX %*% t(betaY)

  # Bias correction: rho_jk * (No/NxNy * sigma_XY)
  noise_scalar <- (overlap_prop / n_y) * covXY_theory
  Bias_mat <- R * noise_scalar

  M_mat <- Cross_prod - Bias_mat
  m_vec <- as.vector(M_mat)

  # --- 3. Construct Weighting Matrix (Omega) ---
  # Need Marginal Sigmas from the Init inputs
  # We assume h_init and v_init are scalars passed in

  # Reconstruct Sigma XX and YY using Init estimates
  # Sigma = h * K + v * R
  Sigma_XX <- res_XX_init$slope * Kappa_X + res_XX_init$incpt * R
  Sigma_YY <- res_YY_init$slope * Kappa_Y + res_YY_init$incpt * R

  # Sigma XY
  Sigma_XY <- h_xy_init * Xi_mat + v_xy_init * R

  # Kronecker Terms
  # Term 1: Sigma_XX (x) Sigma_YY
  Kron_Main <- kronecker( Sigma_YY, Sigma_XX)

  # Term 2: (Sigma_XY (x) Sigma_XY) * K
  Kron_Cross <- kronecker(Sigma_XY, Sigma_XY)
  K_perm <- .commutation_matrix(M_b)

  # Omega = Main + Cross * K
  Omega <- Kron_Main + Kron_Cross[, K_perm]

  # --- 4. Solve ---
  invOmega <- tryCatch({
    solve(Omega)
  }, error = function(e) {
    # Fallback for singular matrices (rare in block GLS but possible)
    MASS::ginv(Omega)
  })

  num_part <- as.numeric(t(a_vec) %*% Omega_inv %*% m_vec)
  den_part <- as.numeric(t(a_vec) %*% Omega_inv %*% a_vec)

  return(list(num = num_part, den = den_part))
}
