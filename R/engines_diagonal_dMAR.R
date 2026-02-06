#' Internal Engine: Estimate Single-Trait Heritability Structure (Optimized)
#'
#' @description
#' Performs Two-Step WLS to estimate heritability from summary statistics.
#' Uses vectorized linear algebra (solving for slope through origin)
#' instead of lm() for maximum computational speed.
#'
#' @param beta Vector of GWAS effect sizes
#' @param se Vector of standard errors
#' @param ld_score Vector of LD scores
#' @param n_snp Total number of SNPs (M)
#' @param n_samp Sample size (N)
#'
#' @return A list containing the slope (h2), intercept, and total expected squared moments.
#' @keywords internal
.estimate_h2_struct <- function(beta, se, ld_score, n_snp, n_samp, idx_subset = NULL) {


  # --- 0. Handle Subsetting ---
  if (!is.null(idx_subset)) {
    beta <- beta[idx_subset]
    se <- se[idx_subset]
    ld_score <- ld_score[idx_subset]
    # Note: n_snp usually stays as the global M for the LD score correction logic,
    # but some methods adjust it. Standard LDSC keeps M fixed.
  }



  # --- 1. Construct Variables ---
  # Constants (Predictor X)
  # kappa = (1 + 1/N)*l_j + M/N
  kappa <- (1 + 1/n_samp) * ld_score + n_snp/n_samp
  A_xx  <- kappa - (n_snp - 1)/n_samp

  # Response (Y)
  M_xx <- beta^2 - se^2

  # --- 2. Step 1: Initial Estimates (A robust starting guess for
  #      WLS with weights as ld_score-squared) ---
  # Weight: 1/l^2 (Heuristic for strong signal)
  w_step1 <- 1 / pmax(ld_score, 1)^2

  # Fast Slope (Step 1)
  # Formula: sum(wxy) / sum(wx^2)
  num_init <- sum(w_step1 * M_xx * A_xx)
  den_init <- sum(w_step1 * A_xx^2)
  h_init   <- num_init / den_init

  # Fast Intercept (Step 1)
  # v_init = Mean_Y_weighted - Slope * Mean_X_weighted ---- This is wrong
  #sum_w    <- sum(w_step1)
  #mean_y_w <- sum(w_step1 * M_xx) / sum_w
  #mean_x_w <- sum(w_step1 * A_xx) / sum_w
  #v_init   <- mean_y_w - h_init * mean_x_w


  # Fast Intercept (Step 1) - STRUCTURAL BACK-CALCULATION
  # The regression intercept is NOT free. It is constrained by the Total Moment.
  # v = Mean(SE^2) - h * (M-1)/N
  mean_se2 <- mean(se^2)
  term_scaling <- (n_snp - 1) / n_samp

  v_init <- mean_se2 - h_init * term_scaling
  # --- 3. Step 2: Final WLS Weights ---
  # We use the Step 1 estimates to predict the variance structure
  # E[gamma^2] = h * kappa + v
  E_sq_moment <- abs(h_init * kappa + v_init)

  # Weight = 1 / (L_j * Predicted_Variance^2)
  w_final <- 1 / (pmax(ld_score, 1) * E_sq_moment^2)

  # --- 4. Step 3: Final Slope (The "WLS" Step) ---
  num_final <- sum(w_final * M_xx * A_xx)
  den_final <- sum(w_final * A_xx^2)
  h_final   <- num_final / den_final

  # --- 5. Back-Calculate Intercept (Constrained) ---
  # Instead of regression intercept, we use the constraint:
  # E[sigma^2] = (M-1)/N * h_xx + v_x
  # Average(sigma^2) is roughly Mean(se^2) across the genome

  # To be precise with the summary stats provided:
  # We assume total phenotypic variance is approx 1 (normalized).
  # If not, Mean(se^2) is the best proxy for E[sigma^2] in summary data.
  # Constraint: Mean(se^2) = ( (M-1)/N ) * h_final + v_final



  incpt_final <- mean_se2 - h_final * term_scaling

  return(list(
    slope = h_final,
    incpt = incpt_final,
    E_total_sq = E_sq_moment, # Needed for cross-trait weighting
    kappa = kappa
  ))
}



#' Internal Engine: Simple OLS Heritability (Unweighted)
#' @description Performs unweighted regression (w=1) for single-trait heritability.
#' Used as a baseline.
#' @keywords internal
.estimate_h2_ols <- function(beta, se, ld_score, n_snp, n_samp, idx_subset = NULL) {

  # --- 0. Handle Subsetting ---
  if (!is.null(idx_subset)) {
    beta <- beta[idx_subset]
    se <- se[idx_subset]
    ld_score <- ld_score[idx_subset]
  }

  # --- 1. Construct Variables (Same as WLS) ---
  kappa <- (1 + 1/n_samp) * ld_score + n_snp/n_samp
  A_xx  <- kappa - (n_snp - 1)/n_samp

  M_xx  <- beta^2 - se^2

  # --- 2. Calculate Slope (Direct OLS) ---
  # Formula: sum(xy) / sum(x^2)  <-- Implicitly w=1
  num <- sum(M_xx * A_xx)
  den <- sum(A_xx^2)
  h_est <- num / den

  # --- 3. Calculate Intercept ---

  # --- 5. Back-Calculate Intercept (Constrained) ---
  # Instead of regression intercept, we use the constraint:
  # E[sigma^2] = (M-1)/N * h_xx + v_x
  # Average(sigma^2) is roughly Mean(se^2) across the genome

  # To be precise with the summary stats provided:
  # We assume total phenotypic variance is approx 1 (normalized).
  # If not, Mean(se^2) is the best proxy for E[sigma^2] in summary data.
  # Constraint: Mean(se^2) = ( (M-1)/N ) * h_final + v_final

  mean_se2 <- mean(se^2)
  term_scaling <- (n_snp - 1) / n_samp

  incpt <- mean_se2 -   h_est * term_scaling

  return(list(
    slope = h_est,
    incpt = incpt
  ))
}

#' Internal Engine: Estimate Cross-Trait Covariance Structure (Optimized)
#'
#' @description
#' Performs Two-Step WLS to estimate genetic covariance.
#' Requires output from the single-trait engine (.estimate_h2_struct)
#' to calculate the correct Isserlis weights.
#'
#' @param betaX Vector of GWAS effect sizes for Trait 1.
#' @param betaY Vector of GWAS effect sizes for Trait 2.
#' @param res_XX List output from .estimate_h2_struct for Trait 1. Must contain element `$E_total_sq`, i.e the expected squared moments
#' @param res_YY List output from .estimate_h2_struct for Trait 2. Contains expected squared moments ($E_Gamma_j_sq) needed for weighting.
#' @param ld_score Vector of LD scores (l_j).
#' @param covXY_theory Numeric scalar. The theoretical/phenotypic covariance between the traits (often estimated via intercept of initial regression or supplied externally).
#' @param n_snp Integer. Total number of SNPs (M) in the reference panel.
#' @param n_x Integer. Sample size for Trait 1 GWAS.
#' @param n_y Integer. Sample size for Trait 2 GWAS.
#' @param overlap_prop Numeric (0-1). Proportion of sample overlap between the two GWAS (N_overlap / N_total).
#'
#' @return A list containing the slope (h_xy), intercept, and fit details.
#' @keywords internal
.estimate_hxy_struct <- function(betaX, betaY,
                                 res_XX, res_YY,
                                 ld_score, covXY_theory,
                                 n_snp, n_x, n_y, overlap_prop,
                                 idx_subset = NULL) {

  # --- 0. Handle Subsetting ---
  if (!is.null(idx_subset)) {
    betaX <- betaX[idx_subset]
    betaY <- betaY[idx_subset]
    ld_score <- ld_score[idx_subset]

    # Crucial: We must also subset the expectations passed from the single-trait runs
    # assuming res_XX$E_total_sq is the full vector
    res_XX$E_total_sq <- res_XX$E_total_sq[idx_subset]
    res_YY$E_total_sq <- res_YY$E_total_sq[idx_subset]
  }
  # --- 1. Construct Variables ---
  # Define the global noise term ONCE
  # This is (No / Ny) * Sigma_XY
  total_noise_term <- (overlap_prop / n_y) * covXY_theory


  term_N <- (n_snp * overlap_prop) / n_y
  xi <- (1 + overlap_prop/n_y) * ld_score + term_N
  A_xy <- xi - term_N # Predictor (X)

  # Construct Response (Y)
  # M_xy = Raw_Product - Total_Noise
  M_xy <- betaX * betaY - total_noise_term

  # --- 2. Step 1: Initial Estimates (The "WLS" with weights as squared-LD Step) ---
  w_step1 <- 1 / pmax(ld_score, 1)^2

  # Fast Slope (Step 1)
  num_init <- sum(w_step1 * M_xy * A_xy)
  den_init <- sum(w_step1 * A_xy^2)
  h_xy_init <- num_init / den_init


  # ---  Back-Calculate Intercept (Constrained Eq 2) ---
  # (No/NxNy) * Sigma_XY = (M * No / NxNy) * h_xy + v_xy
  # Let "Total Noise" = (No/NxNy) * Sigma_XY
  # Let "Scaling"     = (M * No / NxNy)

  # Note: term_N is exactly (M * No / Ny), so we reuse it!
  v_xy_init <- total_noise_term -  h_xy_init * term_N
  # --- 3. Step 2: Final WLS Weights ---
  # Part 1: Product of single expectations (Passed from Single-Trait Engine)
  prod_expectations <- res_XX$E_total_sq * res_YY$E_total_sq

  # Part 2: Squared cross expectation (Current Step 1 Estimate)
  E_cross_sq <- (h_xy_init * xi + v_xy_init)^2

  # Calculate Isserlis Weights
  w_final <- 1 / (pmax(ld_score, 1) * (prod_expectations + E_cross_sq))

  # --- 4. Step 3: Final Estimate (The "WLS" Step) ---
  num_final <- sum(w_final * M_xy * A_xy)
  den_final <- sum(w_final * A_xy^2)
  h_xy_final <- num_final / den_final

  # --- 5. Back-Calculate Intercept (Constrained Eq 2) ---
  # (No/NxNy) * Sigma_XY = (M * No / NxNy) * h_xy + v_xy
  # Let "Total Noise" = (No/NxNy) * Sigma_XY
  # Let "Scaling"     = (M * No / NxNy)

  # Note: term_N is exactly (M * No / Ny), so we reuse it!
  incpt_final <- total_noise_term - h_xy_final * term_N
  return(list(
    slope = h_xy_final,
    incpt = incpt_final
  ))
}


#' Internal Engine: Simple OLS (Unweighted)
#' @description Performs unweighted regression (w=1).
#' Used as a baseline to demonstrate the value of WLS.
#' @keywords internal
#'

.estimate_hxy_ols <- function(betaX, betaY,
                              ld_score, covXY_theory,
                              n_snp, n_x, n_y, overlap_prop) {

  # --- 1. Construct Variables (Same as WLS) ---
  term_N <- (n_snp * overlap_prop) / n_y
  xi <- (1 + overlap_prop/n_y) * ld_score + term_N
  A_xy <- xi - term_N # Predictor (X)

  noise_correction <- covXY_theory * (overlap_prop / n_y)
  M_xy <- betaX * betaY - noise_correction # Response (Y)

  # --- 2. Calculate Slope (Direct OLS) ---
  # Formula: sum(xy) / sum(x^2)  <-- Implicitly w=1
  num <- sum(M_xy * A_xy)
  den <- sum(A_xy^2)
  h_xy_est <- num / den

  # --- 3. Calculate Intercept ---
  incpt_final <- noise_correction - h_xy_est * term_N

  return(list(
    slope = h_xy_est,
    incpt = incpt_final
  ))
}
