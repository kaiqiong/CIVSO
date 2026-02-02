#' Internal Engine: Core CIVSO Estimator (Single Iteration)
#' @keywords internal
.civso_engine <- function(idx_subset = NULL,
                          betaX, betaY, seX, seY, ld_score,
                          n_snp, n_x, n_y, overlap_prop, covXY_theory,
                          method = "diagonal", blocks = NULL) {

  # --- 1. Subset Data ---
  if(!is.null(idx_subset)) {
    bX <- betaX[idx_subset]; bY <- betaY[idx_subset]
    sX <- seX[idx_subset];   sY <- seY[idx_subset]; ld <- ld_score[idx_subset]

    # Note: If subsetting, the global constants (term_N) theoretically change
    # if 'overlap_prop' is defined relative to the specific SNPs.
    # However, standard CIVSO assumes 'overlap_prop' is a global sample property.
    # We keep global constants fixed based on N.
  } else {
    bX <- betaX; bY <- betaY; sX <- seX; sY <- seY; ld <- ld_score
  }


  # --- 2. Calculate Geometric Constants (Shared) ---
  # These are used for Back-Calculation AND Final Correction
  # CRITICAL: Always use GLOBAL n_snp (M)
  kappa_j <- (1 + 1/n_x) * ld + n_snp/n_x
  term_N  <- (n_snp * overlap_prop) / n_y
  xi_j    <- (1 + overlap_prop/n_y) * ld + term_N


  # --- 3. Calibration: Step A (Diagonal WLS - The Baseline) ---
  # We calculate everything here. If method="diagonal", these are final.
  # If method="full_gls", these serve as priors.

  calib_X <- .estimate_h2_struct(bX, sX, ld, n_snp, n_x)
  calib_Y <- .estimate_h2_struct(bY, sY, ld, n_snp, n_y)

  calib_XY <- .estimate_hxy_struct(bX, bY, calib_X, calib_Y, ld,
                                   covXY_theory, n_snp, n_x, n_y, overlap_prop)


  # Initialize the "Active" parameters with Diagonal results
  h_xx_use <- calib_X$slope
  v_x_use  <- calib_X$incpt

  h_yy_use <- calib_Y$slope # Needed for Analytic Variance
  v_y_use  <- calib_Y$incpt

  h_xy_use <- calib_XY$slope
  v_xy_use <- calib_XY$incpt

  # SAVE DIAGONAL ESTIMATES (for reporting)
  params_diag <- list(
    h_xx = calib_X$slope,  v_x  = calib_X$incpt,
    h_yy = calib_Y$slope,  v_y  = calib_Y$incpt,
    h_xy = calib_XY$slope, v_xy = calib_XY$incpt
  )
  # --- 4. Calibration: Step B (Full GLS Refinement) ---
  # If active, we calculate better h values and BACK-CALCULATE consistent v values.
  # We then OVERWRITE the "Active" parameters.
  if (method == "full_gls") {
    if(is.null(blocks)) stop("GLS requires 'blocks' list.")

    # 1. Refine h_xx (Heritability)
    num_gls_xx <- 0; den_gls_xx <- 0
    # 2. Refine h_xy (Genetic Covariance) -> NEW!
    num_gls_xy <- 0; den_gls_xy <- 0
    num_gls_yy <- 0; den_gls_yy <- 0  # <--- NEW ACCUMULATORS
    for(b in 1:length(blocks)) {
      # Refine h_xx
      res_gls_xx <- .estimate_block_gls_h2(
        beta = blocks[[b]]$betaX, se = blocks[[b]]$seX, R = blocks[[b]]$R,
        n_total = n_snp, n_samp = n_x,
        h_init = calib_X$slope, v_init = calib_X$incpt
      )
      num_gls_xx <- num_gls_xx + res_gls_xx$num
      den_gls_xx <- den_gls_xx + res_gls_xx$den

      # 2. Refine h_yy (OPTIONAL - Only affects Analytic SE)
      res_gls_yy <- .estimate_block_gls_h2(
        beta = blocks[[b]]$betaY, se = blocks[[b]]$seY, R = blocks[[b]]$R,
        n_total = n_snp, n_samp = n_y,
        h_init = calib_Y$slope, v_init = calib_Y$incpt
      )
      num_gls_yy <- num_gls_yy + res_gls_yy$num
      den_gls_yy <- den_gls_yy + res_gls_yy$den

      # Refine h_xy (We need to run GLS for Cross-Trait too!)
      # Note: We passed the diagonal priors (calib_X, calib_Y) to build Sigma matrices
      res_gls_xy <- .estimate_block_gls_hxy(
        betaX = blocks[[b]]$betaX, betaY = blocks[[b]]$betaY, R = blocks[[b]]$R,
        n_total = n_snp, n_x = n_x, n_y = n_y, overlap_prop = overlap_prop,
        covXY_theory = covXY_theory,
        res_XX_init = list(slope=calib_X$slope, incpt=calib_X$incpt),
        res_YY_init = list(slope=calib_Y$slope, incpt=calib_Y$incpt),
        h_xy_init = calib_XY$slope,
        v_xy_init = calib_XY$incpt
      )
      num_gls_xy <- num_gls_xy + res_gls_xy$num
      den_gls_xy <- den_gls_xy + res_gls_xy$den
    }

    # --- OVERWRITE WITH REFINED ESTIMATES ---
    # A. Update X (Used for Delta_XX)
    h_xx_use <- num_gls_xx / den_gls_xx
    mean_se2_X <- mean(sX^2)
    v_x_use <- mean_se2_X - h_xx_use * ((n_snp - 1) / n_x)

    # B. Update Y (Used for Analytic Variance only)
    h_yy_use <- num_gls_yy / den_gls_yy
    mean_se2_Y <- mean(sY^2)
    v_y_use <- mean_se2_Y - h_yy_use * ((n_snp - 1) / n_y)

    # C. Update XY (Used for Delta_XY)
    h_xy_use <- num_gls_xy / den_gls_xy
    total_noise <- (overlap_prop / n_y) * covXY_theory
    v_xy_use <- total_noise - h_xy_use * term_N



  }

  # --- 3. Correction & Estimation (Eq 5) ---


  # Delta_XY now uses the Refined, Consistent v_xy
  Delta_XY   <- v_xy_use

  # Delta_XX uses the Refined h_xx
  Delta_XX_j <- (kappa_j - xi_j) * h_xx_use + v_x_use

  num_sum <- sum( (bY * bX) - Delta_XY )
  den_sum <- sum( (bX^2) - Delta_XX_j )

  beta_val <- num_sum / den_sum

  # --- 6. Return Structured Output ---
  return(list(
    beta = num_sum / den_sum,
    denom_value = den_sum,

    # Final (Used) Parameters
    params = list(
      h_xx = h_xx_use, v_x = v_x_use,
      h_yy = h_yy_use, v_y = v_y_use,
      h_xy = h_xy_use, v_xy = v_xy_use
    ),

    # Store BOTH for the print method
    structural_fit = list(
      method = method,
      diagonal = params_diag,
      final = list(
        h_xx = h_xx_use, v_x = v_x_use,
        h_yy = h_yy_use, v_y = v_y_use,
        h_xy = h_xy_use, v_xy = v_xy_use
      )
    )
  ))


}
