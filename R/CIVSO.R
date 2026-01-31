#' Main CIVSO Estimation Function
CIVSO <- function(betaX, betaY, seX, seY, ld_score,
                  n_snp, n_x, n_y, overlap_prop,
                  method = "diagonal_WLS", # Default to the Supp Note method
                  param_alpha = 0.05) {

  # --- 1. Pre-computation (Standardize Inputs) ---
  # Calculate kappa, xi, etc. once here to avoid repetition
  kappa_x <- (1 + 1/n_x) * ld_score + n_snp/n_x
  xi_shared <- (1 + overlap_prop/n_y) * ld_score + (n_snp * overlap_prop)/n_y

  # --- 2. Estimate Nuisance Parameters (Intercepts) ---
  # You might want to break this out if using intercept-based constraints

  # --- 3. Method Dispatch ---
  if (method == "diagonal_OLS") {
    fit <- .civso_engine_OLS(betaX, betaY, seX, seY, kappa_x, xi_shared, ...)

  } else if (method == "diagonal_WLS") {
    # This matches your Supp Note 2.2.1
    fit <- .civso_engine_WLS(betaX, betaY, seX, seY, kappa_x, xi_shared, ...)

  } else if (method == "full_GLS") {
    # The heavy computation
    fit <- .civso_engine_GLS(betaX, betaY, seX, seY, ...)
    # (Optional) Keep your current script here ONLY for legacy comparison
    # But mark it as "experimental" or "not recommended".
    warning("Full-moment OLS is experimental. Use diagonal_wls for best results.")
  }

  # --- 4. Variance Calculation ---
  # Compute SE and P-val using the final point estimate
  stats <- .calc_inference_stats(fit$est, fit$var_est, param_alpha)

  # --- 5. Construct Output ---
  # Essential for the user
  res <- list(
    method = method,
    beta = stats$beta,
    se = stats$se,
    pval = stats$pval,
    ci_lower = stats$ci_lower,
    ci_upper = stats$ci_upper
  )

  # Intermediates for YOU (The Developer/Explorer)
  # Hide these in a "details" list so they don't clutter the screen
  res$details <- list(
    h2_x = fit$h2_x,         # Useful to debug heritability
    h2_y = fit$h2_y,
    genetic_cov = fit$h_xy,
    intercept_x = fit$int_x, # Useful to check for stratification
    intercept_y = fit$int_y,
    weights_used = fit$weights # Great for diagnostic plots
  )

  class(res) <- "CIVSO_result" # S3 class for nice printing
  return(res)
}

# --- Internal Helper for Printing ---
print.CIVSO_result <- function(x) {
  cat(paste0("CIVSO Estimate (", x$method, "):\n"))
  cat(paste0("  Beta: ", round(x$beta, 4), "\n"))
  cat(paste0("  SE:   ", round(x$se, 4), "\n"))
  cat(paste0("  Pval: ", format.pval(x$pval), "\n"))
}
