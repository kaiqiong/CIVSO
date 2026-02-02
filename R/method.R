#' Print method for CIVSO results
#' @export
print.CIVSO <- function(x, ...) {
  cat("\n========================================\n")
  cat("   CIVSO Causal Estimator  \n")
  cat("========================================\n")

  cat(sprintf("Method used   : %s\n",
              ifelse(x$method=="full_gls", "Full-Moment GLS (Refined)", "Diagonal dMAR (Fast)")))
  cat(sprintf("SNPs (Global) : %s\n", format(x$n_snps, big.mark=",")))
  cat("\n")

  # 1. Main Result
  cat("Causal Estimate:\n")
  cat("----------------\n")
  cat(sprintf("Beta          : %.5f\n", x$beta))
  cat("\n")

  # 2. Standard Errors (Show all available)
  cat("Inference:\n")
  cat("----------------\n")

  # Check availability
  has_jack <- !is.null(x$se_jack) && !is.na(x$se_jack)
  has_anal <- !is.null(x$se_analytic) && !is.na(x$se_analytic)

  if (has_jack) {
    cat(sprintf("SE (Jackknife): %.5f  [P = %.3e]\n", x$se_jack, x$p_jack))
  }

  if (has_anal) {
    cat(sprintf("SE (Analytic) : %.5f  [P = %.3e]\n", x$se_analytic, x$p_analytic))
  }

  if (!has_jack && !has_anal) {
    cat("SE            : Not available (Check warnings)\n")
  }
  cat("\n")
  # 2. Print Structural Parameters (Calls print.dMAR)
  print(x$model_fit)

  cat("\n========================================\n")
}

#' Print method for dMAR structural parameters
#' @export
print.dMAR <- function(x, ...) {
  cat("Structural Model Parameters (dMAR/GLS):\n")
  cat("----------------------------------------\n")

  # Helper to format numbers
  fmt <- function(v) sprintf("%.8f", v)

  # Check if we have comparison data (i.e., if method was GLS)
  is_gls <- x$method == "full_gls"

  if (is_gls) {
    # Print Comparison Table
    cat(sprintf("%-10s | %-15s | %-15s\n", "Parameter", "Diagonal (WLS)", "Full-Moment (GLS)"))
    cat(paste0(rep("-", 46), collapse=""), "\n")

    cat(sprintf("%-10s | %-15s | %-15s\n", "h2 (Exp)",
                fmt(x$diagonal$h_xx), fmt(x$final$h_xx)))
    cat(sprintf("%-10s | %-15s | %-15s\n", "h2 (Out)",
                fmt(x$diagonal$h_yy), fmt(x$final$h_yy)))
    cat(sprintf("%-10s | %-15s | %-15s\n", "Cov (Gen)",
                fmt(x$diagonal$h_xy), fmt(x$final$h_xy)))
    cat(paste0(rep("-", 46), collapse=""), "\n")
    cat(sprintf("%-10s | %-15s | %-15s\n", "Int (Exp)",
                fmt(x$diagonal$v_x), fmt(x$final$v_x)))
    cat(sprintf("%-10s | %-15s | %-15s\n", "Int (XY)",
                fmt(x$diagonal$v_xy), fmt(x$final$v_xy)))

  } else {
    # Print Single Column (Diagonal Only)
    cat(sprintf("%-15s : %s\n", "h2 (Exposure)", fmt(x$final$h_xx)))
    cat(sprintf("%-15s : %s\n", "h2 (Outcome)",  fmt(x$final$h_yy)))
    cat(sprintf("%-15s : %s\n", "Cov (Genetic)", fmt(x$final$h_xy)))
    cat(sprintf("%-15s : %s\n", "Int (Overlap)", fmt(x$final$v_xy)))
  }
}
