#' CIVSO: Correlated IVs and Sample Overlap Causal Estimator
#'
#' @description
#' Estimates the causal effect ($\beta_0$) of an exposure on an outcome using GWAS summary statistics.
#' The estimator explicitly corrects for:
#' 1. **Sample Overlap Bias:** Using a phenotypic covariance constraint ($v_{XY}$).
#' 2. **LD-Induced Noise Inflation:** Using structural heritability parameters ($h_{xx}, v_X$).
#'
#' It offers two calibration methods:
#' * **Diagonal (dMAR):** Fast, robust, uses only LD scores. Best for large-scale screens.
#' * **Full-Moment GLS:** High-precision, accounts for local LD correlation structure. Best for final reporting.
#'
#' @param betaX Numeric vector. Effect sizes for the Exposure trait.
#' @param betaY Numeric vector. Effect sizes for the Outcome trait.
#' @param seX Numeric vector. Standard errors for the Exposure trait.
#' @param seY Numeric vector. Standard errors for the Outcome trait.
#' @param ld_score Numeric vector. LD scores ($\ell_j$) for each SNP.
#' @param n_snp Integer. The **global** total number of SNPs ($M$) in the reference panel (e.g., 1,000,000), not just the number of SNPs in the input vectors.
#' @param n_x Integer. Sample size of the Exposure GWAS.
#' @param n_y Integer. Sample size of the Outcome GWAS.
#' @param overlap_prop Numeric. The proportion of sample overlap ($N_{overlap} / N_{total}$). Ranges from 0 (disjoint) to 1 (full overlap).
#' @param covXY_theory Numeric. The theoretical or estimated phenotypic covariance ($\sigma_{XY}$) between the two traits. Often estimated via the intercept of an initial LDSC regression.
#' @param method Character string. Either `"diagonal"` (default) or `"full_gls"`.
#'   * `"diagonal"`: Uses fast Diagonal Moment-Alignment Regression.
#'   * `"full_gls"`: Uses Full-Moment GLS to refine structural parameters ($h_{xx}, h_{xy}$). Requires `blocks`.
#' @param blocks List. A list of LD blocks (required if `method="full_gls"` or for Analytic SE).
#'   Each element of the list must contain:
#'   * `$R`: The LD correlation matrix for the block.
#'   * `$betaX`, `$seX`, `$betaY`, `$seY`: Subset vectors for the block.
#' @param n_jack_blocks Integer. Number of blocks to use for the Block-Jackknife SE estimation (default: 200). Ignored if `blocks` are provided explicitly.
#'
#' @return A list of class `"CIVSO"` containing:
#' \item{beta}{The causal effect point estimate.}
#' \item{se}{The primary standard error (defaults to Jackknife SE, or Analytic SE if Jackknife is impossible).}
#' \item{p}{The P-value associated with the primary SE.}
#' \item{se_jack}{Standard error estimated via Block-Jackknife.}
#' \item{p_jack}{P-value derived from the Jackknife SE.}
#' \item{se_analytic}{Standard error estimated via the asymptotic Delta Method formula (requires `blocks`).}
#' \item{p_analytic}{P-value derived from the Analytic SE.}
#' \item{method}{The method used ("diagonal" or "full_gls").}
#' \item{params}{A list of the estimated structural parameters used in the correction ($h_{xx}, v_x, h_{xy}, v_{xy}, h_{yy}, v_y$).}
#' \item{n_snps}{The global SNP count used in calculations.}
#'
#'@export
CIVSO <- function(betaX, betaY, seX, seY, ld_score, n_snp, n_x, n_y, overlap_prop,
                  covXY_theory, method = "diagonal", blocks = NULL,
                  n_jack_blocks = 200) {

  # --- 1. Point Estimate ---
  est_result <- .civso_engine(
    idx_subset = NULL,
    betaX = betaX, betaY = betaY, seX = seX, seY = seY, ld_score = ld_score,
    n_snp = n_snp, n_x = n_x, n_y = n_y, overlap_prop = overlap_prop,
    covXY_theory = covXY_theory, method = method, blocks = blocks
  )
  beta_point <- est_result$beta

  # --- 2. Analytic Variance (Ideal for Single Block) ---
  se_analytic <- NA
  p_analytic  <- NA

  if (!is.null(blocks)) {
    var_an <- .compute_analytic_variance(
      blocks = blocks,
      beta_hat = beta_point,

      # Structural parameters
      h_xx = est_result$params$h_xx, v_x  = est_result$params$v_x,
      h_yy = est_result$params$h_yy, v_y  = est_result$params$v_y,
      h_xy = est_result$params$h_xy, v_xy = est_result$params$v_xy,

      # Constants
      n_snp = n_snp, n_x = n_x, n_y = n_y, overlap_prop = overlap_prop
    )
    se_analytic <- sqrt(var_an)
    p_analytic  <- 2 * pnorm(abs(beta_point / se_analytic), lower.tail = FALSE)
  }

  # --- 3. Jackknife Variance (Safety Check) ---
  se_jack <- NA
  p_jack  <- NA

  # Check if we rely on 'blocks' (GLS) or automatic blocking (Diagonal)
  # Logic: If blocks are provided explicitly, check length.
  # If doing automatic jackknife (blocks=NULL), check n_jack_blocks vs n_snps.

  can_jackknife <- TRUE
  if (!is.null(blocks) && length(blocks) < 2) {
    can_jackknife <- FALSE
    warning("Only 1 block provided. Skipping Jackknife. Please use Analytic SE.")
  }

  if (can_jackknife) {
    jack_res <- .compute_jackknife_variance(
      estimator_func = .civso_engine,
      n_blocks = n_jack_blocks,
      total_snps = length(betaX),
      # Pass Data
      betaX = betaX, betaY = betaY, seX = seX, seY = seY, ld_score = ld_score,
      n_snp = n_snp, n_x = n_x, n_y = n_y, overlap_prop = overlap_prop,
      covXY_theory = covXY_theory, method = method, blocks = blocks
    )
    se_jack <- jack_res$se
    p_jack  <- 2 * pnorm(abs(beta_point / se_jack), lower.tail = FALSE)
  }

  # --- 4. Return ---
  # Prioritize Analytic SE if Jackknife failed
  report_se <- if(is.na(se_jack)) se_analytic else se_jack
  report_p  <- if(is.na(p_jack)) p_analytic else p_jack


  # --- 4. Return Object ---
  # We return a list with class "CIVSO" so we can write a pretty print method later
  res <- list(
    beta = beta_point,
    se = report_se,       # Smart default
    p = report_p,         # Smart default

    # Standard outputs
    se_jack = se_jack,
    p_jack = p_jack,
    se_analytic = se_analytic,
    p_analytic = p_analytic,

    method = method,
    n_snps = n_snp,

    # The Structural Model (Class: dMAR)
    model_fit = est_result$structural_fit
  )

  # Assign classes
  class(res$model_fit) <- "dMAR"  # <--- MAGIC LINE
  class(res) <- "CIVSO"

  return(res)
}
