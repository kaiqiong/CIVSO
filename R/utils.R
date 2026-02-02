#' Create Block Data for CIVSO
#'
#' @description
#' Helper function to partition genome-wide summary statistics into
#' the block format required by CIVSO's GLS and Analytic SE engines.
#'
#' @param sumstats A data.frame containing global summary statistics.
#'   Must have a column `snp` (ID) and columns for beta/se.
#' @param ld_list A list of LD correlation matrices (R).
#'   Each matrix must have row/col names matching the SNP IDs in `sumstats`.
#' @param col_map Named vector mapping standard names to your column names.
#'   Default: c(snp="snp", betaX="betaX", seX="seX", betaY="betaY", seY="seY")
#'
#' @return A list of blocks formatted for CIVSO.
#' @export
create_civso_blocks <- function(sumstats, ld_list,
                                col_map = c(snp="snp", betaX="betaX", seX="seX",
                                            betaY="betaY", seY="seY")) {

  formatted_blocks <- list()

  # Ensure SNP column is character for matching
  sumstats[[col_map["snp"]]] <- as.character(sumstats[[col_map["snp"]]])

  for (i in seq_along(ld_list)) {
    R <- ld_list[[i]]

    # Get SNPs in this LD block (assuming R has row names)
    if(is.null(rownames(R))) {
      stop(sprintf("LD matrix in block %d is missing row names (SNP IDs).", i))
    }
    block_snps <- rownames(R)

    # Match Summary Stats to this Block
    # We use match() to ensure strict ordering: Data must match R exactly.
    match_idx <- match(block_snps, sumstats[[col_map["snp"]]])

    # Safety Check: Missing SNPs?
    if (any(is.na(match_idx))) {
      missing_count <- sum(is.na(match_idx))
      warning(sprintf("Block %d: %d SNPs in LD matrix not found in summary stats. Dropping missing from R.", i, missing_count))

      # Alignment Strategy: Intersection
      valid_indices <- which(!is.na(match_idx))
      # 1. Update R to keep only found SNPs
      R <- R[valid_indices, valid_indices]
      # 2. Update match index
      match_idx <- match_idx[valid_indices]
    }

    # Extract Data
    subset_data <- sumstats[match_idx, ]

    # Build the List Element
    formatted_blocks[[i]] <- list(
      R = R,
      betaX = subset_data[[col_map["betaX"]]],
      seX   = subset_data[[col_map["seX"]]],
      betaY = subset_data[[col_map["betaY"]]],
      seY   = subset_data[[col_map["seY"]]]
    )
  }

  return(formatted_blocks)
}
