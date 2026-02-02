# ==============================================================================
# TEST SCRIPT: CIVSO Package Validation (Corrected)
# ==============================================================================
message("Loading CIVSO source files...")
# Ensure these match your actual filenames
source("CIVSO.R")
source("engine_civso.R")
source("engines_diagonal_dMAR.R")
source("engines_full_moment_dMAR.R")
source("method.R")
source("utils.R")
source("variance_analytic.R")
source("variance_jackknife.R") # Matches your uploaded filename

if(!require(MASS)) install.packages("MASS")

# ==============================================================================
# 2. Simulate Dummy Data (With Strong Signal & Row Names)
# ==============================================================================
message("\nSimulating dummy GWAS data...")
set.seed(123)

n_snps_block <- 50
n_blocks     <- 4
n_snps_total <- n_snps_block * n_blocks
n_x <- 50000; n_y <- 50000; overlap_prop <- 0.5

# Create SNP IDs first
snp_ids <- paste0("rs", 1:n_snps_total)

betaX <- numeric(n_snps_total)
betaY <- numeric(n_snps_total)
ld_scores <- numeric(n_snps_total)
ld_list <- list()

# Standard Errors
seX <- rep(1/sqrt(n_x), n_snps_total)
seY <- rep(1/sqrt(n_y), n_snps_total)
names(seX) <- names(seY) <- names(betaX) <- names(betaY) <- snp_ids

for (k in 1:n_blocks) {
  # 1. Create LD Matrix
  X_rand <- matrix(rnorm(n_snps_block^2), ncol=n_snps_block)
  R_blk  <- cov2cor(cov(X_rand))
  diag(R_blk) <- 1.05
  R_blk <- cov2cor(R_blk)

  # --- CRITICAL FIX: ASSIGN NAMES ---
  indices <- ((k-1)*n_snps_block + 1) : (k*n_snps_block)
  blk_snps <- snp_ids[indices]
  rownames(R_blk) <- colnames(R_blk) <- blk_snps

  ld_list[[k]] <- R_blk

  # 2. Calculate LD Scores
  scores <- colSums(R_blk^2)
  ld_scores[indices] <- scores

  # 3. Generate Betas CORRELATED with LD (Strong Signal)
  # Signal correlates with sqrt(LD Score)
  signal_X <- rnorm(n_snps_block, sd=0.05) * sqrt(scores)

  betaX[indices] <- signal_X + rnorm(n_snps_block, sd=1/sqrt(n_x))
  betaY[indices] <- 0.3 * signal_X + rnorm(n_snps_block, sd=1/sqrt(n_y))
}

# Combine into dataframe
sumstats <- data.frame(
  snp = snp_ids,
  betaX = betaX, seX = seX,
  betaY = betaY, seY = seY,
  ld = ld_scores
)

# ==============================================================================
# 3. Test Helper & Pipeline
# ==============================================================================
message("\nTesting 'create_civso_blocks' helper...")

civso_blocks <- create_civso_blocks(
  sumstats = sumstats,
  ld_list = ld_list,
  col_map = c(snp="snp", betaX="betaX", seX="seX", betaY="betaY", seY="seY")
)

if(length(civso_blocks) == n_blocks) message("[PASS] Blocks created.")

message("\n------------------------------------------------")
message("Running Test 1: Diagonal dMAR...")
res_diag <- CIVSO(
  betaX=sumstats$betaX, betaY=sumstats$betaY, seX=sumstats$seX, seY=sumstats$seY,
  ld_score=sumstats$ld, n_snp=n_snps_total, n_x=n_x, n_y=n_y, overlap_prop=overlap_prop,
  covXY_theory=0.0, method="diagonal", blocks=NULL, n_jack_blocks=10
)
print(res_diag)

message("\n------------------------------------------------")
message("Running Test 2: Full-Moment GLS...")
res_gls <- CIVSO(
  betaX=sumstats$betaX, betaY=sumstats$betaY, seX=sumstats$seX, seY=sumstats$seY,
  ld_score=sumstats$ld, n_snp=n_snps_total, n_x=n_x, n_y=n_y, overlap_prop=overlap_prop,
  covXY_theory=0.0, method="full_gls", blocks=civso_blocks, n_jack_blocks=4
)
print(res_gls)

message("\n[SUCCESS] Pipeline finished.")
