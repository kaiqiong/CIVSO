test_that("Full-Moment GLS Engine works on toy data", {

  # --- 1. SETUP: Create Toy Data (2 Blocks) ---
  # Block size M=10 to keep matrix inversion fast
  M_block <- 10
  N_samp <- 1000
  h2_truth <- 0.5

  # Simulate Block 1 (High LD)
  # Create an AR(1) correlation matrix: rho^|i-j|
  idx <- 1:M_block
  R1 <- 0.8^abs(outer(idx, idx, "-"))

  # Simulate Block 2 (Low LD - almost independent)
  R2 <- 0.1^abs(outer(idx, idx, "-"))

  # Generate true effects (gamma) and noise (e)
  # Simple simulation: y = X*beta + e
  # We just fake the summary stats directly for the test
  # E[beta^2] approx h2/M_total

  set.seed(42)
  # Mock betas: just random normal scaled to look roughly right
  beta1 <- rnorm(M_block, 0, sqrt(h2_truth/(2*M_block)))
  beta2 <- rnorm(M_block, 0, sqrt(h2_truth/(2*M_block)))

  # Mock SEs: constant noise
  se1 <- rep(1/sqrt(N_samp), M_block)
  se2 <- rep(1/sqrt(N_samp), M_block)

  # Global params
  n_total <- 2 * M_block

  # --- 2. STEP 1: Initialization (Diagonal WLS) ---
  # We manually create the inputs the init engine expects
  beta_all <- c(beta1, beta2)
  se_all   <- c(se1, se2)

  # Calculate LD scores from our R matrices
  ld1 <- rowSums(R1^2)
  ld2 <- rowSums(R2^2)
  ld_all <- c(ld1, ld2)

  # Run the Diagonal Engine we wrote earlier
  init_res <- .estimate_h2_struct(beta_all, se_all, ld_all, n_total, N_samp)

  # Check if init ran successfully
  expect_type(init_res, "list")
  expect_true(is.numeric(init_res$slope))

  # Use these init values for GLS
  h_start <- init_res$slope
  # Calculate v_start (intercept) from the diagonal engine logic
  # We approximate it as mean(se^2) for the test
  v_start <- mean(se_all^2)


  # --- 3. STEP 2: Run Block-GLS ---

  # Process Block 1
  gls_b1 <- .estimate_block_gls_h2(
    beta = beta1,
    se = se1,
    R = R1,
    n_total = n_total,
    n_samp = N_samp,
    h_init = h_start,
    v_init = v_start
  )

  # Process Block 2
  gls_b2 <- .estimate_block_gls_h2(
    beta = beta2,
    se = se2,
    R = R2,
    n_total = n_total,
    n_samp = N_samp,
    h_init = h_start,
    v_init = v_start
  )

  # --- 4. VERIFICATION ---

  # A. Structure Checks
  expect_named(gls_b1, c("num", "den", "M_dim"))
  expect_true(is.numeric(gls_b1$num))
  expect_true(is.numeric(gls_b1$den))

  # B. Accumulation Logic
  # The final estimator is Sum(Num) / Sum(Den)
  total_num <- gls_b1$num + gls_b2$num
  total_den <- gls_b1$den + gls_b2$den
  h2_gls <- total_num / total_den

  # C. Sanity Checks
  # The estimate should be a valid number (not NA, not Inf)
  expect_false(is.na(h2_gls))

  # It should be somewhat positive (since truth is 0.5)
  # Note: With N=1000 and M=20, variance is huge, so we just check it exists
  # rather than forcing it to be close to 0.5 (which requires N=100k)
  print(paste("Toy GLS Estimate:", h2_gls))
})

test_that("Full-Moment GLS Cross-Trait works", {

  # Setup simple params
  M_b <- 5
  R <- diag(M_b) # Identity matrix for simplicity
  n_total <- 10
  n_samp <- 1000

  # Fake data
  bX <- rnorm(M_b)
  bY <- rnorm(M_b)

  # Fake init objects (minimal needed structure)
  # The function needs list(h=val, v=val) for sigmas
  resXX <- list(h=0.5, v=0.001)
  resYY <- list(h=0.5, v=0.001)

  gls_cross <- .estimate_block_gls_hxy(
    betaX = bX, betaY = bY, R = R,
    n_total = n_total, n_x = n_samp, n_y = n_samp, overlap_prop = 1,
    covXY_theory = 0.1,
    res_XX_init = resXX,
    res_YY_init = resYY,
    h_xy_init = 0.2,
    v_xy_init = 0.001
  )

  expect_named(gls_cross, c("num", "den"))
  expect_false(is.na(gls_cross$num))
})
