#---------------------------------------------------------------
# Protein Turnover Model: Single Compartment
#---------------------------------------------------------------

simulate_turnover <- function(half_life = 14, days = 60) {
  # Calculate rate constant (per day)
  k <- log(2) / half_life
  
  # Time vector
  time <- seq(0, days, by = 0.5)
  
  # Fractional synthesis over time
  frac_synth <- 1 - exp(-k * time)
  
  # Plot
  plot(
    time, frac_synth,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = "Time (days)",
    ylab = "Fractional Synthesis",
    main = paste0("Protein Turnover (Half-life = ", half_life, " days)")
  )
  grid()
  
  # Return data frame for downstream use
  return(data.frame(time = time, frac_synth = frac_synth))
}

#---------------------------------------------------------------
# Example: 14-day half-life
#---------------------------------------------------------------
turnover_data <- simulate_turnover(half_life = 14)
#======================================================================#

#===========================================================
# Multi-pool Protein Turnover Model (Miller-style dynamic pool)
#===========================================================

simulate_multipool_turnover <- function(
    half_lives   = c(3, 30),     # vector of half-lives (days) for dynamic pools
    pool_fracs   = c(0.3, 0.6),  # corresponding fractions of total pool
    static_frac  = 0.1,          # fraction that does NOT turn over (on this timescale)
    days         = 60,
    dt           = 0.5,
    normalize    = FALSE,        # if TRUE, normalize by dynamic pool (i.e. / (1 - static_frac))
    main_title   = "Multi-pool Protein Turnover"
) {
  #----- basic checks -----
  if (length(half_lives) != length(pool_fracs)) {
    stop("half_lives and pool_fracs must have the same length.")
  }
  
  total_frac <- sum(pool_fracs) + static_frac
  if (abs(total_frac - 1) > 1e-6) {
    stop("pool_fracs + static_frac must sum to 1.")
  }
  
  #----- parameters -----
  k_vals <- log(2) / half_lives
  time   <- seq(0, days, by = dt)
  
  #----- compute fraction new from each dynamic pool -----
  frac_new_mat <- sapply(seq_along(k_vals), function(i) {
    pool_fracs[i] * (1 - exp(-k_vals[i] * time))
  })
  
  # If only one pool, sapply returns a vector; ensure matrix
  if (is.vector(frac_new_mat)) {
    frac_new_mat <- matrix(frac_new_mat, ncol = 1)
  }
  
  # Total fraction new relative to TOTAL pool
  frac_new_total <- rowSums(frac_new_mat)
  
  # Optionally normalize to dynamic pool only
  if (normalize) {
    dyn_pool <- 1 - static_frac
    frac_new_plot <- frac_new_total / dyn_pool
    ylab_txt <- "Fraction New (of dynamic pool)"
  } else {
    frac_new_plot <- frac_new_total
    ylab_txt <- "Fraction New (of total pool)"
  }
  
  #----- plot -----
  plot(
    time, frac_new_plot,
    type = "l", lwd = 2,
    xlab = "Time (days)",
    ylab = ylab_txt,
    main = main_title
  )
  grid()
  
  #----- return full data frame -----
  out <- data.frame(
    time = time,
    frac_new_total = frac_new_total,
    frac_new_plot  = frac_new_plot
  )
  
  # add columns for each pool's contribution if helpful
  for (i in seq_along(k_vals)) {
    out[[paste0("pool", i, "_new")]] <- frac_new_mat[, i]
  }
  
  return(out)
}

#===========================================================
# Example:
# Fast pool: t1/2 = 3 d, 30% of total
# Slow pool: t1/2 = 30 d, 60% of total
# Static:    10% of total
#===========================================================
turnover_multi <- simulate_multipool_turnover(
  half_lives  = c(3, 30),
  pool_fracs  = c(0.3, 0.6),
  static_frac = 0.1,
  days        = 60,
  main_title  = "Multi-pool Turnover (3d, 30d, 10% static)"
)
#============================================================#
#============================================================#

#===========================================================
# Compare Two Multi-Pool Turnover Curves (Miller-style)
#===========================================================

compare_turnover_curves <- function(
    # Curve A parameters
  half_lives_A   = c(3, 30),
  pool_fracs_A   = c(0.3, 0.6),
  static_frac_A  = 0.1,
  
  # Curve B parameters
  half_lives_B   = c(5, 50),
  pool_fracs_B   = c(0.3, 0.6),
  static_frac_B  = 0.1,
  
  days           = 60,
  dt             = 0.5,
  normalize      = FALSE,
  labels         = c("Curve A", "Curve B"),
  colors         = c("blue", "red"),
  main_title     = "Comparison of Two Protein Turnover Curves"
) {
  #----- Helper function (reuses multi-pool logic) -----
  compute_curve <- function(half_lives, pool_fracs, static_frac) {
    k_vals <- log(2) / half_lives
    time   <- seq(0, days, by = dt)
    
    frac_new_mat <- sapply(seq_along(k_vals), function(i) {
      pool_fracs[i] * (1 - exp(-k_vals[i] * time))
    })
    if (is.vector(frac_new_mat)) frac_new_mat <- matrix(frac_new_mat, ncol = 1)
    frac_new_total <- rowSums(frac_new_mat)
    if (normalize) frac_new_total <- frac_new_total / (1 - static_frac)
    
    data.frame(time = time, frac_new_total = frac_new_total)
  }
  
  #----- Compute both curves -----
  curveA <- compute_curve(half_lives_A, pool_fracs_A, static_frac_A)
  curveB <- compute_curve(half_lives_B, pool_fracs_B, static_frac_B)
  
  #----- Plot -----
  plot(
    curveA$time, curveA$frac_new_total,
    type = "l", lwd = 2, col = colors[1],
    xlab = "Time (days)",
    ylab = ifelse(normalize, "Fraction New (dynamic pool)", "Fraction New (total pool)"),
    main = main_title,
    ylim = range(c(curveA$frac_new_total, curveB$frac_new_total))
  )
  lines(curveB$time, curveB$frac_new_total, lwd = 2, col = colors[2])
  grid()
  legend("bottomright", legend = labels, col = colors, lwd = 2, bty = "n")
  
  #----- Return combined data frame -----
  merged <- data.frame(
    time = curveA$time,
    curveA = curveA$frac_new_total,
    curveB = curveB$frac_new_total
  )
  return(merged)
}

#===========================================================
# Example: Compare "Young" vs "Old" turnover
#===========================================================
compare_turnover_curves(
  half_lives_A  = c(1, 14, 30),  pool_fracs_A = c(0.2,0.2, 0.3), static_frac_A = 0.3,  # young
  half_lives_B  = c(1, 14, 40),  pool_fracs_B = c(0.2,0.2, 0.2), static_frac_B = 0.4,  # old
  main_title    = "Young (blue) vs Old (red) Multi-Pool Turnover"
)
