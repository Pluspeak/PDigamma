# =============================================================================
# simulations/generate_scenarios.R
# Generate simulation scenarios from Section 4 of paper
# =============================================================================

library(PDigamma)  # Your package

#' Generate Simulation Scenario 1: Homogeneous Heavy Tails
#' 
#' Fixed parameters, varying sample size
generate_scenario_1 <- function(n, alpha = 2, delta = 2.5, rho = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y <- rPDigamma(n, alpha, delta, rho)
  
  data.frame(
    y = y,
    x1 = rep(0, n),  # No covariates (intercept only)
    z1 = rep(0, n),
    w1 = rep(0, n),
    true_alpha = alpha,
    true_delta = delta,
    true_rho = rho,
    scenario = "homogeneous"
  )
}

#' Generate Simulation Scenario 2: Heterogeneous Tail Index
#' 
#' Delta varies with covariate: log(delta - 1) = xi0 + xi1 * z1
generate_scenario_2 <- function(n, xi = c(0.5, 0.3), alpha = 2, rho = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Covariate
  z1 <- rnorm(n, mean = 0, sd = 1)
  
  # Varying delta
  delta <- 1 + exp(xi[1] + xi[2] * z1)
  
  # Generate data
  y <- numeric(n)
  for (i in 1:n) {
    y[i] <- rPDigamma(1, alpha, delta[i], rho)
  }
  
  data.frame(
    y = y,
    x1 = rep(0, n),
    z1 = z1,
    w1 = rep(0, n),
    true_alpha = alpha,
    true_delta = delta,
    true_rho = rho,
    scenario = "heterogeneous_delta"
  )
}

#' Generate Scenario 3: Near Sibuya Boundary (delta ≈ 1.2)
generate_scenario_3 <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Very heavy tails
  delta <- 1.2
  
  y <- rPDigamma(n, alpha = 0.5, delta = delta, rho = 0.5)
  
  data.frame(
    y = y,
    x1 = rep(0, n),
    z1 = rep(0, n),
    w1 = rep(0, n),
    true_alpha = 0.5,
    true_delta = delta,
    true_rho = 0.5,
    scenario = "sibuya_zone"
  )
}

# Generate and save all scenarios
set.seed(42)

# Small sample for testing
saveRDS(generate_scenario_1(n = 100, seed = 1), 
        "data/processed/sim_scenario1_n100.rds")
saveRDS(generate_scenario_2(n = 100, seed = 2), 
        "data/processed/sim_scenario2_n100.rds")

# Larger samples for paper
saveRDS(generate_scenario_1(n = 1000, seed = 10), 
        "data/processed/sim_scenario1_n1000.rds")
saveRDS(generate_scenario_2(n = 1000, seed = 20), 
        "data/processed/sim_scenario2_n1000.rds")

message("Simulation datasets generated in data/processed/")