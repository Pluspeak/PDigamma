# inst/examples/example_comparison.R

# =============================================================================
# Comparison with Alternative Methods
# - Poisson regression
# - Negative Binomial regression  
# - Continuous POT (Peaks-Over-Threshold) with discretization
# =============================================================================

library(PDigamma)
library(MASS)  # For glm.nb

#' Fit and evaluate Poisson regression
fit_poisson <- function(y, X) {
  fit <- glm(y ~ X[, -1, drop = FALSE], family = poisson)
  list(fit = fit, loglik = logLik(fit), aic = AIC(fit))
}

#' Fit and evaluate Negative Binomial regression
fit_nb <- function(y, X) {
  fit <- glm.nb(y ~ X[, -1, drop = FALSE])
  list(fit = fit, loglik = logLik(fit), aic = AIC(fit), theta = fit$theta)
}

#' Continuous POT with discretization correction
#' Simplified version - full version would use evd or extRemes package
fit_pot_discrete <- function(y, threshold_quantile = 0.9) {
  # This is a placeholder - real implementation would use proper EVT packages
  threshold <- quantile(y, threshold_quantile)
  exceedances <- y[y > threshold]
  
  # Fit GPD to exceedances (ignoring discreteness for demo)
  # In practice: use discrete GPD correction or continuous approximation
  
  list(
    threshold = threshold,
    n_exceed = length(exceedances),
    mean_excess = mean(exceedances - threshold),
    note = "Full POT implementation requires evd/extRemes packages"
  )
}

#' Run comparison on dataset
compare_methods <- function(data, y_var = "y", 
                            alpha_vars = NULL, 
                            delta_vars = NULL) {
  
  y <- data[[y_var]]
  n <- length(y)
  
  # Build design matrices
  if (is.null(alpha_vars)) {
    X <- matrix(1, n, 1)
  } else {
    X <- model.matrix(as.formula(paste("~", paste(alpha_vars, collapse = "+"))), data)
  }
  
  if (is.null(delta_vars)) {
    Z <- matrix(1, n, 1)
  } else {
    Z <- model.matrix(as.formula(paste("~", paste(delta_vars, collapse = "+"))), data)
  }
  
  cat("=== Method Comparison ===\n")
  cat("Sample size:", n, "\n")
  cat("Mean:", round(mean(y), 2), "| Var:", round(var(y), 2), 
      "| CV:", round(sd(y)/mean(y), 2), "\n\n")
  
  # 1. Poisson
  cat("1. Poisson Regression:\n")
  tryCatch({
    p_fit <- fit_poisson(y, X)
    cat("   Log-likelihood:", round(as.numeric(p_fit$loglik), 2), "\n")
    cat("   AIC:", round(p_fit$aic, 2), "\n")
    cat("   Note: Assumes equidispersion (often inadequate)\n\n")
  }, error = function(e) cat("   FAILED:", conditionMessage(e), "\n\n"))
  
  # 2. Negative Binomial
  cat("2. Negative Binomial Regression:\n")
  tryCatch({
    nb_fit <- fit_nb(y, X)
    cat("   Log-likelihood:", round(as.numeric(nb_fit$loglik), 2), "\n")
    cat("   AIC:", round(nb_fit$aic, 2), "\n")
    cat("   Theta (dispersion):", round(nb_fit$theta, 3), "\n")
    cat("   Note: Exponential tails, may underestimate extremes\n\n")
  }, error = function(e) cat("   FAILED:", conditionMessage(e), "\n\n"))
  
  # 3. PDigamma (our method)
  cat("3. PDigamma Regression:\n")
  tryCatch({
    pd_result <- pdigamma_analysis(
      formula_alpha = if (is.null(alpha_vars)) ~ 1 else as.formula(paste("~", paste(alpha_vars, collapse = "+"))),
      formula_delta = if (is.null(delta_vars)) ~ 1 else as.formula(paste("~", paste(delta_vars, collapse = "+"))),
      formula_rho = ~ 1,
      data = data,
      y_var = y_var,
      B = 200,  # Reduced for speed
      verbose = FALSE
    )
    cat("   Log-likelihood:", round(pd_result$fit_final$loglik, 2), "\n")
    cat("   AIC:", round(-2 * pd_result$fit_final$loglik + 2 * length(unlist(pd_result$coefficients)), 2), "\n")
    cat("   Mean delta:", round(mean(fitted(pd_result, "delta")), 3), "\n")
    cat("   Note: Polynomial tails, captures heavy-tail heterogeneity\n\n")
  }, error = function(e) cat("   FAILED:", conditionMessage(e), "\n\n"))
  
  # 4. POT (placeholder)
  cat("4. Continuous POT (Peaks-Over-Threshold):\n")
  pot_fit <- fit_pot_discrete(y)
  cat("   Threshold:", round(pot_fit$threshold, 2), "\n")
  cat("   Exceedances:", pot_fit$n_exceed, "\n")
  cat("   Mean excess:", round(pot_fit$mean_excess, 2), "\n")
  cat("   Note:", pot_fit$note, "\n")
  cat("   Limitation: Discretization bias, threshold selection\n\n")
  
  invisible(list(poisson = if (exists("p_fit")) p_fit else NULL,
                 nb = if (exists("nb_fit")) nb_fit else NULL,
                 pdigamma = if (exists("pd_result")) pd_result else NULL,
                 pot = pot_fit))
}

# Example usage
if (FALSE) {
  # Generate heavy-tailed data
  set.seed(42)
  sim_data <- generate_scenario_2(n = 500, xi = c(0.5, 0.3))
  
  # Run comparison
  comparison <- compare_methods(
    data = sim_data,
    y_var = "y",
    alpha_vars = NULL,  # Intercept only
    delta_vars = "z1"   # Tail heterogeneity
  )
}