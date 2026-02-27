# =============================================================================
# 2_em_core.R
# Expectation-Maximization algorithm for PDigamma regression (Algorithm 3.1)
# Implements E-step (Proposition 3.1), M-step with block coordinate ascent,
# and convergence monitoring with hard boundary projection
# =============================================================================

#' E-Step: Conditional Expectation of Latent Variables (Proposition 3.1)
#' 
#' Computes E[U_i | Y_i] and E[log V_i | Y_i] for each observation given current parameters.
#' These are the sufficient statistics needed for M-step.
#' 
#' @param y Vector of observed counts (positive integers)
#' @param alpha_vec Vector of individual alpha parameters
#' @param delta_vec Vector of individual delta parameters
#' @param rho_vec Vector of individual rho parameters
#' @return List containing:
#'   \item{u_vec}{E[U_i | Y_i] for each observation}
#'   \item{v_vec}{E[log V_i | Y_i] for each observation}
#'   \item{aux_info}{Additional diagnostics (optional)}
#' @details 
#'   From Proposition 3.1 and latent variable representation (Proposition 2.2):
#'   - U_i | Y_i, V_i ~ Gamma(alpha_i + Y_i, 1 + V_i^(-1/eta_i))
#'   - V_i | Y_i has complex marginal; we use approximation via digamma functions
#'   
#'   Key formulas:
#'   u_i = (alpha_i + y_i) / (1 + E[V_i^(-1/eta_i) | Y_i])
#'   E[log V_i | Y_i] involves psi(gamma_i) - psi(alpha_i + gamma_i) + ...
#' @references Proposition 3.1, Section 3.3.2
#' @seealso M_step, fit_em
#' @export
E_step <- function(y, alpha_vec, delta_vec, rho_vec) {
  n <- length(y)
  
  # Transform to original parameterization
  gamma_vec <- rho_vec * delta_vec / (1 + rho_vec)
  eta_vec <- delta_vec / (1 + rho_vec)
  
  # Initialize output vectors
  u_vec <- numeric(n)
  v_vec <- numeric(n)  # E[log V_i | Y_i]
  
  for (i in 1:n) {
    ai <- alpha_vec[i]
    gi <- gamma_vec[i]
    etai <- eta_vec[i]
    yi <- y[i]
    
    # Compute E[V^(-1/eta) | Y] using numerical integration or approximation
    # Exact form requires integrating over V's posterior
    # We use Laplace approximation: posterior of V given Y is concentrated
    
    # Approximation: treat V as fixed at its mode, then correct
    # Mode of V|Y approximately satisfies: (gi - 1)/V - 1 - (yi/etai)/V^(1/etai + 1) * U_mode = 0
    # This is complex; instead use moment approximation
    
    # Simplified approach: use relationship from Gamma-Poisson hierarchy
    # E[V^(-1/eta) | Y] ≈ exp(-1/eta * E[log V | Y] + (1/(2*eta^2)) * Var(log V | Y))
    
    # First, approximate E[log V | Y] using digamma differences
    # From theory: E[log V | Y] ≈ psi(gi) - psi(ai + gi) + log(etai) + psi(ai + yi) - psi(ai)
    
    elogV <- digamma(gi) - digamma(ai + gi) + log(etai) + digamma(ai + yi) - digamma(ai)
    
    # Approximate variance of log V (using trigamma)
    # Var(log V | Y) ≈ trigamma(gi) - trigamma(ai + gi) + trigamma(ai + yi) - trigamma(ai)
    # But this ignores dependence; use conservative estimate
    varlogV <- max(trigamma(gi) - trigamma(ai + gi), 0.01)
    
    # Compute E[V^(-1/etai)] using log-normal approximation
    # If log V ~ N(elogV, varlogV), then E[V^(-1/etai)] = exp(-elogV/etai + varlogV/(2*etai^2))
    inv_eta <- 1 / etai
    EV_inv <- exp(-elogV * inv_eta + varlogV * inv_eta^2 / 2)
    
    # Compute u_i = E[U | Y]
    # U | Y, V ~ Gamma(ai + yi, 1 + V^(-1/etai))
    # E[U | Y] = E[ E[U | Y, V] | Y ] = E[ (ai + yi) / (1 + V^(-1/etai)) | Y ]
    # Approximate: (ai + yi) / (1 + E[V^(-1/etai) | Y])
    u_vec[i] <- (ai + yi) / (1 + EV_inv)
    
    # Store E[log V | Y]
    v_vec[i] <- elogV
    
    # Numerical stability checks
    if (!is.finite(u_vec[i]) || u_vec[i] <= 0) {
      warning("Non-finite u_i at observation ", i, "; using fallback")
      u_vec[i] <- ai + yi  # Conservative: ignore denominator correction
    }
    if (!is.finite(v_vec[i])) {
      v_vec[i] <- log(gi)  # Fallback: log of prior mean
    }
  }
  
  # Ensure reasonable bounds
  u_vec <- pmax(u_vec, 0.001)
  
  return(list(u_vec = u_vec, v_vec = v_vec))
}

#' M-Step Block (a): Update Beta (Alpha/Proximal Shape Parameters)
#' 
#' Maximizes Q-function for alpha using IRLS (Iteratively Reweighted Least Squares).
#' 
#' @param y Vector of observed counts (not directly used, but kept for interface)
#' @param X Design matrix for alpha
#' @param u_vec Conditional expectations E[U_i | Y_i] from E-step
#' @param beta_old Previous beta values
#' @param alpha_min Minimum alpha for stability (default: 0.01)
#' @param maxit Maximum IRLS iterations (default: 50)
#' @param tol Convergence tolerance (default: 1e-6)
#' @return Updated beta vector
#' @details 
#'   Q_beta = sum_i [ (alpha_i - 1) * E[log U_i | Y_i] - alpha_i - log Gamma(alpha_i) + ... ]
#'   Using u_vec as proxy for E[log U | Y] (since E[U] relates to exp(E[log U]) for Gamma)
#'   
#'   Actually: for U ~ Gamma(alpha, 1), E[log U] = psi(alpha)
#'   So we match psi(alpha_i) to log(u_vec[i]) approximately
#'   
#'   This becomes a Gamma GLM with log link: alpha_i = exp(X_i' beta)
#'   Response: "working" variable derived from u_vec
#' @seealso E_step, update_xi, update_zeta
#' @export
update_beta <- function(y, X, u_vec, beta_old, alpha_min = 0.01, 
                        maxit = 50, tol = 1e-6) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Working response: use log(u_vec) as approximation to E[log U | Y]
  # Adjust for digamma relationship: psi(alpha) ≈ log(alpha) - 1/(2*alpha) for large alpha
  # So log(alpha) ≈ log(u_vec) + 1/(2*alpha) ≈ log(u_vec) (first order)
  z_working <- log(pmax(u_vec, 0.001))
  
  # Weights: approximate from Gamma variance
  # Var(log U) ≈ trigamma(alpha) ≈ 1/alpha for large alpha
  alpha_old <- as.vector(exp(X %*% beta_old))
  weights <- 1 / pmax(trigamma(pmax(alpha_old, 0.1)), 0.001)
  
  # IRLS iterations
  for (iter in 1:maxit) {
    # Weighted least squares
    W_mat <- diag(weights)
    XWX <- crossprod(X, W_mat %*% X)
    
    # Add small ridge for stability
    XWX <- XWX + 0.001 * diag(p)
    
    # Update
    beta_new <- solve(XWX, crossprod(X, W_mat %*% z_working))
    
    # Check convergence
    if (max(abs(beta_new - beta_old)) < tol) {
      break
    }
    
    # Update for next iteration
    beta_old <- beta_new
    alpha_new <- as.vector(exp(X %*% beta_new))
    weights <- 1 / pmax(trigamma(pmax(alpha_new, 0.1)), 0.001)
  }
  
  # Ensure resulting alphas are positive
  alpha_final <- as.vector(exp(X %*% beta_new))
  if (any(alpha_final < alpha_min)) {
    # Adjust intercept if needed
    adjustment <- log(alpha_min / min(alpha_final))
    beta_new[1] <- beta_new[1] + adjustment
  }
  
  return(as.matrix(beta_new))
}

#' M-Step Block (b): Update Xi (Delta/Tail Index Parameters)
#' 
#' Maximizes Q-function for delta using L-BFGS-B with convexity guarantee.
#' 
#' @param y Vector of observed counts
#' @param Z Design matrix for delta
#' @param v_vec Conditional expectations E[log V_i | Y_i] from E-step
#' @param xi_old Previous xi values
#' @param alpha_vec Current alpha values (fixed in this block)
#' @param rho_vec Current rho values (fixed in this block)
#' @param delta_min Minimum delta for hard boundary (default: 1.05)
#' @return Updated xi vector
#' @details 
#'   Q_xi = sum_i [ -delta_i/(1+rho_i) * log(y_i) - rho_i*delta_i/(1+rho_i) * v_vec[i] 
#'                  - log C(alpha_i, delta_i, rho_i) ]
#'   
#'   With delta_i = 1 + exp(Z_i' xi), this is strictly convex in xi (Proposition 3.2)
#' @references Proposition 3.2 (Convexity of Tail Index Subproblem)
#' @seealso E_step, update_beta, update_zeta
#' @export
update_xi <- function(y, Z, v_vec, xi_old, alpha_vec, rho_vec, 
                      delta_min = 1.05) {
  
  n <- length(y)
  
  # Objective function: negative Q (to minimize)
  neg_Q_xi <- function(xi) {
    delta_vec <- 1 + as.vector(exp(Z %*% xi))
    
    # Enforce hard boundary via parameterization (exp ensures delta > 1)
    # But we need delta >= delta_min, so shift: delta = delta_min + exp(...)
    # Actually: use log(delta - 1) = Z %*% xi, so delta = 1 + exp(Z %*% xi)
    # To enforce delta >= delta_min: use log(delta - delta_min) parametrization in optimization
    
    # For now, use penalty method
    penalty <- sum(pmin(delta_vec - delta_min, 0)^2) * 1e6
    
    # Compute Q components
    eta_vec <- delta_vec / (1 + rho_vec)
    gamma_vec <- rho_vec * delta_vec / (1 + rho_vec)
    
    # Log-likelihood contribution (excluding constant)
    ll <- 0
    for (i in 1:n) {
      # -eta * log(y) - gamma * v - log C
      ll <- ll - eta_vec[i] * log(y[i]) - gamma_vec[i] * v_vec[i]
      
      # Subtract log normalizing constant
      log_C <- compute_normalizing_constant(alpha_vec[i], delta_vec[i], rho_vec[i])
      ll <- ll - log_C
    }
    
    return(-ll + penalty)
  }
  
  # Gradient (numerical or analytical)
  # Use numerical for robustness
  grad_Q <- function(xi) {
    grad <- numDeriv::grad(neg_Q_xi, xi)
    return(grad)
  }
  
  # Optimize using L-BFGS-B with box constraints to keep delta reasonable
  # xi bounds: keep delta in [delta_min, 50] roughly
  xi_lower <- rep(-10, length(xi_old))
  xi_upper <- rep(5, length(xi_old))  # exp(5) + 1 ≈ 149, reasonable upper bound
  
  # Try optimization
  opt_result <- tryCatch({
    optim(xi_old, fn = neg_Q_xi, gr = grad_Q, 
          method = "L-BFGS-B",
          lower = xi_lower, upper = xi_upper,
          control = list(maxit = 100, factr = 1e7))
  }, error = function(e) {
    warning("L-BFGS-B failed in update_xi, using fallback: ", conditionMessage(e))
    list(par = xi_old, value = neg_Q_xi(xi_old), convergence = 99)
  })
  
  # If optimization failed to improve, keep old values
  if (opt_result$convergence > 0 && opt_result$value >= neg_Q_xi(xi_old)) {
    return(xi_old)
  }
  
  # Final projection to ensure hard boundary
  xi_new <- opt_result$par
  delta_new <- 1 + as.vector(exp(Z %*% xi_new))
  
  if (any(delta_new < delta_min - 0.01)) {
    # Adjust to satisfy constraint
    adjustment <- log((delta_min - 1) / min(delta_new - 1))
    xi_new[1] <- xi_new[1] + adjustment
  }
  
  return(as.matrix(xi_new))
}

#' M-Step Block (c): Update Zeta (Rho/Coupling Ratio Parameters)
#' 
#' Maximizes Q-function for rho using L-BFGS-B.
#' Only active for Model A (full heterogeneity).
#' 
#' @param y Vector of observed counts
#' @param W Design matrix for rho
#' @param v_vec Conditional expectations E[log V_i | Y_i] from E-step
#' @param zeta_old Previous zeta values
#' @param alpha_vec Current alpha values (fixed)
#' @param delta_vec Current delta values (updated in previous block)
#' @param rho_min Minimum rho for hard boundary (default: 0.05)
#' @return Updated zeta vector
#' @details 
#'   Q_zeta = sum_i [ -delta_i/(1+rho_i) * log(y_i) + gamma_i * v_vec[i] 
#'                   - log Gamma(gamma_i) ]
#'   where gamma_i = rho_i * delta_i / (1 + rho_i)
#'   
#'   With rho_i = exp(W_i' zeta)
#' @seealso E_step, update_beta, update_xi
#' @export
update_zeta <- function(y, W, v_vec, zeta_old, alpha_vec, delta_vec,
                        rho_min = 0.05) {
  
  n <- length(y)
  
  # Objective function
  neg_Q_zeta <- function(zeta) {
    rho_vec <- as.vector(exp(W %*% zeta))
    
    # Penalty for boundary violation
    penalty <- sum(pmin(rho_vec - rho_min, 0)^2) * 1e6
    
    gamma_vec <- rho_vec * delta_vec / (1 + rho_vec)
    eta_vec <- delta_vec / (1 + rho_vec)
    
    ll <- 0
    for (i in 1:n) {
      # Q contribution
      ll <- ll - eta_vec[i] * log(y[i]) + gamma_vec[i] * v_vec[i] - lgamma(gamma_vec[i])
      
      # Subtract log C (depends on rho through gamma and eta)
      log_C <- compute_normalizing_constant(alpha_vec[i], delta_vec[i], rho_vec[i])
      ll <- ll - log_C
    }
    
    return(-ll + penalty)
  }
  
  # Bounds to keep rho reasonable
  zeta_lower <- rep(-10, length(zeta_old))
  zeta_upper <- rep(5, length(zeta_old))  # exp(5) ≈ 148
  
  # Optimize
  opt_result <- tryCatch({
    optim(zeta_old, fn = neg_Q_zeta, 
          method = "L-BFGS-B",
          lower = zeta_lower, upper = zeta_upper,
          control = list(maxit = 100, factr = 1e7))
  }, error = function(e) {
    warning("L-BFGS-B failed in update_zeta: ", conditionMessage(e))
    list(par = zeta_old, value = neg_Q_zeta(zeta_old), convergence = 99)
  })
  
  # Check improvement
  if (opt_result$convergence > 0 && opt_result$value >= neg_Q_zeta(zeta_old)) {
    return(zeta_old)
  }
  
  # Final projection
  zeta_new <- opt_result$par
  rho_new <- as.vector(exp(W %*% zeta_new))
  
  if (any(rho_new < rho_min - 0.01)) {
    adjustment <- log(rho_min / min(rho_new))
    zeta_new[1] <- zeta_new[1] + adjustment
  }
  
  return(as.matrix(zeta_new))
}

#' M-Step: Block Coordinate Ascent
#' 
#' Executes all three blocks of M-step in sequence.
#' 
#' @param y Vector of observed counts
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @param u_vec From E-step
#' @param v_vec From E-step
#' @param theta_old List with beta, xi, zeta from previous iteration
#' @param model_type "A", "B", "C", or "D"
#' @param delta_min Hard boundary for delta (default: 1.05)
#' @param rho_min Hard boundary for rho (default: 0.05)
#' @return List with updated beta, xi, zeta and individual parameters
#' @export
M_step <- function(y, X, Z, W, u_vec, v_vec, theta_old, model_type,
                   delta_min = 1.05, rho_min = 0.05) {
  
  # Extract old values
  beta_old <- theta_old$beta
  xi_old <- theta_old$xi
  zeta_old <- theta_old$zeta
  
  # Get current individual parameters for use in blocks
  alpha_vec <- as.vector(exp(X %*% beta_old))
  delta_vec <- 1 + as.vector(exp(Z %*% xi_old))
  
  if (model_type == "C") {
    rho_vec <- rep(0, length(y))  # Fixed at 0 for Generalized Zeta
  } else {
    rho_vec <- as.vector(exp(W %*% zeta_old))
  }
  
  # Block (a): Update beta (Model A, C, D)
  if (model_type %in% c("A", "C", "D")) {
    beta_new <- update_beta(y, X, u_vec, beta_old)
    alpha_vec <- as.vector(exp(X %*% beta_new))
  } else {
    beta_new <- beta_old  # Model B: fixed alpha
  }
  
  # Block (b): Update xi (all models)
  xi_new <- update_xi(y, Z, v_vec, xi_old, alpha_vec, rho_vec, delta_min)
  delta_vec <- 1 + as.vector(exp(Z %*% xi_new))
  
  # Block (c): Update zeta (Model A only)
  if (model_type == "A") {
    zeta_new <- update_zeta(y, W, v_vec, zeta_old, alpha_vec, delta_vec, rho_min)
    rho_vec <- as.vector(exp(W %*% zeta_new))
  } else if (model_type == "B") {
    # Common rho: update via profile optimization (simplified)
    zeta_new <- zeta_old  # Keep fixed or use simple update
    rho_vec <- rep(exp(zeta_new[1]), length(y))
  } else {
    zeta_new <- zeta_old  # Model C, D: fixed rho
  }
  
  # Hard boundary projection
  delta_vec <- pmax(delta_vec, delta_min)
  rho_vec <- pmax(rho_vec, rho_min)
  alpha_vec <- pmax(alpha_vec, 0.01)
  
  return(list(
    beta = beta_new,
    xi = xi_new,
    zeta = zeta_new,
    alpha_vec = alpha_vec,
    delta_vec = delta_vec,
    rho_vec = rho_vec
  ))
}

#' Main EM Algorithm for PDigamma Regression (Algorithm 3.1)
#' 
#' Fits PDigamma regression model using EM algorithm with hard boundary projection.
#' 
#' @param y Vector of observed counts (positive integers)
#' @param X Design matrix for alpha (NULL for intercept-only)
#' @param Z Design matrix for delta (NULL for intercept-only)
#' @param W Design matrix for rho (NULL for intercept-only)
#' @param theta_init Initial parameter values from initialize_params
#' @param model_type Model type: "A", "B", "C", or "D"
#' @param delta_min Hard boundary for delta (default: 1.05)
#' @param rho_min Hard boundary for rho (default: 0.05)
#' @param max_iter Maximum EM iterations (default: 1000)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param verbose Print progress (default: TRUE)
#' @return List containing:
#'   \item{coefficients}{Named list with beta, xi, zeta estimates}
#'   \item{theta_list}{Individual parameters (alpha_vec, delta_vec, rho_vec)}
#'   \item{loglik}{Final log-likelihood}
#'   \item{iterations}{Number of EM iterations}
#'   \item{converged}{Logical, TRUE if converged}
#'   \item{history}{Iteration history (if verbose)}
#' @references Algorithm 3.1, Theorem 3.2 (Convergence)
#' @seealso initialize_params, E_step, M_step
#' @export
fit_em <- function(y, X = NULL, Z = NULL, W = NULL, theta_init, model_type = "A",
                   delta_min = 1.05, rho_min = 0.05,
                   max_iter = 1000, tol = 1e-6, verbose = TRUE) {
  
  # Setup design matrices
  n <- length(y)
  if (is.null(X)) X <- matrix(1, n, 1)
  if (is.null(Z)) Z <- matrix(1, n, 1)
  if (is.null(W)) W <- matrix(1, n, 1)
  
  # Initialize from starting values
  theta <- list(
    beta = theta_init$beta_init,
    xi = theta_init$xi_init,
    zeta = theta_init$zeta_init
  )
  
  # Compute initial individual parameters
  alpha_vec <- theta_init$theta_list$alpha_vec
  delta_vec <- theta_init$theta_list$delta_vec
  rho_vec <- theta_init$theta_list$rho_vec
  
  # Initial log-likelihood
  ll_old <- log_likelihood(y, alpha_vec, delta_vec, rho_vec)
  
  if (verbose) {
    message("EM Algorithm: model type ", model_type)
    message("Initial log-likelihood: ", round(ll_old, 4))
  }
  
  # Storage for history
  history <- data.frame(
    iter = integer(),
    loglik = numeric(),
    delta_mean = numeric(),
    param_change = numeric()
  )
  
  # EM iterations
  converged <- FALSE
  for (iter in 1:max_iter) {
    
    # E-step
    E_result <- E_step(y, alpha_vec, delta_vec, rho_vec)
    u_vec <- E_result$u_vec
    v_vec <- E_result$v_vec
    
    # M-step
    M_result <- M_step(y, X, Z, W, u_vec, v_vec, theta, model_type,
                       delta_min, rho_min)
    
    # Update parameters
    theta$beta <- M_result$beta
    theta$xi <- M_result$xi
    theta$zeta <- M_result$zeta
    alpha_vec <- M_result$alpha_vec
    delta_vec <- M_result$delta_vec
    rho_vec <- M_result$rho_vec
    
    # Hard boundary projection (Algorithm 3.1 Step 3)
    delta_vec <- pmax(delta_vec, delta_min)
    if (model_type != "C") {
      rho_vec <- pmax(rho_vec, rho_min)
    }
    
    # Compute new log-likelihood
    ll_new <- log_likelihood(y, alpha_vec, delta_vec, rho_vec)
    
    # Check convergence (Algorithm 3.1 Step 4)
    param_change <- max(
      abs(theta$beta - M_result$beta),
      abs(theta$xi - M_result$xi),
      abs(theta$zeta - M_result$zeta)
    )
    ll_change <- abs(ll_new - ll_old) / (abs(ll_old) + 1e-10)
    
    # Store history
    if (verbose && iter %% 10 == 0) {
      message("Iter ", iter, ": loglik = ", round(ll_new, 4), 
              ", delta_mean = ", round(mean(delta_vec), 3),
              ", change = ", round(param_change, 6))
    }
    
    history <- rbind(history, data.frame(
      iter = iter,
      loglik = ll_new,
      delta_mean = mean(delta_vec),
      param_change = param_change
    ))
    
    # Convergence check
    if (param_change < sqrt(tol) && ll_change < tol) {
      converged <- TRUE
      if (verbose) message("Converged at iteration ", iter)
      break
    }
    
    # Update for next iteration
    ll_old <- ll_new
  }
  
  if (!converged && verbose) {
    warning("EM algorithm did not converge within ", max_iter, " iterations")
  }
  
  # Final result assembly
  result <- list(
    coefficients = list(
      beta = theta$beta,
      xi = theta$xi,
      zeta = theta$zeta
    ),
    theta_list = list(
      alpha_vec = alpha_vec,
      delta_vec = delta_vec,
      rho_vec = rho_vec
    ),
    loglik = ll_new,
    iterations = iter,
    converged = converged,
    model_type = model_type,
    history = if (verbose) history else NULL
  )
  
  class(result) <- "pdigamma_fit"
  return(result)
}

#' Print Method for PDigamma Fit
#' 
#' @param x Object of class pdigamma_fit
#' @param ... Additional arguments
#' @export
print.pdigamma_fit <- function(x, ...) {
  cat("PDigamma Regression Fit (Model ", x$model_type, ")\n", sep = "")
  cat("Iterations:", x$iterations, "| Converged:", x$converged, "\n")
  cat("Final log-likelihood:", round(x$loglik, 4), "\n")
  cat("\nMean parameters:\n")
  cat("  Mean alpha:", round(mean(x$theta_list$alpha_vec), 3), "\n")
  cat("  Mean delta:", round(mean(x$theta_list$delta_vec), 3), "\n")
  cat("  Mean rho:", round(mean(x$theta_list$rho_vec), 3), "\n")
  cat("\nTail index range: [", round(min(x$theta_list$delta_vec), 3), ", ", 
      round(max(x$theta_list$delta_vec), 3), "]\n", sep = "")
}

#' Summary Method for PDigamma Fit
#' 
#' @param object Object of class pdigamma_fit
#' @param ... Additional arguments
#' @export
summary.pdigamma_fit <- function(object, ...) {
  structure(object, class = c("summary.pdigamma_fit", "pdigamma_fit"))
}

#' Print Summary Method
#' @export
print.summary.pdigamma_fit <- function(x, ...) {
  print.pdigamma_fit(x)
  cat("\nCoefficient estimates:\n")
  cat("\nAlpha (beta):\n")
  print(x$coefficients$beta)
  cat("\nDelta (xi):\n")
  print(x$coefficients$xi)
  cat("\nRho (zeta):\n")
  print(x$coefficients$zeta)
}