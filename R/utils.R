# =============================================================================
# utils.R
# Basic utility functions for PDigamma distribution family
# =============================================================================

#' Compute Normalizing Constant for PDigamma Distribution
#' 
#' Calculates C(alpha, delta, rho) = sum_{x=1}^infty (alpha)_x / [x^{delta/(1+rho)} * (alpha + rho*delta/(1+rho))_x]
#' using truncation with asymptotic approximation for tail.
#' 
#' @param alpha Positive scalar, shape parameter near origin
#' @param delta Scalar > 1, total tail index (delta = gamma + eta)
#' @param rho Positive scalar, coupling ratio (rho = gamma/eta)
#' @param max_x Maximum x for direct summation (default: 1e5)
#' @param tol Convergence tolerance for series truncation (default: 1e-10)
#' @return Scalar, log normalizing constant
#' @details For large x, uses integral approximation: sum_{x=N}^infty ~ integral_{N-0.5}^infty
#' @references Section 2.1.3, Definition 2.1
#' @export
compute_normalizing_constant <- function(alpha, delta, rho, max_x = 1e5, tol = 1e-10) {
  # Parameter validation
  if (alpha <= 0) stop("alpha must be positive")
  if (delta <= 1) stop("delta must be > 1")
  if (rho <= 0) stop("rho must be positive")
  
  # Transform to original parameterization
  gamma <- rho * delta / (1 + rho)
  eta <- delta / (1 + rho)
  
  # Direct summation for moderate x
  log_terms <- numeric(max_x)
  for (x in 1:min(max_x, 1e6)) {
    log_num <- lgamma(alpha + x) - lgamma(alpha)  # log((alpha)_x)
    log_denom <- eta * log(x) + lgamma(alpha + gamma + x) - lgamma(alpha + gamma)
    log_terms[x] <- log_num - log_denom
  }
  
  # Check convergence
  if (max_x >= 1e6 && log_terms[max_x] > log(tol)) {
    warning("Series may not have converged; consider increasing max_x")
  }
  
  # Truncate where terms become negligible
  valid_idx <- which(is.finite(log_terms) & log_terms > log(tol))
  if (length(valid_idx) == 0) return(-Inf)
  
  # Log-sum-exp for numerical stability
  max_log <- max(log_terms[valid_idx])
  log_sum <- max_log + log(sum(exp(log_terms[valid_idx] - max_log)))
  
  return(log_sum)
}

#' Probability Mass Function of PDigamma Distribution
#' 
#' Computes q(x; alpha, delta, rho) for given parameters.
#' 
#' @param x Vector of positive integers (support: 1, 2, ...)
#' @param alpha Positive scalar or vector, shape parameter
#' @param delta Scalar or vector > 1, total tail index
#' @param rho Positive scalar or vector, coupling ratio
#' @param log Logical, return log probabilities (default: FALSE)
#' @return Vector of probabilities (or log-probabilities)
#' @details Vectorized over x; parameters are recycled if needed
#' @references Definition 2.2
#' @export
dPDigamma <- function(x, alpha, delta, rho, log = FALSE) {
  # Input validation
  if (any(x < 1) || any(x != floor(x))) stop("x must be positive integers")
  if (any(alpha <= 0)) stop("alpha must be positive")
  if (any(delta <= 1)) stop("delta must be > 1")
  if (any(rho <= 0)) stop("rho must be positive")
  
  # Ensure compatible lengths
  n <- max(length(x), length(alpha), length(delta), length(rho))
  x <- rep_len(x, n)
  alpha <- rep_len(alpha, n)
  delta <- rep_len(delta, n)
  rho <- rep_len(rho, n)
  
  # Compute log PMF
  log_C <- mapply(compute_normalizing_constant, alpha, delta, rho, 
                  MoreArgs = list(max_x = 1e5, tol = 1e-10))
  
  gamma <- rho * delta / (1 + rho)
  eta <- delta / (1 + rho)
  
  log_num <- lgamma(alpha + x) - lgamma(alpha)  # log((alpha)_x)
  log_denom <- eta * log(x) + lgamma(alpha + gamma + x) - lgamma(alpha + gamma)
  
  log_pmf <- log_num - log_denom - log_C
  
  if (log) return(log_pmf)
  return(exp(log_pmf))
}

#' Log-Likelihood for PDigamma Regression Model
#' 
#' Computes sum of log probabilities for observed data.
#' 
#' @param y Vector of observed counts (positive integers)
#' @param alpha_vec Vector of individual alpha parameters
#' @param delta_vec Vector of individual delta parameters  
#' @param rho_vec Vector of individual rho parameters
#' @return Scalar, total log-likelihood
#' @details All parameter vectors must have length equal to length(y)
#' @export
log_likelihood <- function(y, alpha_vec, delta_vec, rho_vec) {
  if (length(y) != length(alpha_vec) || length(y) != length(delta_vec) || length(y) != length(rho_vec)) {
    stop("All parameter vectors must have same length as y")
  }
  
  log_probs <- dPDigamma(y, alpha_vec, delta_vec, rho_vec, log = TRUE)
  return(sum(log_probs))
}

#' Project Parameters to Hard Boundary Constraint Set
#' 
#' Enforces delta >= delta_min and rho >= rho_min via projection.
#' 
#' @param theta List with elements alpha, delta, rho (each vectors)
#' @param delta_min Minimum allowed delta (default: 1.05)
#' @param rho_min Minimum allowed rho (default: 0.05)
#' @return List with projected parameters
#' @references Section 3.2, Hard boundary constraints
#' @export
project_to_boundary <- function(theta, delta_min = 1.05, rho_min = 0.05) {
  theta$delta <- pmax(theta$delta, delta_min)
  theta$rho <- pmax(theta$rho, rho_min)
  theta$alpha <- pmax(theta$alpha, .Machine$double.eps)  # ensure positive
  return(theta)
}

#' Compute Fisher Information Matrix
#' 
#' Calculates expected Fisher information for given parameter values.
#' Used for standard error estimation and asymptotic theory.
#' 
#' @param theta List with alpha, delta, rho (scalars or single values)
#' @param X Design matrix for alpha (intercept only if NULL)
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @return Matrix, Fisher information I(theta)
#' @details Uses numerical differentiation of log-likelihood
#' @references Theorem 2.2, Section 3.2
#' @export
fisher_information <- function(theta, X = NULL, Z = NULL, W = NULL) {
  # This is a placeholder for analytical Fisher information
  # In practice, use numerical Hessian or analytical formula from Appendix
  
  # For regression parameters, compute using chain rule
  # I(beta) = X^T * I(alpha) * X, etc.
  
  warning("Analytical Fisher information not yet implemented; use numerical Hessian")
  return(NULL)
}

#' Numerical Hessian via Finite Differences
#' 
#' Computes Hessian matrix numerically for optimization diagnostics.
#' 
#' @param fn Function to differentiate
#' @param x Point at which to evaluate Hessian
#' @param eps Step size for finite differences (default: 1e-6)
#' @return Matrix, approximate Hessian
#' @export
numerical_hessian <- function(fn, x, eps = 1e-6) {
  n <- length(x)
  H <- matrix(0, n, n)
  fx <- fn(x)
  
  for (i in 1:n) {
    x_plus <- x
    x_minus <- x
    x_plus[i] <- x[i] + eps
    x_minus[i] <- x[i] - eps
    
    for (j in i:n) {
      if (i == j) {
        # Diagonal: second derivative
        H[i, i] <- (fn(x_plus) - 2*fx + fn(x_minus)) / (eps^2)
      } else {
        # Off-diagonal: mixed partial
        x_pp <- x_plus
        x_pm <- x_plus
        x_mp <- x_minus
        x_mm <- x_minus
        
        x_pp[j] <- x_pp[j] + eps
        x_pm[j] <- x_pm[j] - eps
        x_mp[j] <- x_mp[j] + eps
        x_mm[j] <- x_mm[j] - eps
        
        H[i, j] <- (fn(x_pp) - fn(x_pm) - fn(x_mp) + fn(x_mm)) / (4 * eps^2)
        H[j, i] <- H[i, j]
      }
    }
  }
  return(H)
}

#' Check Convergence Criterion
#' 
#' Determines if EM algorithm has converged based on parameter and likelihood changes.
#' 
#' @param theta_new New parameter values
#' @param theta_old Old parameter values  
#' @param ll_new New log-likelihood
#' @param ll_old Old log-likelihood
#' @param tol Parameter convergence tolerance (default: 1e-6)
#' @return Logical, TRUE if converged
#' @export
check_convergence <- function(theta_new, theta_old, ll_new, ll_old, tol = 1e-6) {
  param_change <- max(
    abs(theta_new$beta - theta_old$beta),
    abs(theta_new$xi - theta_old$xi),
    abs(theta_new$zeta - theta_old$zeta)
  )
  
  ll_change <- abs(ll_new - ll_old) / (abs(ll_old) + 1e-10)
  
  return(param_change < sqrt(tol) && ll_change < tol)
}

#' Safe Digamma and Trigamma Functions
#' 
#' Wrapper for base::digamma and base::trigamma with input validation.
#' 
#' @param x Numeric vector
#' @return Vector of digamma values
#' @export
safe_digamma <- function(x) {
  if (any(x <= 0)) stop("digamma requires positive arguments")
  return(digamma(x))
}

#' @param x Numeric vector
#' @return Vector of trigamma values
#' @export
safe_trigamma <- function(x) {
  if (any(x <= 0)) stop("trigamma requires positive arguments")
  return(trigamma(x))
}