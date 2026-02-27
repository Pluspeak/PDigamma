# =============================================================================
# distributions.R  
# Distribution functions for PDigamma family: RNG, CDF, quantiles, moments
# =============================================================================

#' Random Number Generation for PDigamma Distribution
#' 
#' Generates random variates using rejection sampling or mixture representation.
#' 
#' @param n Number of observations
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @return Vector of n random positive integers
#' @details Uses latent variable representation (Proposition 2.2):
#'   U ~ Gamma(alpha, 1), V ~ Gamma(gamma, 1), lambda = U * V^(-1/eta),
#'   Y|U,V ~ Poisson(lambda)
#' @references Proposition 2.2 (Latent Variable Representation)
#' @export
rPDigamma <- function(n, alpha, delta, rho) {
  if (n < 1 || n != floor(n)) stop("n must be positive integer")
  if (alpha <= 0) stop("alpha must be positive")
  if (delta <= 1) stop("delta must be > 1")
  if (rho <= 0) stop("rho must be positive")
  
  # Transform parameters
  gamma <- rho * delta / (1 + rho)
  eta <- delta / (1 + rho)
  
  # Generate latent variables
  U <- rgamma(n, shape = alpha, rate = 1)
  V <- rgamma(n, shape = gamma, rate = 1)
  
  # Compute Poisson rate
  lambda <- U * V^(-1/eta)
  
  # Generate Poisson counts
  Y <- rpois(n, lambda)
  
  # Ensure positive support (theoretical model is for x >= 1)
  # If Y=0 generated, resample (zero-truncated)
  zero_idx <- which(Y == 0)
  while (length(zero_idx) > 0) {
    U[zero_idx] <- rgamma(length(zero_idx), shape = alpha, rate = 1)
    V[zero_idx] <- rgamma(length(zero_idx), shape = gamma, rate = 1)
    lambda[zero_idx] <- U[zero_idx] * V[zero_idx]^(-1/eta)
    Y[zero_idx] <- rpois(length(zero_idx), lambda[zero_idx])
    zero_idx <- which(Y == 0)
  }
  
  return(Y)
}

#' Cumulative Distribution Function of PDigamma Distribution
#' 
#' Computes P(X <= q) for given quantiles.
#' 
#' @param q Vector of quantiles (positive real numbers)
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @param lower.tail Logical, if TRUE return P(X <= q), else P(X > q) (default: TRUE)
#' @param log.p Logical, return log probabilities (default: FALSE)
#' @return Vector of probabilities
#' @details Uses direct summation of PMF; for large q uses asymptotic approximation
#' @export
pPDigamma <- function(q, alpha, delta, rho, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0)) stop("q must be non-negative")
  
  n <- length(q)
  result <- numeric(n)
  
  for (i in 1:n) {
    if (q[i] < 1) {
      result[i] <- 0
    } else {
      x_vals <- 1:floor(q[i])
      probs <- dPDigamma(x_vals, alpha, delta, rho)
      result[i] <- sum(probs)
    }
  }
  
  if (!lower.tail) result <- 1 - result
  if (log.p) result <- log(result)
  
  return(result)
}

#' Quantile Function of PDigamma Distribution
#' 
#' Computes smallest x such that P(X <= x) >= p.
#' 
#' @param p Vector of probabilities in [0, 1]
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @param lower.tail Logical, if TRUE p = P(X <= x), else p = P(X > x) (default: TRUE)
#' @param log.p Logical, p is given as log(p) (default: FALSE)
#' @return Vector of quantiles (positive integers)
#' @details Uses binary search for efficiency with heavy tails
#' @export
qPDigamma <- function(p, alpha, delta, rho, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  if (any(p < 0 | p > 1)) stop("p must be in [0, 1]")
  
  n <- length(p)
  result <- numeric(n)
  
  # Find upper bound for search
  # Use asymptotic: P(X > x) ~ C * x^(-delta), so x_max ~ (C/(1-p))^(1/delta)
  log_C <- compute_normalizing_constant(alpha, delta, rho)
  # Rough approximation of constant
  C_approx <- exp(log_C) * gamma(delta) / (gamma(alpha + rho*delta/(1+rho)) / gamma(alpha))
  
  for (i in 1:n) {
    if (p[i] <= 0) {
      result[i] <- 1
    } else if (p[i] >= 1) {
      result[i] <- Inf
    } else {
      # Binary search for quantile
      low <- 1
      high <- max(100, ceiling((C_approx / (1 - p[i]))^(1/delta) * 2))
      
      # Expand high if needed
      while (pPDigamma(high, alpha, delta, rho) < p[i] && high < 1e10) {
        high <- high * 2
      }
      
      # Binary search
      while (high - low > 1) {
        mid <- floor((low + high) / 2)
        if (pPDigamma(mid, alpha, delta, rho) < p[i]) {
          low <- mid
        } else {
          high <- mid
        }
      }
      result[i] <- high
    }
  }
  
  return(result)
}

#' Moment Computation for PDigamma Distribution
#' 
#' Computes theoretical moments E[X^r] if they exist.
#' 
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @param order Integer, moment order r (default: 1)
#' @return Scalar, E[X^r] if delta > r + 1, else Inf
#' @details From Proposition 2.2: E[X^r] = E[U^r] * E[V^(-r/eta)]
#'   where U ~ Gamma(alpha, 1), V ~ Gamma(gamma, 1)
#' @references Corollary 2.1 (Moment Existence Criterion)
#' @export
moment_PDigamma <- function(alpha, delta, rho, order = 1) {
  if (order < 0 || order != floor(order)) stop("order must be non-negative integer")
  if (delta <= order + 1) return(Inf)  # Moment does not exist
  
  gamma <- rho * delta / (1 + rho)
  eta <- delta / (1 + rho)
  
  # E[U^r] = Gamma(alpha + r) / Gamma(alpha)
  moment_U <- gamma(alpha + order) / gamma(alpha)
  
  # E[V^(-r/eta)] = Gamma(gamma - r/eta) / Gamma(gamma), requires gamma > r/eta
  if (gamma <= order / eta) return(Inf)
  moment_V <- gamma(gamma - order / eta) / gamma(gamma)
  
  return(moment_U * moment_V)
}

#' Mean of PDigamma Distribution
#' 
#' @inheritParams moment_PDigamma
#' @return Scalar, E[X] if delta > 2, else Inf
#' @export
mean_PDigamma <- function(alpha, delta, rho) {
  moment_PDigamma(alpha, delta, rho, order = 1)
}

#' Variance of PDigamma Distribution
#' 
#' @inheritParams moment_PDigamma
#' @return Scalar, Var(X) if delta > 3, else Inf
#' @export
var_PDigamma <- function(alpha, delta, rho) {
  if (delta <= 3) return(Inf)
  
  mu <- moment_PDigamma(alpha, delta, rho, 1)
  mu2 <- moment_PDigamma(alpha, delta, rho, 2)
  
  return(mu2 - mu^2)
}

#' Survival Function (Tail Probability)
#' 
#' Computes P(X > x) efficiently using asymptotic approximation for large x.
#' 
#' @param x Vector of positive integers
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @param log Logical, return log survival function (default: FALSE)
#' @return Vector of survival probabilities
#' @details For x > 1000, uses asymptotic: P(X > x) ~ C * x^(-delta) / (delta - 1)
#' @references Theorem 2.1 and Corollary 2.2 (Asymptotic Tail Behavior)
#' @export
survival_PDigamma <- function(x, alpha, delta, rho, log = FALSE) {
  if (any(x < 0)) stop("x must be non-negative")
  
  # For small x, use CDF
  # For large x, use asymptotic approximation
  threshold <- 1000
  
  result <- numeric(length(x))
  
  for (i in seq_along(x)) {
    if (x[i] < threshold) {
      result[i] <- 1 - pPDigamma(x[i], alpha, delta, rho)
    } else {
      # Asymptotic: P(X > x) ~ K * x^(-(delta-1)) where K involves gamma functions
      # More precise: from Corollary 2.2
      gamma <- rho * delta / (1 + rho)
      log_C <- compute_normalizing_constant(alpha, delta, rho)
      
      # Leading constant from asymptotic analysis
      log_K <- lgamma(alpha + gamma) - lgamma(alpha) - log_C - log(delta - 1)
      
      result[i] <- exp(log_K - (delta - 1) * log(x[i]))
    }
  }
  
  if (log) return(log(result))
  return(result)
}

#' Random Deviates via Direct PMF (Alternative Method)
#' 
#' Alternative RNG using inverse transform sampling on truncated support.
#' Useful for validation of mixture method.
#' 
#' @param n Number of observations
#' @param alpha Positive scalar, shape parameter
#' @param delta Scalar > 1, total tail index
#' @param rho Positive scalar, coupling ratio
#' @param max_x Maximum value to consider (default: 1e6)
#' @return Vector of n random positive integers
#' @export
rPDigamma_direct <- function(n, alpha, delta, rho, max_x = 1e6) {
  # Compute PMF on truncated support
  x_vals <- 1:max_x
  pmf <- dPDigamma(x_vals, alpha, delta, rho)
  
  # Normalize (slight truncation error)
  pmf <- pmf / sum(pmf)
  
  # Inverse transform sampling
  u <- runif(n)
  cumprobs <- cumsum(pmf)
  
  result <- findInterval(u, c(0, cumprobs))
  
  # Handle boundary case
  if (any(result == 0)) {
    warning("Some samples hit maximum truncation; increase max_x")
    result[result == 0] <- max_x
  }
  
  return(result)
}