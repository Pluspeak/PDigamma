# =============================================================================
# 0_data_prep.R
# Data preprocessing and preliminary diagnostics for PDigamma regression
# Includes Hill estimator, model pre-selection, and Sibuya warning detection
# =============================================================================

#' Hill Estimator for Tail Index
#' 
#' Computes the Hill estimator of the tail index xi = 1/(delta - 1).
#' Used for initialization of delta parameter and Sibuya boundary detection.
#' 
#' @param y Vector of positive observations
#' @param k Number of upper order statistics to use (default: floor(0.1*n))
#' @return Scalar, Hill estimate of tail index xi
#' @details The Hill estimator is: xi_hat = (1/k) * sum_{j=1}^k log(Y_{(n-j+1)}) - log(Y_{(n-k)})
#'   where Y_{(1)} <= ... <= Y_{(n)} are order statistics.
#' @references Hill (1975), Section 3.4 (Initialization Strategy)
#' @seealso initialize_params
#' @export
hill_estimator <- function(y, k = floor(0.1 * length(y))) {
  n <- length(y)
  
  if (k < 1 || k >= n) stop("k must be between 1 and n-1")
  if (any(y <= 0)) stop("y must be positive")
  
  # Sort in ascending order
  y_sorted <- sort(y)
  
  # Upper order statistics
  y_k <- y_sorted[n - k]  # (n-k)th order statistic
  y_upper <- y_sorted[(n - k + 1):n]  # k largest observations
  
  # Hill estimator
  xi_hat <- mean(log(y_upper)) - log(y_k)
  
  # Ensure positive
  xi_hat <- max(xi_hat, 0.01)
  
  return(xi_hat)
}

#' Adaptive Hill Estimator with Threshold Selection
#' 
#' Selects optimal k using diagnostic plots or stability criterion.
#' 
#' @param y Vector of positive observations
#' @param k_range Vector of k values to evaluate (default: seq(5, floor(n/2), by=5))
#' @param method Selection method: "stable" (most stable region) or "highest" (highest plateau)
#' @return List with elements: xi_hat (estimate), k_opt (optimal k), xi_path (estimates for all k)
#' @details Stability method finds region where successive estimates vary least
#' @export
adaptive_hill <- function(y, k_range = NULL, method = "stable") {
  n <- length(y)
  
  if (is.null(k_range)) {
    k_range <- seq(5, floor(n/2), by = max(1, floor(n/20)))
  }
  
  xi_path <- sapply(k_range, function(k) hill_estimator(y, k))
  
  if (method == "stable") {
    # Find most stable region (minimum variance in moving window)
    window <- max(3, floor(length(k_range) / 10))
    variances <- sapply(1:(length(xi_path) - window + 1), function(i) {
      var(xi_path[i:(i + window - 1)])
    })
    opt_idx <- which.min(variances) + floor(window / 2)
    k_opt <- k_range[min(opt_idx, length(k_range))]
  } else {
    # Use highest plateau (median of upper half)
    k_opt <- k_range[which.max(xi_path)]
  }
  
  xi_hat <- hill_estimator(y, k_opt)
  
  return(list(xi_hat = xi_hat, k_opt = k_opt, xi_path = xi_path, k_range = k_range))
}

#' Pre-Diagnosis and Model Pre-Selection (Stage 0)
#' 
#' Performs preliminary diagnostics on data and selects initial model type.
#' Detects Sibuya boundary cases and determines appropriate model complexity.
#' 
#' @param y Vector of observed counts (positive integers)
#' @param X Design matrix for alpha (can be NULL for intercept-only)
#' @param Z Design matrix for delta (can be NULL for intercept-only)  
#' @param W Design matrix for rho (can be NULL for intercept-only)
#' @param delta_min Minimum delta threshold for Sibuya warning (default: 1.05)
#' @return List containing:
#'   \item{model_type}{Selected model: "A" (full), "B" (tail-focused), "C" (Generalized Zeta), or "D" (proximal-tail separation)}
#'   \item{sibuya_warning}{Logical, TRUE if data suggests infinite mean (delta close to 1)}
#'   \item{xi_hill}{Global Hill estimate of tail index}
#'   \item{delta_global}{Initial global delta estimate}
#'   \item{covariate_richness}{Data frame with p and n for each parameter}
#' @details Selection rules:
#'   - If xi_hill > 10 (delta < 1.1): Sibuya warning, select Model B
#'   - If total covariates p < 5: select Model B (parsimonious)
#'   - Otherwise: select Model A (full heterogeneity)
#' @references Algorithm 3.4 (Stage 0), Section 3.6
#' @seealso initialize_params, fit_em
#' @export
pre_diagnosis <- function(y, X = NULL, Z = NULL, W = NULL, delta_min = 1.05) {
  
  # Validate data
  if (any(y < 1) || any(y != floor(y))) {
    stop("y must be positive integers (1, 2, ...); apply zero-truncation if needed")
  }
  
  n <- length(y)
  
  # Determine dimensions
  p_alpha <- if (is.null(X)) 1 else ncol(X)
  p_delta <- if (is.null(Z)) 1 else ncol(Z)
  p_rho <- if (is.null(W)) 1 else ncol(W)
  p_total <- p_alpha + p_delta + p_rho
  
  # Compute global Hill estimate
  hill_result <- adaptive_hill(y)
  xi_hill <- hill_result$xi_hat
  
  # Convert to delta scale: xi = 1/(delta - 1) => delta = 1 + 1/xi
  delta_global <- 1 + 1 / xi_hill
  
  # Cap at reasonable value for initialization
  delta_global <- min(delta_global, 3)
  
  # Sibuya boundary detection
  sibuya_warning <- (xi_hill > 10) || (delta_global < 1.2)
  
  # Model pre-selection
  if (sibuya_warning) {
    model_type <- "B"  # Tail-focused, avoid complex modeling near boundary
    message("Sibuya warning: data suggests infinite mean (delta ~ 1); selecting Model B")
  } else if (p_total < 5) {
    model_type <- "B"  # Limited covariates, use parsimonious model
    message("Limited covariates (p = ", p_total, "); selecting Model B")
  } else {
    model_type <- "A"  # Full model
  }
  
  # Covariate richness info
  cov_info <- data.frame(
    parameter = c("alpha", "delta", "rho"),
    p = c(p_alpha, p_delta, p_rho),
    n = n,
    ratio = c(p_alpha, p_delta, p_rho) / n
  )
  
  return(list(
    model_type = model_type,
    sibuya_warning = sibuya_warning,
    xi_hill = xi_hill,
    delta_global = delta_global,
    covariate_richness = cov_info,
    hill_details = hill_result
  ))
}

#' Zero-Truncation Handler
#' 
#' Prepares data by handling zero counts if present.
#' PDigamma has support on {1, 2, ...}, so zeros must be truncated or shifted.
#' 
#' @param y Vector of counts (non-negative integers)
#' @param method Handling method: "truncate" (remove zeros) or "shift" (add 1 to all)
#' @param verbose Logical, print information about truncation (default: TRUE)
#' @return List with processed y and truncation info
#' @details For true PDigamma data, zeros should not occur. Their presence suggests
#'   either zero-inflation or model misspecification.
#'   @export
handle_zeros <- function(y, method = "truncate", verbose = TRUE) {
  n_orig <- length(y)
  n_zero <- sum(y == 0)
  
  if (n_zero == 0) {
    return(list(y = y, was_truncated = FALSE, n_removed = 0))
  }
  
  if (verbose) {
    warning(n_zero, " zero values found (", round(100*n_zero/n_orig, 1), 
            "%); applying ", method, " method")
  }
  
  if (method == "truncate") {
    y_new <- y[y > 0]
    if (verbose) message("Removed ", n_zero, " zeros; remaining n = ", length(y_new))
  } else if (method == "shift") {
    y_new <- y + 1
    if (verbose) message("Shifted all values by +1")
  } else {
    stop("method must be 'truncate' or 'shift'")
  }
  
  return(list(
    y = y_new,
    was_truncated = (method == "truncate"),
    n_removed = if (method == "truncate") n_zero else 0,
    original_n = n_orig
  ))
}

#' Design Matrix Constructor
#' 
#' Safely constructs design matrices from formulas or data frames.
#' 
#' @param formula Formula object or NULL (for intercept-only)
#' @param data Data frame containing covariates
#' @param param_name Name of parameter (for error messages)
#' @return Design matrix (including intercept if formula includes one)
#' @details If formula is NULL, returns intercept-only matrix
#' @export
construct_design_matrix <- function(formula, data, param_name = "parameter") {
  if (is.null(formula)) {
    # Intercept only
    return(matrix(1, nrow = nrow(data), ncol = 1, 
                  dimnames = list(NULL, "(Intercept)")))
  }
  
  if (!inherits(formula, "formula")) {
    stop(param_name, " formula must be a formula object or NULL")
  }
  
  # Use model.matrix with safeguards
  mf <- model.frame(formula, data, na.action = na.omit)
  X <- model.matrix(formula, mf)
  
  return(X)
}

#' Data Validation for PDigamma Regression
#' 
#' Comprehensive validation of all inputs before model fitting.
#' 
#' @param y Response vector
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @param delta_min Minimum allowed delta
#' @param rho_min Minimum allowed rho
#' @return Validated and potentially modified list of inputs
#' @details Checks: dimensions, NA handling, collinearity in design matrices,
#'   and warns about potential numerical issues
#' @export
validate_inputs <- function(y, X, Z, W, delta_min = 1.05, rho_min = 0.05) {
  n <- length(y)
  
  # Check dimensions match
  if (nrow(X) != n) stop("X must have ", n, " rows")
  if (nrow(Z) != n) stop("Z must have ", n, " rows")
  if (nrow(W) != n) stop("W must have ", n, " rows")
  
  # Check for NAs
  if (any(is.na(y))) stop("y contains NA values")
  if (any(is.na(X))) stop("X contains NA values")
  if (any(is.na(Z))) stop("Z contains NA values")
  if (any(is.na(W))) stop("W contains NA values")
  
  # Check for collinearity (condition number)
  check_collinearity <- function(M, name) {
    if (ncol(M) > 1) {
      # Remove intercept for collinearity check
      M_no_int <- M[, -1, drop = FALSE]
      if (ncol(M_no_int) > 0) {
        kappa <- kappa(M_no_int, exact = TRUE, norm = TRUE)
        if (kappa > 1e10) {
          warning(name, " design matrix is highly collinear (kappa = ", 
                  format(kappa, digits = 3), ")")
        }
      }
    }
  }
  
  check_collinearity(X, "alpha (X)")
  check_collinearity(Z, "delta (Z)")
  check_collinearity(W, "rho (W)")
  
  # Check response distribution
  y_summary <- summary(y)
  if (y_summary["Mean"] / y_summary["Median"] > 10) {
    warning("Response has high mean/median ratio; check for extreme outliers")
  }
  
  # Check for heavy tail indicators
  y_sorted <- sort(y, decreasing = TRUE)
  top_ratio <- y_sorted[1] / y_sorted[floor(n/10)]
  if (top_ratio > 100) {
    message("Note: Heavy-tailed appearance (top 1% / top 10% ratio = ", 
            round(top_ratio, 1), ")")
  }
  
  return(list(y = y, X = X, Z = Z, W = W, n = n))
}

#' Summary Statistics for Count Data
#' 
#' Computes descriptive statistics useful for model diagnostics.
#' 
#' @param y Vector of counts
#' @return List of summary statistics
#' @export
summary_counts <- function(y) {
  list(
    n = length(y),
    mean = mean(y),
    var = var(y),
    dispersion_index = var(y) / mean(y),  # > 1 indicates overdispersion
    min = min(y),
    max = max(y),
    median = median(y),
    quantiles = quantile(y, probs = c(0.9, 0.95, 0.99, 0.999)),
    zeros = sum(y == 0),
    unique_values = length(unique(y))
  )
}