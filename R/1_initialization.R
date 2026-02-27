# =============================================================================
# 1_initialization.R
# Three-tier initialization strategy for EM algorithm (Algorithm 3.2)
# Includes Hill-based tail index, truncated NB for proximal shape, and 
#' dispersion-based coupling ratio initialization
# =============================================================================

#' Right-Truncated Negative Binomial Regression
#' 
#' Fits NB regression with right truncation to obtain initial alpha estimates.
#' Used for proximal shape parameter initialization.
#' 
#' @param y Vector of counts
#' @param X Design matrix
#' @param truncation_point Upper truncation point (default: floor(max(y)/2))
#' @param maxit Maximum iterations for glm.nb (default: 100)
#' @return List with coefficients (beta) and theta (dispersion parameter)
#' @details Truncation prevents extreme values from distorting initial estimates
#' @seealso initialize_params
#' @export
fit_truncated_nb <- function(y, X, truncation_point = floor(max(y)/2), maxit = 100) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' required for truncated NB initialization")
  }
  
  # Create truncation indicator
  trunc_idx <- y <= truncation_point
  
  if (sum(trunc_idx) < 10) {
    warning("Too few observations below truncation point; using full data")
    trunc_idx <- rep(TRUE, length(y))
  }
  
  y_trunc <- y[trunc_idx]
  X_trunc <- X[trunc_idx, , drop = FALSE]
  
  # Fit NB regression
  fit <- tryCatch({
    MASS::glm.nb(y_trunc ~ X_trunc[, -1, drop = FALSE] - 1,  # Remove intercept, already in X
                 control = glm.control(maxit = maxit))
  }, error = function(e) {
    # Fallback to Poisson if NB fails
    warning("NB initialization failed, using Poisson: ", conditionMessage(e))
    glm(y_trunc ~ X_trunc[, -1, drop = FALSE] - 1, family = poisson)
  })
  
  # Extract coefficients
  beta_init <- coef(fit)
  if (any(is.na(beta_init))) {
    warning("NA coefficients in NB initialization; setting to 0")
    beta_init[is.na(beta_init)] <- 0
  }
  
  # Extract dispersion if available
  theta_disp <- if (inherits(fit, "negbin")) fit$theta else 1
  
  return(list(beta = beta_init, theta = theta_disp, converged = fit$converged))
}

#' Sample Dispersion Index Calculation
#' 
#' Computes local dispersion indices for coupling ratio initialization.
#' 
#' @param y Vector of counts
#' @param W Design matrix for grouping (if NULL, computes global)
#' @return Vector of dispersion indices (sqrt(var/mean))
#' @details Used to initialize rho based on relationship between dispersion and tail thickness
#' @export
compute_dispersion_index <- function(y, W = NULL) {
  if (is.null(W) || ncol(W) == 1) {
    # Global dispersion
    di <- sqrt(var(y) / mean(y))
    return(rep(di, length(y)))
  }
  
  # Group-based dispersion (if W has grouping structure)
  # For continuous covariates, use nearest-neighbor grouping
  n <- length(y)
  di_vec <- numeric(n)
  
  # Simple approach: use quantile-based groups
  n_groups <- min(10, floor(n / 10))
  if (n_groups < 2) {
    di_vec[] <- sqrt(var(y) / mean(y))
    return(di_vec)
  }
  
  # Use first continuous covariate for grouping (if available)
  if (ncol(W) > 1) {
    w_cont <- W[, 2]  # First non-intercept column
    groups <- cut(w_cont, breaks = quantile(w_cont, probs = seq(0, 1, length.out = n_groups + 1)), 
                  include.lowest = TRUE, labels = FALSE)
  } else {
    groups <- rep(1, n)
  }
  
  for (g in unique(groups)) {
    idx <- which(groups == g)
    if (length(idx) > 1) {
      di_vec[idx] <- sqrt(var(y[idx]) / mean(y[idx]))
    } else {
      di_vec[idx] <- 1  # Default for singleton groups
    }
  }
  
  # Ensure finite values
  di_vec[!is.finite(di_vec) | di_vec == 0] <- 1
  
  return(di_vec)
}

#' Three-Tier Initialization Strategy (Algorithm 3.2)
#' 
#' Constructs initial parameter estimates for EM algorithm using:
#' (1) Hill estimator for tail index, (2) Truncated NB for proximal shape,
#' (3) Dispersion index for coupling ratio.
#' 
#' @param y Vector of observed counts
#' @param X Design matrix for alpha (can be NULL)
#' @param Z Design matrix for delta (can be NULL)
#' @param W Design matrix for rho (can be NULL)
#' @param model_type Model type: "A", "B", "C", or "D" (from pre_diagnosis)
#' @param delta_min Minimum delta for hard boundary (default: 1.05)
#' @param rho_min Minimum rho for hard boundary (default: 0.05)
#' @param alpha_min Minimum alpha (default: 0.01)
#' @return List containing:
#'   \item{beta_init}{Initial values for alpha regression coefficients}
#'   \item{xi_init}{Initial values for delta regression coefficients}
#'   \item{zeta_init}{Initial values for rho regression coefficients}
#'   \item{theta_list}{List of individual parameters (alpha_vec, delta_vec, rho_vec)}
#'   \item{delta_global}{Global delta estimate from Hill estimator}
#'   \item{sibuya_flag}{Logical, TRUE if Sibuya boundary detected}
#' @details 
#'   Tier 1: delta_global = min(1 + 1/xi_hill, 3), with Sibuya check
#'   Tier 2: alpha from right-truncated NB (Model A/D) or method of moments (Model B/C)
#'   Tier 3: rho from log(delta_global / max(dispersion_index, 0.1))
#' @references Algorithm 3.2, Section 3.4
#' @seealso pre_diagnosis, fit_em
#' @export
initialize_params <- function(y, X = NULL, Z = NULL, W = NULL, 
                              model_type = "A",
                              delta_min = 1.05, rho_min = 0.05, 
                              alpha_min = 0.01) {
  
  n <- length(y)
  
  # Ensure design matrices exist
  if (is.null(X)) X <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  if (is.null(Z)) Z <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  if (is.null(W)) W <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  
  # ============================================================================
  # TIER 1: Tail Index Initialization (via Hill estimator)
  # ============================================================================
  
  # Compute global Hill estimate
  xi_hill <- hill_estimator(y, k = floor(0.1 * n))
  
  # Convert to delta scale with safeguards
  if (xi_hill > 5) {
    # Sibuya warning: extremely heavy tails
    sibuya_flag <- TRUE
    delta_global <- 1.2  # Safe distance from boundary
    message("Sibuya warning in initialization: xi_hill = ", round(xi_hill, 2), 
            "; setting delta_global = ", delta_global)
  } else {
    sibuya_flag <- FALSE
    delta_global <- 1 + 1 / xi_hill
    # Cap at reasonable value
    delta_global <- min(delta_global, 3)
  }
  
  # Apply hard boundary
  delta_global <- max(delta_global, delta_min)
  
  # Initialize xi (coefficients for delta) based on model type
  if (model_type %in% c("A", "B", "C")) {
    # Free delta: start from global value
    # log(delta - 1) = Z %*% xi => xi = solve(Z'Z) Z' log(delta - 1)
    target <- log(delta_global - 1)
    if (ncol(Z) == 1) {
      xi_init <- matrix(target, nrow = 1)
    } else {
      # Ridge regression for stability
      xi_init <- solve(crossprod(Z) + 0.01 * diag(ncol(Z)), crossprod(Z, rep(target, n)))
    }
  } else {
    # Model D: delta fixed (intercept only, but not free)
    xi_init <- matrix(log(delta_global - 1), nrow = 1)
  }
  
  # ============================================================================
  # TIER 2: Proximal Shape Initialization (alpha)
  # ============================================================================
  
  if (model_type %in% c("A", "D")) {
    # Free alpha: use truncated NB regression
    trunc_point <- floor(max(y) / 2)
    nb_fit <- fit_truncated_nb(y, X, truncation_point = trunc_point)
    beta_init <- matrix(nb_fit$beta, ncol = 1)
    
    # Ensure positive
    alpha_vec <- as.vector(exp(X %*% beta_init))
    alpha_vec <- pmax(alpha_vec, alpha_min)
    beta_init[1] <- beta_init[1] + log(mean(alpha_vec) / alpha_min)  # Adjust intercept
    
  } else {
    # Model B/C: common alpha via method of moments
    # Use relationship: E[Y] ≈ (alpha / gamma) * something (approximate)
    # Simpler: set alpha to make mean roughly match
    mean_y <- mean(y)
    var_y <- var(y)
    
    # Rough initialization: alpha such that Gamma(alpha + 1)/Gamma(alpha) ≈ mean
    # This is just alpha, so set alpha ≈ mean_y / scale_factor
    alpha_global <- max(mean_y / 10, alpha_min)  # Conservative
    
    if (ncol(X) == 1) {
      beta_init <- matrix(log(alpha_global), nrow = 1)
    } else {
      # Start with intercept only
      beta_init <- matrix(0, nrow = ncol(X), ncol = 1)
      beta_init[1] <- log(alpha_global)
    }
  }
  
  # Recompute alpha_vec
  alpha_vec <- as.vector(exp(X %*% beta_init))
  alpha_vec <- pmax(alpha_vec, alpha_min)
  
  # ============================================================================
  # TIER 3: Coupling Ratio Initialization (rho)
  # ============================================================================
  
  if (model_type == "A") {
    # Free rho: based on dispersion index relationship
    di_vec <- compute_dispersion_index(y, W)
    
    # Heuristic: rho proportional to log(delta / dispersion)
    # Higher dispersion -> lower rho (more Zeta-like)
    # Lower dispersion -> higher rho (more Digamma-like)
    rho_target <- delta_global / pmax(di_vec, 0.1)
    rho_target <- pmax(rho_target, rho_min)
    
    log_rho <- log(rho_target)
    
    if (ncol(W) == 1) {
      zeta_init <- matrix(mean(log_rho), nrow = 1)
    } else {
      zeta_init <- solve(crossprod(W) + 0.01 * diag(ncol(W)), crossprod(W, log_rho))
    }
    
  } else if (model_type == "B") {
    # Common rho: preset value
    rho_global <- 1.0  # Balanced between Zeta and Digamma
    zeta_init <- matrix(log(rho_global), nrow = 1)
    
  } else if (model_type == "C") {
    # Generalized Zeta: rho = 0 (fixed)
    zeta_init <- matrix(-Inf, nrow = 1)  # log(0) = -Inf, but will be fixed
    
  } else {
    # Model D: common rho
    rho_global <- 1.0
    zeta_init <- matrix(log(rho_global), nrow = 1)
  }
  
  # Compute initial rho_vec
  if (model_type == "C") {
    rho_vec <- rep(0, n)  # Fixed at boundary
  } else {
    rho_vec <- as.vector(exp(W %*% zeta_init))
    rho_vec <- pmax(rho_vec, rho_min)
  }
  
  # Compute initial delta_vec
  delta_vec <- 1 + exp(as.vector(Z %*% xi_init))
  delta_vec <- pmax(delta_vec, delta_min)
  
  # ============================================================================
  # Final Assembly and Projection
  # ============================================================================
  
  # Hard boundary projection (Algorithm 3.2 Step 4)
  delta_vec <- pmax(delta_vec, delta_min)
  rho_vec <- pmax(rho_vec, rho_min)
  alpha_vec <- pmax(alpha_vec, alpha_min)
  
  # Re-adjust coefficients to match projected values (approximately)
  # This ensures consistency between coefficients and individual parameters
  if (model_type %in% c("A", "B", "C")) {
    xi_init[1] <- xi_init[1] + log(mean(delta_vec) - 1) - mean(log(delta_vec - 1))
  }
  
  theta_list <- list(
    alpha_vec = alpha_vec,
    delta_vec = delta_vec,
    rho_vec = rho_vec
  )
  
  return(list(
    beta_init = beta_init,
    xi_init = xi_init,
    zeta_init = zeta_init,
    theta_list = theta_list,
    delta_global = delta_global,
    sibuya_flag = sibuya_flag,
    model_type = model_type
  ))
}

#' Simple Moment-Based Initialization (Fallback)
#' 
#' Lightweight initialization when NB regression fails or for quick starts.
#' 
#' @param y Vector of counts
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @return List with initial coefficients (same structure as initialize_params)
#' @export
initialize_simple <- function(y, X, Z, W) {
  n <- length(y)
  
  # Simple method of moments
  mean_y <- mean(y)
  var_y <- var(y)
  
  # Delta from coefficient of variation heuristic
  cv2 <- var_y / (mean_y^2)
  delta_init <- 2 + 1 / max(cv2, 0.1)  # Heuristic: higher CV -> lower delta
  
  # Alpha from mean relationship (rough)
  alpha_init <- mean_y / 5
  
  # Rho fixed at 1
  rho_init <- 1
  
  beta_init <- matrix(log(alpha_init), nrow = ncol(X))
  xi_init <- matrix(log(delta_init - 1), nrow = ncol(Z))
  zeta_init <- matrix(log(rho_init), nrow = ncol(W))
  
  list(
    beta_init = beta_init,
    xi_init = xi_init,
    zeta_init = zeta_init,
    theta_list = list(
      alpha_vec = rep(alpha_init, n),
      delta_vec = rep(delta_init, n),
      rho_vec = rep(rho_init, n)
    ),
    delta_global = delta_init,
    sibuya_flag = FALSE,
    model_type = "A"
  )
}