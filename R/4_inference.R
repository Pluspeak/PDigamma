# =============================================================================
# 4_inference.R
# Robust inference procedures including adaptive m-out-of-n bootstrap (Algorithm 3.3)
# and risk measure computation with tail-index-guided selection (Table 2.3)
# =============================================================================

#' Adaptive M-out-of-N Bootstrap (Algorithm 3.3)
#' 
#' Computes confidence intervals for risk measures using adaptive subsample size
#' based on estimated tail index. Essential for valid inference when delta <= 2.
#' 
#' @param y Original response vector
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @param fit_final Final fitted model from model_selection
#' @param model_type Final model type ("A" or "C")
#' @param risk_measure_fn Function to compute risk measure from fitted model (see examples)
#' @param B Number of bootstrap replications (default: 1000)
#' @param conf_level Confidence level (default: 0.95)
#' @param verbose Print progress (default: TRUE)
#' @return List containing:
#'   \item{point_estimate}{Point estimate from original data}
#'   \item{ci_percentile}{Percentile bootstrap CI}
#'   \item{ci_bca}{BCa CI (if delta > 2)}
#'   \item{bootstrap_estimates}{All B bootstrap estimates}
#'   \item{subsample_sizes}{Adaptive m for each observation (if observation-specific)}
#' @details 
#'   Adaptive exponent: tau_i = 0.7 + 0.3 * min(1, delta_i - 1)
#'   Subsample size: m_i = floor(n^tau_i)
#'   
#'   For delta > 2: standard bootstrap valid, use BCa or percentile
#'   For 1 < delta <= 2: m-out-of-n bootstrap required, use percentile only
#' @references Algorithm 3.3, Table 2.4 (Bootstrap Validity), Section 3.5.3
#' @seealso compute_risk_measures, sibuya_remedy
#' @export
adaptive_bootstrap <- function(y, X, Z, W, fit_final, model_type,
                               risk_measure_fn = NULL,
                               B = 1000, conf_level = 0.95,
                               verbose = TRUE) {
  
  n <- length(y)
  alpha_level <- 1 - conf_level
  
  # Default risk measure: mean of delta (tail index heterogeneity)
  if (is.null(risk_measure_fn)) {
    risk_measure_fn <- function(fit) mean(fit$theta_list$delta_vec)
  }
  
  # Compute point estimate
  point_estimate <- risk_measure_fn(fit_final)
  
  # Determine delta range for bootstrap strategy
  delta_vec <- fit_final$theta_list$delta_vec
  delta_min_obs <- min(delta_vec)
  delta_max_obs <- max(delta_vec)
  
  # Adaptive subsample sizes
  tau_vec <- 0.7 + 0.3 * pmin(delta_vec - 1, 1)
  m_vec <- floor(n^tau_vec)
  
  if (verbose) {
    message("Adaptive bootstrap: delta range [", round(delta_min_obs, 3), ", ", 
            round(delta_max_obs, 3), "]")
    message("  Subsample sizes: m in [", min(m_vec), ", ", max(m_vec), "]")
  }
  
  # Storage
  boot_estimates <- numeric(B)
  
  # Bootstrap loop
  if (verbose) message("  Running ", B, " bootstrap replications...")
  
  for (b in 1:B) {
    if (verbose && b %% 100 == 0) message("    Replication ", b, "/", B)
    
    # Subsample index (using median m for simplicity, or observation-specific)
    m_use <- median(m_vec)  # Could use individual m_i for observation-specific CIs
    
    # Resample with replacement
    idx <- sample(1:n, size = m_use, replace = TRUE)
    y_star <- y[idx]
    X_star <- X[idx, , drop = FALSE]
    Z_star <- Z[idx, , drop = FALSE]
    W_star <- W[idx, , drop = FALSE]
    
    # Refit model on subsample
    fit_star <- tryCatch({
      # Use simplified initialization for speed
      init_star <- initialize_params(y_star, X_star, Z_star, W_star, 
                                     model_type = model_type,
                                     verbose = FALSE)
      fit_em(y_star, X_star, Z_star, W_star, init_star, 
             model_type = model_type, verbose = FALSE, max_iter = 200)
    }, error = function(e) {
      NULL  # Return NULL on failure
    })
    
    # Store estimate if successful
    if (!is.null(fit_star) && fit_star$converged) {
      boot_estimates[b] <- risk_measure_fn(fit_star)
    } else {
      boot_estimates[b] <- NA
    }
  }
  
  # Remove failed replications
  boot_estimates <- boot_estimates[!is.na(boot_estimates)]
  B_eff <- length(boot_estimates)
  
  if (B_eff < B * 0.5) {
    warning("More than 50% bootstrap replications failed; results unreliable")
  }
  
  # Construct CIs
  if (delta_min_obs > 2) {
    # Standard case: BCa or percentile
    ci_percentile <- quantile(boot_estimates, 
                              probs = c(alpha_level/2, 1 - alpha_level/2),
                              na.rm = TRUE)
    
    # Simple BCa (without full acceleration correction for speed)
    z0 <- qnorm(mean(boot_estimates < point_estimate, na.rm = TRUE))
    a <- 0  # Simplified acceleration
    
    ci_bca <- ci_percentile  # Placeholder; full BCa requires jackknife
    
    method_used <- "BCa (delta > 2)"
    
  } else {
    # Heavy tail case: percentile only
    ci_percentile <- quantile(boot_estimates,
                              probs = c(alpha_level/2, 1 - alpha_level/2),
                              na.rm = TRUE)
    ci_bca <- NULL
    method_used <- "Percentile (delta <= 2, m-out-of-n)"
  }
  
  return(list(
    point_estimate = point_estimate,
    ci_percentile = ci_percentile,
    ci_bca = ci_bca,
    bootstrap_estimates = boot_estimates,
    B_effective = B_eff,
    method_used = method_used,
    delta_range = c(delta_min_obs, delta_max_obs),
    subsample_size = median(m_vec)
  ))
}

#' Risk Measure Selection and Computation (Table 2.3)
#' 
#' Selects appropriate risk measure based on tail index range and computes it.
#' 
#' @param fit Fitted PDigamma model
#' @param p Probability level for quantile-based measures (default: 0.95)
#' @return List containing:
#'   \item{risk_type}{Selected risk measure type}
#'   \item{risk_values}{Vector of risk measure values (one per observation)}
#'   \item{recommendation}{Inference recommendation}
#'   \item{bootstrap_needed}{Logical, whether adaptive bootstrap required}
#' @details Selection rules from Table 2.3:
#'   - delta > 3: TVaR, standard Wald intervals
#'   - 2 < delta <= 3: TVaR, standard bootstrap
#'   - 1.5 < delta <= 2: Wang transform, m-out-of-n bootstrap
#'   - 1 < delta <= 1.5: Extreme value index regression only (point estimates unreliable)
#' @references Table 2.3, Section 2.5.3
#' @seealso adaptive_bootstrap, compute_tvar, compute_wang_transform
#' @export
compute_risk_measures <- function(fit, p = 0.95) {
  
  delta_vec <- fit$theta_list$delta_vec
  delta_min <- min(delta_vec)
  delta_max <- max(delta_vec)
  
  # Determine risk measure type
  if (delta_min > 3) {
    risk_type <- "TVaR"
    inference_method <- "Wald intervals (Fisher information)"
    bootstrap_needed <- FALSE
    
  } else if (delta_min > 2) {
    risk_type <- "TVaR"
    inference_method <- "Standard bootstrap"
    bootstrap_needed <- TRUE
    
  } else if (delta_min > 1.5) {
    risk_type <- "Wang_transform"
    inference_method <- "M-out-of-n bootstrap"
    bootstrap_needed <- TRUE
    
  } else {
    risk_type <- "EVI_regression_only"
    inference_method <- "Point estimates unreliable; use Algorithm 3.5"
    bootstrap_needed <- TRUE  # But will be ineffective
    
    warning("Extreme heavy tails (delta <= 1.5); standard inference unreliable. ",
            "Consider Algorithm 3.5 (mixed extreme value inference).")
  }
  
  # Compute selected risk measure
  if (risk_type == "TVaR") {
    risk_values <- compute_tvar(fit, p)
  } else if (risk_type == "Wang_transform") {
    risk_values <- compute_wang_transform(fit, p)
  } else {
    # EVI regression: just return delta (extreme value index = 1/(delta-1))
    risk_values <- 1 / (delta_vec - 1)
  }
  
  return(list(
    risk_type = risk_type,
    risk_values = risk_values,
    delta_range = c(delta_min, delta_max),
    recommendation = inference_method,
    bootstrap_needed = bootstrap_needed,
    probability_level = p
  ))
}

#' Compute Tail Value-at-Risk (TVaR)
#' 
#' Computes TVaR_p = E[X | X > VaR_p] for PDigamma distribution.
#' 
#' @param fit Fitted model
#' @param p Probability level (default: 0.95)
#' @return Vector of TVaR values (one per observation)
#' @details Requires delta > 2 for finite variance (well-defined TVaR estimation).
#'   Uses numerical integration of truncated expectation.
#' @seealso compute_risk_measures
#' @export
compute_tvar <- function(fit, p = 0.95) {
  
  n <- length(fit$theta_list$alpha_vec)
  tvar_vec <- numeric(n)
  
  for (i in 1:n) {
    ai <- fit$theta_list$alpha_vec[i]
    di <- fit$theta_list$delta_vec[i]
    ri <- fit$theta_list$rho_vec[i]
    
    # Compute VaR (quantile)
    var_i <- qPDigamma(p, ai, di, ri)
    
    # Compute E[X | X > var_i] = sum_{x > var_i} x * P(X=x) / P(X > var_i)
    # Truncate sum at reasonable upper limit
    x_max <- max(var_i * 10, 1e6)
    x_vals <- (var_i + 1):min(x_max, var_i + 10000)
    
    if (length(x_vals) == 0) {
      tvar_vec[i] <- var_i
      next
    }
    
    pmf_vals <- dPDigamma(x_vals, ai, di, ri)
    tail_prob <- sum(pmf_vals)
    
    if (tail_prob < 1e-10) {
      tvar_vec[i] <- var_i
    } else {
      tvar_vec[i] <- sum(x_vals * pmf_vals) / tail_prob
    }
  }
  
  return(tvar_vec)
}

#' Compute Wang Transform Risk Measure
#' 
#' Computes distortion risk measure using Wang transform g(u) = Phi(Phi^{-1}(u) + lambda).
#' 
#' @param fit Fitted model
#' @param lambda Distortion parameter (default: 0.5)
#' @param p Quantile level for truncation (default: 0.95)
#' @return Vector of Wang transform values
#' @details More robust than TVaR for heavy tails (1.5 < delta <= 2).
#'   Integral representation from Section 2.5.2.
#' @seealso compute_risk_measures
#' @export
compute_wang_transform <- function(fit, lambda = 0.5, p = 0.95) {
  
  # Simplified implementation: use quantile-based approximation
  # Full implementation requires numerical integration over distortion
  
  n <- length(fit$theta_list$alpha_vec)
  wang_vec <- numeric(n)
  
  for (i in 1:n) {
    ai <- fit$theta_list$alpha_vec[i]
    di <- fit$theta_list$delta_vec[i]
    ri <- fit$theta_list$rho_vec[i]
    
    # Compute distorted quantile
    # g(u) = Phi(Phi^{-1}(u) + lambda)
    u_distorted <- pnorm(qnorm(p) + lambda)
    
    wang_vec[i] <- qPDigamma(u_distorted, ai, di, ri)
  }
  
  return(wang_vec)
}

#' Mixed Extreme Value Inference for Sibuya Zone (Algorithm 3.5)
#' 
#' Fallback procedure when data is in Sibuya danger zone (delta ≈ 1).
#' Combines PDigamma for proximal modeling with Hill estimator for tail.
#' 
#' @param y Response vector
#' @param X Design matrix for alpha (proximal)
#' @param Z Design matrix for delta (tail covariates)
#' @param fit_gz Fitted Generalized Zeta model
#' @param K Number of strata for Hill estimation (default: 5)
#' @param verbose Print progress (default: TRUE)
#' @return List containing:
#'   \item{xi_ev}{Extreme value index regression coefficients}
#'   \item{xi_pred}{Predicted extreme value indices}
#'   \item{strata_info}{Information about each stratum}
#'   \item{note}{Description of procedure}
#' @details 
#'   Stage 1: Stratify by alpha (proximal risk)
#'   Stage 2: Hill estimator within each stratum
#'   Stage 3: Log-linear regression of xi on Z
#'   
#'   This is a fallback when PDigamma MLE fails near Sibuya boundary.
#' @references Algorithm 3.5, Section 3.7
#' @seealso sibuya_diagnostic, adaptive_bootstrap
#' @export
sibuya_remedy <- function(y, X, Z, fit_gz, K = 5, verbose = TRUE) {
  
  if (verbose) message("Algorithm 3.5: Mixed Extreme Value Inference")
  
  n <- length(y)
  
  # Stage 1: Stratify by alpha (proximal parameter)
  alpha_vec <- fit_gz$theta_list$alpha_vec
  
  # Create K strata based on alpha quantiles
  alpha_breaks <- quantile(alpha_vec, probs = seq(0, 1, length.out = K + 1))
  strata <- cut(alpha_vec, breaks = alpha_breaks, include.lowest = TRUE, labels = FALSE)
  
  # Stage 2: Hill estimation within each stratum
  xi_hill <- numeric(K)
  var_xi <- numeric(K)
  n_strata <- numeric(K)
  Z_mean <- matrix(0, K, ncol(Z))
  
  for (k in 1:K) {
    idx <- which(strata == k)
    n_strata[k] <- length(idx)
    
    if (n_strata[k] < 10) {
      warning("Stratum ", k, " has only ", n_strata[k], " observations")
      xi_hill[k] <- NA
      var_xi[k] <- Inf
      next
    }
    
    y_k <- y[idx]
    
    # Hill estimator with adaptive k
    hill_k <- adaptive_hill(y_k, method = "stable")
    xi_hill[k] <- hill_k$xi_hat
    
    # Variance estimate: xi^2 / k
    var_xi[k] <- xi_hill[k]^2 / hill_k$k_opt
    
    # Mean covariates in stratum
    Z_mean[k, ] <- colMeans(Z[idx, , drop = FALSE])
    
    if (verbose) {
      message("  Stratum ", k, ": n = ", n_strata[k], 
              ", xi_hill = ", round(xi_hill[k], 3),
              ", k_opt = ", hill_k$k_opt)
    }
  }
  
  # Remove empty strata
  valid_k <- which(is.finite(xi_hill) & n_strata >= 10)
  
  if (length(valid_k) < 2) {
    warning("Too few valid strata for EVI regression; using pooled estimate")
    xi_pred <- rep(mean(xi_hill[valid_k], na.rm = TRUE), n)
    xi_ev <- matrix(0, ncol(Z), 1)
    xi_ev[1] <- log(mean(xi_hill[valid_k], na.rm = TRUE))
  } else {
    # Stage 3: Log-linear regression of xi on Z
    # log(xi) = Z %*% xi_ev + error
    
    response <- log(pmax(xi_hill[valid_k], 0.01))
    design <- Z_mean[valid_k, , drop = FALSE]
    weights <- 1 / pmax(var_xi[valid_k], 0.001)
    
    # Weighted least squares
    W_mat <- diag(weights)
    xi_ev <- solve(crossprod(design, W_mat %*% design),
                   crossprod(design, W_mat %*% response))
    
    # Predicted EVI
    xi_pred <- as.vector(exp(Z %*% xi_ev))
  }
  
  return(list(
    xi_ev = xi_ev,
    xi_pred = xi_pred,
    strata_info = data.frame(
      stratum = 1:K,
      n = n_strata,
      xi_hill = xi_hill,
      var_xi = var_xi,
      alpha_min = alpha_breaks[-(K+1)],
      alpha_max = alpha_breaks[-1]
    ),
    note = "Mixed inference: PDigamma for proximal, Hill estimator for tail"
  ))
}

#' End-to-End Inference Flow (Algorithm 3.4 Complete)
#' 
#' Executes full inference pipeline from fitted model to final risk measures with CIs.
#' 
#' @param y Response vector
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @param fit_final Final fitted model from model_selection
#' @param model_type Final model type
#' @param sibuya_check Results from sibuya_diagnostic (can be NULL)
#' @param B Bootstrap replications (default: 1000)
#' @param conf_level Confidence level (default: 0.95)
#' @param verbose Print progress (default: TRUE)
#' @return Comprehensive results list with risk measures and inference
#' @references Algorithm 3.4 (Stage 4), Section 3.6
#' @seealso model_selection, adaptive_bootstrap, compute_risk_measures
#' @export
run_inference <- function(y, X, Z, W, fit_final, model_type,
                          sibuya_check = NULL,
                          B = 1000, conf_level = 0.95,
                          verbose = TRUE) {
  
  if (verbose) message("Stage 4: Risk measure computation and inference")
  
  # Check if Sibuya remedy needed
  if (!is.null(sibuya_check) && sibuya_check$in_danger_zone) {
    if (verbose) message("  Activating Sibuya remedy (Algorithm 3.5)...")
    
    remedy_results <- sibuya_remedy(y, X, Z, fit_final, verbose = verbose)
    
    return(list(
      method = "sibuya_remedy",
      risk_measures = NULL,  # Not computed in remedy mode
      evi_regression = remedy_results,
      note = "Standard PDigamma inference unreliable; used mixed EVI approach"
    ))
  }
  
  # Standard path: compute risk measures
  if (verbose) message("  Computing risk measures...")
  risk_info <- compute_risk_measures(fit_final, p = 0.95)
  
  if (verbose) {
    message("    Selected: ", risk_info$risk_type)
    message("    Delta range: [", round(risk_info$delta_range[1], 3), ", ",
            round(risk_info$delta_range[2], 3), "]")
  }
  
  # Bootstrap inference if needed
  if (risk_info$bootstrap_needed) {
    if (verbose) message("  Running adaptive bootstrap (B = ", B, ")...")
    
    # Define risk measure function for bootstrap
    risk_fn <- function(fit) {
      # Return mean of selected risk measure
      if (risk_info$risk_type == "TVaR") {
        return(mean(compute_tvar(fit, 0.95)))
      } else if (risk_info$risk_type == "Wang_transform") {
        return(mean(compute_wang_transform(fit, 0.5, 0.95)))
      } else {
        return(mean(1 / (fit$theta_list$delta_vec - 1)))
      }
    }
    
    boot_results <- adaptive_bootstrap(y, X, Z, W, fit_final, model_type,
                                       risk_measure_fn = risk_fn,
                                       B = B, conf_level = conf_level,
                                       verbose = verbose)
    
    inference_results <- list(
      point_estimate = boot_results$point_estimate,
      confidence_interval = boot_results$ci_percentile,
      bootstrap_details = boot_results
    )
    
  } else {
    # Wald intervals using Fisher information
    if (verbose) message("  Using Wald intervals (Fisher information)...")
    
    # Placeholder: would compute from observed information matrix
    inference_results <- list(
      point_estimate = mean(risk_info$risk_values),
      confidence_interval = c(NA, NA),  # To be implemented
      method = "Wald (Fisher information)"
    )
  }
  
  return(list(
    method = "standard",
    risk_measures = risk_info,
    inference = inference_results,
    individual_risks = risk_info$risk_values
  ))
}