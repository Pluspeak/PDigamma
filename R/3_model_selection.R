# =============================================================================
# 3_model_selection.R
# Model selection procedures including boundary likelihood ratio test for 
# Generalized Zeta submodel and Sibuya boundary diagnostics
# =============================================================================

#' Fit Generalized Zeta Submodel (rho = 0)
#' 
#' Fits restricted model with rho fixed at 0, corresponding to limit as rho -> 0+.
#' This is the boundary case of PDigamma family.
#' 
#' @param y Vector of observed counts
#' @param X Design matrix for alpha (NULL for intercept-only)
#' @param Z Design matrix for delta (NULL for intercept-only)
#' @param theta_init Initial values (typically from full model fit or initialization)
#' @param delta_min Hard boundary for delta (default: 1.05)
#' @param max_iter Maximum EM iterations (default: 1000)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param verbose Print progress (default: FALSE)
#' @return Fitted model object of class pdigamma_fit with rho fixed at 0
#' @details 
#'   Generalized Zeta distribution: q(x) = (alpha)_x / [zeta(delta; alpha) * Gamma(alpha+1) * x^delta]
#'   where zeta(delta; alpha) is the Hurwitz zeta function.
#'   
#'   This is the boundary case rho = 0 of PDigamma family (Section 2.4.1).
#' @references Section 2.4.1 (Generalized Zeta Boundary), Proposition 3.5
#' @seealso boundary_lrt, fit_em
#' @export
fit_generalized_zeta <- function(y, X = NULL, Z = NULL, theta_init = NULL,
                                 delta_min = 1.05, max_iter = 1000, 
                                 tol = 1e-6, verbose = FALSE) {
  
  n <- length(y)
  if (is.null(X)) X <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  if (is.null(Z)) Z <- matrix(1, n, 1, dimnames = list(NULL, "(Intercept)"))
  
  # Use simple initialization if not provided
  if (is.null(theta_init)) {
    init <- initialize_simple(y, X, Z, matrix(1, n, 1))
    init$model_type <- "C"
    init$zeta_init <- matrix(-Inf, nrow = 1)  # log(0)
    init$theta_list$rho_vec <- rep(0, n)
    theta_init <- init
  }
  
  # Override to ensure rho = 0
  theta_init$theta_list$rho_vec <- rep(0, n)
  theta_init$model_type <- "C"
  
  # Fit with modified E_step and M_step that fix rho = 0
  # Reuse fit_em but with model_type = "C"
  fit <- fit_em(y, X, Z, W = NULL, theta_init, model_type = "C",
                delta_min = delta_min, rho_min = 0.05,  # rho_min irrelevant
                max_iter = max_iter, tol = tol, verbose = verbose)
  
  fit$model_type <- "C"
  fit$note <- "Generalized Zeta submodel (rho = 0)"
  
  return(fit)
}

#' Boundary Likelihood Ratio Test for Generalized Zeta Submodel (Proposition 3.5)
#' 
#' Tests H0: rho = 0 (Generalized Zeta) vs H1: rho > 0 (full PDigamma).
#' 
#' @param loglik_full Log-likelihood of full model (Model A)
#' @param loglik_gz Log-likelihood of Generalized Zeta submodel (rho = 0)
#' @param n Sample size (for information criteria)
#' @param p_full Number of parameters in full model
#' @param p_gz Number of parameters in GZ model
#' @return List containing:
#'   \item{statistic}{LRT statistic Lambda_n^GZ}
#'   \item{p_value}{P-value based on mixture distribution}
#'   \item{decision}{"Reject H0" or "Fail to reject H0"}
#'   \item{model_selected}{"A" if H0 rejected, "C" otherwise}
#'   \item{aic_full, aic_gz}{AIC for comparison}
#'   \item{bic_full, bic_gz}{BIC for comparison}
#' @details 
#'   Under H0, Lambda ~ 0.5 * chi^2_0 + 0.5 * chi^2_1 (Self & Liang, 1987).
#'   P(Lambda > c) = 0.5 * P(chi^2_1 > c) for c > 0.
#'   Critical value at alpha = 0.05 is 2.71 (not 3.84).
#' @references Proposition 3.5, Self & Liang (1987)
#' @seealso fit_generalized_zeta, model_selection
#' @export
boundary_lrt <- function(loglik_full, loglik_gz, n, p_full, p_gz) {
  
  # LRT statistic
  lambda_stat <- 2 * (loglik_full - loglik_gz)
  lambda_stat <- max(lambda_stat, 0)  # Ensure non-negative
  
  # P-value from mixture: 0.5 * P(chi^2_0 > stat) + 0.5 * P(chi^2_1 > stat)
  # P(chi^2_0 > stat) = 0 for stat > 0, = 1 for stat = 0
  if (lambda_stat == 0) {
    p_value <- 1.0
  } else {
    p_value <- 0.5 * (1 - pchisq(lambda_stat, df = 1))
  }
  
  # Decision at 5% level (critical value = 2.71)
  alpha_level <- 0.05
  critical_value <- qchisq(1 - 2*alpha_level, df = 1)  # Approximately 2.71
  
  if (lambda_stat > critical_value) {
    decision <- "Reject H0"
    model_selected <- "A"  # Full model
  } else {
    decision <- "Fail to reject H0"
    model_selected <- "C"  # Generalized Zeta
  }
  
  # Information criteria
  aic_full <- -2 * loglik_full + 2 * p_full
  aic_gz <- -2 * loglik_gz + 2 * p_gz
  
  bic_full <- -2 * loglik_full + log(n) * p_full
  bic_gz <- -2 * loglik_gz + log(n) * p_gz
  
  return(list(
    statistic = lambda_stat,
    p_value = p_value,
    critical_value = critical_value,
    alpha = alpha_level,
    decision = decision,
    model_selected = model_selected,
    aic_full = aic_full,
    aic_gz = aic_gz,
    bic_full = bic_full,
    bic_gz = bic_gz,
    prefer_full_aic = aic_full < aic_gz,
    prefer_full_bic = bic_full < bic_gz
  ))
}

#' Sibuya Boundary Diagnostic (Algorithm 3.5, Stage 1)
#' 
#' Detects if fitted model is in Sibuya danger zone (delta ≈ 1, infinite mean).
#' 
#' @param fit_gz Fitted Generalized Zeta model (from fit_generalized_zeta)
#' @param alpha_level Significance level for diagnostic (default: 0.05)
#' @return List containing:
#'   \item{d_sibuya}{Diagnostic statistic D_Sibuya = (delta_hat - 1) / SE(delta_hat)}
#'   \item{in_danger_zone}{Logical, TRUE if D_Sibuya < 2 (danger zone)}
#'   \item{delta_hat}{Estimated delta (intercept)}
#'   \item{se_delta}{Standard error of delta estimate}
#'   \item{recommendation}{Action recommendation}
#' @details 
#'   D_Sibuya < 2 indicates proximity to Sibuya boundary where standard MLE theory fails.
#'   This triggers need for Algorithm 3.5 (mixed extreme value inference).
#'   
#'   Note: SE(delta) computed from observed information matrix (may be unreliable near boundary).
#' @references Algorithm 3.5, Section 3.7
#' @seealso sibuya_remedy, fit_generalized_zeta
#' @export
sibuya_diagnostic <- function(fit_gz, alpha_level = 0.05) {
  
  # Extract delta estimate (intercept)
  xi_hat <- fit_gz$coefficients$xi
  delta_hat <- 1 + exp(xi_hat[1])
  
  # Compute standard error (approximate from Fisher information)
  # This is a simplified version; full implementation would use analytical Hessian
  se_delta <- NA
  
  # Try to compute numerical Hessian
  tryCatch({
    # Function to optimize over xi only (fixing beta)
    neg_loglik_xi <- function(xi) {
      delta_vec <- 1 + exp(as.vector(fit_gz$design_matrices$Z %*% xi))
      alpha_vec <- fit_gz$theta_list$alpha_vec  # Fixed
      rho_vec <- rep(0, length(fit_gz$theta_list$rho_vec))
      return(-log_likelihood(fit_gz$y, alpha_vec, delta_vec, rho_vec))
    }
    
    # Numerical Hessian at MLE
    H <- numDeriv::hessian(neg_loglik_xi, xi_hat)
    var_xi <- solve(H)
    se_xi <- sqrt(diag(var_xi)[1])
    
    # Delta = 1 + exp(xi), so by delta method:
    # SE(delta) = exp(xi) * SE(xi) = (delta - 1) * SE(xi)
    se_delta <- (delta_hat - 1) * se_xi
    
  }, error = function(e) {
    warning("Could not compute SE for Sibuya diagnostic; using approximation")
    # Conservative approximation
    se_delta <<- (delta_hat - 1) / sqrt(length(fit_gz$theta_list$alpha_vec))
  })
  
  # Diagnostic statistic
  d_sibuya <- (delta_hat - 1) / se_delta
  
  # Danger zone: D_Sibuya < 2 (roughly 2 SE from boundary)
  in_danger_zone <- (d_sibuya < 2) || (delta_hat < 1.2)
  
  # Recommendation
  if (in_danger_zone) {
    recommendation <- "Activate Algorithm 3.5: Mixed extreme value inference protocol"
  } else {
    recommendation <- "Standard inference valid; proceed with bootstrap"
  }
  
  return(list(
    d_sibuya = d_sibuya,
    in_danger_zone = in_danger_zone,
    delta_hat = delta_hat,
    se_delta = se_delta,
    xi_hat = xi_hat[1],
    se_xi = if (!is.na(se_delta)) se_delta / (delta_hat - 1) else NA,
    recommendation = recommendation
  ))
}

#' Model Selection and Confirmation (Algorithm 3.4, Stage 3)
#' 
#' Executes full model selection procedure: fits GZ submodel, performs boundary LRT,
#' and checks Sibuya condition if needed.
#' 
#' @param y Vector of observed counts
#' @param X Design matrix for alpha
#' @param Z Design matrix for delta
#' @param W Design matrix for rho
#' @param fit_full Fitted full model (Model A) from fit_em
#' @param model_type Initial model type from pre-selection
#' @param alpha_level Significance level for LRT (default: 0.05)
#' @param delta_min Hard boundary for delta (default: 1.05)
#' @param verbose Print progress (default: TRUE)
#' @return List containing:
#'   \item{final_model}{Final selected model type ("A" or "C")}
#'   \item{fit_final}{Fitted model object (full or GZ)}
#'   \item{lrt_results}{Results from boundary LRT}
#'   \item{sibuya_check}{Results from Sibuya diagnostic (if applicable)}
#' @references Algorithm 3.4 (Stage 3), Section 3.5-3.6
#' @seealso boundary_lrt, sibuya_diagnostic, fit_generalized_zeta
#' @export
model_selection <- function(y, X, Z, W, fit_full, model_type = "A",
                            alpha_level = 0.05, delta_min = 1.05,
                            verbose = TRUE) {
  
  # Only do selection if initial model was A
  if (model_type != "A") {
    if (verbose) message("Model type ", model_type, " fixed; skipping selection")
    return(list(
      final_model = model_type,
      fit_final = fit_full,
      lrt_results = NULL,
      sibuya_check = NULL
    ))
  }
  
  if (verbose) message("Stage 3: Model selection (boundary LRT for GZ submodel)")
  
  # Fit Generalized Zeta submodel
  if (verbose) message("  Fitting GZ submodel (rho = 0)...")
  fit_gz <- fit_generalized_zeta(y, X, Z, theta_init = NULL,
                                 delta_min = delta_min, verbose = verbose)
  
  # Count parameters
  p_full <- length(fit_full$coefficients$beta) + 
    length(fit_full$coefficients$xi) + 
    length(fit_full$coefficients$zeta)
  p_gz <- length(fit_gz$coefficients$beta) + 
    length(fit_gz$coefficients$xi)
  
  # Boundary LRT
  lrt_results <- boundary_lrt(
    loglik_full = fit_full$loglik,
    loglik_gz = fit_gz$loglik,
    n = length(y),
    p_full = p_full,
    p_gz = p_gz
  )
  
  if (verbose) {
    message("  LRT statistic: ", round(lrt_results$statistic, 4))
    message("  P-value: ", round(lrt_results$p_value, 4))
    message("  Decision: ", lrt_results$decision)
    message("  AIC: full = ", round(lrt_results$aic_full, 2), 
            ", GZ = ", round(lrt_results$aic_gz, 2))
  }
  
  # Select final model
  if (lrt_results$model_selected == "A") {
    final_model <- "A"
    fit_final <- fit_full
    if (verbose) message("  -> Selected: Full PDigamma model (Model A)")
    
    # No Sibuya check needed for Model A (delta should be away from boundary)
    sibuya_check <- NULL
    
  } else {
    final_model <- "C"
    fit_final <- fit_gz
    if (verbose) message("  -> Selected: Generalized Zeta model (Model C)")
    
    # Sibuya diagnostic for GZ model
    if (verbose) message("  Running Sibuya boundary diagnostic...")
    sibuya_check <- sibuya_diagnostic(fit_gz, alpha_level)
    
    if (verbose) {
      message("    D_Sibuya = ", round(sibuya_check$d_sibuya, 3))
      message("    Danger zone: ", sibuya_check$in_danger_zone)
      message("    Recommendation: ", sibuya_check$recommendation)
    }
  }
  
  return(list(
    final_model = final_model,
    fit_final = fit_final,
    lrt_results = lrt_results,
    sibuya_check = sibuya_check
  ))
}

#' Extract Design Matrices from Fit Object
#' 
#' Helper to store design matrices for later use (e.g., in diagnostics).
#' 
#' @param fit Fit object
#' @param X Original X matrix
#' @param Z Original Z matrix  
#' @param W Original W matrix
#' @param y Response vector
#' @return Modified fit object with design_matrices element
#' @export
store_design_info <- function(fit, X, Z, W, y) {
  fit$design_matrices <- list(X = X, Z = Z, W = W)
  fit$y <- y
  return(fit)
}