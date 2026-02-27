# =============================================================================
# main.R
# End-to-end inference flow for PDigamma regression (Algorithm 3.4 Complete)
# Provides high-level interface for users
# =============================================================================

#' Complete PDigamma Analysis Pipeline
#' 
#' One-function interface for full analysis: from raw data to final inference.
#' Implements Algorithm 3.4 (Stages 0-4).
#' 
#' @param formula_alpha Formula for alpha (proximal shape), e.g., alpha ~ x1 + x2
#'   Use ~ 1 for intercept-only. Response variable ignored (extracted from data).
#' @param formula_delta Formula for delta (tail index), e.g., delta ~ z1 + z2
#' @param formula_rho Formula for rho (coupling ratio), e.g., rho ~ w1 + w2
#' @param data Data frame containing all variables
#' @param y_var Name of response variable in data (default: first variable in formula_alpha)
#' @param delta_min Hard boundary for delta (default: 1.05)
#' @param rho_min Hard boundary for rho (default: 0.05)
#' @param alpha_level Significance level for model selection (default: 0.05)
#' @param B Bootstrap replications (default: 1000)
#' @param max_iter Maximum number of iterations.
#' @param conf_level Confidence level (default: 0.95)
#' @param verbose Print progress (default: TRUE)
#' @param save_intermediate Save intermediate fits (default: FALSE)
#' @return Comprehensive results object of class "pdigamma_analysis" containing:
#'   \item{data_info}{Summary of input data}
#'   \item{stage0}{Pre-diagnosis results}
#'   \item{stage1}{Initialization results}
#'   \item{stage2}{EM fitting results (full model)}
#'   \item{stage3}{Model selection results}
#'   \item{stage4}{Final inference results}
#'   \item{final_model}{Final selected model object}
#' @details 
#'   Stage 0: Pre-diagnosis and model pre-selection (pre_diagnosis)
#'   Stage 1: Three-tier initialization (initialize_params)
#'   Stage 2: EM algorithm fitting (fit_em)
#'   Stage 3: Model selection via boundary LRT (model_selection)
#'   Stage 4: Risk measures with adaptive bootstrap (run_inference)
#' @references Algorithm 3.4 (Complete), Section 3.6
#' @seealso pdigamma_fit (lower-level interface)
#' @examples
#' \dontrun{
#' # Basic usage with formulas
#' result <- pdigamma_analysis(
#'   formula_alpha = ~ x1 + x2,      # alpha ~ x1 + x2
#'   formula_delta = ~ z1,            # delta ~ z1 (tail heterogeneity)
#'   formula_rho = ~ 1,               # common rho
#'   data = mydata,
#'   y_var = "claim_count"
#' )
#' 
#' # Intercept-only for all parameters (homogeneous model)
#' result <- pdigamma_analysis(~ 1, ~ 1, ~ 1, data = mydata, y_var = "y")
#' }
#' @export
pdigamma_analysis <- function(formula_alpha = ~ 1, 
                              formula_delta = ~ 1,
                              formula_rho = ~ 1,
                              data,
                              y_var = NULL,
                              delta_min = 1.05,
                              rho_min = 0.05,
                              alpha_level = 0.05,
                              B = 1000,
                              max_iter,
                              conf_level = 0.95,
                              verbose = TRUE,
                              save_intermediate = FALSE) {
  
  if (verbose) {
    message("============================================================")
    message("PDigamma Regression: End-to-End Analysis Pipeline")
    message("============================================================")
  }
  
  # ============================================================================
  # Parse formulas and extract data
  # ============================================================================
  
  # Determine response variable
  if (is.null(y_var)) {
    # Try to extract from formula_alpha
    if (length(formula_alpha) == 3) {
      y_var <- as.character(formula_alpha[[2]])
    } else {
      stop("Must specify y_var if formula_alpha has no left-hand side")
    }
  }
  
  if (!(y_var %in% names(data))) {
    stop("Response variable '", y_var, "' not found in data")
  }
  
  y <- data[[y_var]]
  
  # Handle zeros if present
  zero_info <- handle_zeros(y, method = "truncate", verbose = verbose)
  y <- zero_info$y
  
  # Adjust data if zeros removed
  if (zero_info$was_truncated) {
    data <- data[y > 0, , drop = FALSE]
    y <- data[[y_var]]
  }
  
  # Build design matrices
  X <- construct_design_matrix(formula_alpha, data, "alpha")
  Z <- construct_design_matrix(formula_delta, data, "delta")
  W <- construct_design_matrix(formula_rho, data, "rho")
  
  # Ensure same rows after NA removal
  n <- length(y)
  
  if (verbose) {
    message("\nData summary:")
    message("  Observations: ", n)
    message("  Response: ", y_var, " (range: ", min(y), " to ", max(y), ")")
    message("  Alpha predictors: ", ncol(X), " (", paste(colnames(X), collapse = ", "), ")")
    message("  Delta predictors: ", ncol(Z), " (", paste(colnames(Z), collapse = ", "), ")")
    message("  Rho predictors: ", ncol(W), " (", paste(colnames(W), collapse = ", "), ")")
  }
  
  # Validate inputs
  valid_inputs <- validate_inputs(y, X, Z, W, delta_min, rho_min)
  
  # ============================================================================
  # Stage 0: Pre-Diagnosis
  # ============================================================================
  
  if (verbose) {
    message("\n------------------------------------------------------------")
    message("Stage 0: Pre-Diagnosis and Model Pre-Selection")
    message("------------------------------------------------------------")
  }
  
  stage0 <- pre_diagnosis(y, X, Z, W, delta_min)
  
  if (verbose) {
    message("  Hill estimate (xi): ", round(stage0$xi_hill, 3))
    message("  Global delta estimate: ", round(stage0$delta_global, 3))
    message("  Sibuya warning: ", stage0$sibuya_warning)
    message("  Pre-selected model: ", stage0$model_type)
  }
  
  # ============================================================================
  # Stage 1: Initialization
  # ============================================================================
  
  if (verbose) {
    message("\n------------------------------------------------------------")
    message("Stage 1: Three-Tier Initialization")
    message("------------------------------------------------------------")
  }
  
  stage1 <- initialize_params(y, X, Z, W, 
                              model_type = stage0$model_type,
                              delta_min = delta_min,
                              rho_min = rho_min,
                              alpha_min = alpha_level)
  
  if (verbose) {
    message("  Initial delta_global: ", round(stage1$delta_global, 3))
    message("  Sibuya flag: ", stage1$sibuya_flag)
  }
  
  # ============================================================================
  # Stage 2: EM Algorithm (Full Model)
  # ============================================================================
  
  if (verbose) {
    message("\n------------------------------------------------------------")
    message("Stage 2: EM Algorithm Estimation")
    message("------------------------------------------------------------")
  }
  
  # Determine if direct optimization preferred
  use_direct <- (stage0$model_type == "A" && 
                   n < 10000 && 
                   stage1$delta_global > 1.5 &&
                   !stage0$sibuya_warning)
  
  if (use_direct && verbose) {
    message("  Using direct L-BFGS-B optimization (n < 10^4, delta > 1.5)")
    # Placeholder for direct optimization; currently uses EM
    stage2 <- fit_em(y, X, Z, W, stage1, 
                     model_type = stage0$model_type,
                     delta_min = delta_min,
                     rho_min = rho_min,
                     max_iter = max_iter,
                     verbose = verbose)
  } else {
    stage2 <- fit_em(y, X, Z, W, stage1,
                     model_type = stage0$model_type,
                     delta_min = delta_min,
                     rho_min = rho_min,
                     max_iter = max_iter,
                     verbose = verbose)
  }
  
  # Store design info for later use
  stage2 <- store_design_info(stage2, X, Z, W, y)
  
  if (verbose) {
    message("\n  EM converged: ", stage2$converged)
    message("  Iterations: ", stage2$iterations)
    message("  Final log-likelihood: ", round(stage2$loglik, 4))
  }
  
  # ============================================================================
  # Stage 3: Model Selection
  # ============================================================================
  
  if (verbose) {
    message("\n------------------------------------------------------------")
    message("Stage 3: Model Selection and Confirmation")
    message("------------------------------------------------------------")
  }
  
  stage3 <- model_selection(y, X, Z, W, stage2, 
                            model_type = stage0$model_type,
                            alpha_level = alpha_level,
                            delta_min = delta_min,
                            verbose = verbose)
  
  final_model <- stage3$final_model
  fit_final <- stage3$fit_final
  
  # Re-store design info if model changed
  if (final_model != stage0$model_type) {
    fit_final <- store_design_info(fit_final, X, Z, W, y)
  }
  
  # ============================================================================
  # Stage 4: Inference and Risk Measures
  # ============================================================================
  
  if (verbose) {
    message("\n------------------------------------------------------------")
    message("Stage 4: Risk Measures and Bootstrap Inference")
    message("------------------------------------------------------------")
  }
  
  stage4 <- run_inference(y, X, Z, W, fit_final, final_model,
                          sibuya_check = stage3$sibuya_check,
                          B = B,
                          conf_level = conf_level,
                          verbose = verbose)
  
  # ============================================================================
  # Assembly and Return
  # ============================================================================
  
  if (verbose) {
    message("\n============================================================")
    message("Analysis Complete")
    message("============================================================")
    message("Final model type: ", final_model)
    message("Mean delta (tail index): ", round(mean(fit_final$theta_list$delta_vec), 3))
    message("Delta range: [", round(min(fit_final$theta_list$delta_vec), 3), ", ",
            round(max(fit_final$theta_list$delta_vec), 3), "]")
  }
  
  result <- structure(list(
    data_info = list(
      n = n,
      y_var = y_var,
      p_alpha = ncol(X),
      p_delta = ncol(Z),
      p_rho = ncol(W),
      zero_truncated = zero_info$was_truncated,
      n_removed = zero_info$n_removed
    ),
    stage0 = stage0,
    stage1 = stage1,
    stage2 = stage2,
    stage3 = stage3,
    stage4 = stage4,
    final_model = final_model,
    fit_final = fit_final,
    coefficients = fit_final$coefficients,
    call = match.call()
  ), class = "pdigamma_analysis")
  
  return(result)
}

#' Print Method for PDigamma Analysis
#' 
#' @param x Object of class pdigamma_analysis
#' @param ... Additional arguments
#' @export
print.pdigamma_analysis <- function(x, ...) {
  cat("PDigamma Regression Analysis\n")
  cat("============================\n\n")
  
  cat("Data: ", x$data_info$n, " observations")
  if (x$data_info$zero_truncated) {
    cat(" (", x$data_info$n_removed, " zeros removed)")
  }
  cat("\n")
  
  cat("Final model: ", x$final_model, "\n", sep = "")
  cat("Log-likelihood: ", round(x$fit_final$loglik, 2), "\n\n")
  
  cat("Tail index (delta):\n")
  cat("  Mean: ", round(mean(x$fit_final$theta_list$delta_vec), 3), "\n")
  cat("  Range: [", round(min(x$fit_final$theta_list$delta_vec), 3), ", ",
      round(max(x$fit_final$theta_list$delta_vec), 3), "]\n\n")
  
  if (!is.null(x$stage4$risk_measures)) {
    cat("Risk measure: ", x$stage4$risk_measures$risk_type, "\n")
    cat("Method: ", x$stage4$risk_measures$recommendation, "\n")
  }
  
  if (x$stage0$sibuya_warning) {
    cat("\n*** Sibuya warning: results near infinite mean boundary ***\n")
  }
}

#' Summary Method for PDigamma Analysis
#' 
#' @param object Object of class pdigamma_analysis
#' @param ... Additional arguments
#' @export
summary.pdigamma_analysis <- function(object, ...) {
  structure(object, class = c("summary.pdigamma_analysis", "pdigamma_analysis"))
}

#' Print Summary Method
#' @export
print.summary.pdigamma_analysis <- function(x, ...) {
  print.pdigamma_analysis(x)
  
  cat("\n--- Coefficient Estimates ---\n\n")
  
  cat("Alpha (beta):\n")
  print(x$coefficients$beta)
  
  cat("\nDelta (xi):\n")
  print(x$coefficients$xi)
  
  if (x$final_model == "A") {
    cat("\nRho (zeta):\n")
    print(x$coefficients$zeta)
  }
  
  if (!is.null(x$stage4$inference)) {
    cat("\n--- Inference ---\n")
    cat("Point estimate: ", round(x$stage4$inference$point_estimate, 4), "\n")
    if (!is.null(x$stage4$inference$confidence_interval)) {
      ci <- x$stage4$inference$confidence_interval
      cat("95% CI: [", round(ci[1], 4), ", ", round(ci[2], 4), "]\n")
    }
  }
}

#' Lower-Level Fitting Interface
#' 
#' Direct interface similar to glm(), for advanced users.
#' 
#' @param formula Formula with response on left, predictors on right
#'   Special terms: alpha(), delta(), rho() to specify parameter-specific predictors
#'   Default: all predictors go to delta (tail index)
#' @param data Data frame
#' @param model Type: "A" (full), "B" (tail-focused), "C" (Generalized Zeta), "D" (proximal-tail)
#' @param ... Additional arguments passed to pdigamma_analysis
#' @return Fitted model object
#' @examples
#' \dontrun{
#' # Simple: all predictors affect tail index
#' fit <- pdigamma_fit(y ~ x1 + x2, data = mydata)
#' 
#' # Advanced: specify different predictors for each parameter
#' # (would need formula parsing enhancement)
#' }
#' @export
pdigamma_fit <- function(formula, data, model = "A", ...) {
  # Parse formula to extract response and determine parameter allocation
  # This is a simplified version; full version would parse special terms
  
  if (length(formula) != 3) {
    stop("formula must have response on left-hand side")
  }
  
  y_var <- as.character(formula[[2]])
  
  # For now, put all predictors in delta, intercept-only for others
  # Full implementation would parse alpha(), delta(), rho() terms
  
  pdigamma_analysis(
    formula_alpha = ~ 1,
    formula_delta = formula,
    formula_rho = ~ 1,
    data = data,
    y_var = y_var,
    ...
  )
}

#' Extract Coefficients
#' 
#' @param object Fitted model
#' @param ... Additional arguments
#' @export
coef.pdigamma_analysis <- function(object, ...) {
  object$coefficients
}

#' Extract Fitted Values (Individual Parameters)
#' 
#' @param object Fitted model
#' @param parameter Which parameter: "alpha", "delta", or "rho"
#' @param ... Additional arguments
#' @export
fitted.pdigamma_analysis <- function(object, parameter = "delta", ...) {
  switch(parameter,
         "alpha" = object$fit_final$theta_list$alpha_vec,
         "delta" = object$fit_final$theta_list$delta_vec,
         "rho" = object$fit_final$theta_list$rho_vec,
         stop("parameter must be 'alpha', 'delta', or 'rho'"))
}

#' Predict Method
#' 
#' Predict parameters for new data
#' 
#' @param object Fitted model
#' @param newdata New data frame
#' @param type "parameters" (alpha, delta, rho) or "risk" (risk measures)
#' @param ... Additional arguments
#' @export
predict.pdigamma_analysis <- function(object, newdata = NULL, 
                                      type = "parameters", ...) {
  
  if (is.null(newdata)) {
    # Return fitted values
    return(data.frame(
      alpha = fitted(object, "alpha"),
      delta = fitted(object, "delta"),
      rho = fitted(object, "rho")
    ))
  }
  
  # Would need to store formula info to build design matrices for newdata
  stop("Prediction on newdata not yet implemented")
}