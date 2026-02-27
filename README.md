# PDigamma: Heavy-Tailed Count Regression with Tail Index Modeling

[![R-CMD-check](https://github.com/Pluspeak/PDigamma/workflows/R-CMD-check/badge.svg)](https://github.com/Pluspeak/PDigamma/actions)

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Power-Denominator Digamma Family for Heavy-Tailed Count Data**

## Overview

The `PDigamma` R package implements a novel distribution family for regression modeling of heavy-tailed count data, addressing a critical gap in statistical methodology. Unlike classical models that assume exponential tails (Poisson, Negative Binomial) or fixed power-law tails (Zeta), PDigamma allows **covariate-dependent tail heaviness** through explicit tail index regression.

### Key Features

- **Pure Discrete Framework**: No continuous approximation or discretization bias
- **Tail Index Regression**: Model how tail heaviness varies with covariates

  log(δᵢ − 1) = Zᵢᵀξ

  where δ > 1 is the tail index (smaller = heavier tails)
- **Hard Boundary Constraints**: Valid inference even with infinite variance (1 < δ ≤ 2)
- **Adaptive Bootstrap**: m-out-of-n subsampling for heavy-tailed regimes
- **Sibuya Diagnostics**: Automatic detection and remediation for extreme heavy tails

### When to Use PDigamma

| Your Data | Classical Models | PDigamma |
|:---|:---|:---|
| Catastrophe insurance claims | Underestimate extremes | ✓ Captures tail heterogeneity |
| Network attack frequencies | Assume light tails | ✓ Adapts to varying threat levels |
| Epidemic case counts | Fixed dispersion | ✓ Time-varying tail behavior |
| Financial transaction counts | Exponential decay | ✓ Power-law with covariates |

## Installation

### Quick Install (Recommended)

Fast installation without vignettes (vignettes are computationally intensive):

```r
devtools::install_github("yourusername/PDigamma")
```

### With Vignettes

To build vignettes locally (requires significant computation time):

```r
devtools::install_github("yourusername/PDigamma", build_vignettes = TRUE)
```

**Note**: Vignette compilation involves EM algorithm fitting and bootstrap inference which are computationally intensive. The quick install provides all functionality. Subsequent releases are optimized for this issue.

## Quick Start

### One-Line Analysis

```r
library(PDigamma)

# Complete pipeline: data → model → inference
result <- pdigamma_analysis(
  formula_alpha = ~ coverage_amount,      # affects claim frequency
  formula_delta = ~ risk_region,           # affects tail heaviness  
  formula_rho   = ~ 1,                    # common shape
  data          = insurance_data,
  y_var         = "claim_count"
)

summary(result)
```

### Step-by-Step Control

```r
# Stage 0: Diagnose tail behavior
diag <- pre_diagnosis(y, X, Z, W)

# Stage 1: Initialize parameters  
init <- initialize_params(y, X, Z, W, model_type = "A")

# Stage 2: Fit via EM algorithm
fit <- fit_em(y, X, Z, W, init, model_type = "A")

# Stage 3: Select between full model and Generalized Zeta submodel
selected <- model_selection(y, X, Z, W, fit)

# Stage 4: Bootstrap inference
inference <- run_inference(y, X, Z, W, selected$fit_final, 
                           selected$final_model, B = 1000)
```

## Example: Catastrophe Insurance

```r
# Load example data (included)
data(insurance_claims)

# Model: tail index varies by region risk level
fit <- pdigamma_analysis(
  formula_alpha = ~ log(coverage),
  formula_delta = ~ region_risk + earthquake_exposure,
  formula_rho   = ~ 1,
  data          = insurance_claims,
  y_var         = "annual_claims",
  B             = 1000
)

# Visualize tail heterogeneity
plot(fit, type = "tail_index")

# Compare with Negative Binomial (underestimates extremes)
comparison <- compare_methods(insurance_claims, "annual_claims")
```

## Methodology

The PDigamma distribution is defined by the probability mass function:

$$q(x; \alpha, \delta, \rho) = \frac{1}{C(\alpha,\delta,\rho)} \cdot \frac{(\alpha)_x}{x^{\delta/(1+\rho)} \left(\alpha + \frac{\rho\delta}{1+\rho}\right)_x}, \quad x = 1, 2, \ldots$$

where:

- **α > 0**: proximal shape (near-origin behavior)
- **δ > 1**: total tail index (tail decay rate: P(X > x) ~ x^{-(δ-1)})
- **ρ > 0**: coupling ratio (balance between Zeta-like and Digamma-like tails)

### Model Hierarchy

| Model                    | ρ     | Heterogeneity    | Use Case                          |
| :----------------------- | :---- | :--------------- | :-------------------------------- |
| **A** (Full)             | Free  | α(X), δ(Z), ρ(W) | Rich covariates, full flexibility |
| **B** (Tail-focused)     | Fixed | α fixed, δ(Z)    | Limited data, focus on tail       |
| **C** (Generalized Zeta) | 0     | α(X), δ(Z)       | ρ ≈ 0 supported by data           |
| **D** (Proximal-tail)    | Fixed | α(X), δ fixed    | Homogeneous tail, varying mean    |

## Citation

If you use this package, please cite:

```bibtex
@article{...}
```

## Documentation

- [Package vignette](vignettes/PDigamma-introduction.html): Comprehensive tutorial
- [Function reference](https://Pluspeak.github.io/PDigamma/reference): Full API documentation
- [Paper PDF](https://arxiv.org/abs/XXXX.XXXXX): Theoretical foundations

## Comparison with Alternatives

| Method            | Tail Type               | Covariate Effect | Discrete        | Inference      |
| :---------------- | :---------------------- | :--------------- | :-------------- | :------------- |
| Poisson           | Exponential             | Mean only        | ✓               | Standard       |
| Negative Binomial | Exponential             | Mean/dispersion  | ✓               | Standard       |
| Zeta              | Power-law (fixed)       | None             | ✓               | Standard       |
| **PDigamma**      | **Power-law (varying)** | **Tail index**   | **✓**           | **Adaptive**   |
| Continuous POT    | Power-law               | Scale/shape      | ✗ (discretized) | Threshold-dep. |

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.
