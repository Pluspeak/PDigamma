# PDigamma 0.1.0

## Initial Release

### New Features
- Core PDigamma distribution with `dPDigamma()`, `rPDigamma()`, `pPDigamma()`, `qPDigamma()`
- Tail index regression via `pdigamma_analysis()` with heterogeneous parameters
- EM algorithm with hard boundary projection (`fit_em()`)
- Three-tier initialization strategy (`initialize_params()`)
- Model selection via boundary likelihood ratio test (`model_selection()`)
- Adaptive m-out-of-n bootstrap for heavy tails (`adaptive_bootstrap()`)
- Sibuya boundary diagnostics and remedy protocol (`sibuya_diagnostic()`, `sibuya_remedy()`)

### Documentation
- Vignette: "Introduction to PDigamma" with 5 examples
- Example scripts: simulation study, real data analysis, method comparison
- Complete function documentation via roxygen2

### Methods
- Comparison with Poisson, Negative Binomial, and continuous POT
- Tail index visualization and interpretation guides