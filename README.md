
<!-- README.md is generated from README.Rmd. Please edit that file -->

# faircfint

<!-- badges: start -->
<!-- badges: end -->

The goal of faircfint is to facilitate estimation and inference for
intersectional, counterfactual unfairness metrics. The metrics and
techniques are described in [“An intersectional framework for
counterfactual fairness in risk
prediction”](https://doi.org/10.1093/biostatistics/kxad021) and
[“Counterfactual fairness for small
subgroups”](https://arxiv.org/abs/2310.19988).

## Installation

You can install the development version of faircfint from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("swastvedt/faircfint")
```

## Usage

The main function, analysis_estimation, calculates the intersectional,
counterfactual unfairness metrics from data (which must include risk
predictions).

The ‘normal’ option for estimator_type uses the methods from [“An
intersectional framework for counterfactual fairness in risk
prediction”](https://doi.org/10.1093/biostatistics/kxad021). Options
‘small_internal’ and ‘small_borrow’ use the methods from
[“Counterfactual fairness for small
subgroups”](https://arxiv.org/abs/2310.19988).

``` r
library(faircfint)

# Build propensity score model
pi_model_ex <- glm(as.factor(as.character(D)) ~ A1*A2 + A1 + A2 + S_prob + X.1 + X.2 + X.3 + X.4, 
                   data = faircfint::ex_data_estimation, family = "binomial")

# Estimate unfairness metrics, including bootstrap and null distribution
results <- analysis_estimation(data = faircfint::ex_data_estimation, cutoff = 0.5, estimator_type = "normal",
                               pi_model = pi_model_ex, pi_model_type = "glm", pi_xvars = c("X.1", "X.2", "X.3", "X.4"),
                               gen_null = T, bootstrap = 'rescaled')

# Plot results
results_plots <- get_plots(results, sampsize = 5000, alpha = 0.05, m_factor = 0.75,
                           plot_labels = c("Group 1", "Group 2", "Group 3", "Group 4"),
                           plot_values = c(15,16,17,18), plot_colors = c("black", "gray", "red", "dodgerblue"),
                           delta_uval = 0.1)

# Example plot of negative (cFNR) versions of metrics
#results_plots$cfnr
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
