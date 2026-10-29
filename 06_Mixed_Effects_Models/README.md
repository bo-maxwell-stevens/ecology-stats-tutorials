# Linear regression

This version fits an ordinary multiple linear regression because the supplied metadata have one independent environmental record per `site.id`. It estimates the association between a continuous response and the listed environmental predictors while holding other model predictors constant.

Run `python 06_Mixed_Effects_Models/mixed_effects.py`. Inspect residual-versus-fitted and Q-Q plots for approximate constant variance and residual normality. Predictor effects are associations, not causal estimates. If a future dataset has genuinely repeated observations within fields or sites, a random-intercept model may then be appropriate. Count or binary disease outcomes generally need a different model family. Documentation: [statsmodels OLS](https://www.statsmodels.org/stable/regression.html).
