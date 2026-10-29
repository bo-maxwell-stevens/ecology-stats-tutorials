# Field regression and mixed effects

The script fits ordinary regression using one row per independent field, then a mixed-effects model using the two soil samples nested in each field. A field random intercept represents shared field conditions and prevents related samples being treated as independent fields.

Run `python 06_Mixed_Effects_Models/mixed_effects.py`. The OLS model asks how field black-dot severity varies with pH and crop rotation. The MixedLM model asks how a synthetic sample-level fungal response varies with pH and management while accounting for field. Inspect residual plots; effects are associations, not causal estimates. Documentation: [statsmodels](https://www.statsmodels.org/stable/mixed_linear.html).
