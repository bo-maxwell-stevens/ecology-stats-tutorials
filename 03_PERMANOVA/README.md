# PERMANOVA

Tests a multivariable Bray-Curtis model with `vegan::adonis2`; pseudo-F, R-squared, p-value, and sample size are written by R. The null hypothesis is no association between composition and each predictor after the other terms. This tutorial uses **marginal** tests; sequential tests depend on formula order. Correlated environmental variables and spatial structure complicate individual-term interpretation.

Run `python 03_PERMANOVA/permanova.py`. The included data have no replicated categorical management variable, so dispersion testing is only activated after setting `categorical_predictor` to a valid column. Differences in dispersion can be mistaken for differences in centroids. Documentation: [vegan](https://vegandevs.github.io/vegan/reference/adonis.html).
