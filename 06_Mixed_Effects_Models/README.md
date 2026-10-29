# Mixed-effects models

Fixed effects estimate population-average predictor associations. A random intercept allows measurements from the same site to have a shared baseline. Gaussian mixed models assume approximately normal, constant-variance residuals, conditional on the model.

Run `python 06_Mixed_Effects_Models/mixed_effects.py`. It correctly stops on the supplied data because each environmental `site.id` has one record; a site random effect cannot be estimated. Use a replicated field/site grouping in your own study. Count or binary disease outcomes generally need a different model family. Documentation: [statsmodels MixedLM](https://www.statsmodels.org/stable/mixed_linear.html).
