"""
TOPIC: Linear mixed-effects model

WHAT DOES THIS ANALYSIS DO? Relates a numeric response to predictors while allowing site-specific intercepts.
WHEN SHOULD I USE IT? When several independent observations occur within each site/field.
INPUT: Metadata with response, predictors, and replicated group. OUTPUT: model table and diagnostics.
MAIN STATISTIC: Fixed-effect coefficients estimate average associations; random effects model group variation.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import statsmodels.api as sm

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data"); output_folder = Path("Outputs/Mixed_Models")
response_variable = "alpha"  # Replace with fungal diversity or a suitable continuous disease response.
fixed_effects = ["pH_KCl", "P_Mehlich3_mg_100g"]
group_column = "site.id"  # Replace with field/site only when it genuinely has repeated observations.

# =====================================================
# STEP 1: Check whether the data support a random intercept
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv")
model_data = data[[response_variable, group_column] + fixed_effects].dropna()
group_sizes = model_data[group_column].value_counts()
if len(group_sizes) < 2 or group_sizes.max() < 2:
    print("Mixed model not run: no site has repeated observations in this data. A random intercept is not identifiable.")
    print("For potato data, use a real replicated field/site column; do not treat sequence replicates as independent fields.")
    raise SystemExit

# =====================================================
# STEP 2: Fit fixed and random effects
# =====================================================
formula = response_variable + " ~ " + " + ".join(fixed_effects)
model = smf.mixedlm(formula, model_data, groups=model_data[group_column])
result = model.fit(reml=True)
confidence_intervals = result.conf_int()
pd.DataFrame({"coefficient": result.params, "standard_error": result.bse, "p_value": result.pvalues, "CI_low": confidence_intervals[0], "CI_high": confidence_intervals[1]}).to_csv(output_folder / "mixed_model_coefficients.csv")

# =====================================================
# STEP 3: Check model assumptions
# =====================================================
figure, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].scatter(result.fittedvalues, result.resid, alpha=.7); axes[0].axhline(0, color="firebrick")
axes[0].set(xlabel="Fitted values", ylabel="Residuals", title="Residuals versus fitted")
sm.qqplot(result.resid, line="45", ax=axes[1]); axes[1].set_title("Normal Q-Q plot")
figure.tight_layout(); figure.savefig(output_folder / "mixed_model_diagnostics.png", dpi=300); figure.savefig(output_folder / "mixed_model_diagnostics.pdf")
print(result.summary())
