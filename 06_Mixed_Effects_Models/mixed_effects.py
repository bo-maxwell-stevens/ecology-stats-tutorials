"""
TOPIC: Multiple linear regression

WHAT DOES THIS ANALYSIS DO? Relates a numeric response to multiple environmental predictors.
WHEN SHOULD I USE IT? When observations are independent and the response is approximately Gaussian.
INPUT: Metadata with a response and predictors. OUTPUT: coefficient table and diagnostics.
MAIN STATISTIC: Regression coefficients estimate associations while holding the other predictors constant.
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

# =====================================================
# STEP 1: Load independent observations
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv")
model_data = data[[response_variable] + fixed_effects].dropna()

# =====================================================
# STEP 2: Fit the multiple linear regression
# =====================================================
formula = response_variable + " ~ " + " + ".join(fixed_effects)
model = smf.ols(formula, data=model_data)
result = model.fit()
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
