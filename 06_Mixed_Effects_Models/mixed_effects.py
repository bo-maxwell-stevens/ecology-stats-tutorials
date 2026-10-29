"""TOPIC: Field regression and sample-level mixed effects.

OLS uses one independent row per field. MixedLM uses the two soil samples within
each field, with a field random intercept so replicates are not treated as fields.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf

# USER SETTINGS
data_folder = Path("Data/Synthetic_Potato_Example")
output_folder = Path("Outputs/Synthetic/Mixed_Models")
metadata_file = "Final_data_with_diversity_prefixed.csv"

# STEP 1: Load field and sample data
output_folder.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(data_folder / metadata_file)
fields = samples.groupby("site.id", as_index=False).first()

# STEP 2: Ordinary regression: one observation per independent field
ols = smf.ols("Severity_Har_Cc ~ pH_KCl + C(crop_rotation)", data=fields).fit()
ols_ci = ols.conf_int()
pd.DataFrame({"coefficient": ols.params, "standard_error": ols.bse, "p_value": ols.pvalues, "CI_low": ols_ci[0], "CI_high": ols_ci[1]}).to_csv(output_folder / "field_ols_coefficients.csv")

# STEP 3: Mixed model: a documented sample-level response, nested within fields
mixed = smf.mixedlm("Sample_fungal_response ~ pH_KCl + C(management)", samples, groups=samples["site.id"])
mixed_result = mixed.fit()
pd.DataFrame({"coefficient": mixed_result.params, "standard_error": mixed_result.bse, "p_value": mixed_result.pvalues}).to_csv(output_folder / "sample_mixedlm_coefficients.csv")

# STEP 4: OLS diagnostic plot
figure, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].scatter(ols.fittedvalues, ols.resid); axes[0].axhline(0, color="firebrick")
axes[0].set(xlabel="Fitted severity", ylabel="Residual", title="Field OLS residuals")
sm.qqplot(ols.resid, line="45", ax=axes[1]); axes[1].set_title("OLS Q-Q plot")
figure.tight_layout(); figure.savefig(output_folder / "field_ols_diagnostics.png", dpi=300); figure.savefig(output_folder / "field_ols_diagnostics.pdf")
print(f"OLS: {len(fields)} independent fields. MixedLM: {len(samples)} non-empty samples nested in {samples['site.id'].nunique()} fields.")
