"""
TOPIC: Leave-one-site-out Random Forest validation
WHAT DOES THIS ANALYSIS DO? Predicts each held-out site from all other sites.
WHEN SHOULD I USE IT? When a site/field grouping contains independent prediction units.
INPUT: Same metadata/settings as random_forest.py. OUTPUT: pooled held-out predictions.
MAIN STATISTICS: Pooled R-squared, RMSE, and MAE across all held-out sites.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data"); output_folder = Path("Outputs/Random_Forest")
response_variable = "alpha"; predictor_columns = ["pH_KCl", "P_Mehlich3_mg_100g", "K_Mehlich3_mg_100g", "N_pct", "C_pct", "bio12now.100"]
group_column = "site.id"; random_seed = 42

# =====================================================
# STEP 1: Hold out one group at a time
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").dropna(subset=[response_variable, group_column])
X, y, groups = data[predictor_columns], data[response_variable], data[group_column]
predictions = []
for train_index, test_index in LeaveOneGroupOut().split(X, y, groups):
    numeric = X.select_dtypes(include="number").columns.tolist()
    categorical = [column for column in predictor_columns if column not in numeric]
    preprocess = ColumnTransformer([("numeric", SimpleImputer(strategy="median"), numeric), ("categorical", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("one_hot", OneHotEncoder(handle_unknown="ignore"))]), categorical)])
    forest = Pipeline([("preprocess", preprocess), ("forest", RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=random_seed, n_jobs=-1))])
    forest.fit(X.iloc[train_index], y.iloc[train_index])
    predictions.append(pd.DataFrame({"site": groups.iloc[test_index].values, "observed": y.iloc[test_index].values, "predicted": forest.predict(X.iloc[test_index]), "baseline_prediction": y.iloc[train_index].mean()}))
all_predictions = pd.concat(predictions, ignore_index=True)

# =====================================================
# STEP 2: Report pooled performance and plot it
# =====================================================
metrics = pd.DataFrame([{"model": "random_forest", "n_sites": groups.nunique(), "n_observations": len(data), "pooled_R2": r2_score(all_predictions.observed, all_predictions.predicted), "pooled_RMSE": mean_squared_error(all_predictions.observed, all_predictions.predicted) ** .5, "pooled_MAE": mean_absolute_error(all_predictions.observed, all_predictions.predicted)}, {"model": "training_mean_baseline", "n_sites": groups.nunique(), "n_observations": len(data), "pooled_R2": r2_score(all_predictions.observed, all_predictions.baseline_prediction), "pooled_RMSE": mean_squared_error(all_predictions.observed, all_predictions.baseline_prediction) ** .5, "pooled_MAE": mean_absolute_error(all_predictions.observed, all_predictions.baseline_prediction)}])
all_predictions.to_csv(output_folder / "leave_one_site_out_predictions.csv", index=False); metrics.to_csv(output_folder / "leave_one_site_out_metrics.csv", index=False)
figure, axis = plt.subplots(figsize=(6, 6)); axis.scatter(all_predictions.observed, all_predictions.predicted)
limits = [min(all_predictions.observed.min(), all_predictions.predicted.min()), max(all_predictions.observed.max(), all_predictions.predicted.max())]; axis.plot(limits, limits, "--", color="black")
axis.set(xlabel="Observed response", ylabel="Held-out prediction", title=f"Leave-one-site-out: pooled R2={metrics.iloc[0].pooled_R2:.2f}")
figure.tight_layout(); figure.savefig(output_folder / "leave_one_site_out.png", dpi=300); figure.savefig(output_folder / "leave_one_site_out.pdf")
