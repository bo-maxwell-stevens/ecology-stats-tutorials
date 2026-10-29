"""
TOPIC: Random Forest regression

WHAT DOES THIS ANALYSIS DO? Predicts a continuous ecological response from predictors.
WHEN SHOULD I USE IT? After defining independent train/test groups and a prediction question.
INPUT: Metadata. OUTPUT: test metrics, predictions, and predictor importance.
MAIN STATISTICS: R-squared, RMSE, and MAE quantify held-out prediction error.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data"); output_folder = Path("Outputs/Random_Forest")
response_variable = "alpha"
predictor_columns = ["pH_KCl", "P_Mehlich3_mg_100g", "K_Mehlich3_mg_100g", "N_pct", "C_pct", "bio12now.100"]
group_column = "site.id"  # Use field ID when samples are nested within fields.
random_seed = 42

# =====================================================
# STEP 1: Split independent groups before preprocessing
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
data = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").dropna(subset=[response_variable, group_column])
X, y, groups = data[predictor_columns], data[response_variable], data[group_column]
splitter = GroupShuffleSplit(n_splits=1, test_size=.20, random_state=random_seed)
train_index, test_index = next(splitter.split(X, y, groups))
X_train, X_test, y_train, y_test = X.iloc[train_index], X.iloc[test_index], y.iloc[train_index], y.iloc[test_index]

# =====================================================
# STEP 2: Build a training-only preprocessing/model pipeline
# =====================================================
numeric_columns = X.select_dtypes(include="number").columns.tolist()
categorical_columns = [column for column in predictor_columns if column not in numeric_columns]
preprocessing = ColumnTransformer([
    ("numeric", SimpleImputer(strategy="median"), numeric_columns),
    ("categorical", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("one_hot", OneHotEncoder(handle_unknown="ignore"))]), categorical_columns),
])
model = Pipeline([("preprocessing", preprocessing), ("forest", RandomForestRegressor(n_estimators=500, min_samples_leaf=2, random_state=random_seed, n_jobs=-1))])
model.fit(X_train, y_train)

# =====================================================
# STEP 3: Evaluate held-out samples and save results
# =====================================================
predicted = model.predict(X_test)
baseline = np.repeat(y_train.mean(), len(y_test))  # Baseline uses training response only.
metrics = pd.DataFrame([{"model": "random_forest", "R2": r2_score(y_test, predicted), "RMSE": mean_squared_error(y_test, predicted) ** .5, "MAE": mean_absolute_error(y_test, predicted)}, {"model": "training_mean_baseline", "R2": r2_score(y_test, baseline), "RMSE": mean_squared_error(y_test, baseline) ** .5, "MAE": mean_absolute_error(y_test, baseline)}])
metrics.to_csv(output_folder / "random_forest_metrics.csv", index=False)
pd.DataFrame({"observed": y_test, "predicted": predicted, "site": groups.iloc[test_index]}).to_csv(output_folder / "random_forest_test_predictions.csv", index=False)

# =====================================================
# STEP 4: Plot predictions and importance
# =====================================================
figure, axis = plt.subplots(figsize=(6, 6)); axis.scatter(y_test, predicted)
limits = [min(y_test.min(), predicted.min()), max(y_test.max(), predicted.max())]; axis.plot(limits, limits, "--", color="black")
axis.set(xlabel="Observed response", ylabel="Predicted response", title=f"Held-out prediction: R2={metrics.iloc[0].R2:.2f}")
figure.tight_layout(); figure.savefig(output_folder / "random_forest_observed_predicted.png", dpi=300); figure.savefig(output_folder / "random_forest_observed_predicted.pdf")
importance = pd.Series(model.named_steps["forest"].feature_importances_, index=model.named_steps["preprocessing"].get_feature_names_out()).sort_values()
importance.to_csv(output_folder / "random_forest_importance.csv", header=["importance"])
figure, axis = plt.subplots(figsize=(7, 4)); importance.plot.barh(ax=axis); axis.set(xlabel="Impurity-based importance", title="Predictor importance")
figure.tight_layout(); figure.savefig(output_folder / "random_forest_importance.png", dpi=300); figure.savefig(output_folder / "random_forest_importance.pdf")
