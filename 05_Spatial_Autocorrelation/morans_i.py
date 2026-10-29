"""
TOPIC: Moran's I spatial autocorrelation

WHAT DOES THIS ANALYSIS DO? Tests whether similar values occur near each other geographically.
WHEN SHOULD I USE IT? When independent sites have coordinates and a numeric response.
INPUT: Metadata, ITS and AMF count tables. OUTPUT: Moran's I tables and scatterplots.
MAIN STATISTIC: Moran's I is positive for local similarity and negative for local dissimilarity.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from libpysal.weights import KNN
from esda import Moran

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data"); output_folder = Path("Outputs/Spatial")
variables_to_test = ["ITS_richness", "AMF_richness", "pH_KCl"]
number_of_neighbours = 5
permutations = 999
random_seed = 42

# =====================================================
# STEP 1: Load metadata and calculate richness
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
np.random.seed(random_seed)
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").set_index("site.id")
for community, filename in {"ITS": "ITS_OTU_table_final.tsv", "AMF": "AMF_OTU_table_final.tsv"}.items():
    counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    counts["canonical"] = counts.index.str.replace(r"\.N\d+$", "", regex=True)
    counts = counts.groupby("canonical").sum()
    metadata[f"{community}_richness"] = (counts > 0).sum(axis=1)

# =====================================================
# STEP 2: Build geographic k-nearest-neighbour weights
# =====================================================
# Convert longitude/latitude to approximate kilometre coordinates before KNN.
coordinates = metadata[["lat", "lon"]].copy()
coordinates["x_km"] = coordinates["lon"] * 111 * np.cos(np.deg2rad(coordinates["lat"]))
coordinates["y_km"] = coordinates["lat"] * 111
results = []
for variable in variables_to_test:
    table = metadata[[variable]].join(coordinates[["x_km", "y_km"]]).dropna()
    # Repeated positions have zero distance. Aggregate them rather than inventing distance with jitter.
    table = table.groupby(["x_km", "y_km"], as_index=False)[variable].mean()
    if len(table) <= number_of_neighbours:
        print(f"Skipped {variable}: too few unique locations for k={number_of_neighbours}."); continue
    weights = KNN.from_array(table[["x_km", "y_km"]].to_numpy(), k=number_of_neighbours)
    weights.transform = "r"

    # =====================================================
    # STEP 3: Calculate Moran's I and plot it
    # =====================================================
    moran = Moran(table[variable].to_numpy(), weights, permutations=permutations)
    results.append({"variable": variable, "n_unique_locations": len(table), "Moran_I": moran.I, "expected_I": moran.EI, "permutation_p_value": moran.p_sim})
    standardized_value = (table[variable] - table[variable].mean()) / table[variable].std()
    spatial_lag = weights.sparse @ standardized_value
    figure, axis = plt.subplots(figsize=(7, 5))
    axis.scatter(standardized_value, spatial_lag, alpha=.75)
    slope = np.polyfit(standardized_value, spatial_lag, 1)[0]
    axis.axline((0, 0), slope=slope, color="firebrick")
    axis.set(xlabel=f"Standardized {variable}", ylabel="Spatial lag", title=f"Moran scatterplot: I={moran.I:.3f}, p={moran.p_sim:.3f}")
    figure.tight_layout(); figure.savefig(output_folder / f"{variable}_moran.png", dpi=300); figure.savefig(output_folder / f"{variable}_moran.pdf"); plt.close(figure)
pd.DataFrame(results).to_csv(output_folder / "morans_i_results.csv", index=False)
