"""
TOPIC: NMDS and Principal Coordinates Analysis (PCoA)

WHAT DOES THIS ANALYSIS DO? Visualizes differences in community composition.
WHEN SHOULD I USE IT? When I have a community abundance table and sample metadata.
INPUT: Community abundance table and sample metadata.
OUTPUT: NMDS/PCoA plots and coordinates.
MAIN STATISTIC: Bray-Curtis dissimilarity describes compositional differences.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import MDS
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data")
output_folder = Path("Outputs/Ordination")
community_files = {"AMF": "AMF_OTU_table_final.tsv", "ITS": "ITS_OTU_table_final.tsv"}
colour_variable = "pH_KCl"  # Other real numeric choice: P_Mehlich3_mg_100g
random_seed = 42
aggregate_to_site = True  # Required for supplied site-level metadata; do not assume .N samples are technical replicates.

# =====================================================
# STEP 1: Load metadata and align samples
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="talk")
# `site.id` (for example D018) matches a sequencing ID after its .N# suffix is removed.
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").set_index("site.id")
for community_name, filename in community_files.items():
    counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    if aggregate_to_site:
        counts["site.id"] = counts.index.str.replace(r"\.N\d+$", "", regex=True)
        counts = counts.groupby("site.id").sum()
    else:
        raise ValueError("The supplied metadata are site-level. Use sample-level metadata before setting aggregate_to_site=False.")
    shared_samples = counts.index.intersection(metadata.index)
    counts, matched_metadata = counts.loc[shared_samples], metadata.loc[shared_samples]
    matched_metadata = matched_metadata.dropna(subset=[colour_variable])
    counts = counts.loc[matched_metadata.index]
    zero_reads = counts.sum(axis=1) == 0
    print(f"{community_name}: {len(counts)} sites matched; {zero_reads.sum()} zero-read sites excluded.")
    counts, matched_metadata = counts.loc[~zero_reads], matched_metadata.loc[~zero_reads]
    if len(counts) < 3:
        raise ValueError("At least three non-empty matched communities are required for ordination.")

    # =====================================================
    # STEP 2: Prepare the community table
    # =====================================================
    # Relative abundance removes library-size scale; square roots reduce dominance by abundant taxa.
    relative_abundance = counts.div(counts.sum(axis=1), axis=0).fillna(0)
    transformed = np.sqrt(relative_abundance)

    # =====================================================
    # STEP 3: Calculate Bray-Curtis dissimilarities
    # =====================================================
    distance_matrix = beta_diversity("braycurtis", transformed.values, ids=transformed.index)

    # =====================================================
    # STEP 4: Run NMDS and PCoA
    # =====================================================
    nmds = MDS(n_components=2, metric=False, dissimilarity="precomputed", n_init=10, max_iter=500, random_state=random_seed, normalized_stress="auto")
    nmds_coordinates = pd.DataFrame(nmds.fit_transform(distance_matrix.data), index=counts.index, columns=["NMDS1", "NMDS2"])
    pcoa_result = pcoa(distance_matrix)
    pcoa_coordinates = pcoa_result.samples.iloc[:, :2].rename(columns={pcoa_result.samples.columns[0]: "PCoA1", pcoa_result.samples.columns[1]: "PCoA2"})
    explained = pcoa_result.proportion_explained.iloc[:2] * 100
    # Bray-Curtis can yield negative PCoA eigenvalues. Calculate them from the
    # doubly centred squared-distance matrix because some software truncates them to zero.
    sample_count = len(distance_matrix.ids)
    centring_matrix = np.eye(sample_count) - np.ones((sample_count, sample_count)) / sample_count
    raw_eigenvalues = np.linalg.eigvalsh(-0.5 * centring_matrix @ (distance_matrix.data ** 2) @ centring_matrix)
    negative_eigenvalues = (raw_eigenvalues < -1e-10).sum()
    nmds_coordinates.join(matched_metadata[[colour_variable]]).to_csv(output_folder / f"{community_name}_nmds_coordinates.csv")
    pcoa_coordinates.join(matched_metadata[[colour_variable]]).to_csv(output_folder / f"{community_name}_pcoa_coordinates.csv")

    # =====================================================
    # STEP 5: Plot and save results
    # =====================================================
    for method, coordinates, x_column, y_column, subtitle in [
        ("NMDS", nmds_coordinates, "NMDS1", "NMDS2", f"Stress = {nmds.stress_:.3f}"),
        ("PCoA", pcoa_coordinates, "PCoA1", "PCoA2", f"Axes explain {explained.iloc[0]:.1f}% and {explained.iloc[1]:.1f}%; negative eigenvalues: {negative_eigenvalues}"),
    ]:
        plot_data = coordinates.join(matched_metadata[[colour_variable]])
        figure, axis = plt.subplots(figsize=(7, 6))
        if pd.api.types.is_numeric_dtype(plot_data[colour_variable]):
            scatter = axis.scatter(plot_data[x_column], plot_data[y_column], c=plot_data[colour_variable], cmap="viridis", s=70, edgecolor="black", linewidth=.3)
            figure.colorbar(scatter, ax=axis, label=colour_variable)
        else:
            for group, group_data in plot_data.groupby(colour_variable):
                axis.scatter(group_data[x_column], group_data[y_column], label=group, s=70, edgecolor="black", linewidth=.3)
            axis.legend(title=colour_variable, bbox_to_anchor=(1.02, 1), loc="upper left")
        axis.axhline(0, color="grey", linewidth=.7); axis.axvline(0, color="grey", linewidth=.7)
        axis.set(title=f"{community_name} {method}: coloured by {colour_variable}\n{subtitle}")
        figure.tight_layout()
        figure.savefig(output_folder / f"{community_name}_{method.lower()}.png", dpi=300)
        figure.savefig(output_folder / f"{community_name}_{method.lower()}.pdf")
        plt.close(figure)
print("Saved ordinations. NMDS stress is a lack-of-fit measure: lower is generally better.")
