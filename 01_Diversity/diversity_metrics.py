"""
TOPIC: Alpha diversity

WHAT DOES THIS ANALYSIS DO?
Calculates within-sample diversity from microbial count tables.

WHEN SHOULD I USE IT?
When each row is an independently collected community sample.

INPUT: AMF and ITS abundance tables. OUTPUT: diversity tables and histograms.
MAIN STATISTIC: Shannon diversity combines richness and evenness.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data")
output_folder = Path("Outputs/Diversity")
community_files = {"AMF": "AMF_OTU_table_final.tsv", "ITS": "ITS_OTU_table_final.tsv"}
random_seed = 42

# =====================================================
# STEP 1: Load the data
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
sns.set_theme(style="whitegrid", context="talk")

all_summaries = []
for community_name, filename in community_files.items():
    counts = pd.read_csv(data_folder / filename, sep="\t", index_col=0).fillna(0)
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    # =====================================================
    # STEP 2: Calculate diversity metrics
    # =====================================================
    total_reads = counts.sum(axis=1)
    richness = (counts > 0).sum(axis=1)
    proportions = counts.div(total_reads.replace(0, np.nan), axis=0).fillna(0)
    shannon = -(proportions.where(proportions > 0, 1) * np.log(proportions.where(proportions > 0, 1))).sum(axis=1)
    simpson = 1 - (proportions ** 2).sum(axis=1)
    inverse_simpson = 1 / (proportions ** 2).sum(axis=1).replace(0, np.nan)
    evenness = shannon / np.log(richness.replace(0, np.nan))

    summary = pd.DataFrame({
        "sample_id": counts.index, "community": community_name,
        "total_reads": total_reads, "observed_richness": richness,
        "shannon": shannon, "pielou_evenness": evenness,
        "simpson": simpson, "inverse_simpson": inverse_simpson,
        "effective_species_number": np.exp(shannon),
    })
    summary.to_csv(output_folder / f"{community_name}_diversity_summary.csv", index=False)
    all_summaries.append(summary)

    # =====================================================
    # STEP 3: Plot and save the results
    # =====================================================
    figure, axis = plt.subplots(figsize=(8, 5))
    sns.histplot(summary["observed_richness"], bins=20, ax=axis, color="#237a57")
    axis.set(xlabel="Observed richness (number of taxa)", ylabel="Number of samples", title=f"{community_name}: richness distribution")
    figure.tight_layout()
    figure.savefig(output_folder / f"{community_name}_richness_histogram.png", dpi=300)
    figure.savefig(output_folder / f"{community_name}_richness_histogram.pdf")
    plt.close(figure)

# A single table makes later tutorials easier to use.
pd.concat(all_summaries, ignore_index=True).to_csv(output_folder / "all_diversity_metrics.csv", index=False)
print(f"Saved diversity tables and figures to {output_folder}")
