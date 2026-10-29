"""
TOPIC: ITS fungal functional groups

WHAT DOES THIS ANALYSIS DO? Summarizes ITS reads by annotated ecological lifestyle.
WHEN SHOULD I USE IT? When ITS feature identifiers match a functional annotation table.
INPUT: ITS count table, ITS metadata, and environmental metadata.
OUTPUT: Annotation coverage, group relative abundance, barplot, and correlations.
MAIN STATISTIC: Spearman rho describes a rank-based association with an environmental predictor.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# =====================================================
# USER SETTINGS
# =====================================================
data_folder = Path("Data"); output_folder = Path("Outputs/Functional_Groups")
environmental_predictor = "pH_KCl"
functional_groups = {
    "saprotroph": ["soil_saprotroph", "litter_saprotroph", "wood_saprotroph", "dung_saprotroph", "pollen_saprotroph", "nectar/tap_saprotroph"],
    "symbiotroph": ["arbuscular_mycorrhizal", "ectomycorrhizal", "root_endophyte", "foliar_endophyte", "lichenized", "moss_symbiont"],
    "pathotroph_or_parasite": ["plant_pathogen", "unspecified_pathotroph", "animal_parasite", "lichen_parasite", "mycoparasite"],
}

# =====================================================
# STEP 1: Load data and match SH identifiers
# =====================================================
output_folder.mkdir(parents=True, exist_ok=True)
all_counts = pd.read_csv(data_folder / "ITS_OTU_table_final.tsv", sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
annotations = pd.read_csv(data_folder / "ITS_feature_metadata.tsv", sep="\t").set_index("SH")
shared_taxa = all_counts.columns.intersection(annotations.index)
counts, annotations = all_counts[shared_taxa], annotations.loc[shared_taxa]
zero_reads = all_counts.sum(axis=1) == 0
all_counts, counts = all_counts.loc[~zero_reads], counts.loc[~zero_reads]
print(f"{len(all_counts)} ITS libraries used; {zero_reads.sum()} zero-read libraries excluded.")

# =====================================================
# STEP 2: Transparently map primary lifestyles to broad groups
# =====================================================
# We use only primary_lifestyle. Secondary annotations are retained in the coverage table
# but not double-counted, because summing ambiguous assignments can inflate total abundance.
def broad_group(lifestyle):
    for group, labels in functional_groups.items():
        if lifestyle in labels:
            return group
    return "unclassified"
annotations["broad_group"] = annotations["primary_lifestyle"].map(broad_group)
coverage = pd.DataFrame({
    "total_SH_in_count_table": [len(counts.columns)],
    "primary_lifestyle_annotated_SH": [annotations["primary_lifestyle"].notna().sum()],
    "mapped_to_broad_group_SH": [(annotations["broad_group"] != "unclassified").sum()],
    "reads_in_mapped_SH": [counts.loc[:, annotations["broad_group"] != "unclassified"].sum().sum()],
    "all_ITS_reads": [all_counts.sum().sum()],
})
coverage["fraction_reads_mapped"] = coverage.reads_in_mapped_SH / coverage.all_ITS_reads
coverage.to_csv(output_folder / "functional_annotation_coverage.csv", index=False)

# =====================================================
# STEP 3: Calculate relative functional-group abundance
# =====================================================
group_counts = pd.DataFrame(index=counts.index)
for group in functional_groups:
    group_counts[group] = counts.loc[:, annotations["broad_group"] == group].sum(axis=1)
# The denominator is all ITS reads, including unmatched and unclassified taxa.
group_counts["unclassified_or_unmatched"] = all_counts.sum(axis=1) - group_counts.sum(axis=1)
group_counts["site.id"] = group_counts.index.str.replace(r"\.N\d+$", "", regex=True)
site_group_counts = group_counts.groupby("site.id").sum()
site_total_reads = all_counts.copy(); site_total_reads["site.id"] = site_total_reads.index.str.replace(r"\.N\d+$", "", regex=True)
site_total_reads = site_total_reads.groupby("site.id").sum().sum(axis=1)
relative_abundance = site_group_counts.div(site_total_reads, axis=0)
relative_abundance.index.name = "site.id"
relative_abundance.to_csv(output_folder / "functional_group_relative_abundance.csv")
relative_abundance.sum(axis=1).rename("proportion_sum").to_csv(output_folder / "functional_group_proportion_check.csv")

# =====================================================
# STEP 4: Join environment, plot, and calculate correlations
# =====================================================
metadata = pd.read_csv(data_folder / "Final_data_with_diversity_prefixed.csv").set_index("site.id")
plot_data = relative_abundance.join(metadata[[environmental_predictor]]).dropna()
plot_data.sort_values(environmental_predictor)[list(functional_groups)].plot(kind="bar", stacked=True, figsize=(12, 5), colormap="viridis", width=.9)
plt.ylabel("Relative read abundance (of all ITS reads)"); plt.xlabel("Sample, ordered by " + environmental_predictor); plt.title("ITS functional groups")
plt.tight_layout(); plt.savefig(output_folder / "functional_groups_stacked.png", dpi=300); plt.savefig(output_folder / "functional_groups_stacked.pdf"); plt.close()
correlations = []
for group in functional_groups:
    rho, p_value = spearmanr(plot_data[environmental_predictor], plot_data[group])
    correlations.append({"functional_group": group, "predictor": environmental_predictor, "n_samples": len(plot_data), "spearman_rho": rho, "p_value": p_value})
correlation_table = pd.DataFrame(correlations)
correlation_table["BH_FDR"] = multipletests(correlation_table["p_value"], method="fdr_bh")[1]
correlation_table.to_csv(output_folder / "functional_group_spearman.csv", index=False)
