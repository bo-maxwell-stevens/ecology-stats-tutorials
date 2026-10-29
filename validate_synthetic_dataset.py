"""Readable checks for the completely synthetic potato teaching dataset."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import braycurtis
import matplotlib.pyplot as plt

data = Path("Data/Synthetic_Potato_Example")
out = Path("Outputs/Synthetic/Validation")
out.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(data / "Final_data_with_diversity_prefixed.csv").set_index("sample_id")
amf = pd.read_csv(data / "AMF_OTU_table_final.tsv", sep="\t", index_col=0)
its = pd.read_csv(data / "ITS_OTU_table_final.tsv", sep="\t", index_col=0)
features = pd.read_csv(data / "ITS_feature_metadata.tsv", sep="\t")

assert len(metadata) == 120 and metadata.index.is_unique and metadata["site.id"].nunique() == 60
assert amf.shape == (120, 80) and its.shape == (120, 200) and len(features) == 200
assert amf.index.equals(metadata.index) and its.index.equals(metadata.index)
assert (amf.to_numpy() >= 0).all() and (its.to_numpy() >= 0).all()
assert np.equal(amf.to_numpy(), np.floor(amf.to_numpy())).all()
assert metadata.filter(like="Severity").apply(lambda x: x.between(0, 100).all()).all()
assert features["SH"].is_unique and features["SH"].isin(its.columns).all()

# Bray-Curtis ignores zero-read libraries because their composition is undefined.
nonempty = its.loc[its.sum(axis=1) > 0]
within, between = [], []
for i, first_id in enumerate(nonempty.index):
    for second_id in nonempty.index[i + 1:]:
        value = braycurtis(nonempty.loc[first_id], nonempty.loc[second_id])
        (within if metadata.loc[first_id, "site.id"] == metadata.loc[second_id, "site.id"] else between).append(value)
summary = pd.DataFrame({"comparison": ["within_field", "between_field"], "n_pairs": [len(within), len(between)], "median_bray_curtis": [np.median(within), np.median(between)]})
assert np.median(within) < np.median(between), "Replicates should be more similar than different fields."
summary.to_csv(out / "synthetic_validation_summary.csv", index=False)
plt.boxplot([within, between], tick_labels=["Within field", "Between fields"])
plt.ylabel("ITS Bray-Curtis dissimilarity")
plt.tight_layout(); plt.savefig(out / "within_vs_between_bray_curtis.png", dpi=300); plt.close()
print(summary.to_string(index=False))
