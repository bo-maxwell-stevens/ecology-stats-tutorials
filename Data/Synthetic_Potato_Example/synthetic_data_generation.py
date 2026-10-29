"""Create a completely fictional potato-soil teaching dataset.

The simulation has 60 fields, two biological soil samples per field, field effects,
within-field variation, uneven sequencing depth, and observed multinomial reads.
It does not use or modify any real research records.
"""
from pathlib import Path
import numpy as np
import pandas as pd

SEED = 20261029
N_FIELDS, REPLICATES, N_AMF, N_ITS = 60, 2, 80, 200
OUT = Path(__file__).parent
rng = np.random.default_rng(SEED)

regions = ["North_Fjord", "East_Plains", "South_Moraine", "West_Heath"]
region_centres = [(56.9, 9.1), (55.7, 11.2), (55.2, 9.8), (56.2, 8.3)]  # Fictional teaching coordinates.
rotations = ["short_rotation", "diverse_rotation", "cereal_dominant"]
cultivars = ["Fictional_A", "Fictional_B", "Fictional_C", "Fictional_D"]

fields = []
for i in range(N_FIELDS):
    region_index = i % 4
    latitude, longitude = region_centres[region_index]
    pH = rng.normal(6.1 + .18 * region_index, .45)
    phosphorus = max(.2, 2.0 + .55 * pH + rng.normal(0, .75))
    organic_carbon = max(.5, rng.normal(2.4 - .12 * pH, .35))
    rotation = rotations[i % len(rotations)]
    management = "organic" if i % 5 < 2 else "conventional"
    rotation_effect = {"short_rotation": 8, "diverse_rotation": -5, "cereal_dominant": 2}[rotation]
    pathogen = rng.normal(0, 1)
    black_dot = np.clip(38 - 3.0 * pH + rotation_effect + 8 * pathogen + rng.normal(0, 10), 0, 100)
    silver_scurf = np.clip(25 + 5 * pathogen + rng.normal(0, 11), 0, 100)
    fields.append(dict(field_id=f"Field_{i + 1:03d}", region=regions[region_index], lat=latitude + rng.normal(0, .10), lon=longitude + rng.normal(0, .12), pH_KCl=pH, P_Mehlich3_mg_100g=phosphorus, Mg_mg_100g=max(2, rng.normal(16 + phosphorus, 4)), CN_ratio=max(7, rng.normal(12 - .3 * pH, 1.1)), SOC_pct=organic_carbon, soil_moisture_pct=np.clip(rng.normal(24 + 1.5 * organic_carbon, 4), 8, 45), crop_rotation=rotation, previous_crop=["cereal", "legume", "oilseed"][i % 3], cultivar=cultivars[(i * 3) % 4], management=management, irrigation_per_month=int(rng.integers(0, 7)), cover_crop="yes" if i % 3 else "no", fertilizer_category=["low", "medium", "high"][i % 3], field_effect=rng.normal(0, .55), Severity_Har_Cc=black_dot, Severity_Har_Hs=silver_scurf))

field_table = pd.DataFrame(fields)
field_table["Severity_60_Cc"] = np.clip(field_table.Severity_Har_Cc + rng.normal(8, 7, N_FIELDS), 0, 100)
field_table["Severity_60_Hs"] = np.clip(field_table.Severity_Har_Hs + rng.normal(12, 8, N_FIELDS), 0, 100)
for pathogen, severity in [("Cc", "Severity_Har_Cc"), ("Hs", "Severity_Har_Hs")]:
    field_table[f"qPCR_Har_{pathogen}"] = np.exp(-1 + .045 * field_table[severity] + rng.normal(0, .45, N_FIELDS))
    field_table[f"qPCR_60_{pathogen}"] = np.exp(-1 + .045 * field_table[f"Severity_60_{pathogen}"] + rng.normal(0, .45, N_FIELDS))

metadata_rows = []
for _, field in field_table.iterrows():
    for replicate in range(1, REPLICATES + 1):
        row = field.drop(labels=["field_effect"]).to_dict()
        row.update(sample_id=f"{field.field_id}.N{replicate}", canonical=f"{field.field_id}.N{replicate}", **{"site.id": field.field_id}, biological_replicate=replicate)
        row["pH_KCl"] += rng.normal(0, .10); row["soil_moisture_pct"] += rng.normal(0, 1.5)
        row["alpha"] = np.clip(42 + 5 * row["pH_KCl"] - 1.5 * row["P_Mehlich3_mg_100g"] + rng.normal(0, 8), 5, 100)
        metadata_rows.append(row)
metadata = pd.DataFrame(metadata_rows)

def depths(n):
    categories = rng.choice(4, n, p=[.10, .20, .40, .30])
    ranges = [(50, 300), (300, 1500), (1500, 5000), (5000, 15000)]
    result = np.array([rng.integers(*ranges[c]) for c in categories])
    result[rng.choice(n, 3, replace=False)] = 0
    return result

def community(n_taxa, prefix, depth_values, disease_direction):
    base = rng.dirichlet(np.linspace(3, .15, n_taxa))
    taxon_pH = rng.normal(0, .22, n_taxa); taxon_disease = rng.normal(0, .12, n_taxa)
    taxon_disease[:5] = disease_direction * .55; taxon_disease[5:10] = -disease_direction * .45
    # Each field has a taxon-specific composition shift shared by its two soil samples.
    # Unlike one common constant, this does not cancel when relative abundances are normalized.
    field_taxon_effects = rng.normal(0, .85, (N_FIELDS, n_taxa))
    table = []
    for sample_index, sample in metadata.iterrows():
        field_index = int(sample["site.id"].split("_")[1]) - 1
        log_abundance = np.log(base + 1e-12) + taxon_pH * (sample.pH_KCl - 6) + taxon_disease * ((sample.Severity_Har_Cc - 35) / 20) + field_taxon_effects[field_index] + rng.normal(0, .10, n_taxa)
        probabilities = np.exp(log_abundance - log_abundance.max()); probabilities /= probabilities.sum()
        table.append(rng.multinomial(depth_values[sample_index], probabilities))
    return pd.DataFrame(table, index=metadata.sample_id, columns=[f"{prefix}{i:03d}" for i in range(1, n_taxa + 1)])

amf = community(N_AMF, "VTX_SYN", depths(len(metadata)), -1)
its = community(N_ITS, "SH_SYN", depths(len(metadata)), 1)
amf.to_csv(OUT / "AMF_OTU_table_final.tsv", sep="\t")
its.to_csv(OUT / "ITS_OTU_table_final.tsv", sep="\t")
metadata.to_csv(OUT / "Final_data_with_diversity_prefixed.csv", index=False)

lifestyles = ["soil_saprotroph", "plant_pathogen", "ectomycorrhizal", "root_endophyte", "arbuscular_mycorrhizal", "unspecified"]
feature_rows = []
for i, taxon in enumerate(its.columns):
    lifestyle = lifestyles[i % len(lifestyles)] if i % 11 else ""
    feature_rows.append({"SH": taxon, "taxonomy": f"k__Fungi;p__SyntheticPhylum;c__SyntheticClass;o__SyntheticOrder;f__SyntheticFamily;g__SyntheticGenus_{i % 25};s__Synthetic_taxon_{i + 1}", "primary_lifestyle": lifestyle, "Secondary_lifestyle": "foliar_endophyte" if i % 17 == 0 else "", "Specific_hosts": "potato-associated" if lifestyle == "plant_pathogen" else "", "Growth_form_template": "filamentous_mycelium", "Aquatic_habitat_template": "terrestrial", "annotation_confidence": ["high", "medium", "low"][i % 3]})
pd.DataFrame(feature_rows).to_csv(OUT / "ITS_feature_metadata.tsv", sep="\t", index=False)
print(f"Wrote {len(metadata)} samples from {N_FIELDS} fictional fields with seed {SEED}.")
