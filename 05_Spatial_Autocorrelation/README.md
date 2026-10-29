# Spatial autocorrelation

Moran's I is positive when neighbouring sites are similar, negative when they are dissimilar, and near zero when no spatial pattern is detected. The script uses approximate kilometre coordinates and row-standardised 5-nearest-neighbour weights. Repeated coordinates are aggregated before analysis.

Run `python 05_Spatial_Autocorrelation/morans_i.py`. Add disease severity to `variables_to_test` when available. Also examine residual spatial autocorrelation after fitting models. Documentation: [PySAL](https://pysal.org/), [esda Moran](https://pysal.org/esda/generated/esda.Moran.html).
