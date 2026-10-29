# Spatial autocorrelation

Moran's I is positive when neighbouring sites are similar, negative when they are dissimilar, and near zero when no spatial pattern is detected. The script projects coordinates before using row-standardised k-nearest-neighbour weights and reports a sensitivity check at k = 4, 6, and 8. Treat these multiple exploratory tests cautiously rather than selecting the most significant result.

Run `python 05_Spatial_Autocorrelation/morans_i.py`. Add disease severity to `variables_to_test` when available. Also examine residual spatial autocorrelation after fitting models. Documentation: [PySAL](https://pysal.org/), [esda Moran](https://pysal.org/esda/generated/esda.Moran.html).
