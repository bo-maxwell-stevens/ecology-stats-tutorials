# Ordination

NMDS preserves ranked dissimilarities; inspect its stress, because a high value means a poor two-dimensional representation. PCoA represents distances directly; its axes report explained variation and negative eigenvalues warn that Bray-Curtis is not perfectly Euclidean.

The script uses square-root relative abundance to reduce the influence of highly abundant taxa. Run `python 02_Ordination/nmds_pcoa.py`. Change `colour_variable` to a real metadata variable, or a disease/management column in your own data. Documentation: [scikit-bio](https://scikit.bio/docs/latest/), [scikit-learn MDS](https://scikit-learn.org/stable/modules/manifold.html#multidimensional-scaling).
