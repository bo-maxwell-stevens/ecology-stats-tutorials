args <- commandArgs(trailingOnly = TRUE)
community_name <- args[1]; output_folder <- args[2]
predictors <- strsplit(args[3], ",")[[1]]; permutations <- as.integer(args[4]); categorical <- args[5]
library(vegan)
community <- read.csv(file.path(output_folder, paste0(community_name, "_community_for_R.csv")), row.names = 1, check.names = FALSE)
metadata <- read.csv(file.path(output_folder, paste0(community_name, "_metadata_for_R.csv")), row.names = 1, check.names = FALSE)
# Square-root relative abundance reduces the effect of dominant taxa before Bray-Curtis.
community <- sqrt(community / rowSums(community))
formula <- as.formula(paste("community ~", paste(predictors, collapse = " + ")))
# by="margin" asks for marginal terms. Omit it for sequential (order-dependent) tests.
result <- adonis2(formula, data = metadata, method = "bray", permutations = permutations, by = "margin")
write.csv(as.data.frame(result), file.path(output_folder, paste0(community_name, "_permanova_marginal.csv")))
if (categorical != "NONE") {
  group <- metadata[[categorical]]
  if (all(table(group) >= 2)) {
    dispersion <- betadisper(vegdist(community, method = "bray"), group)
    write.csv(as.data.frame(permutest(dispersion, permutations = permutations)$tab), file.path(output_folder, paste0(community_name, "_dispersion_", categorical, ".csv")))
  } else message("Dispersion test skipped: every category needs at least two samples.")
}
