args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) stop("Expected community, output folder, predictors, permutations, categorical predictor, blocks, and seed.")
community_name <- args[1]; output_folder <- args[2]
predictors <- strsplit(args[3], ",")[[1]]; permutations <- as.integer(args[4]); categorical <- args[5]; blocks_file <- args[6]; set.seed(as.integer(args[7]))
if (!requireNamespace("vegan", quietly = TRUE)) stop("Install R package vegan before running PERMANOVA.")
library(vegan)
community <- read.csv(file.path(output_folder, paste0(community_name, "_community_for_R.csv")), row.names = 1, check.names = FALSE)
metadata <- read.csv(file.path(output_folder, paste0(community_name, "_metadata_for_R.csv")), row.names = 1, check.names = FALSE)
# Square-root relative abundance reduces the effect of dominant taxa before Bray-Curtis.
if (!identical(rownames(community), rownames(metadata))) stop("Community and metadata sample IDs are not aligned.")
if (any(rowSums(community) <= 0)) stop("Zero-read communities must be excluded before PERMANOVA.")
community <- sqrt(community / rowSums(community))
formula <- as.formula(paste("community ~", paste(predictors, collapse = " + ")))
# by="margin" asks for marginal terms. Omit it for sequential (order-dependent) tests.
permutation_control <- permutations
if (blocks_file != "NONE") permutation_control <- how(nperm = permutations, blocks = read.csv(blocks_file, row.names = 1)[, 1])
result <- adonis2(formula, data = metadata, method = "bray", permutations = permutation_control, by = "margin")
result_table <- as.data.frame(result); result_table$n_samples <- nrow(community); result_table$permutations <- permutations
write.csv(result_table, file.path(output_folder, paste0(community_name, "_permanova_marginal.csv")))
if (categorical != "NONE") {
  group <- metadata[[categorical]]
  if (all(table(group) >= 2)) {
    dispersion <- betadisper(vegdist(community, method = "bray"), group)
    write.csv(as.data.frame(permutest(dispersion, permutations = permutations)$tab), file.path(output_folder, paste0(community_name, "_dispersion_", categorical, ".csv")))
  } else message("Dispersion test skipped: every category needs at least two samples.")
}
