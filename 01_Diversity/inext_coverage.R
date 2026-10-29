# iNEXT is an R package. This script is intentionally separate from Python so
# that it uses the published iNEXT implementation rather than an approximation.
library(iNEXT)
library(readr)
library(ggplot2)

input_file <- "Data/Synthetic_Potato_Example/ITS_OTU_table_final.tsv"  # Change to Data/ITS_OTU_table_final.tsv for permitted DarkDivNet data.
output_file <- "Outputs/Diversity/ITS_inext_coverage.csv"
target_coverage <- 0.90
n_bootstrap <- 0  # Point estimates only: much faster for many high-dimensional community samples.

counts <- read_tsv(input_file, show_col_types = FALSE)
sample_ids <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
storage.mode(count_matrix) <- "numeric"
if (any(!is.finite(count_matrix)) || any(count_matrix < 0)) stop("Count table must contain finite, non-negative values.")

# iNEXT expects a list of species-abundance vectors, one vector per sample.
sample_list <- lapply(seq_len(nrow(count_matrix)), function(i) count_matrix[i, count_matrix[i, ] > 0])
names(sample_list) <- sample_ids
nonempty <- lengths(sample_list) > 0
if (any(!nonempty)) write_csv(data.frame(sample_id = sample_ids[!nonempty], reason = "zero_reads"), "Outputs/Diversity/ITS_inext_zero_read_excluded.csv")
sample_list <- sample_list[nonempty]

# The endpoint is twice each sample's observed read depth. Samples that do not
# reach 90% coverage are recorded as missing rather than silently extrapolated.
results <- lapply(names(sample_list), function(id) {
  reads <- sum(sample_list[[id]])
  estimate <- iNEXT(sample_list[id], q = 0, datatype = "abundance", endpoint = 2 * reads, nboot = n_bootstrap)
  # iNEXT stores size- and coverage-based estimates in a named list.
  # Size-based estimates let us enforce the twice-observed-read endpoint.
  table <- estimate$iNextEst$size_based
  eligible <- table[table$SC >= target_coverage, ]
  if (nrow(eligible) == 0) return(data.frame(sample_id = id, richness_at_90_coverage = NA, reaches_target = FALSE))
  first <- eligible[which.min(eligible$m), ]
  data.frame(sample_id = id, richness_at_90_coverage = first$qD, reaches_target = TRUE)
})
coverage_results <- do.call(rbind, results)
observed_richness <- data.frame(sample_id = names(sample_list), observed_richness = lengths(sample_list))
coverage_results <- merge(coverage_results, observed_richness, by = "sample_id", all.x = TRUE)
write_csv(coverage_results, output_file)

# Only samples that reached the target within twice their read depth can be
# compared on the same coverage basis.
plot_data <- coverage_results[coverage_results$reaches_target, ]
if (nrow(plot_data) > 0) {
  plot <- ggplot(plot_data, aes(x = observed_richness, y = richness_at_90_coverage)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Observed richness", y = "Richness at 90% coverage", title = "Observed and coverage-standardized richness") +
    theme_bw()
  ggsave("Outputs/Diversity/ITS_observed_vs_coverage_richness.png", plot, width = 7, height = 5, dpi = 300)
  ggsave("Outputs/Diversity/ITS_observed_vs_coverage_richness.pdf", plot, width = 7, height = 5)
}
