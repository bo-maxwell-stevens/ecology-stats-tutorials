# iNEXT is an R package. This script is intentionally separate from Python so
# that it uses the published iNEXT implementation rather than an approximation.
library(iNEXT)
library(readr)

input_file <- "Data/ITS_OTU_table_final.tsv"  # Change to AMF_OTU_table_final.tsv if wanted.
output_file <- "Outputs/Diversity/ITS_inext_coverage.csv"
target_coverage <- 0.90

counts <- read_tsv(input_file, show_col_types = FALSE)
sample_ids <- counts[[1]]
count_matrix <- as.matrix(counts[, -1])
storage.mode(count_matrix) <- "numeric"

# iNEXT expects a list of species-abundance vectors, one vector per sample.
sample_list <- lapply(seq_len(nrow(count_matrix)), function(i) count_matrix[i, count_matrix[i, ] > 0])
names(sample_list) <- sample_ids

# The endpoint is twice each sample's observed read depth. Samples that do not
# reach 90% coverage are recorded as missing rather than silently extrapolated.
results <- lapply(names(sample_list), function(id) {
  reads <- sum(sample_list[[id]])
  estimate <- iNEXT(sample_list[id], q = 0, datatype = "abundance", endpoint = 2 * reads)
  table <- estimate$iNextEst
  eligible <- table[table$SC >= target_coverage, ]
  if (nrow(eligible) == 0) return(data.frame(sample_id = id, richness_at_90_coverage = NA, reaches_target = FALSE))
  first <- eligible[which.min(eligible$m), ]
  data.frame(sample_id = id, richness_at_90_coverage = first$qD, reaches_target = TRUE)
})
write_csv(do.call(rbind, results), output_file)
