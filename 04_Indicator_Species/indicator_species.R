args <- commandArgs(trailingOnly = TRUE); output_folder <- args[1]; permutations <- as.integer(args[2])
library(indicspecies); library(ggplot2)
community <- read.csv(file.path(output_folder, "community_for_R.csv"), row.names = 1, check.names = FALSE)
groups <- read.csv(file.path(output_folder, "groups_for_R.csv"), row.names = 1)$analysis_group
result <- multipatt(community, groups, func = "IndVal.g", duleg = TRUE, control = how(nperm = permutations))
table <- result$sign
table$taxon <- rownames(table); table$FDR <- p.adjust(table$p.value, method = "BH")
table$associated_group <- apply(table[, grep("^s\\.", names(table)), drop = FALSE], 1, function(x) paste(names(x)[x > 0], collapse = ";"))
write.csv(table, file.path(output_folder, "indicator_species_results.csv"), row.names = FALSE)
significant <- table[table$FDR < 0.05, ]
if (nrow(significant) > 0) {
  significant <- significant[order(significant$stat), ]; significant <- tail(significant, 20)
  plot <- ggplot(significant, aes(x = reorder(taxon, stat), y = stat, fill = associated_group)) + geom_col() + coord_flip() + labs(x = "Taxon", y = "IndVal.g", title = "FDR-significant indicator taxa") + theme_bw()
  ggsave(file.path(output_folder, "significant_indicator_taxa.png"), plot, width = 9, height = 6, dpi = 300)
  ggsave(file.path(output_folder, "significant_indicator_taxa.pdf"), plot, width = 9, height = 6)
} else message("No FDR-significant taxa; no figure created.")
