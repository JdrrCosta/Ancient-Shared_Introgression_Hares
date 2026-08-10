#code with cytoscape plot
library(dplyr)
library(signet)
library(graphite)

# ------------------------------
# Load mouse KEGG pathways and convert to SYMBOLs
# ------------------------------
paths <- graphite::pathways("mmusculus", "kegg")
paths_symbol <- lapply(paths, convertIdentifiers, to = "symbol")
kegg_mice <- lapply(paths_symbol, graphite::pathwayGraph)

# ------------------------------
# Load gene score table (symbols + scores)
# ------------------------------
scores <- read.table("gene_scores_matrix_4_SYMBOL.txt", header = TRUE, stringsAsFactors = FALSE)

# Quick diagnostic: overlap with pathways
empty_paths <- sapply(kegg_mice, function(p) length(intersect(nodes(p), scores$gene)) == 0)
cat("Fraction of empty pathways:", sum(empty_paths) / length(kegg_mice), "\n")

# ------------------------------
# Run simulated annealing search
# ------------------------------
HSS <- searchSubnet(kegg_mice, scores)

# Generate null distribution (can take time)
null <- nullDist(kegg_mice, scores, n = 10000)

# Test subnetworks
HSS <- testSubnet(HSS, null)

# ------------------------------
# Summarize and save results
# ------------------------------
tab <- summary(HSS)
write.table(tab, "Cactus_full_genes_Signet_output_kegg.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------
# CLEAN UP SUBNETWORKS BEFORE EXPORT
# ------------------------------
for (i in seq_along(HSS@results)) {
  sigObj <- HSS@results[[i]]
  
  if (!is.null(sigObj@subnet_genes) && !is.null(sigObj@connected_comp)) {
    # Remove NA values and nodes not present in the graph
    valid_nodes <- na.omit(sigObj@subnet_genes)
    valid_nodes <- intersect(valid_nodes, nodes(sigObj@connected_comp))
    
    # Assign as factor to match Signet object requirement
    sigObj@subnet_genes <- factor(valid_nodes)
    
    # Put it back into the results list
    HSS@results[[i]] <- sigObj
  }
}

# Check for any remaining NA nodes
has_na <- sapply(HSS@results, function(x) anyNA(x@subnet_genes))
cat("Subnets with NA nodes remaining:", sum(has_na), "\n")

# ------------------------------
# Export significant subnetworks for Cytoscape
# ------------------------------
writeXGMML(HSS, filename = "Cactus_full_genes_Cytoscape_input_kegg.xgmml", threshold = 0.01)
