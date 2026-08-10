analysed_genes_file <- "all_analysed_genes.txt"
cactus12_file <- "genes_in_cactus4_no_Loc.txt"

# Read gene lists
analysed_genes <- readLines(analysed_genes_file)
cactus12_genes <- readLines(cactus12_file)

# Create data frame
df <- data.frame(
  gene = analysed_genes,
  score = 0
)

# Assign 1 to genes present in cactus
df$score[df$gene %in% cactus12_genes] <- 1

# Format gene names: capitalize only the first letter
df$gene <- paste0(
  toupper(substr(df$gene, 1, 1)),
  tolower(substr(df$gene, 2, nchar(df$gene)))
)

# Write output
write.table(
  df,
  "gene_scores_matrix.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Preview
head(df)
