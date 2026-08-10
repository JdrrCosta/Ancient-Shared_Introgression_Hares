# Load your 2000 random percentage overlaps
random_perc <- scan("cactus4_bootstrap_gene_number_conservative.txt")  # Make sure it's just one number per line

# Your observed empirical overlap percentage
empirical_perc <- 12916  # Replace with your actual empirical percentage

svg("cactus4_gene_conservative.svg", width = 8, height = 6)
# Plot distribution with the empirical line
hist_data <- hist(random_perc, breaks = 50, plot = FALSE)  # Get the data for the histogram without plotting
hist(random_perc, breaks = 50, col = "lightgray",
     main = "Number of genes: Null Distribution vs. Empirical",
     xlab = "Gene count", xlim = range(c(random_perc, empirical_perc)),
     ylab = "Bootstrap Counts")
abline(v = empirical_perc, col = "red", lwd = 2)
legend("topright", legend = "Empirical", col = "red", lwd = 2)

# Calculate p-value (two-sided test)
mean_random <- mean(random_perc)
p_value_two_sided <- sum(abs(random_perc - mean_random) >= abs(empirical_perc - mean_random)) / length(random_perc)
cat("Two-sided empirical p-value:", p_value_two_sided, "\n")

# Print p-value on the plot
text(x = empirical_perc + 45, 
     y = max(hist_data$counts) * 0.9,  # Set y position near the top of the histogram
     labels = paste("p =", round(p_value_two_sided, 3)), 
     col = "black", pos = 4, cex = 1.2)

dev.off()
