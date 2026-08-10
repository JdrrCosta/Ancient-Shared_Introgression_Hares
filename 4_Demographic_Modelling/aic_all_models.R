# Define model results manually (replace with your values)
# logLik = marginal log-likelihood from BPP
# k = number of free parameters in each model (count theta, tau, M, etc.)
models <- data.frame(
  model = c("MSC", "MSC-I_proto", "MSC-I_proto_cas", "MSC-I_proto_cor","MSC-I_proto_cascor","MSC-I_only_cas","MSC-I_only_cor",
            "MSC-I_only_cascor","MSC-M","MSC-M_proto_cas","MSC-M_proto_cascor","MSC-I_proto_cascor_more","MSC-I_proto_cascor_flipped"),
  logLik = c(-1595106.060000, -1595069.441000, -1594621.709500, -1594754.294000,-1594533.766500,-1594831.321000,-1594867.830500,
             -1594716.236500,-1594701.163500,-1594637.187500,-1594668.706500,-1594395.269000,-1594361.415000),
  k = c(7, 12, 17, 17,22, 12,12,17,9,11,13,22,22)   # <-- fill in parameter counts for each model
)
# Calculate AIC
models$AIC <- 2 * models$k - 2 * models$logLik

# Compare models by ΔAIC
models$deltaAIC <- models$AIC - min(models$AIC)

# Print results
print(models)
write.csv(models, file = "Models_aic_comparison_final.csv")
