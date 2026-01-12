
source('UPDATEDRMEPMTV.R')

set.seed(123)

# this has p = 1, k=2, with different pi.
# Generate parameters using MixSim
Params <- MixSim(BarOmega = 0.01, K = 2, p = 1, hom = TRUE)

# Extract means (Mu) and covariance matrix (S)
Mu <- Params$Mu
CovarianceMatrix <- Params$S


Sw = 0.1*Params$S[,,1]

m <- 200
n_i <- 10
B <- 1


# SIM 1

pi <- c(0.7, 0.2, 0.1)


sim_results <- Simulate_compare_general(pi = pi, Mu = Mu, S = CovarianceMatrix, Sw = Sw, m = m, n_i = n_i, B = B)

# Summarize the results with means and standard deviations 
summary_results <- summarize(sim_results)


write(t(sim_results$RMEP_Results), file = 'Results/results_p1_k2_M1_A100RMEP.txt', ncol = 8)
write(t(sim_results$RMED_Results), file = 'Results/results_p1_k2_M1_A100RMED.txt', ncol = 5)
write(t(summary_results$mean_RMEP), file = 'Results/mean_results_p1_k2_A100M1_RMEP.txt', ncol = 8)
write(t(summary_results$sd_RMEP), file = 'Results/sd_results_p1_k2_M1_A100RMEP.txt', ncol = 8)
write(t(summary_results$mean_RMED), file = 'Results/mean_results_p1_k2_A100M1_RMED.txt', ncol = 5)
write(t(summary_results$sd_RMED), file = 'Results/sd_results_p1_k2_M1_A100RMED.txt', ncol = 5)



# SIM 2

pi <- c(0.5, 0.4, 0.1)


sim_results <- Simulate_compare_general(pi = pi, Mu = Mu, S = CovarianceMatrix, Sw = Sw, m = m, n_i = n_i, B = B)

# Summarize the results with means and standard deviations 
summary_results <- summarize(sim_results)


write(t(sim_results$RMEP_Results), file = 'Results/results_p1_k2_M2_A100RMEP.txt', ncol = 8)
write(t(sim_results$RMED_Results), file = 'Results/results_p1_k2_M2_A100RMED.txt', ncol = 5)
write(t(summary_results$mean_RMEP), file = 'Results/mean_results_p1_k2_A100M2_RMEP.txt', ncol = 8)
write(t(summary_results$sd_RMEP), file = 'Results/sd_results_p1_k2_M2_A100RMEP.txt', ncol = 8)
write(t(summary_results$mean_RMED), file = 'Results/mean_results_p1_k2_M2_A100RMED.txt', ncol = 5)
write(t(summary_results$sd_RMED), file = 'Results/sd_results_p1_k2_M2_A100RMED.txt', ncol = 5)


# SIM 3

pi <- c(0.6, 0.3, 0.1)


sim_results <- Simulate_compare_general(pi = pi, Mu = Mu, S = CovarianceMatrix, Sw = Sw, m = m, n_i = n_i, B = B)

# Summarize the results with means and standard deviations 
summary_results <- summarize(sim_results)


write(t(sim_results$RMEP_Results), file = 'Results/results_p1_k2_M3_A100RMEP.txt', ncol = 8)
write(t(sim_results$RMED_Results), file = 'Results/results_p1_k2_M3_A100RMED.txt', ncol = 5)
write(t(summary_results$mean_RMEP), file = 'Results/mean_results_p1_A100k2_M3_RMEP.txt', ncol = 8)
write(t(summary_results$sd_RMEP), file = 'Results/sd_results_p1_k2_M3_A100RMEP.txt', ncol = 8)
write(t(summary_results$mean_RMED), file = 'Results/mean_results_p1_k2_M3_A100RMED.txt', ncol = 5)
write(t(summary_results$sd_RMED), file = 'Results/sd_results_p1_k2_M3_A100RMED.txt', ncol = 5)


# SIM 4

pi <- c(0.6, 0.4, 0.2)


sim_results <- Simulate_compare_general(pi = pi, Mu = Mu, S = CovarianceMatrix, Sw = Sw, m = m, n_i = n_i, B = B)

# Summarize the results with means and standard deviations 
summary_results <- summarize(sim_results)


write(t(sim_results$RMEP_Results), file = 'Results/results_p1_k2_M4_A100RMEP.txt', ncol = 8)
write(t(sim_results$RMED_Results), file = 'Results/results_p1_k2_M4_A100RMED.txt', ncol = 5)
write(t(summary_results$mean_RMEP), file = 'Results/mean_results_p1_k2_M4_A100RMEP.txt', ncol = 8)
write(t(summary_results$sd_RMEP), file = 'Results/sd_results_p1_k2_M4_A100RMEP.txt', ncol = 8)
write(t(summary_results$mean_RMED), file = 'Results/mean_results_p1_k2_M4_A100RMED.txt', ncol = 5)
write(t(summary_results$sd_RMED), file = 'Results/sd_results_p1_k2_M4_A100RMED.txt', ncol = 5)


# SIM 5

pi <- c(0.33, 0.33)


sim_results <- Simulate_compare_general(pi = pi, Mu = Mu, S = CovarianceMatrix, Sw = Sw, m = m, n_i = n_i, B = B)

# Summarize the results with means and standard deviations 
summary_results <- summarize(sim_results)


write(t(sim_results$RMEP_Results), file = 'Results/results_p1_k2_M5_A100RMEP.txt', ncol = 8)
write(t(sim_results$RMED_Results), file = 'Results/results_p1_k2_M5_A100RMED.txt', ncol = 5)
write(t(summary_results$mean_RMEP), file = 'Results/mean_results_p1_k2_M5_A100RMEP.txt', ncol = 8)
write(t(summary_results$sd_RMEP), file = 'Results/sd_results_p1_k2_M5_A100RMEP.txt', ncol = 8)
write(t(summary_results$mean_RMED), file = 'Results/mean_results_p1_k2_M5_A100RMED.txt', ncol = 5)
write(t(summary_results$sd_RMED), file = 'Results/sd_results_p1_k2_M5_A100RMED.txt', ncol = 5)




