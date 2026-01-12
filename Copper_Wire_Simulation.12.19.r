#######################
# Author:  Isaac Gbene & Dylan
# Generating_Copper_Data.R
# Generating heirarchical synthetic Chemometric data based on
# copper wire models from:
# Dettman, Joshua R. and Cassabaum, Alyssa A. and Saunders, Christopher P. and Snyder, Deanna L. and Buscaglia, JoAnn
# Forensic Discrimination of Copper Wire Using Trace Element Concentrations
# Analytical Chemistry (2014)
#last updated 12/12/2025
#######################

library(MixSim)
library(mvtnorm)
library(dplyr)
#install.packages("comparison")
library(comparison)
library(ggplot2)
set.seed(12345)  # For overall reproducibility

# order of trace elements in vectors
labels <- c('Ag', 'Sb', 'Pb', 'Bi', 'Co', 'Ni', 'As', 'Se')
#Ag (Silver)
#Sb (Antimony)
#Pb (Lead)
#Bi (Bismuth)
#Co (Cobalt)
#Ni (Nickel)

# average concentration vectors for groups 1 and 2
mu_hat_1 <- c(7.429, 0.761, 0.774, 0.070, 0.030, 2.033, 0.917, 0.756)
mu_hat_2 <- c(4.524, 0.465, 0.687, 0.037, 0.28, 1.065, 0.512, 0.536)

mu_hat <- rbind(mu_hat_1, mu_hat_2)

# Between Source Covariance Matrix
S_b_vals <- c(5.83e-1, 4.33e-2, -4.17e-2, 3.61e-3, -1.15e-3, 8.30e-2, 3.47e-2, 4.89e-2,
              4.33e-2, 5.54e-3, -1.01e-2, 4.17e-4, -2.52e-4, 8.44e-3, 4.04e-3, 4.34e-3,
              -4.17e-2, -1.01e-2, 1.53e-1, -8.82e-4, 2.09e-3, -1.49e-2, -8.78e-3, -1.65e-2,
              3.61e-3, 4.17e-4, -8.82e-4, 7.80e-5, -1.66e-5, 7.84e-4, 3.09e-4, 2.53e-4,
              -1.15e-3, -2.52e-4, 2.09e-3, -1.66e-5, 7.04e-5, -4.66e-4, -2.85e-4, -2.43e-4,
              8.30e-2, 8.44e-3, -1.49e-2, 7.84e-4, -4.66e-4, 5.57e-2, 1.47e-2, 2.91e-3,
              3.47e-2, 4.04e-3, -8.78e-3, 3.09e-4, -2.85e-4, 1.47e-2, 9.25e-3, 1.53e-3,
              4.89e-2, 4.34e-3, -1.65e-2, 2.53e-4, -2.43e-4, 2.91e-3, 1.53e-3, 1.93e-2)
S_b <- matrix(S_b_vals, nrow=8, ncol=8)

# Within Source Covariance Matrix
S_w_vals <- c(4.33e-3, 4.05e-4, 3.89e-4, 1.66e-5, -2.45e-5, 9.05e-4, 6.79e-4, 4.28e-4,
              4.05e-4, 3.95e-4, 1.34e-5, -1.18e-6, 2.39e-6, 3.76e-4, 2.36e-4, 1.57e-4,
              3.89e-4, 1.34e-5, 1.71e-4, 9.30e-6, 4.89e-6, 7.65e-5, 4.46e-5, -1.37e-4,
              1.66e-5, -1.18e-6, 9.30e-6, 1.67e-6, 7.29e-7, 6.64e-6, -5.34e-6, -1.08e-5,
              -2.45e-5, 2.39e-6, 4.89e-6, 7.29e-7, 1.69e-6, 1.52e-5, -8.06e-6, -1.24e-5,
              9.05e-4, 3.76e-4, 7.65e-5, 6.64e-6, 1.52e-5, 3.47e-3, 1.49e-4, 2.20e-4,
              6.79e-4, 2.36e-4, 4.46e-5, -5.34e-6, -8.06e-6, 1.49e-4, 3.61e-4, 1.81e-4,
              4.28e-4, 1.57e-4, -1.37e-4, -1.08e-5, -1.24e-5, 2.20e-4, 1.81e-4, 2.27e-3)
S_w <- matrix(S_w_vals, nrow=8, ncol=8)

# Estimated Mixing proportions
pi_hat_1 <- 0.172
pi_hat_2 <- 1-pi_hat_1

pi_hat <- c(pi_hat_1, pi_hat_2)

# between covariance array
copper_S <- array(c(S_b_vals), dim=c(8,8,2))

sim_dat_fun2 <- function(n_a, n_w, pi, mu, S, sig_w){
  # sample n_a sources from the between source distribution
  mixture = simdataset(n_a, pi, mu, S) # S is the between source covariance matrix
  means = mixture$X
  ids = mixture$id
  dat_tmp <- NULL
  for(obj in 1:n_a){
    # draw samples from each source
    obj_samps = rmvnorm(n_w, means[obj,], sig_w)
    dat_tmp = data.frame(rbind(dat_tmp, cbind(obj, id = ids[obj], obj_samps)))
  }
  return(dat_tmp)
}


# number of sources
N_a <- 350
# number of samples within each source
N_w <- 5
set.seed(2345)
my_data <- sim_dat_fun2(n_a = N_a, 
                        n_w = N_w, 
                        pi = pi_hat, 
                        mu = mu_hat, 
                        S = copper_S, 
                        sig_w = S_w)
  
colnames(my_data) <- c('obj', 'subpop', labels)

head(my_data)
print(my_data)
dim(my_data) # Should be 1000 rows (100 sources × 10 samples) × 10 columns
# Check group distribution
table(my_data$subpop)  # Should be approximately 17.2% group 1, 82.8% group 2

table(my_data$obj)




## Real Data Algorithm and heatmap plots 

#glass_data <- read.csv('glass.txt', sep = ' ')
#my_data <- glass_data[, c(6:9,5)]
# Keep features + group + window
my_data <- my_data[, c("Ag","Sb","Pb","Bi","Co","Ni","As","Se","subpop","obj")]


train_data <- my_data

# Ensure 'Window' and 'Group' columns are correctly identified in the copper wire data
uni_sources <- unique(my_data$obj)
source_ids <- integer(length(uni_sources))
names(source_ids) <- uni_sources

for (i in seq_along(uni_sources)) {
  # Find the first index where each unique window appears and take the corresponding group
  index <- which(my_data$obj == uni_sources[i])[1]
  source_ids[i] <- my_data$subpop[index]
}




Realdatanalysis <- function(my_data, B) {
  set.seed(1235)### seed for the function 
  p <- 8  # eight concentrations 
  RMEP_Results <- NULL
  RMED_Results <- NULL
  
  
  
  for (b in 1:B) {
    K = 2# 
    # Group_sources <-list(K)
    # train_indices <- NULL
    # 
    # 
    # for ( k in 1:K){ 
    #   
    #   groupk_sources <-which(source_ids==k) #  in the little kth position store that, all the sources that came from group k 
    #   
    #   train_indices<-c(train_indices, sample(groupk_sources,length(groupk_sources)/2))
    #   
    #   Group_sources[[k]] <-groupk_sources
    # }
    #train_indices=c(sample(groupk_sources,length(groupk_sources)/2))
    
   # test_indices <- setdiff(1:62, train_indices)  # Remaining sources for testing
    # test_indices <- setdiff(1:length(uni_sources), train_indices)
    # 
    # train_data<-my_data[my_data[,(p+2)]%in%train_indices,]
    # 
    # test_data<-my_data[my_data[,(p+2)]%in%test_indices,]
    set.seed(12345 + b)  #seed for each data generation
    
    my_data <- sim_dat_fun2(n_a = N_a, 
                            n_w = N_w, 
                            pi = pi_hat, 
                            mu = mu_hat, 
                            S = copper_S, 
                            sig_w = S_w)
    
    colnames(my_data) <- c('obj', 'subpop', labels)
    
    #glass_data <- read.csv('glass.txt', sep = ' ')
    #my_data <- glass_data[, c(6:9,5)]
    # Keep features + group + window
    my_data <- my_data[, c("Ag","Sb","Pb","Bi","Co","Ni","As","Se","subpop","obj")]
    
    
    train_data <- my_data
    
    test_data <- sim_dat_fun2(n_a = 100, 
                            n_w = 5, 
                            pi = c(0.5, 0.5), 
                            mu = mu_hat, 
                            S = copper_S, 
                            sig_w = S_w)
    
    colnames(test_data) <- c('obj', 'subpop', labels)
    
    # Keep features + group + window now subpop and obj
    test_data <- test_data[, c("Ag","Sb","Pb","Bi","Co","Ni","As","Se","subpop","obj")]
    
    
    
    
    
    
    ## Step 5: Label observations in the test set
    trace_indices <- 1:2      # Then observations from each source so we choose one 1-2 for the trace indces and 3-5 for control
    control_indices <- 3:5
    
    Results <- NULL # Initialize empty data frame for results
    # Obtain the e_a_format from training data
    e_a_format <- two.level.components(train_data, data.columns = 1:p, item.column = p+2)
    
    # Running comparisons
    test_indices <- unique(test_data$obj)
    for (i in seq_along(test_indices)) {
      for (j in seq(i, length(test_indices))) {
        # Extract data for the current pair of indices
        test_i <- test_data[test_data$obj == test_indices[i], ]
        test_j <- test_data[test_data$obj == test_indices[j], ]
        
        # Assuming 'control_indices' and 'trace_indices' are predefined
        control_i <- test_i[control_indices, , drop = FALSE]
        control_j <- test_j[control_indices, , drop = FALSE]
        trace_i <- test_i[trace_indices, , drop = FALSE]
        trace_j <- test_j[trace_indices, , drop = FALSE]
        
        # Format trace and control data
        trace_i_format <- two.level.comparison.items(trace_i, 1:p)
        trace_j_format <- two.level.comparison.items(trace_j, 1:p)
        control_i_format <- two.level.comparison.items(control_i, 1:p)
        control_j_format <- two.level.comparison.items(control_j, 1:p)
        
        # Compute likelihood ratios
        lind_LR_ij <- two.level.normal.LR(trace_i_format, control_j_format, e_a_format)
        lind_LR_ji <- two.level.normal.LR(trace_j_format, control_i_format, e_a_format)
        density_LR_ij <- two.level.density.LR(trace_i_format, control_j_format, e_a_format)
        density_LR_ji <- two.level.density.LR(trace_j_format, control_i_format, e_a_format)
        
        # Group determination for i and j
        
        
        #i_group <- my_data$subpop[my_data$obj == test_indices[i]][1]
        #j_group <- my_data$subpop[my_data$obj == test_indices[j]][1]
        
        i_group <- test_i$subpop[1]
        j_group <- test_j$subpop[2]
        
        # Store results
        Results <- rbind(Results, data.frame(trace_id = test_indices[i], control_id = test_indices[j], trace_group = i_group, control_group = j_group, LR1 = lind_LR_ij, LR2 = density_LR_ij))
        Results <- rbind(Results, data.frame(trace_id = test_indices[j], control_id = test_indices[i], trace_group = j_group, control_group = i_group, LR1 = lind_LR_ji, LR2 = density_LR_ji))
      }
    }
    
    c <- 1  # According to the report/threshold 
    
    Results <-as.data.frame(Results)
    
    Temp_RME<- RME(Results = Results, c = 1,  RMEP_Results, RMED_Results, B)
    RMEP_Results <- Temp_RME$RMEP_Results
    RMED_Results<- Temp_RME$RMED_Results
    
    
    
    print(b)
    flush.console()
  }
  
  return(list(RMEP_Results = RMEP_Results, RMED_Results = RMED_Results, B = B))
}

RME <- function(Results, c = 1, RMEP_Results, RMED_Results, B){
  
  
  
  K  <- length(unique(Results$trace_group))
  RMEP_total1<-mean(Results$LR1[Results$trace_id != Results$control_id] > c) # finding the proportions of LRs that were greater than c
  RMEP_total2<-mean(Results$LR2[Results$trace_id != Results$control_id] > c) # finding the proportions of LRs that were greater than c
  RMED_total1<-mean(Results$LR1[Results$trace_id == Results$control_id] < c) # finding the proportions of LRs that were Less than c within the same group
  RMED_total2<-mean(Results$LR2[Results$trace_id == Results$control_id] < c) # finding the proportions of LRs that were Less than c within the same group
  
  RMEP_k_ALL1<- rep(NA, K)
  
  RMEP_k_ALL2<- rep(NA, K)
  
  RMED_k_Sub1<- rep(NA, K)
  
  RMED_k_Sub2<- rep(NA, K)
  
  RMEP_k_sub1<- rep(NA, K)
  
  RMEP_k_sub2<- rep(NA, K)
  for(k in 1:K){
    RMEP_k_ALL1[k] <- mean(Results$LR1[Results$trace_id != Results$control_id & Results$trace_group == k] > c)
    #The trace is coming from groupk and the control could come from group k or k + 1 finding the proportions of LRs that were greater than c
    RMEP_k_sub1[k]  <- mean(Results$LR1[Results$trace_id != Results$control_id & 
                                          Results$trace_group==k & 
                                          Results$control_group==k] > c) #Different source comparism where the trace came from group 1 and control came from group 1
    
    
    RMEP_k_ALL2[k] <- mean(Results$LR2[Results$trace_id != Results$control_id & Results$trace_group == k] > c)
    #The trace is coming from groupk and the control could come from group k or any of the subpopulations finding the proportions of LRs that were greater than c
    RMEP_k_sub2[k]  <- mean(Results$LR2[Results$trace_id != Results$control_id & 
                                          Results$trace_group==k & 
                                          Results$control_group==k] > c) #Different source comparism where the trace came from group 1 and control came from group 1
    
    RMED_k_Sub1[k] <- mean(Results$LR1[Results$trace_id == Results$control_id&Results$trace_group==k] < c)
    #The trace is coming from groupk and the control comes from group k or finding the proportions of LRs that were less than c
    
    RMED_k_Sub2[k] <- mean(Results$LR2[Results$trace_id == Results$control_id&Results$trace_group==k] < c)
    #The trace is coming from groupk and the control comes from group k or finding the proportions of LRs that were less than c
    
  }
  RMEP_Results<-rbind(RMEP_Results, c(RMEP_total1, RMEP_k_ALL1, RMEP_k_sub1, 1))
  
  RMEP_Results<-rbind(RMEP_Results, c(RMEP_total2, RMEP_k_ALL2, RMEP_k_sub2, 2))
  RMED_Results<-rbind(RMED_Results, c(RMED_total1, RMED_k_Sub1, 1))
  RMED_Results<-rbind(RMED_Results, c(RMED_total2, RMED_k_Sub2, 2))
  
  colnames(RMEP_Results) = c('total', paste0(1:K, '_all'), paste0(1:K, '_subpopulation'), 'LR')
  
  colnames(RMED_Results) = c('total', paste0(1:K, '_subpopulation'), 'LR')
  
  return(list(RMEP_Results =RMEP_Results, RMED_Results = RMED_Results, B = B))
  
}

Analysis_results<- Realdatanalysis(my_data, 100)
saveRDS(Analysis_results, file = paste0('copper_sim_res_na_', N_a,'.rds'))


summarize <- function(Analysis_results) {
  RMEP <- Analysis_results$RMEP_Results
  RMED <- Analysis_results$RMED_Results
  
  RMEP_half1 <- RMEP[seq(from =1, to =nrow(RMEP), by =2), ]
  RMEP_half2 <- RMEP[seq(from =2, to=nrow(RMEP), by =2),  ]
  
  RMED_half1 <- RMED[seq(from =1, to =nrow(RMED), by =2), ]
  RMED_half2 <- RMED[seq(from =2, to=nrow(RMED), by =2),  ]
  
  mean_RMEP <- matrix(0, nrow = 2, ncol = ncol(RMEP))
  sd_RMEP <- matrix(0, nrow = 2, ncol = ncol(RMEP))
  
  mean_RMED <- matrix(0, nrow = 2, ncol = ncol(RMED))
  sd_RMED <- matrix(0, nrow = 2, ncol = ncol(RMED))
  
  for (i in 1:(ncol(RMEP)-1)) {
    mean_RMEP[, i] <- c(mean(RMEP_half1[, i]), mean(RMEP_half2[, i]))
    sd_RMEP[, i] <- c(sd(RMEP_half1[, i]), sd(RMEP_half2[, i]))
  }
  
  for (i in 1:(ncol(RMED))-1) {
    mean_RMED[, i] <- c(mean(RMED_half1[, i]), mean(RMED_half2[, i]))
    sd_RMED[, i] <- c(sd(RMED_half1[, i]), sd(RMED_half2[, i]))
  }
  
  mean_RMEP[, ncol(RMEP)] <- c(1,2)
  sd_RMEP[, ncol(RMEP)] <- c(1,2)
  mean_RMED[, ncol(RMED)] <- c(1,2)
  sd_RMED[, ncol(RMED)] <- c(1,2)
  
  colnames(mean_RMEP) <- colnames(Analysis_results$RMEP_Results)
  colnames(sd_RMEP) <- colnames(Analysis_results$RMEP_Results)
  
  colnames(mean_RMED) <- colnames(Analysis_results$RMED_Results)
  colnames(sd_RMED) <- colnames(Analysis_results$RMED_Results)
  
  return(list(
    mean_RMEP = mean_RMEP,
    sd_RMEP = sd_RMEP,
    mean_RMED = mean_RMED,
    sd_RMED = sd_RMED
  ))
}





summary_results <- summarize(Analysis_results)


#Ealy calculation of the relative change
rel_change_all1 <- (summary_results$mean_RMEP[1,3] -  summary_results$mean_RMEP[1,2])/summary_results$mean_RMEP[1,3]
rel_cahnge_all2 <- (summary_results$mean_RMEP[2,3] -  summary_results$mean_RMEP[2,2])/summary_results$mean_RMEP[2,3]

rel_change_sub1 <-  (summary_results$mean_RMEP[1,5] -  summary_results$mean_RMEP[1,4])/summary_results$mean_RMEP[1,5]
rel_change_sub2 <- (summary_results$mean_RMEP[2,5] -  summary_results$mean_RMEP[2,4])/summary_results$mean_RMEP[2,5]

RMEP_res <- rbind(c("change in LR 1 subpop" = rel_change_sub1, 
                    "change in LR 2 subpop" = rel_change_sub2, 
                    "change in LR 1 all" = rel_change_all1, 
                    "change in LR 2 all" = rel_cahnge_all2))



results <- readRDS('copper_sim_res_na_200.rds')
results_RMEP <- results$RMEP_Results


# Calculating the standard errors associated to the relative change in RMEP per iterations
# Safe relative change helper: (a - b) / a, with zero-protection
set.seed(124)
# Safe relative change helper with Gaussian noise
safe_rel <- function(a, b, eps = 1e-8) {
  #a_noisy <- a + rnorm(length(a), mean = 0, sd = noise_sd)
  #b_noisy <- b + rnorm(length(b), mean = 0, sd = noise_sd)
  ifelse(abs(a) < eps, NA_real_, (a - b) / a)
}

# Pull all per-iteration rows
R <- as.data.frame(results_RMEP)

# Compute relative changes per row (per iteration & LR)
R$rel_change_all   = safe_rel(R$`2_all`,           R$`1_all`)
R$rel_change_sub   = safe_rel(R$`2_subpopulation`, R$`1_subpopulation`)

# Split by LR and summarize
summ_one <- function(x) {
  c(Mean = mean(x, na.rm = TRUE),
    SD   = sd(x, na.rm = TRUE), # describes how much relative change varies per iteration (model variability)
    SE   = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))) # describes how precise my mean estimate is (statistical reliability)
}

RMEP_rel_summary <- rbind(
  `LR1_all`    = summ_one(R$rel_change_all[R$LR == 1]),
  `LR2_all`    = summ_one(R$rel_change_all[R$LR == 2]),
  `LR1_subpop` = summ_one(R$rel_change_sub[R$LR == 1]),
  `LR2_subpop` = summ_one(R$rel_change_sub[R$LR == 2])
)

print(RMEP_rel_summary)


colSums(!is.na(R[, c("rel_change_all", "rel_change_sub")]))  # count valid values
colSums(is.na(R[, c("rel_change_all", "rel_change_sub")]))   # count missing values




-------------------------------------------------------
### Heatmap generation
"C:/Users/igben/Desktop/RESEARCH PAPER/Results"
library(pheatmap)
type <- 'RMEPReal'
p <- 8

mean_mat <- summary_results$mean_RMEP[1, -8, drop=FALSE]
mean_mat2 <- summary_results$mean_RMEP[2, -8, drop=FALSE]
sd_mat <- summary_results$sd_RMEP[1, -8, drop=FALSE]
sd_mat2 <- summary_results$sd_RMEP[2, -8, drop=FALSE]

# Create row and column names

#unique_M<- c(0.1, 0.2, 0.3, 0.4, 0.5)
#unique_M<- c("A", "B", "C", "D", "E")
selected_columns<-c("T","1B","2B","3B","1W","2W", "3W")
#rownames(sd_mat) <- (unique_M)
colnames(sd_mat) <- selected_columns


#rownames(sd_mat2) <- (unique_M)
colnames(sd_mat2) <- selected_columns

#rownames(mean_mat) <- (unique_M)
colnames(mean_mat) <- selected_columns

#rownames(mean_mat2) <- (unique_M)
colnames(mean_mat2) <- selected_columns

breaks <- seq(min(c(mean_mat,mean_mat2, sd_mat, sd_mat2)), max(c(mean_mat,mean_mat2, sd_mat, sd_mat2)), length.out = 100)

x11()
pheatmap((mean_mat), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/mean', type, 'P', p, 'K_3LR1.pdf'), width = 10, height = 10)
dev.off()

x11()
pheatmap((sd_mat), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30, legend=FALSE, labels_row = c('','','','',''))
dev.copy2pdf(file = paste0('Images 7/sd', type, 'P', p, 'K_3LR1.pdf'), width = 5, height = 10)
dev.off()


x11()
pheatmap((mean_mat2), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/mean', type, 'P', p, 'K_3LR2.pdf'), width = 10, height = 10)
dev.off()

x11()
pheatmap((sd_mat2), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/sd', type, 'P', p, 'K_3LR2.pdf'), width = 5, height = 10)
dev.off()




#### RMED Heatmap
"C:/Users/igben/Desktop/RESEARCH PAPER/Results"

type <- 'RMEDReal'
p <- 8

mean_mat <- summary_results$mean_RMED[1, -5, drop=FALSE]
mean_mat2 <- summary_results$mean_RMED[2, -5, drop=FALSE]
sd_mat <- summary_results$sd_RMED[1, -5, drop=FALSE]
sd_mat2 <- summary_results$sd_RMED[2, -5, drop=FALSE]

# Create row and column names

#unique_M<- c(0.1, 0.2, 0.3, 0.4, 0.5)
#unique_M<- c("A", "B", "C", "D", "E")
selected_columns<-c("T","1W","2W", "3W")
#rownames(sd_mat) <- (unique_M)
colnames(sd_mat) <- selected_columns


#rownames(sd_mat2) <- (unique_M)
colnames(sd_mat2) <- selected_columns

#rownames(mean_mat) <- (unique_M)
colnames(mean_mat) <- selected_columns

#rownames(mean_mat2) <- (unique_M)
colnames(mean_mat2) <- selected_columns

breaks <- seq(min(c(mean_mat,mean_mat2, sd_mat, sd_mat2)), max(c(mean_mat,mean_mat2, sd_mat, sd_mat2)), length.out = 100)

x11()
pheatmap((mean_mat), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/mean', type, 'P', p, 'K_3LR1.pdf'), width = 10, height = 10)
dev.off()

x11()
pheatmap((sd_mat), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30, legend=FALSE, labels_row = c('','','','',''))
dev.copy2pdf(file = paste0('Images 7/sd', type, 'P', p, 'K_3LR1.pdf'), width = 5, height = 10)
dev.off()


x11()
pheatmap((mean_mat2), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/mean', type, 'P', p, 'K_3LR2.pdf'), width = 10, height = 10)
dev.off()

x11()
pheatmap((sd_mat2), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 7/sd', type, 'P', p, 'K_3LR2.pdf'), width = 5, height = 10)
dev.off()



