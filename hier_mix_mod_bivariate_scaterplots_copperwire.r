# code to visualize the heirarchical mixture models
#The bivariate scater plots
# Author: Isaac Gbene 
# Last updated: 1/12/2026

library(MixSim)
library(mvtnorm)

# Parameters

# order of trace elements in vectors
labels <- c('Ag', 'Sb', 'Pb', 'Bi', 'Co', 'Ni', 'As', 'Se')

# average concentration vectors for groups 1 and 2
mu_hat_1 <- c(7.429, 0.761, 0.774, 0.070, 0.030, 2.033, 0.917, 0.756)
mu_hat_2 <- c(4.524, 0.465, 0.687, 0.037, 0.28, 1.065, 0.512, 0.536)

mu <- rbind(mu_hat_1, mu_hat_2)

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

pi <- c(pi_hat_1, pi_hat_2)

# between covariance array
copper_S <- array(c(S_b_vals), dim=c(8,8,2))

copper_params <- list(mu = mu_hat_2, Sigma_a = S_b, Sigma_epsilon = S_w)


# Corrected sim_dat_fun2 function
sim_dat_fun2 <- function(n_a, n_w, pi, mu, cov, S_w) {
  library(MASS)
  
  k <- nrow(mu)  # number of groups
  p <- ncol(mu)  # number of variables
  
  # Initialize data frame with proper column names
  the_data <- data.frame()
  
  # For each group
  for(g in 1:k) {
    # Number of objects for this group (ensure at least 1)
    n_obj <- max(1, round(pi[g] * n_a))
    
    # For each object in this group
    for(obj in 1:n_obj) {
      # Generate object mean from between-source distribution
      obj_mean <- mvrnorm(1, mu = mu[g,], Sigma = cov[,,g])
      
      # Generate n_w measurements from within-source distribution
      measurements <- mvrnorm(n_w, mu = obj_mean, Sigma = S_w)
      
      # Create data frame for this object with proper column names
      obj_data <- data.frame(
        obj = rep(paste0("obj_", g, "_", obj), n_w),
        measurements,
        group = rep(g, n_w)
      )
      
      # Set column names for the trace elements
      colnames(obj_data)[2:(p+1)] <- paste0("V", 1:p)
      
      # Combine with main data
      the_data <- rbind(the_data, obj_data)
    }
  }
  
  return(the_data)
}

# point plot

plotting_points <- function(cov, mu, S_w, pi, n_w, n_a, labels, bivariate = FALSE){
  library(ellipse)
  if(bivariate == FALSE){
    
    the_data = sim_dat_fun2(n_a, n_w, pi, mu, cov, S_w)
    
    uni_obj = unique(the_data$obj)
    
    par(mar = c(0.5,0.5,0.5,0.5), oma = c(4,4,4,4))
    col_vec <- c('red','blue', 'darkgreen')
    pch_vec <- c(15,16,17)
    lty_vec <- c(1,2,3)
    
    k <- nrow(mu)
    p <- ncol(mu)
    
    par(mfrow = c(p,p))
    
    for(j in 1:p){
      for(i in 1:p){
        if(j == i){
          
          plot(x = 0, y = 0, type = 'n', yaxt = 'n', xaxt = 'n', ylab = ' ', xlab = ' ')
          text(x = 0, y = 0, labels[j], cex = 2)
          
        }else{
          
          if (j < i){
            
            tmp_ellipses <- list(NULL)
            tmp_within_ellipses <- list(NULL)
            for(g in 1:k){
              tmp_ellipses[[g]] <- ellipse(cov[c(i,j),c(i,j), g], centre = mu[g, c(i,j)], npoints = 1000)
              tmp_within_ellipses[[g]] <- ellipse(S_w[c(i,j),c(i,j)], centre = mu[g, c(i,j)])
            }
            
            min_y <- min(unlist(lapply(tmp_ellipses, function(x){min(x[,2])})),
                         unlist(lapply(tmp_within_ellipses, function(x){min(x[,2])})))
            max_y <- max(unlist(lapply(tmp_ellipses, function(x){max(x[,2])})),
                         unlist(lapply(tmp_within_ellipses, function(x){max(x[,2])})))
            
            min_x <- min(unlist(lapply(tmp_ellipses, function(x){min(x[,1])})),
                         unlist(lapply(tmp_within_ellipses, function(x){min(x[,1])})))
            max_x <- max(unlist(lapply(tmp_ellipses, function(x){max(x[,1])})),
                         unlist(lapply(tmp_within_ellipses, function(x){max(x[,1])})))
            
            plot(tmp_ellipses[[1]], type = 'n', col = 'red', lty = 2,
                 xlim = c(min_x, max_x), 
                 ylim = c(min_y, max_y),
                 yaxt = 'n',
                 xaxt = 'n',
                 ylab = ' ',
                 xlab = ' ')
            
            for(the_obj in uni_obj){
              # Get measurements for this object (exclude obj and group columns)
              this_obj_data <- the_data[the_data$obj == the_obj, ]
              this_obj_measurements <- as.matrix(this_obj_data[, 2:(ncol(this_obj_data)-1)])
              
              # Calculate object means
              this_obj_means <- colMeans(this_obj_measurements)
              
              # Get group (last column)
              this_obj_group <- unique(this_obj_data$group)
              
              # Plot segments and points
              for(samps in 1:nrow(this_obj_measurements)){
                segments(this_obj_means[i], this_obj_means[j],
                         this_obj_measurements[samps,i], this_obj_measurements[samps,j],
                         col = 'darkgrey', lwd = 1.5)  # Reduced lwd for better visibility
              }
              points(this_obj_means[i], this_obj_means[j], 
                     pch = pch_vec[this_obj_group], 
                     cex = 0.5,  # Increased size
                     col = col_vec[this_obj_group])
            }
            for(g in 1:k){
              
              lines(tmp_ellipses[[g]], col = col_vec[g], lty = lty_vec[g], lwd = 1.5)
              
              # points(mu[g, i], mu[g, j], col = col_vec[g], pch = pch_vec[g], cex = 0.5)
              
              
              
              if(j == 1){
                axis(3, cex.axis = 1.5, cex.lab = 1.5, las = 2)
              }
              
              if(i == p){
                axis(4, cex.axis = 1.5, cex.lab = 1.5, las = 2)
              }
              
            } 
          }
          
          if (j > i){
            
            plot(0,0, type = 'n', bty = 'n', axes = FALSE)
            
          }
        }
      }
    }
  }else{
    
    ########################################################################
    
    
    the_data = sim_dat_fun2(n_a, n_w, pi, mu, cov, S_w)
    
    uni_obj = unique(the_data$obj)
    
    par(mar = c(3,3,3,3))
    col_vec <- c('red','blue', 'darkgreen', 'black', 'magenta')
    pch_vec <- c(15,16,17, 18, 19)
    lty_vec <- c(1,2,3,4,5)
    
    k <- nrow(mu)
    p <- ncol(mu)
    
    # used to be in loop
    
    tmp_ellipses <- list(NULL)
    tmp_within_ellipses <- list(NULL)
    for(g in 1:k){
      tmp_ellipses[[g]] <- ellipse(cov[,,g], centre = mu[g,])
      tmp_within_ellipses[[g]] <- ellipse(S_w, centre = mu[g,])
    }
    
    min_y <- min(unlist(lapply(tmp_ellipses, function(x){min(x[,2])})),
                 unlist(lapply(tmp_within_ellipses, function(x){min(x[,2])})))
    max_y <- max(unlist(lapply(tmp_ellipses, function(x){max(x[,2])})),
                 unlist(lapply(tmp_within_ellipses, function(x){max(x[,2])})))
    
    min_x <- min(unlist(lapply(tmp_ellipses, function(x){min(x[,1])})),
                 unlist(lapply(tmp_within_ellipses, function(x){min(x[,1])})))
    max_x <- max(unlist(lapply(tmp_ellipses, function(x){max(x[,1])})),
                 unlist(lapply(tmp_within_ellipses, function(x){max(x[,1])})))
    
    
    plot(tmp_ellipses[[1]], type = 'n', col = 'red', lty = 2,
         xlim = c(min_x, max_x), 
         ylim = c(min_y, max_y),
         ylab = ' ',
         xlab = ' ',
         cex.lab = 1.4,
         cex.axis = 1.4)
    
    for(the_obj in uni_obj){
      ncols = ncol(the_data)
      this_obj_measurements = as.matrix(the_data[the_data$obj == the_obj, -c(1,ncols)])
      
      this_obj_means = colMeans(this_obj_measurements)
      this_obj_group = unique(the_data[the_data$obj == the_obj, ncols])
      
      
      for(samps in 1:n_w){
        segments(this_obj_means[1], this_obj_means[2],
                 this_obj_measurements[samps,1], this_obj_measurements[samps,2],
                 col = 'grey', lwd = 2)
        # points(this_obj_measurements[samps,i], this_obj_measurements[samps,j], 
        #        pch = 16, cex = 0.5,
        #        col = col_vec[this_obj_group])
      }
      points(this_obj_means[1], this_obj_means[2], pch = pch_vec[this_obj_group], 
             cex = 0.5, col = col_vec[this_obj_group])
    }
    
    for(g in 1:k){
      
      lines(tmp_ellipses[[g]], col = col_vec[g], lty = lty_vec[g], lwd = 2)
      
      # points(mu[g, i], mu[g, j], col = col_vec[g], pch = pch_vec[g], cex = 0.5)
    } 
    
    
    
    
    ####################################################################
  }
}

plotting_points(cov = copper_S, # plotting the bivariate scater plots
                mu = mu, 
                S_w = S_w, 
                pi = pi,
                n_w = 5, 
                n_a = 15, 
                labels = labels)
par(new = TRUE, mfrow = c(1,1))
plot(0,0, type = 'n', bty = 'n', axes = FALSE)
legend('bottomleft', legend = c('Subpop. 1', 'Subpop. 2'), 
       pch = c(15,16), lty = c(1,2), col = c('red', 'blue'))





