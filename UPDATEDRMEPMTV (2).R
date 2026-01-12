library(MixSim)
library(comparison)
library(mvtnorm)
#Last Updated: 1/12/2026
#Function for calculating RMEP/RMED multivariate version

generate_observations <- function(source_mixture, n_i, p) {
  
  if(p == 1){
    observations <- NULL
    for (i in 1:length(source_mixture$id)) {
      source_means <- source_mixture$X[i,]
      temp_observations <- rnorm(n_i, mean = source_means, sd = sqrt(Sw))
      observations <- rbind(observations, cbind(temp_observations, rep(i, n_i)))
    }
    #hist(observations)
  }else{
    observations <- NULL
    for (i in 1:length(source_mixture$id)) {
      source_means <- source_mixture$X[i,]
      temp_observations <- rmvnorm(n_i, mean = source_means, sigma = Sw)
      observations <- rbind(observations, cbind(temp_observations, rep(i, n_i)))
    }
  }
  return(observations)
}
Simulate_compare_general <- function (pi, Mu, S,Sw, m, n_i, B){
  
  p <- ncol(Mu)
  
  RMEP_Results <- NULL
  RMED_Results <- NULL
  # Loop over replicates of B
  for (b in 1:B) {   
    
    source_mixture <- simdataset(m, Pi = pi, Mu = Mu, S=S)
    observations <- generate_observations(source_mixture, n_i, p)  
    
    
    source_ids <- source_mixture$id # source IDs
    
    
    K=length(pi)
    
    Group_sources <-list(K)
    train_indices <- NULL
    
    
    for ( k in 1:K){ 
      
      groupk_sources <-which(source_ids==k) #  in the little kth position store that, all the sources that came from group k 
      
      train_indices<-c(train_indices, sample(groupk_sources,length(groupk_sources)/2))
      
      Group_sources[[k]] <-groupk_sources
    }
    
    
    
    #train_indices=c(sample(groupk_sources,length(groupk_sources)/2))
    
    test_indices <- setdiff(1:m, train_indices)  # Remaining sources for testing
    train_data<-observations[observations[,(p+1)]%in%train_indices,]
    
    test_data<-observations[observations[,(p+1)]%in%test_indices,]
    
    
    
    ## Step 5: Label observations in the test set
    trace_indices <- 1:(n_i/2)
    control_indices <- (n_i/2 + 1):n_i
    
    
    Results <- NULL
    
    
    for (i in 1:(length(test_indices) )) {
      for (j in  (i): length(test_indices)){
        #extract data from the current pair of indices
        test_i <- test_data[test_data[, (p+1)] == test_indices[i], ]
        
        test_j <- test_data[test_data[, (p+1)] == test_indices[j], ]
        
        control_i <-test_i[control_indices, ]
        
        control_j <- test_j[control_indices, ]
        
        trace_i <- test_i[trace_indices, ]
        
        trace_j <- test_j[trace_indices,]
        
        source_i <- test_indices[i]
        source_j <- test_indices[j]
        
        
        trace_j_format <- two.level.comparison.items(trace_j, 1:p)
        
        trace_i_format <- two.level.comparison.items(trace_i, 1:p)
        
        control_j_format <- two.level.comparison.items(control_j, 1:p)
        
        control_i_format <- two.level.comparison.items(control_i, 1:p)
        
        e_a_format <- two.level.components(as.data.frame(train_data), data.columns = 1:p, item.column= (p+1))
        
        if (p==1){
          lind_LR_ij <- two.level.lindley.LR(trace_i_format, control_j_format, e_a_format )
        
          lind_LR_ji <- two.level.lindley.LR(trace_j_format, control_i_format, e_a_format )
        
          
        }else{
          lind_LR_ij <- two.level.normal.LR(trace_i_format, control_j_format, e_a_format )
          
          lind_LR_ji <- two.level.normal.LR(trace_j_format, control_i_format, e_a_format )
          
        }
        
      
        
        density_LR_ij <- two.level.density.LR(trace_i_format, control_j_format, e_a_format)
        
        density_LR_ji <- two.level.density.LR(trace_j_format, control_i_format, e_a_format)
        
        
        
        for ( g in 1:K){ 
          if(source_i%in%Group_sources[[g]]){
            i_group <- g# looping through one to 3, see if the ith object is in group 1, if not jumpt to group if not check group 3
          }
          if(source_j%in%Group_sources[[g]]){
            j_group <- g# looping through one to 3, see if the j is the the group.
          }
        }
        #i_group<- ifelse(source_i %in% groupk_sources, 1 ,2) ## if the ith object is in group1_source then result 1 if not 2
        #j_group <- ifelse(source_j %in% groupk_sources, 1 ,2) ## if the jth object is in group1_source then result 1 if not 2
        Results <-rbind(Results, c(trace_id = source_i, control_id = source_j,trace_group = i_group , control_group = j_group, LR1= lind_LR_ij, LR2= density_LR_ij)  )
        Results <-rbind(Results, c(trace_id = source_j, control_id = source_i,trace_group = j_group , control_group = i_group, LR1= lind_LR_ji, LR2= density_LR_ji)  ) # swapping i and j to match 
        
        
        
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
  
  
  #RMEP[[b]] <- RMEP_Results
  #RMED[[b]] <- RMED_Results
  
  return(list(RMEP_Results = RMEP_Results, RMED_Results = RMED_Results, B = B))
}




RME <- function(Results, c = 1, RMEP_Results, RMED_Results, B){
  #RMEP-Results <- NULL
  #RMED_Results < - NULL
  
  
  
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

summarize <- function(sim_results) {
  RMEP <- sim_results$RMEP_Results
  RMED <- sim_results$RMED_Results
  
  RMEP_half1 <- RMEP[seq(from =1, to =nrow(RMEP), by =2), ]
  RMEP_half2 <- RMEP[seq(from =2, to=nrow(RMEP), by =2),  ]
  
  RMED_half1 <- RMED[seq(from =1, to =nrow(RMED), by =2), ]
  RMED_half2 <- RMED[seq(from =2, to=nrow(RMED), by =2),  ]
  
  mean_RMEP <- matrix(0, nrow = 2, ncol = ncol(RMEP))
  sd_RMEP <- matrix(0, nrow = 2, ncol = ncol(RMEP))
  
  mean_RMED <- matrix(0, nrow = 2, ncol = ncol(RMED))
  sd_RMED <- matrix(0, nrow = 2, ncol = ncol(RMED))
  
  for (i in 1:ncol(RMEP)) {
    mean_RMEP[, i] <- c(mean(RMEP_half1[, i]), mean(RMEP_half2[, i]))
    sd_RMEP[, i] <- c(sd(RMEP_half1[, i]), sd(RMEP_half2[, i]))
  }
  
  for (i in 1:ncol(RMED)) {
    mean_RMED[, i] <- c(mean(RMED_half1[, i]), mean(RMED_half2[, i]))
    sd_RMED[, i] <- c(sd(RMED_half1[, i]), sd(RMED_half2[, i]))
  }
  
  colnames(mean_RMEP) <- colnames(sim_results$RMEP_Results)
  colnames(sd_RMEP) <- colnames(sim_results$RMEP_Results)
  
  colnames(mean_RMED) <- colnames(sim_results$RMED_Results)
  colnames(sd_RMED) <- colnames(sim_results$RMED_Results)
  
  return(list(
    mean_RMEP = mean_RMEP,
    sd_RMEP = sd_RMEP,
    mean_RMED = mean_RMED,
    sd_RMED = sd_RMED
  ))
}





















