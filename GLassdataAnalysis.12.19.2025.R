



### Real data exploratory analysis
library(dplyr)
library(MixSim)
library(comparison)
library(mvtnorm)
ibrary(ggplot2)
# Last Updated: 10/29/2025

glass_data <- read.csv('glass.txt', sep = ' ')
glass_data <- glass_data[, c(6:9,5)]


### Check for missing values

# Check for missing values
sum(is.na(glass_data))

# Plot to identify outliers
boxplot(glass_data$logCaK, main="Boxplot for logCaK")
boxplot(glass_data$logCaSi, main="Boxplot for logCaSi")
boxplot(glass_data$logCaFe, main="Boxplot for logCaFe")




# Calculate the frequency of each group
group_counts <- table(glass_data$Group)

# Calculate the proportions of each group
group_proportions <- prop.table(group_counts)

# Print the proportions
group_proportions


### Visualizing data distribution


# Histograms of chemical properties
truehist(glass_data$logCaK, main="Histogram of logCaK", ylab= "density", xlab="logCaK", col="blue")
truehist(glass_data$logCaSi, main="Histogram of logCaSi", ylab= "density", xlab="logCaSi", col="red")
truehist(glass_data$logCaFe, main="Histogram of logCaFe",ylab= "density", xlab="logCaFe", col="green")

# Using ggplot2 for better visualization
library(ggplot2)
ggplot(glass_data, aes(x=logCaK)) + geom_histogram(bins=20, fill="blue", color="black") + ggtitle("Distribution of logCaK")
ggplot(glass_data, aes(x=logCaSi)) + geom_histogram(bins=20, fill="red", color="black") + ggtitle("Distribution of logCaSi")
ggplot(glass_data, aes(x=logCaFe)) + geom_histogram(bins=20, fill="green", color="black") + ggtitle("Distribution of logCaFe")

## Advanced visualization 

# Scatter plots to examine relationships between variables
ggplot(glass_data, aes(x=logCaSi, y=logCaK, color=factor(Group))) + geom_point() + ggtitle("logCaSi vs logCaK by Group")
ggplot(glass_data, aes(x=logCaSi, y=logCaFe, color=factor(Group))) + geom_point() + ggtitle("logCaSi vs logCaFe by Group")


### Analysing Relationship that exist between variables
# Scatter plots to observe relationships between chemical properties
pairs(glass_data[,1:3], pch= 16,  main="Scatterplot Matrix of Chemical Properties",
      col=rainbow(3,alpha=0.5)[glass_data$Group], lower.panel = NULL)
par(new=TRUE, mar = c(1.5,1.5,0.5,0.5))
plot(0,0,type='n', axes = FALSE, xlab= '', ylab='')
legend('bottomleft', legend=sort(unique(glass_data$Group)), pch = 18, col=rainbow(3)[sort(unique(glass_data$Group))])
par(mar = c(5.1, 4.1, 4.1, 2.1))


# Correlation matrix
cor(glass_data[,1:3])

# Boxplots to see distributions of chemical properties by Group and Window
ggplot(glass_data, aes(x=factor(Group), y=logCaK)) + geom_boxplot() + ggtitle("logCaK by Group")
ggplot(glass_data, aes(x=factor(Window), y=logCaK)) + geom_boxplot() + ggtitle("logCaK by Window")

# Using facet_wrap to compare across groups
ggplot(glass_data, aes(x=logCaK, y=logCaFe)) + geom_point() + facet_wrap(~Group) + ggtitle("logCaK vs logCaFe by Group")








## Real Data Algorithm and heatmap plots 

glass_data <- read.csv('glass.txt', sep = ' ')

glass_data <- glass_data[, c(6:9,5)]
summarise(glass_data)



# Ensure 'Window' and 'Group' columns are correctly identified
uni_sources <- unique(glass_data$Window)
source_ids <- integer(length(uni_sources))
names(source_ids) <- uni_sources

for (i in seq_along(uni_sources)) {
  # Find the first index where each unique window appears and take the corresponding group
  index <- which(glass_data$Window == uni_sources[i])[1]
  source_ids[i] <- glass_data$Group[index]
}


Realdatanalysis <- function(glass_data, B) {
  
  p <- 3  # three properties: logCaK, logCaSi, logCaFe
  RMEP_Results <- NULL
  RMED_Results <- NULL
  
  
  
  for (b in 1:B) {
    K = 3# 
    Group_sources <-list(K)
    train_indices <- NULL
    
    
    for ( k in 1:K){ 
      
      groupk_sources <-which(source_ids==k) #  in the little kth position store that, all the sources that came from group k 
      
      train_indices<-c(train_indices, sample(groupk_sources,length(groupk_sources)/2))
      
      Group_sources[[k]] <-groupk_sources
    }
    #train_indices=c(sample(groupk_sources,length(groupk_sources)/2))
    
    test_indices <- setdiff(1:62, train_indices)  # Remaining sources for testing
    train_data<-glass_data[glass_data[,(p+1)]%in%train_indices,]
    
    test_data<-glass_data[glass_data[,(p+1)]%in%test_indices,]
    
    
    
    ## Step 5: Label observations in the test set
    trace_indices <- 1:2       # Five observations from each source so we choose one 1-2 for the trace indces 
    control_indices <- 3:5
    
    Results <- NULL # Initialize empty data frame for results
    # Obtain the e_a_format from training data
    e_a_format <- two.level.components(train_data, data.columns = 1:p, item.column = p+1)
    
    # Running comparisons
    test_indices <- unique(test_data$Window)
    for (i in seq_along(test_indices)) {
      for (j in seq(i, length(test_indices))) {
        # Extract data for the current pair of indices
        test_i <- test_data[test_data$Window == test_indices[i], ]
        test_j <- test_data[test_data$Window == test_indices[j], ]
        
        # control_indices' and 'trace_indices' are predefined
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
        i_group <- glass_data$Group[glass_data$Window == test_indices[i]][1]
        j_group <- glass_data$Group[glass_data$Window == test_indices[j]][1]
        
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
Analysis_results<- Realdatanalysis(glass_data, 100)



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
  
  for (i in 1:ncol(RMEP)) {
    mean_RMEP[, i] <- c(mean(RMEP_half1[, i]), mean(RMEP_half2[, i]))
    sd_RMEP[, i] <- c(sd(RMEP_half1[, i]), sd(RMEP_half2[, i]))
  }
  
  for (i in 1:ncol(RMED)) {
    mean_RMED[, i] <- c(mean(RMED_half1[, i]), mean(RMED_half2[, i]))
    sd_RMED[, i] <- c(sd(RMED_half1[, i]), sd(RMED_half2[, i]))
  }
  
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


### Heatmap generation
"C:/Users/igben/Desktop/RESEARCH PAPER/Results"
library(pheatmap)
type <- 'RMEPReal'
p <- 3

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
dev.copy2pdf(file = paste0('Images 8/mean', type, 'P', p, 'K_3LR1.pdf'), width = 6, height = 6)
dev.off()

x11()
pheatmap((sd_mat), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30, legend=FALSE, labels_row = c('','','','',''))
dev.copy2pdf(file = paste0('Images 8/sd', type, 'P', p, 'K_3LR1.pdf'), width = 6, height = 6)
dev.off()


x11()
pheatmap((mean_mat2), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 8/mean', type, 'P', p, 'K_3LR2.pdf'), width = 6, height = 6)
dev.off()

x11()
pheatmap((sd_mat2), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 8/sd', type, 'P', p, 'K_3LR2.pdf'), width = 6, height = 6)
dev.off()




#### RMED Heatmap
"C:/Users/igben/Desktop/RESEARCH PAPER/Results"

type <- 'RMEDReal'
p <- 3

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
dev.copy2pdf(file = paste0('Images 8/mean', type, 'P', p, 'K_3LR1.pdf'), width = 6, height = 6)
dev.off()

x11()
pheatmap((sd_mat), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30, legend=FALSE, labels_row = c('','','','',''))
dev.copy2pdf(file = paste0('Images 8/sd', type, 'P', p, 'K_3LR1.pdf'), width = 6, height = 6)
dev.off()


x11()
pheatmap((mean_mat2), cluster_rows = FALSE, cluster_cols = FALSE,xlab =  'RMED',
         ylab = 'Mixture proportion', breaks = breaks, angle_col = "0", legend=FALSE, labels_row = c('','','','',''),
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 8/mean', type, 'P', p, 'K_3LR2.pdf'), width = 6, height = 6)
dev.off()

x11()
pheatmap((sd_mat2), cluster_rows = FALSE, cluster_cols = FALSE, xlab =  'RMED',ylab = 'Mixture proportion',breaks = breaks, angle_col = "0",
         fontsize = 30)
dev.copy2pdf(file = paste0('Images 8/sd', type, 'P', p, 'K_3LR2.pdf'), width = 6, height = 6)
dev.off()



