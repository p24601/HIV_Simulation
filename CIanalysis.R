# t_n-1_alpha/2
# At 95% CI, and n = 10, n-1 = 9, alpha/2 = .25
computed_mean_sd = read.csv("~/Desktop/HIV/results_broad_sweep/Combo_1/results69_06_03_01/computed_mean_sd6.csv", sep = ",")        
computed_mean_sd = read.csv("~/Desktop/HIV/results_broad_sweep/Combo_1/results54_06_03_01/computed_mean_sd6.csv", sep = ",")        
computed_mean_sd = read.csv("~/Desktop/HIV/results_broad_sweep/Combo_2/results69_06_03_01/computed_mean_sd6.csv", sep = ",")        
computed_mean_sd = read.csv("~/Desktop/HIV/results_broad_sweep/Combo_2/results54_06_03_01/computed_mean_sd6.csv", sep = ",")        
t <- 2.262

# Error margin in % about the mean (10 cells)
error_margin <- seq(5, 20, by = 5)

sample_size_results <- data.frame()
for (e in error_margin){
  
  sample_size_R1 <- c((100*computed_mean_sd$R1_sd*t)/(e*computed_mean_sd$R1_mean))^2
  sample_size_R3 <- c((100*computed_mean_sd$R3_sd*t)/(e*computed_mean_sd$R3_mean))^2
  sample_size_R2 <- c((100*computed_mean_sd$R2_sd*t)/(e*computed_mean_sd$R2_mean))^2
  df <- data.frame(e, computed_mean_sd$r, sample_size_R1, sample_size_R2, sample_size_R3)
 
  sample_size_results <- rbind(sample_size_results, df) 
}
colnames(sample_size_results)[2]<-'r'

