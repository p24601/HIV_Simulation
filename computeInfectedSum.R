setwd("~/Desktop/HIV/CIed_results/Combo 2")
library(plyr)

Rvalues <- c(150, 175, 187, 191, 192, 193,
             194, 195, 200, 210, 220, 230, 260,
             270, 280, 290, 300, 310, 320, 330) 

values <- seq(1,1500, by=1)
timestep <- 200
resultsData <- data.frame()
missed = 0
base_dest_states = "/states_count_raw.csv"

sims = list("results69_06_03_01", "results54_06_03_01")


next_sim = sims[[2]]

for (r in Rvalues){
  for (i in values){
    f = paste("./", next_sim, "/", r, "/", i, "", base_dest_states, sep = "")
    if (file.exists(f)){
      temp<- read.csv(file=f, header = TRUE, sep = ",")
      row <- temp[timestep,]
      df <- data.frame(r, i, row$Infected, row$Neutralized)
      
      resultsData <- rbind(resultsData, df)
    }else{
      missed = missed +1      
    }
  }
}

resultsData[5] = resultsData[3] + resultsData[4]
names(resultsData)[5] = "Infected_sum"
resultsData <- ddply(resultsData, 'r', summarise, I_mean=mean(Infected_sum),I_sd=sd(Infected_sum))
names(resultsData)[2] = "Infected_sum_mean"
names(resultsData)[3] = "Infected_sum_sd"

nf = paste("./", next_sim, "/computed_mean_sd6.csv", sep = "")
new_df = read.csv(nf, sep = ",")
new_df[9] = resultsData[2]
new_df[10] = resultsData[3]
new_df = new_df[,-1]
write.csv(new_df, nf)

resultsData <- data.frame()
missed = 0
