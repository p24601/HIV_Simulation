setwd("/home/pcervell/hiv_simulation/results_1")

Rvalues <- c(100, 150, 175, 187, 191, 192, 192.5, 193, 193.5,
             194, 194.5, 195, 195.5, 200, 210, 220, 230, 260,
             270, 280, 290, 300, 310, 320, 330) 
values <- seq(1,1500, by=1)


base_dest_states = "/states_count_raw.csv"
base_dest_resistance = "/resistance_count_raw.csv"

timestep <- 200
resultsData <- data.frame()
missingData = as.list(rep(0, 25))
names(missingData) <- c('100', '150', '175', '187', '191', '192', '192.5', '193', '193.5',
                        '194', '194.5', '195', '195.5', '200', '210', '220', '230', '260',
                        '270', '280', '290', '300', '310', '320', '330')

for (r in Rvalues){
  for (i in values){
    f = paste(r, "/", i, base_dest_states, sep = "")
    if (file.exists(f)){
      temp<- read.csv(file=f, header = TRUE, sep = ",")
      row <- temp[timestep,]
      df <- data.frame(r, i, row$Infected, row$Neutralized)
      
      temp<- read.csv(file=paste(r, "/", i, base_dest_resistance, sep = ""), header = TRUE, sep = ",")
      row <- temp[timestep,]
      df$R1 <- row$X1
      df$R2 <- row$X2
      df$R3 <- row$X3
      df$Overall <- row$Overall
      
      resultsData <- rbind(resultsData, df)
    }else{
      missingData[as.character(r)] = as.numeric(missingData[as.character(r)]) +1
    }
  }
}

library(plyr)
computed_mean_sd <- ddply(resultsData, 'r', summarise, R1_mean=mean(R1),R1_sd=sd(R1))
computed_mean_sd <- cbind(computed_mean_sd, ddply(resultsData, 'r', summarise, R2_mean=mean(R2),R2_sd=sd(R2)))
computed_mean_sd <- cbind(computed_mean_sd, ddply(resultsData, 'r', summarise, R3_mean=mean(R3),R3_sd=sd(R3)))
computed_mean_sd <- computed_mean_sd[,-4]
computed_mean_sd <- computed_mean_sd[,-6]

write.csv(computed_mean_sd, file = "computed_mean_sd6.csv")
write.csv(resultsData, file = "resultsData.csv")

library(erer)
missingData.df = as.data.frame(missingData)
write.csv(missingData.df, file = "missingData6.csv")
