setwd("~/Desktop/HIV/results_broad_sweep")

Rvalues <- c(100, 150, 175, 187, 191, 192, 192.5, 193, 193.5,
             194, 194.5, 195, 195.5, 200, 210, 220, 230, 260,
             270, 280, 290, 300, 310, 320, 330) 
values <- seq(1,1500, by=1)

base_dest_states = "/states_count_raw.csv"
base_dest_resistance = "/resistance_count_raw.csv"
base_dest_epochs = "/infectedEpochs_count_raw.csv"
base_dest_genotypes = "/genotypes_count_raw.csv"

missed = 0
timestep <- 200

for (r in Rvalues){
  for (i in values){
    f = paste(r, "/", i, base_dest_states, sep = "")
    if (file.exists(f)){
      f_res = paste(r, "/", i, base_dest_resistance, sep = "")
      f_epo = paste(r, "/", i, base_dest_epochs, sep = "")
      f_gen = paste(r, "/", i, base_dest_genotypes, sep = "")
      
      res = read.csv(file=f_res, header = TRUE, sep = ",")
      epo = read.csv(file=f_epo, header = TRUE, sep = ",")
      gen = read.csv(file=f_gen, header = TRUE, sep = ",")
      
      resultsData<- read.csv(file=f, header = TRUE, sep = ",")
      resultsData[6]<- res[2]
      resultsData[7]<- res[3]
      resultsData[8]<- res[4]
      resultsData[9] <- epo[2]
      resultsData[10]<- epo[3]
      resultsData[11]<- epo[4]
      resultsData[12]<- gen[2]
      resultsData = resultsData[1:200,]
      
      names(resultsData) = c("Timestep", "Healthy", "Infected", "Dead", "Neutralized", 
                            "Resistant to 1", "Resistant to 2", "Resistant to 3", 
                            "Alive for 1", "Alive for 2", "Alive for 3", "Number of genotypes")
      
      f = paste(r, "/r_70_69.0_0.6_0.3_0.1_", i, sep = "")
      write.csv(resultsData, f)
    }else{
      missed = missed +1
    }
  }
}