# Set working directory to results folder
setwd("~/Desktop/HIV/results_broad_sweep/Combo 3")
setwd("~/Desktop/HIV/results_broad_sweep/Combo 1")

Rvalues <- c(150, 175, 187, 191, 192, 193,
             194, 195, 200, 210, 220, 230, 260,
             270, 280, 290, 300, 310, 320, 330)

####################### Broad Sweep Plot mean and Sd #######################
r_100_99.0_0.6_0.3_0.1      = read.csv("./results99_06_03_01/computed_mean_sd6.csv", sep = ",")
r_100_70.0_25.0_4.0_1.0     = read.csv("./results70_25_4_1/computed_mean_sd6.csv", sep = ",")
r_100_70.0_15.0_13.0_2.0    = read.csv("./results70_15_13_2/computed_mean_sd6.csv", sep = ",")
r_100_50.0_25.0_12.5_12.5   = read.csv("./results50_25_125_125/computed_mean_sd6.csv", sep = ",")      

r_85_84.0_0.6_0.3_0.1       = read.csv("./results84_06_03_01/computed_mean_sd6.csv", sep = ",")
r_85_59.5_21.3_3.5_0.8      = read.csv("./results595_2125_35_075/computed_mean_sd6.csv", sep = ",")    
r_85_59.5_12.8_11.1_1.6     = read.csv("./results595_1275_1105_16/computed_mean_sd6.csv", sep = ",")   
r_85_42.5_21.3_10.6_10.6    = read.csv("./results425_2125_1063_1063/computed_mean_sd6.csv", sep = ",") 

r_70_69.0_0.6_0.3_0.1       = read.csv("./results69_06_03_01/computed_mean_sd6.csv", sep = ",")        
r_70_49.0_17.1_2.5_0.7      = read.csv("./results49_171_25_07/computed_mean_sd6.csv", sep = ",")       
r_70_49.0_10.5_9.10_1.4     = read.csv("./results49_105_910_14/computed_mean_sd6.csv", sep = ",")      
r_70_35.0_17.5_8.8_8.8      = read.csv("./results35_175_875_875/computed_mean_sd6.csv", sep = ",")

r_55_54.0_0.6_0.3_0.1       = read.csv("./results54_06_03_01/computed_mean_sd6.csv", sep = ",")        
r_55_38.4_29.8_2.2_0.6      = read.csv("./results384_2975_22_055/computed_mean_sd6.csv", sep = ",")    
r_55_38.4_8.3_7.2_1.1       = read.csv("./results384_825_715_11/computed_mean_sd6.csv", sep = ",")     
r_55_27.5_13.8_6.9_6.9      = read.csv("./results275_1375_688_688/computed_mean_sd6.csv", sep = ",")

r_85_85.0_0_0_0             = read.csv("./results85_0_0_0/computed_mean_sd6.csv", sep = ",")       
r_85_59_21.3_3.5_0.9        = read.csv("./results59_213_35_9/computed_mean_sd6.csv", sep = ",")    
r_70_70.0_0_0_0             = read.csv("./results70_0_0_0/computed_mean_sd6.csv", sep = ",")      
r_55_55.0_0_0_0             = read.csv("./results55_0_0_0/computed_mean_sd6.csv", sep = ",")          
r_55_38.4_13.8_2.2_0.6      = read.csv("./results384_138_22_6/computed_mean_sd6.csv", sep = ",")

sims3 = list(r_100_99.0_0.6_0.3_0.1, r_100_70.0_25.0_4.0_1.0, r_100_70.0_15.0_13.0_2.0, r_100_50.0_25.0_12.5_12.5,  
             r_85_85.0_0_0_0,        r_85_59_21.3_3.5_0.9,    r_85_59.5_12.8_11.1_1.6,  r_85_42.5_21.3_10.6_10.6,   
             r_70_70.0_0_0_0,        r_70_49.0_17.1_2.5_0.7,  r_70_49.0_10.5_9.10_1.4,  r_70_35.0_17.5_8.8_8.8,     
             r_55_55.0_0_0_0,        r_55_38.4_13.8_2.2_0.6,  r_55_38.4_8.3_7.2_1.1,    r_55_27.5_13.8_6.9_6.9)

sims12 = list(r_100_99.0_0.6_0.3_0.1, r_100_70.0_25.0_4.0_1.0, r_100_70.0_15.0_13.0_2.0, r_100_50.0_25.0_12.5_12.5,   
              r_85_84.0_0.6_0.3_0.1,  r_85_59.5_21.3_3.5_0.8,  r_85_59.5_12.8_11.1_1.6,  r_85_42.5_21.3_10.6_10.6,    
              r_70_69.0_0.6_0.3_0.1,  r_70_49.0_17.1_2.5_0.7,  r_70_49.0_10.5_9.10_1.4,  r_70_35.0_17.5_8.8_8.8,      
              r_55_54.0_0.6_0.3_0.1,  r_55_38.4_29.8_2.2_0.6,  r_55_38.4_8.3_7.2_1.1,    r_55_27.5_13.8_6.9_6.9)      

sims3_names = list("r_100_99.0_0.6_0.3_0.1", "r_100_70.0_25.0_4.0_1.0", "r_100_70.0_15.0_13.0_2.0", "r_100_50.0_25.0_12.5_12.5",  
                   "r_85_85.0_0_0_0",        "r_85_59_21.3_3.5_0.9",    "r_85_59.5_12.8_11.1_1.6",  "r_85_42.5_21.3_10.6_10.6",   
                   "r_70_70.0_0_0_0",        "r_70_49.0_17.1_2.5_0.7",  "r_70_49.0_10.5_9.10_1.4",  "r_70_35.0_17.5_8.8_8.8",     
                   "r_55_55.0_0_0_0",        "r_55_38.4_13.8_2.2_0.6",  "r_55_38.4_8.3_7.2_1.1",    "r_55_27.5_13.8_6.9_6.9")

sims12_names = list("r_100_99.0_0.6_0.3_0.1", "r_100_70.0_25.0_4.0_1.0", "r_100_70.0_15.0_13.0_2.0", "r_100_50.0_25.0_12.5_12.5",   
                    "r_85_84.0_0.6_0.3_0.1",  "r_85_59.5_21.3_3.5_0.8",  "r_85_59.5_12.8_11.1_1.6",  "r_85_42.5_21.3_10.6_10.6",    
                    "r_70_69.0_0.6_0.3_0.1",  "r_70_49.0_17.1_2.5_0.7",  "r_70_49.0_10.5_9.10_1.4",  "r_70_35.0_17.5_8.8_8.8",      
                    "r_55_54.0_0.6_0.3_0.1",  "r_55_38.4_29.8_2.2_0.6",  "r_55_38.4_8.3_7.2_1.1",    "r_55_27.5_13.8_6.9_6.9")

####################### CIed Plot mean and Sd #######################
library(reshape)
library(ggplot2)

# Combo 1
setwd("~/Desktop/HIV/CIed_results/Combo 1")
setwd("~/Desktop/HIV/CIed_results/Combo 2")
setwd("~/Desktop/HIV/CIed_results/Combo 3")
r_70_69.0_0.6_0.3_0.1       = read.csv("./results69_06_03_01/computed_mean_sd6.csv", sep = ",")        
r_55_54.0_0.6_0.3_0.1       = read.csv("./results54_06_03_01/computed_mean_sd6.csv", sep = ",")        

####################### Graph Setup and Lables #######################
next_sim = 11
x = 1:20
next_graph = "r_70_69.0_0.6_0.3_0.1"
computed_mean_sd = r_70_69.0_0.6_0.3_0.1
computed_mean_sd$total = computed_mean_sd$R1_mean+computed_mean_sd$R2_mean+computed_mean_sd$R3_mean
computed_mean_sd$total_sd = sqrt(computed_mean_sd$R1_sd^2 +computed_mean_sd$R2_sd^2+computed_mean_sd$R3_sd^2)


# Joint Plots
title_parts = strsplit(next_graph, "_")
title = paste('Total Infection Pr (%): ', title_parts[[1]][2],
              ' 
              d1 =', title_parts[[1]][3],
              ', d2 =', title_parts[[1]][4],
              ', d3 =', title_parts[[1]][5],
              ', beyond =', title_parts[[1]][6]
              , sep = " ")

####################### Stacked graphs #######################
reshaped_data = computed_mean_sd
reshaped_data = reshaped_data[-1]
reshaped_data = reshaped_data[-3]
reshaped_data = reshaped_data[-4]
reshaped_data = reshaped_data[-5]
reshaped_data = reshaped_data[-6]
reshaped_data = reshaped_data[-7]
reshaped_data[2] = reshaped_data[5] - reshaped_data[6]
reshaped_data = reshaped_data[-3]
reshaped_data[1] = reshaped_data[1]*3/10
names(reshaped_data) = c("r", "Infected, not Resistant ", "Resistant to 3", "Total Infected", "Resistant to at least 1")

reshaped_data <- melt(reshaped_data, id="r")
reshaped_data$variable <- factor(reshaped_data$variable, unique(reshaped_data[order(reshaped_data$value, decreasing = T),"variable"]) )

p <- ggplot(reshaped_data, aes(r, value))
p <- p + labs(x = "Effectiveness of triple regimen(%)", y = "# of cells with resistance") + ggtitle(title)  
p <- p + geom_area(aes(colour = variable, fill= variable), position = 'identity') 
p <- p + theme(legend.position="none")
p <- p + ylim(0, 5200)
p  

####################### Line graph with error bars #######################

plot(computed_mean_sd$R1_mean, cex=1,xaxt='n', ylim=c(0, 1800), 
     xlab='Effectiveness of triple regimen(%)',ylab='# of cells with resistance', 
     main= title, col='olivedrab3', pch=16, type = 'l')
legend("topleft", inset = 0.03, title = "Resistance to:",
       legend = c("1 drug","2 drugs","3 drugs","At least 1"),
       lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5), pch = c(16, 15, 17, 18),
       col=c("olivedrab3","orange", "purple", "deepskyblue2"),
       cex = 0.8)


points(computed_mean_sd$R1_mean, cex=1,xaxt='n',col='olivedrab3',pch=16)
axis(1, at=x, labels=Rvalues*0.3)
se.up = computed_mean_sd$R1_mean+(computed_mean_sd$R1_sd/2)
se.dn = computed_mean_sd$R1_mean-(computed_mean_sd$R1_sd/2)
arrows(x,se.dn,x,se.up,code=3,length=0.2,angle=90,col='olivedrab3')

lines(computed_mean_sd$R2_mean, cex=1,xaxt='n',col='orange',type = 'l')
points(computed_mean_sd$R2_mean, cex=1,xaxt='n',col='orange',pch=15)
se.up = computed_mean_sd$R2_mean+(computed_mean_sd$R2_sd/2)
se.dn = computed_mean_sd$R2_mean-(computed_mean_sd$R2_sd/2)
arrows(x,se.dn,x,se.up,code=3,length=0.2,angle=90,col='orange')

se.up = computed_mean_sd$R3_mean+(computed_mean_sd$R3_sd/2)
se.dn = computed_mean_sd$R3_mean-(computed_mean_sd$R3_sd/2)
lines(computed_mean_sd$R3_mean, cex=1,xaxt='n',col='purple',type = 'l')
points(computed_mean_sd$R3_mean, cex=1,xaxt='n',col='purple',pch=17)
arrows(x,se.dn,x,se.up,code=3,length=0.2,angle=90,col='purple')

se.up = computed_mean_sd$total+(computed_mean_sd$total_sd/2)
se.dn = computed_mean_sd$total-(computed_mean_sd$total_sd/2)
lines(computed_mean_sd$total, cex=1,xaxt='n',col='deepskyblue2',type = 'l')
points(computed_mean_sd$total, cex=1,xaxt='n',col='deepskyblue2',pch=18)
arrows(x,se.dn,x,se.up,code=3,length=0.2,angle=90,col='deepskyblue2')
