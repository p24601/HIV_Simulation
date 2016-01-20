############# To Do List ############# 
# Programming:
# - Make Grid an object so multiple grids can be run at the same time
# - Implement multithreading support. One grid per core.
# - Split code by classes(cell and grid)
#
# Model:
# - When is drug therapy started? Currently starts at week 20 as in the
#   original paper.
# - Can a cell that acquire a mutation in an epoch, transfer it to a newly 
#   infected cell in the same epoch?
# - Improve decision logic as to which infected cell is passing on mutations
#   when multiple are available in a neighbourhood
# - Add an Immune System Resistance parameter to improve infection phase 
#   transitions, and the curve shapes for latent and final phases

############# Install Packages ############# 
### install.packages("plotly")
### install.packages("hash")

############# Install Packages ############# 
#install.packages(devtools)
#install.packages("viridis", type="source")
#devtools::install_github("hadley/ggplot2")
#devtools::install_github("ropensci/plotly")
#install.packages("hash")

############# Libraries ############# 
library(plotly)
library(ggplot2)
library(hash)
library(plyr)

############# Drug resistance-conferring mutations ############# 
#Drugs
abacavir_list = list("65" = c("R", "E", "N"), "74" = "V", "115" = "F", "184" = "V")

# Locations are what is on the paper+230. 230 is the esitmate length of the previous protein in the sequence
efavirenz_list = list("330" = "I", "331" = "P", "333" = c("N", "S"), "336" = "M", 
                      "338" = "I", "411" = c("C", "I"), "418" = "L", "420" = c("S", "A"), 
                      "455" = "H", "460" = "L")

# Locations are what is on the paper + (230 + 230). 230 is the esitmate length of the previous protein in the sequence
darunavir_list = list("471" = "I", "492" = "I", "493" = "F", "507" = "V", "510" = "V", 
                      "514" = c("M", "L"), "534" = "P", "536" = "V", "544" = "V", 
                      "549" = "V")

# Build a hashmap of resistance sites
resistanceSites_list = append(abacavir_list, efavirenz_list)
resistanceSites_list = append(resistanceSites_list, darunavir_list)
resistanceSites_env = list2env(resistanceSites_list)

# Remove unused lists
remove(abacavir_list)
remove(efavirenz_list)
remove(darunavir_list)
remove(resistanceSites_list)

############# Cell Class Definition and Methods ###############
cell <- setClass(
  # Set the name for the class
  "cell",
  
  # Define the slots
  slots = c(
    state = "numeric",
    infected_epochs = "numeric",
    mutations = "environment",
    resistance = "numeric"
  ),
  
  # Set the default values for the slots.
  prototype=list(
    state = 1,
    infected_epochs = 1,
    resistance = 0
  )
)

# Set method for the state of a cell
setGeneric(name="setState",
           def= function(aCell, aState){
           standardGeneric("setState")
           })

setMethod(f = "setState", signature = "cell",
          definition = function(aCell, aState) {
          aCell@state = aState
          return (aCell)
          })

# Get method for the state of a cell
setGeneric(name="getState",
           def= function(aCell){
             standardGeneric("getState")
           })

setMethod(f = "getState", signature = "cell",
          definition = function(aCell) {
            return (aCell@state)
          })

# Get method for the infected_epochs of a cell
setGeneric(name="getInfected_epochs",
           def= function(aCell){
             standardGeneric("getInfected_epochs")
           })

setMethod(f = "getInfected_epochs", signature = "cell",
          definition = function(aCell) {
            return (aCell@infected_epochs)
          })

# Get method for the mutation list of a cell
setGeneric(name="getMutations",
           def= function(aCell){
             standardGeneric("getMutations")
           })

setMethod(f = "getMutations", signature = "cell",
          definition = function(aCell) {
            return (aCell@mutations)
          })

# Get method for the resistance levels of a cell
setGeneric(name="getResistance",
           def= function(aCell){
             standardGeneric("getResistance")
           })

setMethod(f = "getResistance", signature = "cell",
          definition = function(aCell) {
            return (aCell@resistance)
          })

############# CA States and Parameters ###############
# States
# State 1: H:   Healthy         
# State 3: I_1: Infected  		  
# State 2: D:   Dead            

# Parameters
n = 100                         # grid dimensions n x n
P_HIV = 0.05                    # initial grid will have P_hiv acute infected cells
P_i = 0.997/7             	      # probability of infection by neighbors
P_v = 0.00001/7                   # probability of infection by random viral contact
P_rep = 0.99/7                    # probability of dead cell being replaced by healthy
tau = 4                         # time delay for an I cell to become D  
hiv_total_aa = 2876             # Estimated total number of amino acids of the HIV-1 proteinome
base_drug_efficiency = 0.30     # base probability that the triple cocktail will kill and infected cell
start_of_therapy = 20           # epoch at which to start drug therapy
totalsteps = 100                 # total number of weeks of simulation to be performed
resiliance = 10                 # number of epochs before the start of phase 2 of infection

# Immune system parameters
is_capacity = 100 + (resiliance/totalsteps)
fatigue_ir = (is_capacity/totalsteps)

# Vector of possible amino-acids 
aa = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", 
       "S", "T", "W", "Y", "V")

# Timesteps for which we want to save simulation status
savesteps = seq(6, 60, by = 5)

# Initialize lists to track the state of cell attribute at each savestep.
stateGrid_list = list()
k = 2

############# CA Grid Definition and Preparation ###############
# Create initial n x n grid and populate it with cell S4 Objects. Default parameters for each 
# cell are set: initial state: a random value is generante and if equal or less than P_HIV, 
# state = 3 (Infected), otherwise state = 1 (Healthy). Each cell is also assigned an empty 
# environment to track mutations. 
grid <- matrix( sapply(1:(n*n), function(x){ new("cell",
                                             state = ifelse(runif(1)<= P_HIV, 3, 1),
                                             mutations = new.env(hash = TRUE))
                                             }), n, n)

# Set the state of the cells at the edges of the grid to Healthy state
grid[,c(1, n)] = lapply(grid[,c(1,n)], function(x) setState(x, 1))
grid[c(1, n),] = lapply(grid[c(1,n),], function(x) setState(x, 1))

# 		  NOTE: Our CA only simulates from rows 2 to n-1, and columns 2 to n-1.
#       This is to prevent the edge row and column cells from having an
#       out-of-bounds error when checking the neighbors around them.
#       The edge values are all set to H state so that it does not affect
#       Rule 1 of the CA for cells next to them

# Show initial status of the grid
stateGrid = matrix(0, nrow = n, ncol = n)
stateGrid[,] = sapply(grid[,], function(x) getState(x))
plot_ly(z = stateGrid, type = "heatmap")
stateGrid_list[[1]] = stateGrid

# Initialize a grid to subset a cell's negihbourhood
current_neighbourhood = matrix(nrow = 3, ncol = 3)

# Initialize a grid to keep track of cell attributes during simulation
states_count = matrix(0, nrow = totalsteps, ncol = 3)
colnames(states_count) = c("Healthy", "Infected", "Dead")

infectedEpochsGrid = matrix(0, nrow = n, ncol = n)
infectedEpochs_count = matrix(0, nrow = totalsteps, ncol = 4)
colnames(infectedEpochs_count) = c("0", "1", "2", "3")

resistanceGrid = matrix(0, nrow = n, ncol = n)
resistance_count = matrix(0, nrow = totalsteps, ncol = 4)
colnames(resistance_count) = c("1", "2", "3", "Overall")

genotypes_count = matrix(0, nrow = totalsteps, ncol = 1)
colnames(genotypes_count)[1] = "Number of Genomes"

############# Simulation ###############
timestep = 1; 
while(timestep <= totalsteps){
  	nextgrid = matrix( sapply(1:(n*n), function(x){ new("cell", state = 1,
                                            	      mutations = new.env(hash = TRUE))
                                            	      }), n, n)
  		for (x in 2:(n-1)){
 	 		  for (y in 2:(n-1)){
          
  	  		# Rule 1
  	  		# If the cell is in H state and at least one of its neighbors
  	  		# is in I state then the cell becomes I with a probability of 
   	 		  # P_i (Conditions c1 and c2). The cell may also becomes I1 by 
   	 		  # randomly coming in contact with a virus from outside its 
   	 		  # neighborhood with a probability of P_v (Condition c3).
 	 		    if(grid[[x,y]]@state == 1){
 	 		      
 	 		      # Subset the main grid for the Moore neighbourhood around the current cell
 	 		      current_neighbourhood[,] = sapply(grid[c(x-1, x, x+1),c(y-1, y, y+1)], function(x) getState(x))
 	 		      
 	 		      # Extract the position of any available infected cells
 	 		      position = which(current_neighbourhood == 3, arr.ind = TRUE) - c(2,2)
 	 		      
 	 		      # Evaluate infectivity conditions 
 	 		      c1 = runif(1)<=P_i;
 	 		      c2 = ifelse(length(position) == 0, FALSE, TRUE)
 	 		      c3 = runif(1)<=P_v
 	 		      
 	 		      # Evaluate Rule 1
 	 		      if((c1 && c2) || c3){  
 	 		        nextgrid[[x,y]]@state = 3
 	 		        
 	 		        # If the infecting cell is from the neighbourhood, calculate correct 
 	 		        # indices and pass mutation map to the infected cell. If infection is
 	 		        # coming from a remote cell, pick a random infected cell from the CA. 
 	 		        if (c2){
 	 		          r = sample(1:nrow(position),1)
 	 		          x_c = x + position[[r,1]]
 	 		          y_c = y + position[[r,2]]
 	 		        }else{
 	 		          stateGrid[,] = sapply(grid[,], function(x) getState(x))
 	 		          position = which(stateGrid[,] == 3, arr.ind = TRUE)
 	 		          r = sample(1:nrow(position),1)
 	 		          x_c = position[[r,1]]
 	 		          y_c = position[[r,2]]
 	 		        }
 	 		        list2env(as.list.environment(grid[[x_c,y_c]]@mutations, all.names = TRUE),nextgrid[[x,y]]@mutations)
 	 		        
 	 		        # Generate mutation and save it in the cell mutations hashmap.
 	 		        mutation_site = sample(1:hiv_total_aa,1)
 	 		        mutation_aa = aa[sample(1:20, 1)]
 	 		        nextgrid[[x,y]]@mutations[[as.character(mutation_site)]] = mutation_aa
 	 		        
 	 		        # If the mutation occurs at a drug resistance-conferring site, check if the
 	 		        # amino acid it is mutating to grants resistance. If yes, increase cell 
 	 		        # resistance count, otherwise, decrease it.
 	 		        potential_resistance = resistanceSites_env[[as.character(mutation_site)]]
 	 		        if(!is.null(potential_resistance)){
 	 		          present = mutation_aa %in% potential_resistance
 	 		          
 	 		          if (present){
 	 		            nextgrid[[x,y]]@resistance = nextgrid[[x,y]]@resistance+1
 	 		            
 	 		          }else if (!present && (nextgrid[[x,y]]@resistance > 0)){
 	 		            nextgrid[[x,y]]@resistance = nextgrid[[x,y]]@resistance-1
 	 		          }
 	 		        }
 	 		      }
 	 		      next
 	 		    }#End Rule 1
  
    		  # Rule 2 a and b
          # If the cell is in I state, the viral genome is subject to 
 	 		    # mutations. Mutations occur at a rate of 1/day or 7/epoch (2.a)
 	 		    # If a cell has been in this state for tau timesteps or does not 
 	 		    # have sufficient drug resistance, the cell becomes D state (dead). 
 	 		    # In this case the accumulated mutations are lost (2.b) 
          if((grid[[x,y]]@state == 3)){
              nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1
              list2env(as.list.environment(grid[[x,y]]@mutations, all.names = TRUE),nextgrid[[x,y]]@mutations)
              
              # Our model simplifies drug resistance by assuming all drugs are equally efficient,
              # and letting each single succesful mutation represent resistance to a single drug. 
              # Therefore, once the number of succesful mutations has reached 3, we consider this 
              # cell to be resistant to the 3 drugs begin employed in the triple cocktail. Further 
              # succesful mutations are accounted for, but not used to increase cell fitness further.
              cell_resistance = ifelse(grid[[x,y]]@resistance <= 2, grid[[x,y]]@resistance, 3)
              drug_efficiency = ifelse(timestep >= start_of_therapy, (base_drug_efficiency * (3 - cell_resistance)), 0) 
              
              # Kill cell if it has lived enough epochs or drug takes over
              if((grid[[x,y]]@infected_epochs >= tau) || (runif(1) <= drug_efficiency)){
                nextgrid[[x,y]]@state = 2
                nextgrid[[x,y]]@infected_epochs = 0
                
                # If a cell dies, its mutations are lost
                nextgrid[[x,y]]@mutations = new.env(hash = TRUE)
                
              }else{
                nextgrid[[x,y]]@state = 3
              }#End else
              next
          }#End Rule 2

		      # Rule 3
          # If the cell is in D state, then the cell will become H state
          # with probability P_rep 
 	 		    if(grid[[x,y]]@state == 2){ 
 	 		        strain = (is_capacity - (fatigue_ir*timestep))/100
   	 		      if(runif(1) <= (P_rep*strain)){
   	 		        nextgrid[[x,y]]@state = 1
   	 		      }else{
   	 		        nextgrid[[x,y]]@state = 2
   	 		      }  
   	 		      next
          }#End Rule 3
 	 		    
 	 		 } #End inner for loop
  	} #End outer for loop
  
  	# Assign the updates of this timestep in nextgrid back to our grid
  	grid = nextgrid
  	
  	# Update aggregate counts
  	stateGrid[,] = sapply(grid[,], function(x) getState(x))
  	states_count[timestep,1] = sum(stateGrid[] == 1)
  	states_count[timestep,2] = sum(stateGrid[] == 3)
  	states_count[timestep,3] = sum(stateGrid[] == 2)
  	
  	infectedEpochsGrid[,] = sapply(grid[,], function(x) getInfected_epochs(x))
  	infectedEpochs_count[timestep,0] = sum(infectedEpochsGrid[] == 0)
  	infectedEpochs_count[timestep,1] = sum(infectedEpochsGrid[] == 1)
  	infectedEpochs_count[timestep,2] = sum(infectedEpochsGrid[] == 2)
  	infectedEpochs_count[timestep,3] = sum(infectedEpochsGrid[] >= 3)
  	
  	resistanceGrid[,] = sapply(grid[,], function(x) getResistance(x))
  	resistance_count[timestep,1] = sum(resistanceGrid[] == 1)
  	resistance_count[timestep,2] = sum(resistanceGrid[] == 2)
  	resistance_count[timestep,3] = sum(resistanceGrid[] >= 3)
  	
  	genotypes = sapply(grid[,], function(x) getMutations(x))
  	genotypes_count[timestep,1] = sum(sapply(genotypes, function (x) length(x)>0) == TRUE)
  	
  	# Save plot of current state of the grid to variable
  	if(timestep %in% savesteps){
  	  stateGrid_list[[k]] = stateGrid
  	  k = k + 1
  	}  	
  	
  	# Move to the next timestep
  	timestep = timestep+1
  	
} #End while loop

############# Clean-Up: CA and Simulation parameters ###############
remove(P_HIV)
remove(P_i)
remove(P_v)
remove(P_rep)
remove(P_repI)
remove(tau)
remove(hiv_total_aa)
remove(totalsteps)
remove(base_drug_efficiency)
remove(start_of_therapy)
remove(aa)
remove(savesteps)
remove(timestep)
remove(x)
remove(x_c)
remove(y)
remove(y_c)
remove(c1)
remove(c2)
remove(c3)
remove(mutation_aa)
remove(mutation_site)
remove(position)
remove(potential_resistance)
remove(present)
remove(cell_resistance)
remove(drug_efficiency)
remove(resistanceSites_env)
remove(i)
remove(grid)
remove(nextgrid)
remove(cell)
remove(getInfected_epochs)
remove(getMutations)
remove(getResistance)
remove(getState)
remove(setState)
remove(r)

############# Analysis: Simulation Overview  ###############
# Plot the state of infection while in progress. 
state_sp1 = subplot(
  plot_ly(z = stateGrid_list[[1]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[2]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[3]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[4]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[5]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[6]], type = "heatmap"),
  nrows = 2,
  margin =  0.05
)

state_sp2 = subplot(
  plot_ly(z = stateGrid_list[[7]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[8]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[9]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[10]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[11]], type = "heatmap"),
  plot_ly(z = stateGrid_list[[12]], type = "heatmap"),
  nrows = 2,
  margin =  0.05
)
state_sp1
state_sp2

# Graphing of cell counts per epoch
number_of_cells = n*n
Healthy = states_count[,1]/number_of_cells
Infected = states_count[,2]/number_of_cells
Dead = states_count[,3]/number_of_cells

simulation_overview = plot(Healthy, 
                           type="l", lty=1, lwd = 2,
                           col="green", 
                           xlab = "Weeks", 
                           ylab = "Cell Count", 
                           main = "Simulation Overview")

lines(Infected, type="l", lty=1, col="red")
lines(Dead, type="l", lty=1, col="blue")

legend("topright", 
       c("Healthy", "Infected", "Dead"),
       lty= c(1,1,1), lwd = c(2,2,2), 
       col = c("green", "red", "blue"),
       cex = .75)

remove(Healthy)
remove(Infected)
remove(Dead)

############# Analysis: Cell Infected Epochs ###############
# Line plot the infected epochs per cell over time. 
Zero_eps = infectedEpochs_count[,1]
One_eps = infectedEpochs_count[,2]
Two_eps = infectedEpochs_count[,3]
Three_eps = infectedEpochs_count[,4]

plot(Zero_eps, 
     type="l", lty=1, lwd = 2,
     col="green", 
     xlab = "Weeks", 
     ylab = "Cell Count", 
     main = "Cell Infected Epochs")

lines(One_eps, type="l", lty=1, col="red")
lines(Two_eps, type="l", lty=1, col="blue")
lines(Three_eps, type="l", lty=1, col="orange")

legend("bottomright", 
       c("Cells alive for 0 epoch","Cells alive for 1 epoch", 
         "Cells alive for 2 epoch", "Cells alive for 3 epoch"),
       lty= c(1,1,1,1), lwd = c(2,2,2,2), 
       col = c("green", "red", "blue", "orange"),
       cex = .75)

remove(Zero_eps)
remove(One_eps)
remove(Two_eps)
remove(Three_eps)

############# Analysis: Cells with Resistance ###############
# Plot the number of cells with resistance over time.
resistance_count[,4] = rowSums(resistance_count)

One_r = resistance_count[,1]
Two_r = resistance_count[,2]
Three_r = resistance_count[,3]
Overall_r = resistance_count[,4]

plot(One_r, 
     type="l", lty=1, lwd = 2,
     col="green", 
     xlab = "Weeks", 
     ylab = "Cell Count", 
     main = "Cells with Resistance")

lines(Two_r, type="l", lty=1, col="blue")
lines(Three_r, type="l", lty=1, col="orange")
lines(Overall_r, type="l", lty=1, col="red")

legend("topleft", 
       c("Cells resistant to 1 drug","Cells resistant to 2 drug", 
         "Cells resistant to 3 drug", "Overall count of cells resistant to drugs"), 
       lty= c(1,1,1,1), lwd = c(2,2,2,2), 
       col = c("green", "blue", "orange", "red"),
       cex = .75)
       

remove(One_r)
remove(Two_r)
remove(Three_r)
remove(Overall_r)

############# Analysis: Genotypes ###############
# Plot the number of genotypes present in the grid
plot(genotypes_count[,1], 
     type="l", lty=1, lwd =2, 
     col="red", 
     xlab = "Weeks", 
     ylab = "Genotypes", 
     main = "Number of Genotypes")

legend("topright", "# of Genotypes", lty= 1, lwd =2, col = "red" )
