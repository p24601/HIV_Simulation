############# To Do List ############# 
# Programming:
# - Implement drugs
# - Implement drug resistance 
# - Add saving of plots programmatically
# - Implement deep copy for grids. Importatnt for tranfer or mutations. 
# - Make Grid an object so multiple grids can be run at the same time
# - Implement multithreading support. One grid per core.
# - Generate 7 mutations per epoch. Currently generating one.
# - Split code by classes(cell and grid)
#
# Model:
# - Can a cell that acquire a mutation in an epoch, transfer it to a newly 
#   infected cell in the same epoch?
# - Improve decision logic as to which infected cell is passing on mutations
#   when multiple are available in a neighbourhood
# - Add an Immune System Resistance parameter to improve infection phase 
#   transitions, and the curve shapes for latent and final phases

############# Install Packages ############# 
### install.packages("plotly")
### install.packages("hash")
###########################################

############# Install Packages ############# 
#install.packages(devtools)
#install.packages("viridis", type="source")
#devtools::install_github("hadley/ggplot2")
#devtools::install_github("ropensci/plotly")
#install.packages("hash")
###########################################


############# Libraries ############# 
library(plotly)
library(ggplot2)
library(hash)
library(plyr)

############# CA States and Parameters ###############
# States
# State 1: H:   Healthy          (Color- Green)
# State 3: I_1: Infected  		   (Color- Cyan)
# State 2: D:   Dead             (Color- Black)

# Parameters
n = 100;            # grid dimensions n x n
P_HIV = 0.05;       # initial grid will have P_hiv acute infected cells
P_i = 0.997;     	  # probability of infection by neighbors
P_v = 0.00001;      # probability of infection by random viral contact
P_rep = 0.99;       # probability of dead cell being replaced by healthy
P_repI = 0.00001;   # probability of dead cell being replaced by infected
tau = 4;            # time delay for an I cell to become D  
totalsteps = 40;   # total number of weeks of simulation to be performed

# Vector of possible amino-acids 
aa = c("A", "R", "N", "D", "C", "Q", "E", 
       "G", "H", "I", "L", "K", "M", "F",
       "P", "S", "T", "W", "Y", "V")

# Timesteps for which we want to save simulation image
savesteps = c(3, 7, 11, 15, 20, 25, 50, 100, 150, 
              200, 250, 300, 350, 400, 450, 500)

############# Cell Class Definition and Methods ###############
cell <- setClass(
  # Set the name for the class
  "cell",
  
  # Define the slots
  slots = c(
    state = "numeric",
    infected_epochs = "numeric",
    mutations = "hash"
  ),
  
  # Set the default values for the slots.
  prototype=list(
    state = 1,
    infected_epochs = 0
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


# Set method for the infected_epochs of a cell
setGeneric(name="setInfected_epochs",
           def= function(aCell, ninfected_epochs){
             standardGeneric("setInfected_epochs")
           })

setMethod(f = "setInfected_epochs", signature = "cell",
          definition = function(aCell, ninfected_epochs) {
            aCell@infected_epochs = ninfected_epochs
            return (aCell)
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

############# CA Grid Definition and Preparation ###############
# Create initial n x n grid and populate it with cell S4 Objects. Default parameters for each 
# cell are set: initial state: a random value is generante and if equal or less than P_HIV, 
# state = 3 (Infected), otherwise state = 1 (Healthy). Each cell is also assigned an empty 
# environment to track mutations. 
grid <- matrix( sapply(1:(n*n), function(x){
                                            new("cell",
                                            state = ifelse(runif(1)<= P_HIV, 3, 1),
                                            infected_epochs = 0,
                                            mutations = hash())
                                            }), n, n)

#Set the state of the cells at the edges of the grid to Healthy state
grid[,c(1, n)] = lapply(grid[,c(1,n)], function(x) setState(x, 1))
grid[c(1, n),] = lapply(grid[c(1,n),], function(x) setState(x, 1))

# 		  NOTE: Our CA only simulates from rows 2 to n-1, and columns 2 to n-1.
#       This is to prevent the edge row and column cells from having an
#       out-of-bounds error when checking the neighbors around them.
#       The edge values are all set to H state so that it does not affect
#       Rule 1 of the CA for cells next to them

# Initialize a grid for to visualize cell states all at once
stateGrid = matrix(ncol = n, nrow = n)
stateGrid[,] = sapply(grid[,], function(x) getState(x))
gridPlot_0 = plot_ly(z = stateGrid, type = "heatmap")
gridPlot_0

# Initialize a grid to subset a cell's negihbourhood
current_neighbourhood = matrix(nrow = 3, ncol = 3)

# Initialize a grid to keep count of cells in every state
states_count = matrix(0, nrow = totalsteps, ncol = 3)
colnames(states_count) <- c("Healthy", "Infected", "Dead")

############# Simulation ###############
timestep = 1; 
while(timestep <= totalsteps){
  	nextgrid = grid;
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
 	 		        # indeces and pass mutation map to the infected cell. If infection is
 	 		        # coming from a remote cell, pick a random infected cell from the CA. 
 	 		        if (c2){
 	 		          x_c = x + position[[1,1]]
 	 		          y_c = y + position[[1,2]]
 	 		        }else{
 	 		          stateGrid[,] = sapply(grid[,], function(x) getState(x))
 	 		          position = which(stateGrid[,] == 3, arr.ind = TRUE)
 	 		          r = sample(1:nrow(position),1)
 	 		          x_c = position[[r,1]]
 	 		          y_c = position[[r,2]]
 	 		        }
 	 		        nextgrid[[x,y]]@mutations = grid[[x_c, y_c]]@mutations
 	 		      }
 	 		      next
 	 		    }
  
 	 		    
    		  # Rule 2 a and b
          # If the cell is in I state, the viral genome is subject to 
 	 		    # mutations. Mutations occur at a rate of 1/day or 7/epoch (2.a)
 	 		    # If a cell has been in this state for tau timesteps, the cell
 	 		    # becomes D state (dead). In this case the accumulated mutations
 	 		    # are lost (2.b) 
          if((grid[[x,y]]@state == 3)){
              nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1
              
              # Randomly genearate a mutation site, and a mutation amino acid. 
              mutation_site = sample(1:3333,1)
              mutation_aa = aa[sample(1:20, 1)]
              
              # Assign the mutation to a hashmap.
              nextgrid[[x,y]]@mutations[[as.character(mutation_site)]] = mutation_aa
              
              # Reset epochs to 0 and clear mutations once a cell dies.
              if(grid[[x,y]]@infected_epochs == tau){
                  nextgrid[[x,y]]@state = 2
                  nextgrid[[x,y]]@infected_epochs = 0
                  clear(nextgrid[[x,y]]@mutations)
              }
              next
          }

		      # Rule 3 a and b
          # If the cell is in D state, then the cell will become H state
          # with probability P_rep (3.a) or I state with probability
          # P_repI (3.b)
          if(grid[[x,y]]@state == 2 && runif(1) <= P_rep){
              nextgrid[[x,y]]@state = ifelse(runif(1) <= P_repI, 3, 1)
              next
          }
 	 		    
 	 		 } #End inner for loop
  	} #End outer for loop
  	
  	# Assign the updates of this timestep in nextgrid back to our grid
  	grid = nextgrid
  
  	# Update aggregate counts
  	stateGrid[,] = sapply(grid[,], function(x) getState(x))
  	states_count[timestep,1] = sum(stateGrid[] == 1)
  	states_count[timestep,2] = sum(stateGrid[] == 3)
  	states_count[timestep,3] = sum(stateGrid[] == 2)
  	
  	# Save plot of current state of the grid to variable
  	if(timestep %in% savesteps){
  	  assign(paste("gridPlot_", timestep, sep = ""), plot_ly(z = stateGrid, type = "heatmap"))
  	}  	
  	
  	# Move to the next timestep
  	timestep = timestep+1

} #End while loop

############# Aggregate Data and Analysis ###############
# Graphing of cell counts per epoch
series1 <- states_count[,1]
series2 <- states_count[,2]
series3 <- states_count[,3]
plot(series1, type="l", lty=1, col="green", xlab = "Weeks", ylab = "Cell Count")
lines(series2, type="l", lty=1, col="red")
lines(series3, type="l", lty=1, col="blue")

# Display mutation list per cell, at the end of the simulation. 
# Caution: Generates a very large element.
cell_mutations = matrix(0, nrow = n*n, ncol = 1)
colnames(cell_mutations) <- c("Mutations List")
cell_mutations[,] = sapply(grid[,], function(x) as.list(getMutations(x)))


