#############  Libraries ############# 
library(plotly)
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
totalsteps = 50;   # total number of weeks of simulation to be performed
savesteps = c(3, 7, 11, 15, 20, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500); 
#timesteps for which we want to save simulation image

############# Cell Class Definition and Methods ###############
cell <- setClass(
  # Set the name for the class
  "cell",
  
  # Define the slots
  slots = c(
    state = "numeric",
    infected_epochs = "numeric",
    mutations = "environment"
  ),
  
  # Set the default values for the slots.
  prototype=list(
    state = 1,
    infected_epochs = 0,
    mutations = new.env()
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

############# CA Grid Definition and Preparation ###############
# Create initial n x n grid and populate it with cell S4 Objects. Default parameters for each 
# cell are set: initial state: a random value is generante and if equal or less than P_HIV, 
# state = 3 (Infected), otherwise state = 1 (Healthy). Each cell is also assigned an empty 
# environment to track mutations. 
grid <- matrix( sapply(1:(n*n), function(x){
                new("cell",
                    state = ifelse(runif(1)<= P_HIV, 3, 1),
                    infected_epochs = 0,
                    mutations = new.env()) }
                ), n, n)


#Set the state of the cells at the edges of the grid to Healthy state
grid[,c(1, n)] = lapply(grid[,c(1,n)], function(x) setState(x, 1))
grid[c(1, n),] = lapply(grid[c(1,n),], function(x) setState(x, 1))

# 		  NOTE: Our CA only simulates from rows 2 to n-1, and columns 2 to n-1.
#       This is to prevent the edge row and column cells from having an
#       out-of-bounds error when checking the neighbors around them.
#       The edge values are all set to H state so that it does not affect
#       Rule 1 of the CA for cells next to them

stateGrid = matrix(ncol = n, nrow = n)
stateGrid[,] = sapply(grid[,], function(x) getState(x))
plot_ly(z = stateGrid, type = "heatmap")

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
 	 		      
 	 		      c1 = runif(1)<=P_i;
 	 		      c2 = grid[[x-1,y-1]]@state == 3 || grid[[x-1,y]]@state == 3 || grid[[x-1,y+1]]@state == 3 || 
 	 		           + grid[[x,y-1]]@state == 3 || grid[[x,y+1]]@state == 3 || grid[[x+1,y-1]]@state == 3 || 
 	 		           + grid[[x+1,y]]@state == 3 || grid[[x+1,y+1]]@state == 3
 	 		      c3 = runif(1)<=P_v
 	 		      
 	 		      if((c1 && c2) || c3){                    
 	 		        nextgrid[[x,y]]@state = 3
 	 		      }
 	 		    }
  
    		  # Rule 2
          # If the cell is in I state, the viral genome is subject to 
 	 		    # mutations. Mutations occur at a rate of 1/day or 7/epoch,
 	 		    # with probability P_m.
 	 		    
 	 		    # Also and has been in this state for
          # tau timesteps, the cell becomes D state (dead)
          if((grid[[x,y]]@state == 3)){
              nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1
              
              if(grid[[x,y]]@infected_epochs == tau){
                  nextgrid[[x,y]]@state = 2
                  nextgrid[[x,y]]@infected_epochs = 0
              }
          }

		      # Rule 3 a and b
          # If the cell is in D state, then the cell will become H state
          # with probability P_rep (3.a) or I state with probability
          # P_repI (3.b)
          if(grid[[x,y]]@state == 2 && runif(1) <= P_rep){
              nextgrid[[x,y]]@state = ifelse(runif(1) <= P_repI, 3, 1)
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
  	
  	if(timestep %in% savesteps){
  	    plot_ly(z = grid, type = "heatmap")
  	}  	
  	# Move to the next timestep
  	timestep = timestep+1

} #End while loop


series1 <- states_count[,1]
series2 <- states_count[,2]
series3 <- states_count[,3]
plot(series1, type="l", lty=1, col="green", xlab = "Weeks", ylab = "Cell Count")
lines(series2, type="l", lty=1, col="red")
lines(series3, type="l", lty=1, col="blue")








