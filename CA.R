# l: number of simulations produced by a single run of the code
for (l in 1:10){

        # base_dest: folder directory for depositing simulation results. Must match
        # the base_drug_efficiency for the simulation.
        base_dest = "/home/pcervell/hiv_simulation/results_1/300/"

base_drug_efficiency = 0.300     # base probability that the triple cocktail will kill and infected cell
Pi_n1 = 0.54                     # probablity of infection in a cell's immediate negihbourhood
Pi_n2 = 0.546                    # probablity of infection in a cell's immediate negihbourhood+1
Pi_n3 = 0.549                    # probablity of infection in a cell's immediate negihbourhood+2
P_v = -1                        # probablity of infection from elsewhere in the grid

# Generate simulation results subfolder for this run
dir.create(paste(base_dest, l, sep = ""))

############# Libraries #############
library(ggplot2)
library(hash)
library(plyr)

############# Drug resistance-conferring mutations #############
# Drugs are represented as key value pairs where the key is the position of the AA
# prone to conferring resistance, and the value represents the mutant AA that
# confers resistance. (Ex: "74" = "V" means that a V at position 74 confers resistance,
# and "65" = c("R", "E", "N") means that either an R, E, or N at position 65 confer resistance.)
# Key value pairs are entered as reported in the "2014 IAS update on hiv drug resistance mutations"
abacavir_list   = list("65" = c("R", "E", "N"), "74" = "V", "115" = "F", "184" = "V")

lamivudine_list = list("65" = c("R", "E", "N"), "184" = c("V", "I"))

tenofovir_list  = list("65" = c("R", "E", "N"), "70" = "E")

zidovudine_list = list("41" = "L", "67" = "N", "70" = "R", "210" = "W", "215" = c("Y", "F"), "219" = c("Q", "E"))

# Locations are what is on the paper+230. 230 is the esitmated length of the previous protein in the sequence
nevirapine_list = list("330" = "I", "331" = "P", "333" = c("N", "S"), "336" = c("M", "A"),
                      "338" = "I", "411" = c("C", "I"), "418" = c("L", "C", "H"), "420" = "A",
                      "460" = "L")

efavirenz_list  = list("330" = "I", "331" = "P", "333" = c("N", "S"), "336" = "M",
                      "338" = "I", "411" = c("C", "I"), "418" = "L", "420" = c("S", "A"),
                      "455" = "H", "460" = "L")

# Locations are what is on the paper + (230 + 230). 230 is the esitmated length of the previous protein in the sequence
darunavir_list  = list("471" = "I", "492" = "I", "493" = "F", "507" = "V", "510" = "V",
                      "514" = c("M", "L"), "534" = "P", "536" = "V", "544" = "V",
                      "549" = "V")




# Build a complete hashmap of resistance sites (resistanceSites_env)
# Drug combinations simulated:
# 1)
# 2)
# 3)
resistanceSites_drug1  = list2env(lamivudine_list)
resistanceSites_drug2  = list2env(zidovudine_list)
resistanceSites_drug3  = list2env(efavirenz_list)


############# Cell Class Definition and Methods ###############
cell <- setClass(
  # Set the name for the class
  "cell",

  # Define the slots
  # state: one of Healthy(1), Infected(3), Neutralized(4), dead(2)
  # infected_epochs: number of timesteps spent in state 3
  # mutations: hash map of cell's mutations
  # resistance: number of drugs cell is resistant to (0 to 3)
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
Healthy     = 1
Dead        = 2
Infected    = 3
Neutralized = 4

# Parameters
n                = 100       # grid dimensions n x n
P_rep            = 0.99      # probability of dead cell being replaced by healthy
tau              = 4         # time delay for an I cell to become D
hiv_total_aa     = 3239      # Estimated total number of amino acids of the HIV-1 proteinome
start_of_therapy = 20*7      # epoch at which to start drug therapy
totalsteps       = 200       # total number of weeks of simulation to be performed

logging = FALSE
if (logging){logFile = paste(base_dest, "log_file.txt", sep = "")}

# Vector of possible amino-acids
aa = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
       "S", "T", "W", "Y", "V")

# Timesteps for which we want to save simulation status
savesteps = seq(1, totalsteps, by = 1)

# Initialize lists to track the state of cell attribute at each savestep.
stateGrid_list = list()
resistanceGrid_list = list()
k = 2

############# CA Grid Definition and Preparation ###############
# Create initial n x n grid and populate it with cell S4 Objects.
# Each cell is assigned an empty environment to track mutations.
# All cell but the one in the center of the grid are started in a
# healthy state

# Simulating infection with a single entry point.
grid <- matrix( sapply(1:(n*n), function(x){ new("cell",
                                                 state = 1,
                                                 mutations = new.env(hash = TRUE))
}), n, n)

grid[[n/2,n/2]]@state = 3


# Show initial status of the grid
stateGrid = matrix(0, nrow = n, ncol = n)
stateGrid[,] = sapply(grid[,], function(x) getState(x))
stateGrid_list[[1]] = stateGrid

mutationGrid_list = list()
mutation_list = matrix(0, nrow = n, ncol = n)

# Initialize a grid to subset expanding cell's negihbourhood
neighbourhood_1 = matrix(nrow = 3, ncol = 3)
neighbourhood_2 = matrix(nrow = 5, ncol = 5)
neighbourhood_3 = matrix(nrow = 7, ncol = 7)

# Initialize matrices to collect simulation results.
states_count = matrix(0, nrow = totalsteps, ncol = 4)
colnames(states_count) = c("Healthy", "Infected", "Dead", "Neutralized")

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

    # NOTE: Our CA only simulates from rows 4 to n-3, and columns 4 to n-3.
    # This is to prevent the edge row and column cells from having an
    # out-of-bounds error when checking the neighbors around them.
    # The edge values are all set to H state so that it does not affect
    # Rule 1 of the CA for cells next to them

    # Each timestep, the CA grid is updated based on 4 rules.
    # Rules
    # Rule 1: H -> I
    # Rule 2: D -> H
    # Rule 3: I -> D
    # Rule 4: I -> N
    for (x in 4:(n-3)){
        for (y in 4:(n-3)){

            # Rule 2
            # If the cell is in I state, the viral genome is subject to
            # mutations. Mutations occur at a rate of 1/day or 7/epoch (2.a)
            # If a cell has been in this state for tau timesteps or does not
            # have sufficient drug resistance, the cell becomes D state (dead).
            # In this case the accumulated mutations are lost (2.b)
            if((grid[[x,y]]@state == Infected)){

                # Drug resistance is simplified by assuming all drugs are equally efficient,
                # and letting each single succesful mutation represent resistance to a single drug.
                # Therefore, once the number of succesful mutations has reached 3, we consider this
                # cell to be resistant to the 3 drugs begin employed in the triple cocktail.

                # Extract cell drug resistance attribute
                cell_resistance = grid[[x,y]]@resistance

                # If a cell has for tau timesteps, it is at the end of its life cycle.
                # The cell state is set to dead, its attributes are reset to 0, and its
                # mutation map to empty.
                if(grid[[x,y]]@infected_epochs >= tau){
                    # Grid is updated immediately to allow Rules 1 and 3 to complete immediately after
                    grid[[x,y]]@state = 2

                    # NextGrid is updated to set up for the next timestep
                    nextgrid[[x,y]]@state = 2
                    nextgrid[[x,y]]@infected_epochs = 0

                    # If a cell dies, its mutations are lost
                    nextgrid[[x,y]]@mutations = new.env(hash = TRUE)

                # After the start of drug therapy, infected cells can be neutralized with Pr
                # dependent on the amount of drug resistance accumulated. The more resistance,
                # the less likely to be neutralized. A neutralized infected cell, cannot spread infection.
                # The probablity of becoming neutralized is calculated as base_drug_efficiency * (3 - cell_resistance)
                }else if((timestep >= start_of_therapy) & (runif(1) <= (base_drug_efficiency * (3 - cell_resistance)))){
                    # Grid is updated immediately to allow Rules 1 and 3 to complete immediately after
                    grid[[x,y]]@state = 4

                    # NextGrid is updated to set up for the next timestep
                    nextgrid[[x,y]]@state = 4
                    nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1

                # If an infected cell is not killed, or neutralized, it remains in infected state
                # No new mutations are assigned. This only happens when a healthy cell becomes
                # infected.
                }else{
                    # NextGrid is updated to set up for the next timestep
                    nextgrid[[x,y]]@state = 3
                    nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1

                    # Cell's mutation map is copied over to the same cell in nextgrid
                    list2env(as.list.environment(grid[[x,y]]@mutations, all.names = TRUE),nextgrid[[x,y]]@mutations)

                }#End else
                next
            }#End Rule 2


            # Rule 4
            # A cell that enters neutralized state will remain in neutralized
            # state until it has lived tau timesteps and finally transition to
            # state dead.
            if((grid[[x,y]]@state == Neutralized)){
                # Kill neutralized cell if it has lived enough epochs
                if(grid[[x,y]]@infected_epochs >= tau){
                    # Grid is updated immediately to allow Rules 1 and 3 to complete immediately after
                    grid[[x,y]]@state = 2

                    # NextGrid is updated to set up for the next timestep
                    nextgrid[[x,y]]@state = 2
                    nextgrid[[x,y]]@infected_epochs = 0

                    # If a cell dies, its mutations are lost
                    nextgrid[[x,y]]@mutations = new.env(hash = TRUE)

                }else{
                    # NextGrid is updated to set up for the next timestep
                    nextgrid[[x,y]]@state = 4
                    nextgrid[[x,y]]@infected_epochs = grid[[x,y]]@infected_epochs+1

                    # Cell's mutation map is copied over to the same cell in nextgrid
                    list2env(as.list.environment(grid[[x,y]]@mutations, all.names = TRUE),nextgrid[[x,y]]@mutations)
                }
                next
            }

            # Rule 1
            # In the following section  of code we evaluate the conditions for a
            # healthy cell to becombe infected. c1, c2, c3 represent those necessary
            # conditions:
            # c1: Probablity of infection has been met
            # c2: There is an infected cell in the neighbourhood
            # c3: Probabiliy of infection from outside the neighbourhood has been met
            if(grid[[x,y]]@state == Healthy){
                # Initialize all conditions as false
                c1 = FALSE
                c2 = FALSE
                c3 = FALSE

                # Initialize helper variables.
                # position: list of coordinates of infected cells in a given neighbourhood
                # negihbourhood: neighbourhood being searched
                position = list()
                neighbourhood = ''

                # Evaluate condition 2
                # Extract submatrix of states of a cell's immediate negihbourhood
                neighbourhood_1[,] = sapply(grid[c(x-1, x, x+1),
                                                 c(y-1, y, y+1)], function(x) getState(x))

                # Create a list of all infected cells from extracted submatrix
                position = which(neighbourhood_1 == 3, arr.ind = TRUE) - c(2,2)

                # If there are no infected cells in the list, set condition 2 to false,
                # otherwise to true
                c2 = ifelse(length(position) == 0, FALSE, TRUE)

                # If there are no infected cells in the immediate neighbourhood, consider
                # a negihbourhood expanded by one cell
                if(c2 == FALSE){
                    # Check neighbourhood expanded by 1
                    neighbourhood_2[,] = sapply(grid[c(x-2, x-1, x, x+1, x+2),
                                                     c(y-2, y-1, y, y+1, y+2)], function(x) getState(x))

                    # Create a list of all infected cells from extracted submatrix
                    position = which(neighbourhood_2 == 3, arr.ind = TRUE) - c(3,3)

                    # Create a list of all infected cells from extracted submatrix
                    c2 = ifelse(length(position) == 0, FALSE, TRUE)

                    # If there are no infected cells in this expanded neighbourhood, consider
                    # a negihbourhood expanded by one more cell
                    if(c2 == FALSE){
                        # Check neighbourhood expanded by 2
                        neighbourhood_3[,] = sapply(grid[c(x-3, x-2, x-1, x, x+1, x+2, x+3),
                                                         c(y-3, y-2, y-1, y, y+1, y+2, y+3)], function(x) getState(x))

                        # Create a list of all infected cells from extracted submatrix
                        position = which(neighbourhood_3 == 3, arr.ind = TRUE) - c(4,4)

                        # Create a list of all infected cells from extracted submatrix
                        c2 = ifelse(length(position) == 0, FALSE, TRUE)

                        # If there are no infected cells in this expanded neighbourhood, consider
                        # infection coming from outside the neighbourhood. At this point the first
                        # condition for infection (c1 && c2) is false because c2 is flase.
                        # Infection may only happen if c3 is true.
                        if(c2 == FALSE){
                            c3 = runif(1)<=P_v

                        # If an infected cell that can pass on the infection is found, evaluate
                        # condition 1:
                        # The probability of infection changes with the neighbourhood that the
                        # infected cells is found in. Those probabilities are expressed as
                        # Pi_n1, Pi_n2, and Pi_n3. We generate an random number Pi, and check if,
                        # given a neighbourhood type, the infection is succesful.
                        # For example if Pi_n1 = 0.6 (60% chance of infection), and Pi = 0.5,
                        # infection is succesful, if Pi = 0.7, infection is not succesful.
                        }else{c1 = runif(1) <= Pi_n3}
                    }else{c1 = runif(1) <= Pi_n2}
                }else{c1 = runif(1) <= Pi_n1}

                # Rule 1
                # If the cell is in H state and at least one of its neighbors
                # is in I state then the cell becomes I with a probability of
                # P_n1, P_n2, or P_n3 depending on which neighborhood the cell is
                # in (Conditions c1 and c2). The cell may also becomes I by
                # randomly coming in contact with a virus from outside its
                # neighborhood with a probability of P_v (Condition c3).
                if((c1 || c3){
                    nextgrid[[x,y]]@state = 3

                    # If the infecting cell is from the neighbourhood, calculate correct
                    # indices and pass mutation map to the infected cell. If infection is
                    # coming from a remote cell, pick a random infected cell from the CA.
                    if (c2){
                        # From the infected cells in the neighbourhood, select one randomly
                        # to pass on the infection to the healthy cell.
                        r = sample(1:nrow(position),1)

                        # Translate the coordinates from neighbourhood to grid reference system
                        x_c = x + position[[r,1]]
                        y_c = y + position[[r,2]]
                    }
                    #else{
                        # If infection comes from outside neighbourhood, list all infected cells,
                        # then select one randomly.
                    #    stateGrid[,] = sapply(grid[,], function(x) getState(x))
                    #    position = which(stateGrid[,] == 3, arr.ind = TRUE)
                    #    r = sample(1:nrow(position),1)

                        # Translate the coordinates from neighbourhood to grid reference system
                    #    x_c = position[[r,1]]
                    #    y_c = position[[r,2]]
                    #}
                    # When a cell gets infected, it acquires the mutations carried by the cell
                    # it is being infected by.
                    list2env(as.list.environment(grid[[x_c,y_c]]@mutations, all.names = TRUE),nextgrid[[x,y]]@mutations)

                    # We also make its resistance attribute = to that of the infecting cell.
                    nextgrid[[x,y]]@resistance = grid[[x_c,y_c]]@resistance

                    # Upon infection of a healthy cell a new set of mutations is also generated
                    # and appended to the list passed by the infecting cell.

                    # randomly select an amino acid, and a location.
                    mutation_site = sample(1:hiv_total_aa,1)
                    mutation_aa = aa[sample(1:20, 1)]

                    # Add the the new key value pair in the mutation hashmap
                    nextgrid[[x,y]]@mutations[[as.character(mutation_site)]] = mutation_aa

                    # Extract possible drug resistance conferring mutation at generated site
                    potential_resistance = c(resistanceSites_drug1[[as.character(mutation_site)]],
                                             resistanceSites_drug2[[as.character(mutation_site)]],
                                             resistanceSites_drug3[[as.character(mutation_site)]])

                    # If there are possible canditates, check that the generated amino acids
                    # appears among them.
                    if(!is.null(potential_resistance)){
                        # Count the how many of the drugs, this mutation provides resistance to 1-3
                        present = sum(potential_resistance == mutation_aa)

                        # If present = 0, the mutation changed a drug resistance conferring mutation
                        # into a non drug resistance conferring mutation. Therefore, we decrease the
                        # resistance attribute.
                        if ((present == 0) && (nextgrid[[x,y]]@resistance > 0)){
                            nextgrid[[x,y]]@resistance = nextgrid[[x,y]]@resistance-1

                        # Otherwise we increase the resistance parameter by the number of drugs affected
                        # by this muation. The resistance parameter can increase to a maximum of 3.
                        }else{
                            if(logging){
                                cat("Resistance Acquired: ", file=logFile, append=TRUE)
                                cat("(", file=logFile, append=TRUE)
                                cat(timestep, file=logFile, append=TRUE)
                                cat(",", file=logFile, append=TRUE)
                                cat(x, file=logFile, append=TRUE)
                                cat(",", file=logFile, append=TRUE)
                                cat(y, file=logFile, append=TRUE)
                                cat(") ", file=logFile, append=TRUE, sep = "\t" )
                                cat(mutation_site, file=logFile, append=TRUE)
                                cat(":", file=logFile, append=TRUE)
                                cat(mutation_aa, file=logFile, append=TRUE, sep = "\n")
                            }
                            new_resistance = nextgrid[[x,y]]@resistance+present
                            nextgrid[[x,y]]@resistance = ifelse(new_resistance > 3, 3, new_resistance)
                        }
                    }


                    if(logging){
                        cat("(", file=logFile, append=TRUE)
                        cat(timestep-1, file=logFile, append=TRUE)
                        cat(",", file=logFile, append=TRUE)
                        cat(x_c, file=logFile, append=TRUE)
                        cat(",", file=logFile, append=TRUE)
                        cat(y_c, file=logFile, append=TRUE)
                        cat(") ", file=logFile, append=TRUE, sep = "\t")

                        cat("(", file=logFile, append=TRUE)
                        cat(timestep, file=logFile, append=TRUE)
                        cat(",", file=logFile, append=TRUE)
                        cat(x, file=logFile, append=TRUE)
                        cat(",", file=logFile, append=TRUE)
                        cat(y, file=logFile, append=TRUE)
                        cat(") ", file=logFile, append=TRUE)
                        g = 1
                        for (v in ls(nextgrid[[x,y]]@mutations)) {
                          cat((ls(nextgrid[[x,y]]@mutations)[g]), file=logFile, append=TRUE)
                          cat(":", file=logFile, append=TRUE)
                          cat(nextgrid[[x,y]]@mutations[[v]], file=logFile, append=TRUE, sep = "\t")
                          cat("   ", file=logFile, append=TRUE)
                          g = g +1
                        }
                        cat(" ", file=logFile, append=TRUE, sep = "\n")
                    }
                }else{
                    # If conditions c1, c2 and c3 cannot be adequately satisfied, the cell
                    # remains in healthy state.
                    nextgrid[[x,y]]@state = 1
                }
                next
            }#End Rule 1


            # Rule 3
            # If the cell is in D state, then the cell will become H state
            # with probability P_rep
            if(grid[[x,y]]@state == Dead){
                if(runif(1) <= (P_rep)){
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

    #Update aggregate counts
    stateGrid[,] = sapply(grid[,], function(x) getState(x))
    states_count[timestep,1] = sum(stateGrid[] == 1)
    states_count[timestep,2] = sum(stateGrid[] == 3)
    states_count[timestep,3] = sum(stateGrid[] == 2)
    states_count[timestep,4] = sum(stateGrid[] == 4)

    infectedEpochsGrid[,] = sapply(grid[,], function(x) getInfected_epochs(x))
    infectedEpochs_count[timestep,0] = sum(infectedEpochsGrid[] == 0)
    infectedEpochs_count[timestep,1] = sum(infectedEpochsGrid[] == 1)
    infectedEpochs_count[timestep,2] = sum(infectedEpochsGrid[] == 2)
    infectedEpochs_count[timestep,3] = sum(infectedEpochsGrid[] >= 3)

    resistanceGrid[,] = sapply(grid[,], function(x) getResistance(x))
    resistance_count[timestep,1] = sum(resistanceGrid[] == 1)
    resistance_count[timestep,2] = sum(resistanceGrid[] == 2)
    resistance_count[timestep,3] = sum(resistanceGrid[] == 3)

    genotypes = sapply(grid[,], function(x) getMutations(x))
    genotypes_count[timestep,1] = sum(sapply(genotypes, function (x) length(x)>0) == TRUE)

    # Save plot of current state of the grid to variable
    if(timestep %in% savesteps){
        stateGrid_list[[k]] = stateGrid
        resistanceGrid_list[[k]] = resistanceGrid
        k = k + 1
    }

    # Move to the next timestep
    timestep = timestep+1

} #End while loop

############# Clean-Up: CA and Simulation parameters ###############
remove(P_v)
remove(P_rep)
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
remove(grid)
remove(nextgrid)
remove(cell)
remove(getInfected_epochs)
remove(getMutations)
remove(getResistance)
remove(getState)
remove(setState)
remove(r)

############# Write all simulation data to file  ###############
write.csv(states_count, file = paste(base_dest, l, "/states_count_raw.csv", sep = ""))
write.csv(infectedEpochs_count, file = paste(base_dest, l, "/infectedEpochs_count_raw.csv", sep = ""))
write.csv(resistance_count, file = paste(base_dest, l, "/resistance_count_raw.csv", sep = ""))
write.csv(genotypes_count, file = paste(base_dest, l, "/genotypes_count_raw.csv", sep = ""))

remove(Zero_eps)
remove(One_eps)
remove(Two_eps)
remove(Three_eps)
remove(One_r)
remove(Two_r)
remove(Three_r)
rm(list = ls(all = TRUE))
}
