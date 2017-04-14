### Author: Eric Kalosa-Kenyon
### STA 243 Hw2 - Genetic algorithms and simulated annealing

##### BEGIN: Genetic algorithm

### Purpose of this algorithm
#   The goal of this program is to find an optimal solution for a piecewise
#   constant regression problem. The method used is genetic algorithms.
#   Chromosomes are binary representations of the presense or absense of
#   discontinuities in the fitted piecewise constant regression.

### Pseudocode overview
# Generate observations used for regression
# Parameterize algorithm
# Initialize random population
# Until convergence
#   Calculate fitness in current generation
#   Mate individuals proportional to fitness
# Plot result

### BEGIN: Code
rm(list=ls())

## Imports
library("ggplot2")

## Subroutines
random_chromosome = function(num_genes, nucleobases=c(0,1)){
    # Generates a chromosome
    #   lenght of chromosome is (num_genes)
    #   each element of the chromosome is an element of (nucleobases)
    x = sample(nucleobases, num_genes, replace=TRUE)
    x = paste(x, collapse='')
    return(x)
}

mate_chromosomes = function(parents, xr=0.1, mr=0.1){
    # Produce a single offspring from parent chromosomes (p1, p2) with a
    # crossover rate (xr) and a mutation rate (mr)

    p1 = parents[1]
    p2 = parents[2]
    stopifnot(nchar(p1) == nchar(p2))

    child = paste(rep("0", nchar(p1)), collapse='') # initialize child "000..."

    w = 1 # crossover parameter
    for(i in 1:nchar(p1)){

        # Give child parents' genes for each gene (i)
        #   from p1 if w == 1
        #   from p2 if w == 0
        if(w){
            child[i] = p1[i]
        }else{
            child[i] = p2[i]
        }

        # Mutate gene
        if(runif(1)<mr){
            child[i] = (child[i] + 1) %% 2
        }

        # Crossover if required
        if(runif(1)<xr){
            w = (w + 1) %% 2
        }
    }
    return(child)
}

fitness = function(population, raw){
    #TODO:
}

## Generate observations
# TODO: replace simple raw data with the generator given in the homework
raw = c(0,0,0,1,1,1,1,5,5,5,6,6,2,2,2,2,2,2,3,1,1,1)
sd = 0.3
N = length(raw)
raw = raw + rnorm(N, sd=sd)

# Plot the generated data
raw_df = data.frame(y=raw, x=1:length(raw))
plt1 = ggplot(raw_df, aes(x=x, y=y)) +
    geom_point() +
    xlab("Index") + ylab("Value") +
    ggtitle("Raw observations")

## Parameterize algoritm and initialize major objects
pop_col_names = c("Chromosome", "Fitness", "Generation")
pop = data.frame(matrix(ncol = 3, nrow = 0))   # Holds all individuals
colnames(pop) = pop_col_names
K = 20  # Number of chromosomes in each generation
G = 70  # Maximum number of generations

## Initialize first generation
g = 1
for(i in 1:K){
    x = random_chromosome(num_genes=N-1) # returns string of 1s and 0s
    individual = c(x, NA, g)
    pop[dim(pop)[1]+1,] = individual
}

## Calculate fitness and mate each generation
for(g in 2:G){
    parent_generation = pop[pop$Generation == g-1,]
    parent_fitness = fitness(parent_generation, raw)
    # Put parent fitness into main dataframe
    pop[rownames(parent_generation),]$Fitness = parent_fitness

    for(k in 1:K){
        parents = sample(parent_generation,
                         2,
                         replace=TRUE,
                         prob=parent_fitness)
        child_chromosome = mate_chromosomes(parents)
        child = c(child_chromosome, NA, g)
        pop[dim(pop)[1]+1,] = child
    }
}
#TODO: plot max_fitness as a function of generation

## Recover and plot max fitness chromosome
max_fitness_chromosome = pop$Chromosome[pop$Fitness == max(pop$Fitness)]
#TODO: ensure there is only one max_fitness_chromosome
#TODO: plot this chromosome with the raw data

##### END: Genetic algorithm

### NOTES:

# stopifnot(a==b) # R's version of assert

# # Given function in hw handout, used for evaluating fitness
# truefunction<-function(x){
#     t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
#     h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
#     temp <- 0
#     for(i in 1:11) {
#         temp <- temp + h[i]/2 * (1 + sign(x - t[i]))
#     }
#     return(temp)
# }
# n<-512
# x<-(0:(n-1))/n
# f<-truefunction(x)
# set.seed(0401)
# y<-f+rnorm(f)/3
# plot(x,y)
# lines(x,f)
