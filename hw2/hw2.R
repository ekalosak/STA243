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
pop_col_names = c("Chromosome", "Loss", "Generation")
pop = data.frame(matrix(ncol = 3, nrow = 0))   # Holds all individuals
colnames(pop) = pop_col_names
K = 50  # Number of chromosomes in each generation
G = 10  # Maximum number of generations

## Initialize first generation
g = 1
for(i in 1:K){
    x = random_chromosome(num_genes=N-1) # returns string of 1s and 0s
    individual = c(x, NA, g)
    pop[dim(pop)[1]+1,] = individual
}

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

# enc_chromosome = function(x){
#     # Encode the string of 1's and 0's into an integer
#     return(strtoi(x, base=2))
# }

# dec_chromosome = function(x){
#     # Decode the integer into a string of 1's and 0's
#     return(paste(as.integer(intToBits(x)), collapse=''))
# }

