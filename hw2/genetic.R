### Author: Eric Kalosa-Kenyon
### STA 243 Hw2 - Genetic algorithms

##### BEGIN

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
library("dplyr")

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
    browser()
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

loss_fxn = function(segment_observns){
    # Sum of square distance from mean
    return(
        sum((segment_observns - mean(segment_observns))^2)
    )
}

fitness = function(population, raw, complexity="AIC"){
    # Calculate fitness for each individual in the population relative to the
    #   raw data.
    #   Complexity is either "AIC" or "MDL". Raw is the raw data to fit.
    # TODO: Incorporate complexity loss as option of MDL, AIC
    # Return a vector of floats between 0 and 1.
    #   The AIC and MDL are normalized to [0,1]

    # Check that input is ok
    stopifnot(nchar(population$Chromosomes[1]) + 1 == length(raw))
    stopifnot(complexity == "AIC" || complexity == "MDL")

    # Initialize parameters
    losses = c()
    K = dim(population)[1]  # (K) is number of individuals in population
    n = length(raw)

    for(k in 1:K){ # for each individual in the population
        x = population$Chromosome[k] # the individual's chromosome
        segment_observations = c(raw[1])
        loss = 0 # initialize base loss to 0

        for(i in 1:nchar(x)){ # for each gene in the chromosome
            g = as.integer(substr(x, i, i)) # extract the gene from the string

            if(g){ # if the gene indicates a breakpoint
                # calculate and increment overall loss by segment loss
                loss = loss + loss_fxn(segment_observations)
                segment_observations = c(raw[i+1]) # reset to next segment

            }else{ # if the gene does not indicate a breakpoint
                # extend the segment
                segment_observations = c(segment_observations, raw[i+1])
            }
        }

        # Calculate loss for final segment
        loss = loss + loss_fxn(segment_observations)

        losses = c(losses, loss) # record that individual's loss
    }

    ## Penalize for complexity
    ns_params = sapply( # vector of integers summing number of breakpoints
                    lapply(
                        strsplit(population$Chromosome,""),
                        as.integer),
                    sum) + 1 # number of constants is breaks + 1

    if(complexity=="AIC"){
        # AIC is k-ln(sum sq losses)
        #   losses is vector of sum(square losses) so
        AIC = 2*log(n)*ns_params + n*log(losses/n)
        # Note that this AIC can be alternatively parameterized
        fitnesses = -AIC + max(AIC) + 0.001 # low AIC is high fitness

    }else if(complexity=="MDL"){
        # calculate sum of log(\hat{n_j})
        # sum_n_hat = for each chromosome, length of 0s between 1s plus 1
        sum_n_hat = sapply(
                        lapply( # for each chromosome, get \hat{n_j}'s
                            strsplit( # get strs of 0's (whose lengths = n_j-1)
                                population$Chromosome,
                                "1"
                            ),
                            nchar),
                        (function (x) sum(log(x+1))) # sum(log(\hat{n_j}))
                    )
        MDL = ns_params*log(n) +
            1/2*sum_n_hat +
            n/2*log(losses/n)

        fitnesses = -MDL + max(MDL) + 0.001 # low MDL is high fitness
    }

    fitnesses = fitnesses/max(fitnesses) # normalize to (0,1)
    return(fitnesses)
}

## Generate observations
# TODO: replace simple raw data with the generator given in the homework
raw = c(0,0,0,1,1,1,1,5,5,5,6,6,2,2,2,2,2,2,3,1,1,1)
truth = "001000100101000001100"
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
which_penalty = "AIC" # or "MDL", both are calculated
for(g in 2:G){
    parent_generation = pop[pop$Generation == g-1,]
    parent_aic_mdl = fitness(parent_generation, raw, "MDL")
    # Put parent fitness into main dataframe
    pop[rownames(parent_generation),]$Fitness = parent_aic_mdl
    parent_generation$Fitness = parent_aic_mdl # inelegant but effective

    for(k in 1:K){
        parents = sample_n(parent_generation,
                         2,
                         replace=FALSE,
                         weight=Fitness)
        parent_chromosomes = parents$Chromosome
        child_chromosome = mate_chromosomes(parent_chromosomes)
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
