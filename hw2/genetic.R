### Author: Eric Kalosa-Kenyon
### STA 243 Hw2 - Genetic algorithms

##### BEGIN: Introduction

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

##### END: Introduction

##### BEGIN: Code
rm(list=ls())

### Imports
library("ggplot2")
library("dplyr")

### Parameterize script
K = 20  # Number of chromosomes in each generation
G = 30  # Maximum number of generations
Pcross = 0.1 # Crossover rate
Pmutate = 0.05 # Mutation rate
which_penalty = "AIC" # "AIC" or "MDL"

### BEGIN: Subroutines
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

    child = ""
    w = 1 # crossover parameter
    for(i in 1:nchar(p1)){

        # Give child parents' genes for each gene (i)
        #   from p1 if w == 1
        #   from p2 if w == 0
        if(w){
            chr = substr(p1, i, i) # get i'th character from parent 1
        }else{
            chr = substr(p2, i, i) # get i'th character from parent 2
        }

        # Mutate gene
        if(runif(1)<mr){
            chr = as.character((as.integer(chr) + 1) %% 2)
        }

        # Extend child with the i'th parental character
        child = paste(child, chr, sep="")

        # Crossover if required
        if(runif(1)<xr){
            w = (w + 1) %% 2
        }
    }
    return(child)
}

regression_loss = function(reg, raw){
    # Sum of square error betweeen regression and raw
    # len(reg) == len(raw)
    return(sum((reg-raw)^2))
}

regress_chromosome = function(chrom, raw){
    # return c(fhat(x))
    # chrom is a string of "0001110101.."
    # raw is a c(num), len(chrom) == len(raw) - 1

    # Initialize working objects
    segment_observations = c(raw[1])
    regression = c()

    for(i in 1:nchar(chrom)){ # for each gene in the chromosome
        g = as.integer(substr(chrom, i, i)) # extract the gene from the string

        if(g){ # if the gene indicates a breakpoint
            # calculate the regression for that segment

            regression = c(regression,
                           rep(mean(segment_observations),
                               length(segment_observations)
                               )
                           )
            segment_observations = c(raw[i+1]) # reset to next segment

        }else{ # if the gene does not indicate a breakpoint
            # extend the segment
            segment_observations = c(segment_observations, raw[i+1])
        }
    }

    # calculate final segment regression
    regression = c(regression,
                   rep(mean(segment_observations),
                       length(segment_observations)
                       )
                   )

    return(regression)
}

calc_AIC = function(n, k, loss){
    # AIC is k-ln(sum sq losses)
    #   n is positive integer
    #   k is num params
    #   loss is losses
    #   k and loss can be c() as long as they have the same length
    # Note that this AIC can be alternatively parameterized

    r = 2*log(n)*k + n*log(loss/n)
    return(r)
}

count_params = function(vec_of_chromosomes){
    # given a vector of chromosomes (each is "00111010111...")
    #   count the number of segments this represents i.e. the number of
    #   parameters that are implicitly contained in this chromosome.

    r = sapply( # vector of integers summing number of breakpoints
            lapply(
                strsplit(vec_of_chromosomes, ""),
                as.integer),
            sum) + 1 # number of constants is breaks + 1
    return(r)
}

fitness = function(population, raw){
    # Calculate fitness for each individual in the population relative to the
    #   raw data.
    #   Complexity is either "AIC" or "MDL". Raw is the raw data to fit.
    # Return a vector of floats between 0 and 1.
    #   The AIC and MDL are normalized to [0,1]

    # Check that input is ok
    stopifnot(nchar(population$Chromosomes[1]) + 1 == length(raw))

    # Initialize parameters
    losses = c()
    K = dim(population)[1]  # (K) is number of individuals in population
    n = length(raw) # (n) is the number of data observations

    # Calculate losses for each individual in the generation
    for(k in 1:K){
        x = population$Chromosome[k] # extract a chromosome
        regr = regress_chromosome(x, raw) # calculate the regression
        loss = regression_loss(regr, raw) # and calculate the loss
        losses = c(losses, loss) # record that individual's loss
    }

    ## Penalize for complexity
    ns_params = count_params(population$Chromosome)
    AIC = calc_AIC(n, ns_params, losses)

    # calculate sum of log(\hat{n_j}) for MDL calculation
    # sum_n_hat = for each chromosome, length of 0s between 1s plus 1
    # TODO: refactor into a separate subroutine
    #   MDL = calc_MDL(n, ns_params, losses)
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

    return(list(AIC, MDL))
}

plot_chromosome = function(chrom, raw_df, color="blue"){
    # return a ggplot object with piecewise constants plotted on the scatter

    fhat = regress_chromosome(chrom, raw_df$y)
    fhat_df = data.frame(x=raw_df$x, y=fhat)

    plt = ggplot(raw_df, aes(x=x, y=y)) +
        geom_point() +
        geom_line(data=fhat_df, color=color) +
        xlab("Index") + ylab("Value")

    return(plt)
}

plot_generation = function(gen, raw_df, wh_score="AIC"){
    # return a ggplot object with piecewise constants plotted on the scatter
    #   for each chromosome in the generation
    # Currently, colorized lines based on AIC/MDL are unimplemented

    wh_gen = gen$Generation[1]

    # prepare base plot
    plt = ggplot(raw_df, aes(x=x, y=y)) +
        geom_point() +
        xlab("Index") + ylab("Value") +
        ggtitle(paste("Generation", wh_gen, wh_score))

    # prepare regression data for plotting
    fhat_df = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(fhat_df) = c("x","y","Group","AIC","MDL")
    for(i in 1:dim(gen)[1]){
        chrom = gen$Chromosome[i]
        fhat = regress_chromosome(chrom, raw_df$y)
        fhat_df_i = data.frame(x=raw_df$x, y=fhat)
        fhat_df_i$Group = i
        fhat_df_i$AIC = as.numeric(gen$AIC[i])
        fhat_df_i$MDL = as.numeric(gen$MDL[i])
        fhat_df = rbind(fhat_df, fhat_df_i)
    }

    # add geomlines for aic or mdl depending on which was passed into this fxn
    if(wh_score=="AIC"){
        plt = plt + geom_line(data=fhat_df, aes(color=AIC, group=Group))
    }else if(wh_score=="MDL"){
        plt = plt + geom_line(data=fhat_df, aes(color=MDL, group=Group))
    }

    # plot AIC or MDL in a color gradient
    plt = plt + scale_colour_gradientn(colours=rainbow(4))

    return(plt)
}

### END: Subroutines

### BEGIN: Simulation

## Generate observations

# # Small testing dataset in this paragraph
# raw = c(0,0,0,1,1,1,1,5,5,5,6,6,2,2,2,2,2,2,3,1,1,1)
# truth = "001000100101000001100"
# sd = 0.3 # used for generating observations from (raw)
# N = length(raw)
# raw = raw + rnorm(N, sd=sd)

# True data given in homework
truefunction<-function(x){
    t <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    h <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
    temp <- 0
    for(i in 1:11) {
        temp <- temp + h[i]/2 * (1 + sign(x - t[i]))
    }
    return(temp)
}
ndata = 512
xdata = (0:(ndata-1))/ndata
ydata = truefunction(xdata)
set.seed(0401)
N = length(ydata)
raw = ydata + rnorm(N)/3

## Initialize major objects
pop_col_names = c("Chromosome", "AIC", "MDL", "Generation")
pop = data.frame(matrix(ncol = 4, nrow = 0))   # Holds all individuals
colnames(pop) = pop_col_names

## Initialize first generation
g = 1
for(i in 1:K){
    x = random_chromosome(num_genes=N-1) # returns string of 1s and 0s
    individual = c(x, NA, NA, g)
    pop[dim(pop)[1]+1,] = individual
}

## Calculate fitness and mate each generation
for(g in 2:G){
    print(paste("starting on generation", g))
    parent_generation = pop[pop$Generation == g-1,]
    parent_aic_mdl = fitness(parent_generation, raw)
    # Put parent fitness into main dataframe
    pop[rownames(parent_generation),]$AIC = parent_aic_mdl[[1]]
    pop[rownames(parent_generation),]$MDL = parent_aic_mdl[[2]]

    # Update current working object
    parent_generation$AIC = parent_aic_mdl[[1]]
    parent_generation$MDL = parent_aic_mdl[[2]]

    # invert AIC and MDL for sampling
    parent_generation$nAIC =
        -(parent_generation$AIC - min(parent_generation$AIC))
    parent_generation$nAIC =
        parent_generation$nAIC - min(parent_generation$nAIC) + 0.001
    parent_generation$nMDL =
        -(parent_generation$MDL - min(parent_generation$MDL))
    parent_generation$nMDL =
        parent_generation$nMDL - min(parent_generation$nMDL) + 0.001

    for(k in 1:K){
        # Sample parents according to AIC or MDL
        if(which_penalty == "AIC"){
            parents = sample_n(parent_generation,
                             2,
                             replace=TRUE, # do not require sexual reproduction
                             weight=nAIC)
        } else if(which_penalty == "MDL"){
            parents = sample_n(parent_generation,
                             2,
                             replace=TRUE, # require sexual reproduction
                             weight=nMDL)
        }

        parent_chromosomes = parents$Chromosome
        child_chromosome = mate_chromosomes(
                                parent_chromosomes,
                                Pcross,
                                Pmutate
                            )
        child = c(child_chromosome, NA, NA, g)
        pop[dim(pop)[1]+1,] = child
    }
}

# Calculate final AIC/MDL
last_gen = pop[pop$Generation == G,]
aic_mdl = fitness(last_gen, raw)
# Put last fitness into main dataframe
pop[rownames(last_gen),]$AIC = aic_mdl[[1]]
pop[rownames(last_gen),]$MDL = aic_mdl[[2]]

### END: Simulation

### BEGIN: Plotting

# Plot observations
raw_df = data.frame(y=raw, x=1:length(raw))
plt_raw = ggplot(raw_df, aes(x=x, y=y)) +
    geom_point() +
    xlab("Index") + ylab("Value") +
    ggtitle("Raw observations")

# # Plot true chromosome
# plt_tru = plot_chromosome(truth, raw_df, color="green") +
#     ggtitle("Ground truth")

# Plot first and last generation colored by AIC or MDL
fgen = pop[pop$Generation==1,]
lgen = pop[pop$Generation==G,]
plt_fgen = plot_generation(fgen, raw_df, which_penalty)
plt_lgen = plot_generation(lgen, raw_df, which_penalty)

# Plot best AIC or MDL
if(which_penalty=="AIC"){
    best_score = min(as.numeric(pop$AIC))
    which_organism_is_best = which(as.numeric(pop$AIC) == best_score)
}else if(which_penalty=="MDL"){
    best_score = min(as.numeric(pop$MDL))
    which_organism_is_best = which(as.numeric(pop$MDL) == best_score)
}

best_chromosome = pop$Chromosome[which_organism_is_best]
plt_best = plot_chromosome(best_chromosome, raw_df, color="green") +
    ggtitle(
        paste(
            "Chromosome ", substr(best_chromosome, 1, 12), "... ",
            "with ", which_penalty, "=", best_score,
            sep=""
        ))

### END: Plotting

##### END: Code

### NOTES:
