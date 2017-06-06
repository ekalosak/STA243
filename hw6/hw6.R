## Imports
library(ggplot2)
library(stats)
library(reshape)

## Parameterize
n = 50
t = 3
B = 5000

## Generate dataset
xs = runif(n, 0, t)

## Generate bootstrap samples and calculate var(\hat{\theta})

# generate bootstrap
bs_raw = sample(xs, size=n*B, replace=TRUE)
bs = matrix(data=bs_raw, nrow=B, ncol=n)

# calculate var

### PR 2
rm(list=ls())

library("ggplot2")
library("dplyr")
source("genetic_supplement.R")

n = 512
x1s = (0:(n-1))/n
fx1s = model1(x1s)
y1s = fx1s + noise1(x1s)

G = 10 # generations
K = 10 # individuals per generation
Pcross = 0.9 # Crossover rate
Pmutate = 0.10 # Mutation rate
which_penalty = "MDL" # "AIC" or "MDL"

pop_col_names = c("Chromosome", "AIC", "MDL", "Generation")
pop = data.frame(matrix(ncol = 4, nrow = G*K))   # Holds all individuals
colnames(pop) = pop_col_names

g = 1
for(i in 1:K){
    x = random_chromosome(num_genes=n-1) # returns string of 1s and 0s
    pop$Generation[i] = g
    pop$Chromosome[i] = x
}

i = K + 1
for(g in 2:G){
    print(paste("starting on generation", g))
    parent_generation = pop[pop$Generation == g-1 & !is.na(pop$Generation),]
    parent_aic_mdl = fitness(parent_generation, y1s)
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
                             replace=TRUE,
                             weight=nMDL)
        }

        parent_chromosomes = parents$Chromosome
        child_chromosome = mate_chromosomes(
                                parent_chromosomes,
                                Pcross,
                                Pmutate
                            )

        pop$Generation[i] = g
        pop$Chromosome[i] = child_chromosome
        i = i + 1
    }
}

# Calculate final AIC/MDL
last_gen = pop[pop$Generation == G,]
aic_mdl = fitness(last_gen, y1s)
# Put last fitness into main dataframe
pop[rownames(last_gen),]$AIC = aic_mdl[[1]]
pop[rownames(last_gen),]$MDL = aic_mdl[[2]]

## Plot some results
genetic_df = data.frame(y=y1s, x=1:n)
fgen = pop[pop$Generation==1,]
lgen = pop[pop$Generation==G,]
plt_fgen = plot_generation(fgen, genetic_df, which_penalty)
plt_lgen = plot_generation(lgen, genetic_df, which_penalty)

if(which_penalty=="AIC"){
    best_score = min(as.numeric(pop$AIC))
    which_organism_is_best = which(as.numeric(pop$AIC) == best_score)
}else if(which_penalty=="MDL"){
    best_score = min(as.numeric(pop$MDL))
    which_organism_is_best = which(as.numeric(pop$MDL) == best_score)
}
best_chromosome = pop$Chromosome[which_organism_is_best]
plt_best = plot_chromosome(best_chromosome, genetic_df, color="green") +
    ggtitle(
        paste(
            "Chromosome ", substr(best_chromosome, 1, 12), "... ",
            "with ", which_penalty, "=", best_score,
            sep=""
        ))
