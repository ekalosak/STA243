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
source("genetic_supplement.R")

n = 512
x1s = (0:(n-1))/n
fx1s = model1(x1s)
y1s = fx1s + noise1(x1s)

G = 40 # generations
K = 20 # individuals per generation
Pcross = 0.8 # Crossover rate
Pmutate = 0.05 # Mutation rate
Pimmigrate = 0.2
Psex = 0.7
which_penalty = "AIC" # "AIC" or "MDL"

pop = genetic_piecewise_regression(y1s,
    G, K, Pcross, Pmutate, Pimmigrate, Psex, which_penalty)

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
plt_best = plot_chromosome(best_chromosome, genetic_df, color="coral") +
    ggtitle(
        paste(
            "Chromosome ", substr(best_chromosome, 1, 12), "... ",
            "with ", which_penalty, "=", best_score,
            sep=""
        ))
