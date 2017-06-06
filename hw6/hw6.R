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

G = 10 # generations
K = 10 # individuals per generation
Pcross = 0.9 # Crossover rate
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

best_organism = get_best_org(pop, which_penalty)
best_score = best_organism$AIC # NOTE: use an if/else block here
best_chromosome = best_organism$Chromosome

plt_best = plot_chromosome(best_chromosome, genetic_df, color="coral") +
    ggtitle(
        paste(
            "Chromosome ", substr(best_chromosome, 1, 12), "... ",
            "with ", which_penalty, "=", best_score,
            sep=""
        ))

## Resid boot
# get residuals
# for each boot
#   sample n residuals w repl
#   add these to regressed y
#   use this new sample to compute another regression
#   record regressed values at each point
#   these serve as the emperical distribution for each point

# B = 50
# y1hats = regress_chromosome(best_chromosome, y1s)
# res1 = y1s - y1hats
# point_ests1 = data.frame(matrix(nrow=B, ncol=n))
# for(b in 1:B){
#     print(paste(b, "bootstrap"))
#     boot_resid1 = base::sample(res1, n)
#     boot_ys1 = y1hats + boot_resid1
#     boot_pop = genetic_piecewise_regression(boot_ys1,
#         G, K, Pcross, Pmutate, Pimmigrate, Psex, which_penalty)
#     boot_best_org = get_best_org(boot_pop, which_penalty)
#     boot_best_chrom = boot_best_org$Chromosome
#     boot_y1hats = regress_chromosome(boot_best_chrom, boot_ys1)
#     point_ests1[b,] = boot_y1hats
# }
# # Plot a band around the original bootstrap using the outer quantiles of the
# # point_ests1
# boot_quantiles1 = sapply(
#         point_ests1,
#         (function(x) quantile(x, probs=c(0.025, 0.975)))
#     )

# df_95c_resid1 = data.frame(cbind(t(boot_quantiles1), y1hats, x1s))
# colnames(df_95c_resid1) = c("low", "up", "y", "x")
# plt_boot_resid1 = ggplot(data=df_95c_resid1) +
#     geom_ribbon(aes(x=x, ymin=low, ymax=up))

## Bootstrapping pairs
# for each boot
#   resample x1s y1s
#   fit the model
#   yhat = regress against x1s
#   use these yhat as pointwise quantil

B = 10
y1hats = regress_chromosome(best_chromosome, y1s)
point_ests1 = data.frame(matrix(nrow=B, ncol=n))
for(b in 1:B){
    print(paste(b, "bootstrap"))
    boot_ixs = sample(x=1:n, size=n, replace=T)
    boot_ys1 = y1s[boot_ixs]
    boot_pop = genetic_piecewise_regression(boot_ys1,
        G, K, Pcross, Pmutate, Pimmigrate, Psex, which_penalty)
    boot_best_org = get_best_org(boot_pop, which_penalty)
    boot_best_chrom = boot_best_org$Chromosome
    boot_y1hats = regress_chromosome(boot_best_chrom, boot_ys1)
    point_ests1[b,] = boot_y1hats
}
# Plot a band around the original bootstrap using the outer quantiles of the
# point_ests1
boot_quantiles1 = sapply(
        point_ests1,
        (function(x) quantile(x, probs=c(0.025, 0.975)))
    )
