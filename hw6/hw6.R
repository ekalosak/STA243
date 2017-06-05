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

