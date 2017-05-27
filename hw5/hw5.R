## Imports
library(ggplot2)
library(stats)

### NOISE LEVEL SIMULATION

## Subroutines
f = function(x){
    # prob density function defined in problem specification
    1.5*dnorm((x-0.35)/0.15) - dnorm((x-0.8)/0.04)
}

noise_sd = function(j){
    return(0.02+0.04*(j-1)^2)
}

fhl = function(x, y, f, l){
    n = length(x)
    rss = 1/n*(sum((y-f(x))^2))
    penalty = l*
}

## Parameterize
J = 6
k = 30
n = 200
is = 1:n
xis = (is-0.5)/n

## Generate observations
fxs = f(xis)
noise_df = data.frame(x=xis, yt=fxs)
for(j in 1:J){
    yis = fxs + rnorm(n, mean=0, sd=noise_sd(j))
    colname = paste("y", j, sep="")
    noise_df$colname = yis
}

## Select lambda for each method
# TODO
lam0 = 1
lam1 = 1.5
lam2 = 2
lam3 = 2.5

## Fit the spline
model_cv = smooth.spline(   # regular CV
    x=xis, y=yis, nknots=k, lambda=lam0)
model_gcv = smooth.spline(  # generalized CV
    x=xis, y=yis, nknots=k, lambda=lam1)
model_el = smooth.spline(   # expected loss
    x=xis, y=yis, nknots=k, lambda=lam2)
model_aic = smooth.spline(   # Akaike information criterion
    x=xis, y=yis, nknots=k, lambda=lam3)

## Plot
plt_pdf = ggplot(data=noise_df) +
    geom_line(aes(x=x, y=yt)) +
    geom_point(aes(x=x, y=y1)) + #NOTE only showing one noise schedule
    labs(x="X", y="Y", title="Raw data and true f(x)")
