### PR 1.a

### BEGIN: Code

## Preamble
rm(list=ls())
set.seed(10)

## Imports
library(ggplot2)

## Subroutines
f = function(x){
    return(x^2)
}
p = function(x){
    return(1)
}
f = Vectorize(f)
p = Vectorize(p)

## Parameterize
N = 5000 # Samples to average for the estimation of the integral of f(x)
K = 1000 # Resolution on the interval over which to sample w.r.t p(x)
v = 1 # int_0^1 dx = 1

## Simulate
pxs = seq(0, 1, length.out=K)
xs = sample(pxs, N, prob=p(pxs), replace=T)
i1 = v*mean(f(xs))

### PR 1.c

## Imports
library(reshape)

## Parameterize
N = 5000 # Upper limit of samples to average for the MC integral estimation
m = 50 # Sample N/m, 2N/m, .., N times to observe convergence rate
K = 1000 # Resolution on the interval over which to sample w.r.t p(x)
v = 1 # int_0^inf dp(x)
mx = 4 # point after which p(mx) is practically 0

## Subroutines
k = 2^(2/3)*gamma(4/3)
f = function(x){
    return(k*3/4*x^4)
}
p = function(x){
    return(1/k*exp(-x^3/4))
}
f = Vectorize(f)
p = Vectorize(p)

## Precompute the CDF for inverse sampling
d = mx/K
xs = seq(0, mx, by=d)
pdf = p(xs)
cdf = cumsum(pdf)*d

# must calculate inverse cdf carefully so that it is uniformly sampled over 0,1
us = seq(0, 1, length.out=K) # uniform sampling grid
invcdf = rep(NA, K)
for(i in 1:length(us)){
    ix = sum(cdf<us[i]) + 1
    invcdf[i] = cdf[ix]
}

## Simulate
i3s = rep(NA, m)
ns = (1:m)*N/m

for(i in 1:m){ # m is small so this is fast enough
    xs = sample(xs, ns[i], prob=p(xs), replace=T) # sample uniformly
    i3 = v*mean(f(xs)) # sample mean of f(x)
    i3s[i] = i3
}

## Plot
# pdf and cdf
pc_df = data.frame(x=xs, CDF=cdf, PDF=pdf)
pc_df_m = melt(pc_df, id="x") # longform for plotting
plt3 = ggplot(pc_df_m, aes(x=x, y=value, color=variable)) +
    geom_line() +
    labs(x="X", y="Probability", title="PDF and CDF")
plt3

# simulation results
tru = 2^(4/3)*gamma(5/3)
df_plt4 = data.frame(x = ns, y = i3s)
plt4 = ggplot(data=df_plt4, aes(x=x, y=y)) +
    geom_line() +
    geom_hline(yintercept = tru) + # true evaluation
    labs(
         x="Samples",
         y="MC Estimation",
         title=latex2exp(paste(
            "MC Estimation of",
            "\\int_0^{\\infty} \\frac{3}{4} x^{4} e^{-x^3/4}"
            ))
         )
