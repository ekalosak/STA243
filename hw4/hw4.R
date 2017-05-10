### PR 1.a

### BEGIN: Code

## Preamble
rm(list=ls())
set.seed(10)

## Imports
library(ggplot2, reshape)

# ## Subroutines
# f = function(x){
#     return(x^2)
# }
# p = function(x){
#     return(1)
# }
# f = Vectorize(f)
# p = Vectorize(p)

# ## Parameterize
# N = 5000 # Samples to average for the estimation of the integral of f(x)
# K = 1000 # Resolution on the interval over which to sample w.r.t p(x)
# v = 1 # int_0^1 dx = 1

# ## Simulate
# pxs = seq(0, 1, length.out=K)
# xs = sample(pxs, N, prob=p(pxs), replace=T)
# i1 = v*mean(f(xs))

# ### PR 1.c

# ## Imports
# library(reshape)

# ## Parameterize
# N = 5000 # Upper limit of samples to average for the MC integral estimation
# m = 50 # Sample N/m, 2N/m, .., N times to observe convergence rate
# K = 1000 # Resolution on the interval over which to sample w.r.t p(x)
# v = 1 # int_0^inf dp(x)
# mx = 4 # point after which p(mx) is practically 0

# ## Subroutines
# k = 2^(2/3)*gamma(4/3)
# f = function(x){
#     return(k*3/4*x^4)
# }
# p = function(x){
#     return(1/k*exp(-x^3/4))
# }
# f = Vectorize(f)
# p = Vectorize(p)

# ## Precompute the CDF for inverse sampling
# d = mx/K
# xs = seq(0, mx, by=d)
# pdf = p(xs)
# cdf = cumsum(pdf)*d

# # must calculate inverse cdf carefully so that it is uniformly sampled over 0,1
# us = seq(0, 1, length.out=K) # uniform sampling grid
# invcdf = rep(NA, K)
# for(i in 1:length(us)){
#     ix = sum(cdf<us[i]) + 1
#     invcdf[i] = cdf[ix]
# }

# ## Simulate
# i3s = rep(NA, m)
# ns = (1:m)*N/m

# for(i in 1:m){ # m is small so this is fast enough
#     xs = sample(xs, ns[i], prob=p(xs), replace=T) # sample uniformly
#     i3 = v*mean(f(xs)) # sample mean of f(x)
#     i3s[i] = i3
# }

# ## Plot
# # pdf and cdf
# pc_df = data.frame(x=xs, CDF=cdf, PDF=pdf)
# pc_df_m = melt(pc_df, id="x") # longform for plotting
# plt3 = ggplot(pc_df_m, aes(x=x, y=value, color=variable)) +
#     geom_line() +
#     labs(x="X", y="Probability", title="PDF and CDF")
# plt3

# # simulation results
# tru = 2^(4/3)*gamma(5/3)
# df_plt4 = data.frame(x = ns, y = i3s)
# plt4 = ggplot(data=df_plt4, aes(x=x, y=y)) +
#     geom_line() +
#     geom_hline(yintercept = tru) + # true evaluation
#     labs(
#          x="Samples",
#          y="MC Estimation",
#          title=latex2exp(paste(
#             "MC Estimation of",
#             "\\int_0^{\\infty} \\frac{3}{4} x^{4} e^{-x^3/4}"
#             ))
#          )


# ### PR2
# g = function(n, u, sd){
#     # pdf of importance distribution
#     return(dnorm(n, mean=u, sd=sd))
# }
# p = function(x){
#     # pdf of nominal distribution
#     if(x<1){return(0)}
#     else if(x>2){return(0)}
#     else{return(1)}
# }
# p = Vectorize(p)
# f = function(x){
#     # pdf of numerator of likelihood ration
#     r = exp(-x^2/2)/sqrt(2*pi)
#     return(r)
# }

# m = 5000 # number of samples
# k = 10000 # resolution of sample domain
# u = 1.5 # mean of importance sampling
# vs = c(0.1, 1, 10) # standard deviations

# samples = list()
# summands = list()
# i = 1
# for(v in vs){
#     d = c(-1, 1)*(4*v) + u         # sample domain
#     domain = seq(d[1], d[2], length.out=k)
#     ps = g(domain, u=u, sd=v)
#     xs = sample(domain, m, prob=ps, replace=T) # sample from g(x)
#     samples[[i]] = xs
#     summands[[i]] = (f(xs)*p(xs))/(g(xs, u=u, sd=v))
#     i = i + 1
# }

# library(ggplot2, reshape)
# df_plt = data.frame(
#                     sum1 = summands[[1]],
#                     sum2 = summands[[2]],
#                     sum3 = summands[[3]]
#                 )
# names(df_plt) = c(
#                   paste("v =", vs[1]),
#                   paste("v =", vs[2]),
#                   paste("v =", vs[3])
#                 )
# df_plt_m = melt(df_plt)
# plt = ggplot(df_plt_m) +
#     geom_histogram(aes(x=value, color=variable), position="dodge")


### PR 3
## a
h = function(x){
    if(x<0){return(0)}
    else if(x>1){return(0)}
    else{
        return(1/(1+x))
    }
}
h = Vectorize(h)
N = 1500
us = runif(N, 0, 1)
i1 = mean(h(us))
itru = log(2)

library(ggplot2)
library(latex2exp)

xs = seq(0,1,length.out=N)
df_3.a = data.frame(x=xs, y=h(xs))
df_3.a.2 = data.frame(x=xs, y=cumsum((1/N)*h(xs)))
plt_3.a = ggplot() +
    geom_line(data=df_3.a, aes(x=x, y=y), color="coral") +
    geom_line(data=df_3.a.2, aes(x=x, y=y), color="steelblue") +
    # geom_hline(yintercept=i1, color="green") +
    # geom_hline(yintercept=itru, color="blue") +
    labs(title=latex2exp("\\frac{1}{1+x}"))

## b

cc = function(x){
    if(x<0){return(0)}
    else if(x>1){return(0)}
    else{
        return(1+x)
    }
}
cc = Vectorize(cc)

b = -cov(h(xs), cc(xs)) / var(cc(xs))
i2 = mean(h(us)) + b*(mean(cc(us)) - 1.5)

## Parameterize
M = 300 # number of samples for variance estimation

## Simulate
i1s = rep(NA, M)
i2s = rep(NA, M)

for(i in 1:M){
    us = runif(N)
    i1 = mean(h(us))
    i2 = mean(h(us)) + b*(mean(cc(us)) - 1.5)
    i1s[i] = i1
    i2s[i] = i2
}

v1 = var(i1)
v2 = var(i2)
