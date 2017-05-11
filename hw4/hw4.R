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


# ### PR 3
# ## a
# h = function(x){
#     if(x<0){return(0)}
#     else if(x>1){return(0)}
#     else{
#         return(1/(1+x))
#     }
# }
# h = Vectorize(h)
# N = 1500
# us = runif(N, 0, 1)
# i1 = mean(h(us))
# itru = log(2)

library(ggplot2)
library(latex2exp)

# xs = seq(0,1,length.out=N)
# df_3.a = data.frame(x=xs, y=h(xs))
# df_3.a.2 = data.frame(x=xs, y=cumsum((1/N)*h(xs)))
# plt_3.a = ggplot() +
#     geom_line(data=df_3.a, aes(x=x, y=y), color="coral") +
#     geom_line(data=df_3.a.2, aes(x=x, y=y), color="steelblue") +
#     # geom_hline(yintercept=i1, color="green") +
#     # geom_hline(yintercept=itru, color="blue") +
#     labs(title=latex2exp("\\frac{1}{1+x}"))

# ## b

# cc = function(x){
#     if(x<0){return(0)}
#     else if(x>1){return(0)}
#     else{
#         return(1+x)
#     }
# }
# cc = Vectorize(cc)

# b = -cov(h(xs), cc(xs)) / var(cc(xs))
# i2 = mean(h(us)) + b*(mean(cc(us)) - 1.5)

# ## Parameterize
# M = 300 # number of samples for variance estimation

# ## Simulate
# i1s = rep(NA, M)
# i2s = rep(NA, M)

# for(i in 1:M){
#     us = runif(N)
#     i1 = mean(h(us))
#     i2 = mean(h(us)) + b*(mean(cc(us)) - 1.5)
#     i1s[i] = i1
#     i2s[i] = i2
# }

# v1 = var(i1)
# v2 = var(i2)

# ### Pr 5

# ## Parameterize
# n = 100
# ltrue = 2
# ptrue = 0.3
# N = 400 # samples from the Gibbs sampler
# A = seq(1/2, 2, by=1/2)
# B = A

# ## Generate x's
# rxs = runif(n) < ptrue # r Bernoulli(p)
# xs = rxs * rpois(n, ltrue)

# ## Instantiate main data frame
# df.main = data.frame(matrix(NA,
#                             length(A)*length(B)*N, # N samples for each (a,b)
#                             (5+n) # i, p, lambda, r1, ..., rn
#                         ))
# names(df.main) = c("n", "p", "l", "a", "b", paste("r", 1:n, sep=""))

# z = 0
# for(a in A){for(b in B){

# ## Sample from prior
# r = runif(n) < ptrue # r Bernoulli(p)

# ## Initialize lambda, r, p randomly
# r0 = runif(n) > rbeta(n, a, b) # uninformative prior
# l0 = ceiling(runif(1)*5) # could envelope sample from p(l) = 1/N *l^(-1/2)
# p0 = rbeta(1, a, b)
# r1 = r0
# l1 = l0
# p1 = p0

# Rs = matrix(NA, N, n)
# ls = rep(NA, N)
# ps = rep(NA, N)
# Rs[1,] = r1
# ls[1] = l1
# ps[1] = p1

# ## Run the Gibbs sampler

# for(i in 2:N){
#     # Update according to conditional distributions
#     r1 = runif(n) < (p1*exp(-l0))/(p1*exp(-l0) + (1-p1)*(xs==0))
#     l1 = rgamma(1, shape = a + sum(xs), rate = b + sum(r1))
#     p1 = rbeta(1, 1 + sum(r1), n + 1 - sum(r1))

#     # Record those samples
#     Rs[i,] = r1
#     ls[i] = l1
#     ps[i] = p1
# }

# # Update the main dataframe (n, p, l, r1, ..., rn)
# ixs = (1:N)+z*N
# df.main[ixs,][1] = 1:n  # iteration number
# df.main[ixs,][2] = ps   # estimate of p
# df.main[ixs,][3] = ls   # estimate of lambda
# df.main[ixs,][4] = a    # prior a
# df.main[ixs,][5] = b    # prior b
# df.main[ixs,][6:dim(df.main)[2]] = as.integer(Rs)   # hidden \mathfb{r}_j
# z = z + 1

# }} # end for A for B

# ## Plot
# df = data.frame(1:n, ps, ls, Rs)
# names(df) = c("i", "p", "l", paste("r", 1:n, sep=""))

# plt_lambda = ggplot(data=df.main) +
#     geom_line(aes(x=n, y=l), color="steelblue") +
#     geom_hline(aes(yintercept=ltrue), color="coral") +
#     facet_grid(a ~ b,
#             labeller=labeller(
#                 a = (function (x) paste("a=", x, sep="")),
#                 b = (function (x) paste("b=", x, sep=""))
#             )) +
#     labs(x="Gibbs sample iteration",
#          y=TeX("$\\lambda_i$"),
#          title=TeX(paste(
#             "Gibbs sampler convergence for $\\lambda$",
#             "using prior $f(p)~\\beta(a,b)(p)$"
#             ))
#         )

# plt_rs = ggplot() +
#     geom_raster(
#                 data=melt(ceiling(Rs)),
#                 mapping=aes(X1, X2, fill=value)
#             ) +
#     labs(
#          x="Gibbs sample iteration",
#          y=TeX("$\\mathbf{r}_{j}$"),
#          title=TeX("Gibbs sampler convergence for $\\mathbf{r}$")
#         )

# plt_xs = ggplot(melt(matrix(xs)), aes(X1, value)) +
#     geom_point() +
#     labs(
#          x=TeX("Index $i$"),
#          y=TeX("$X_i$"),
#          title=TeX("True observations $\\mathbf{X}$")
#         )

# plt_rn = ggplot(melt(matrix(ceiling(Rs[N,]))), aes(X1, value)) +
#     geom_point() +
#     labs(
#          x=TeX("Index $i$"),
#          y=TeX(paste("$r_{", N, ",i}$", sep="")),
#          title=TeX(paste("True observations $\\mathbf{r_{", N,"}}$", sep=""))
#         )

### PR 6

N = 1000    # Number of samples
ran = seq(0, 10, by=0.01)   # Sampling domain
t1 = 1.5    # Parameters for f_t1,t2(z)
t2 = 2

f = function(z){
    # Target density f(x)

    # Domain enforcement
    if(z<=0){return(0)}
    else if(t1<=0){return(0)}
    else if(t2<=0){return(0)}

    # Density
    r = z^(-3/2)*exp(-t1*z-t2/z+2*sqrt(t1*t2)+log(sqrt(2*t2)))
    return(r)
}
f = Vectorize(f)

g = function(z){
    # Envelope density g(z)
    # return(dgamma(z, t1*t2, t1*t2))
    return(dgamma(z, 1, 1))
}
g = Vectorize(g)

r = function(x,y){
    return(min(f(y)/f(x)*g(x)/g(y),1))
}

dists_df = data.frame(x=ran, f=f(ran), g=g(ran))
plt_dists = ggplot(data=melt(dists_df, id="x")) +
    geom_line(aes(x=x, y=value, color=variable))

samp = rep(NA, N) # want sample of size N from x(z)
envelope = g(ran) # proposition envelope
z0 = sample(size=1, x=ran, prob=envelope) # initial sample from g(z)
samp[1] = z0
for(i in 2:N){
    z1 = sample(size=1, x=ran, prob=envelope) # proposal from the envelope
    a = r(z0,z1)
    if(a == 1){
        # accept
        z0 = z1
    }else{
        # accept with pr(a)
        if(a > runif(1)){
            # accept
            z0 = z1
        } # else, reject and do nothing this time around the loop
    }
    samp[i] = z0
}

plt_methast = ggplot() +
    geom_histogram(data=melt(samp), aes(x=value, y=..density..), alpha=0.8) +
    geom_line(data=melt(dists_df, id="x"), aes(x=x, y=value, color=variable))
