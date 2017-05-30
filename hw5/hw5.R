# E(||e^\top e||^2 + ||(H_\lambda e)^\top H_\lambda e||^2 -
#||e(H_\lambda e)^\top||^2 - ||H_\lambda e e^\top||^2) =
rm(list=ls())

## Imports
library(ggplot2)
library(stats)
library(reshape)

## Subroutines
f = function(x){
    # prob density function defined in problem specification
    1.5*dnorm((x-0.35)/0.15) - dnorm((x-0.8)/0.04)
}

## Parameterize
n = 200
is = 1:n
xis = (is-0.5)/n

## Subroutines
noise_sd = function(j){
    # varaible noise schedule
    return(0.02+0.04*(j-1)^2)
}

## Parameterize
J = 6 # for generating noisy observations
k = 30 # number of knots
p = 3 # cubic spline
a = 0
b = 1 # bounds of inegration for lambda*int_a^b f^''(x)^2 dx penalty
knots = seq(a, b, length.out=k)

fxs = f(xis)
noise_df = data.frame(x=xis, yt=fxs)
for(j in 1:J){
    yis = fxs + rnorm(n, mean=0, sd=noise_sd(j))
    colname = paste("y", j, sep="")
    noise_df[colname] = yis
}

# Perform cubic penalized spline regression for multiple lambda values
ex_xs = noise_df$x
ex_ys = noise_df$y3 # arbitrary noise profile for example
lam_df = data.frame(x=xis, y=fxs, yn=ex_ys)
lams = 10^(-2:2)

D = diag(c(rep(0, p+1), rep(1, k)))
X1 = outer(ex_xs, 0:p, "^")
X2 = outer(ex_xs, knots, ">")*outer(ex_xs, knots, "-")^p # source [2]
X = cbind(X1, X2)
Y = ex_ys

for(lam in lams){
    H = X %*% solve(t(X) %*% X + lam*D) %*% t(X)
    fhys = H %*% Y
    colname = paste("lambda", lam, sep="=")
    lam_df[colname] = fhys
}

aicc = function(X=X, Y=Y, k=k, l){
    # k is num knots
    # l is lambda, penalty
    # X is design matrix
    # Y is response

    n = dim(X)[1]
    D = diag(c(rep(0, dim(X)[2]-k), rep(1, k)))
    H = X %*% solve(t(X) %*% X + l*D) %*% t(X)

    fhys = H %*% Y
    tr = sum(diag(H))

    r1 = log(mean((Y-fhys)^2))
    r2n = 2*(tr+1)
    r2d = n - tr - 2
    r = r1 + r2n/r2d
    return(r)
}

aiccv = Vectorize(aicc, vectorize.args = c("l"))
lambdas = 10^(seq(-10, -3, length.out=100))
aiccs = aiccv(X=X, Y=Y, l=lambdas, k=k)
