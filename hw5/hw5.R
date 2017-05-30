# E(||e^\top e||^2 + ||(H_\lambda e)^\top H_\lambda e||^2 -
#||e(H_\lambda e)^\top||^2 - ||H_\lambda e e^\top||^2) =
# rm(list=ls())

## Imports
library(ggplot2)
library(stats)
library(reshape)

## Subroutines
f = function(x){
    # prob density function defined in problem specification
    1.5*dnorm((x-0.35)/0.15) - dnorm((x-0.8)/0.04)
}

# ## Parameterize
# n = 200
# is = 1:n
# xis = (is-0.5)/n

## Subroutines
noise_sd = function(j){
    # varaible noise schedule
    return(0.02+0.04*(j-1)^2)
}

# ## Parameterize
# J = 6 # for generating noisy observations
# k = 30 # number of knots
# p = 3 # cubic spline
# a = 0
# b = 1 # bounds of inegration for lambda*int_a^b f^''(x)^2 dx penalty
# knots = seq(a, b, length.out=k)

# fxs = f(xis)
# noise_df = data.frame(x=xis, yt=fxs)
# for(j in 1:J){
#     yis = fxs + rnorm(n, mean=0, sd=noise_sd(j))
#     colname = paste("y", j, sep="")
#     noise_df[colname] = yis
# }

# # Perform cubic penalized spline regression for multiple lambda values
# ex_xs = noise_df$x
# ex_ys = noise_df$y3 # arbitrary noise profile for example
# lam_df = data.frame(x=xis, y=fxs, yn=ex_ys)
# lams = 10^(-2:2)

# D = diag(c(rep(0, p+1), rep(1, k)))
# X1 = outer(ex_xs, 0:p, "^")
# X2 = outer(ex_xs, knots, ">")*outer(ex_xs, knots, "-")^p # source [2]
# X = cbind(X1, X2)
# Y = ex_ys

# for(lam in lams){
#     H = X %*% solve(t(X) %*% X + lam*D) %*% t(X)
#     fhys = H %*% Y
#     colname = paste("lambda", lam, sep="=")
#     lam_df[colname] = fhys
# }

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

gcv = function(X=X, Y=Y, k=k, l){
    # k is num knots
    # l is lambda, penalty
    # X is design matrix
    # Y is response

    n = dim(X)[1]
    D = diag(c(rep(0, dim(X)[2]-k), rep(1, k)))
    H = X %*% solve(t(X) %*% X + l*D) %*% t(X)

    fhys = H %*% Y
    tr = diag(H)

    numer = mean((Y-fhys)^2)
    denom = (1-mean(tr))^2
    r = numer/denom
    return(r)
}

gcvv = Vectorize(gcv, vectorize.args = c("l"))

er = function(X=X, Y=Y, k=k, l, lp){
    # k is num knots
    # l is lambda, penalty
    # X is design matrix
    # Y is response

    n = dim(X)[1]
    D = diag(c(rep(0, dim(X)[2]-k), rep(1, k)))
    H = X %*% solve(t(X) %*% X + l*D) %*% t(X)
    Hcv = X %*% solve(t(X) %*% X + lp*D) %*% t(X)

    fcvhys = Hcv %*% Y
    fhys = H %*% Y
    sig = sighat(X=X, Y=Y, k=k, l=lp)

    I = diag(dim(H)[1])
    r = (1/n)*(sum(((H-I)%*%fcvhys)^2) + sig*sum(diag(H%*%t(H))))
    return(r)
}

erv = Vectorize(er, vectorize.args = c("l"))

cv = function(l, k=k, X=X, Y=Y){
    # k is num knots
    # l is lambda, penalty
    # X is design matrix
    # Y is response

    n = dim(X)[1]
    D = diag(c(rep(0, dim(X)[2]-k), rep(1, k)))
    H = X %*% solve(t(X) %*% X + l*D) %*% t(X)
    fhys = H %*% Y
    hii = diag(H)
    r = mean(((Y-fhys)/(1-hii))^2)
    return(r)
}

cvv = Vectorize(cv, vectorize.args = c("l"))

sighat = function(Y=Y, X=X, k=k, l){
    D = diag(c(rep(0, dim(X)[2]-k), rep(1, k)))
    H = X %*% solve(t(X) %*% X + l*D) %*% t(X)
    fhys = H %*% Y
    r = mean((Y-fhys)^2)
    return(r)
}

### SIMULATE

J = 2 # number of noise regimes
M = 5 # number of simulations
p = 3 # dimension of spline
k = 30 # number of knots
n = 200 # number of observations
a = 0
b = 1
knots = seq(a, b, length.out=k)
xs = ((1:n)-0.5)/n
fxs = f(xs)
lambdas = 10^seq(-11,-3, length.out=50)

## Setup
# Data storage objects
noise_raw_df = data.frame(x=xs, yt=fxs) # observations for plotting
noise_lam_df = data.frame( # best lambda for each simulation 1:M
    cv=rep(NA, J*M), gcv=rep(NA, J*M), aic=rep(NA, J*M), er=rep(NA, J*M),
    J=rep(NA, J*M))
noise_sco_df = data.frame( # score of the corresponding best lambda
    cv=rep(NA, J*M), gcv=rep(NA, J*M), aic=rep(NA, J*M), er=rep(NA, J*M),
    J=rep(NA, J*M))

D = diag(c(rep(0, p+1), rep(1, k)))
X1 = outer(xs, 0:p, "^")
X2 = outer(xs, knots, ">")*outer(xs, knots, "-")^p # source [2]
X = cbind(X1, X2)

for(j in 1:J){

    for(i in 1:M){
        # generate observatons
        Y = fxs + rnorm(n, mean=0, sd=noise_sd(j))

        # calculate lambda using each method
        cvs = cvv(k=k, l=lambdas, X=X, Y=Y)
        gcvs = gcvv(k=k, l=lambdas, X=X, Y=Y)
        aiccs = aiccv(k=k, l=lambdas, X=X, Y=Y)

        bcv = min(cvs) # best cv score is lowest
        cv_lam = lambdas[which(cvs == bcv)]

        ers = erv(k=k, l=lambdas, X=X, Y=Y, lp=cv_lam) # need cv_lam for pilot

        bgcv = min(gcvs)
        gcv_lam = lambdas[which(gcvs == bgcv)]
        baicc = min(aiccs)
        aicc_lam = lambdas[which(aiccs == baicc)]
        ber = min(ers)
        er_lam = lambdas[which(ers == ber)]

        print(j, i)
        ix = (j-1)*M + i
        noise_lam_df$J[ix] = j
        noise_lam_df$cv[ix] = cv_lam
        noise_lam_df$gcv[ix] = gcv_lam
        noise_lam_df$aic[ix] = aicc_lam
        noise_lam_df$er[ix] = er_lam

        noise_sco_df$cv[ix] = bcv
        noise_sco_df$gcv[ix] = bgcv
        noise_sco_df$aic[ix] = baicc
        noise_sco_df$er[ix] = ber

    }

    # record an example of the observations generated in this noise regime
    colname = paste("J=", j, sep="")
    noise_raw_df[colname] = Y
}

plt_noise_obv = ggplot(data=melt(noise_raw_df, id=c("x", "yt"))) +
    # geom_line(aes(x=x, y=yt)) +
    geom_point(aes(x=x, y=value)) +
    facet_wrap(~variable) +
    labs(x="X", y="Y", title="Observations with different noise")

noise_lam_df_m = melt(noise_lam_df, id=c("J"))
noise_lam_df_m$J = as.factor(noise_lam_df_m$J)

plt_noise_lams = ggplot(data=noise_lam_df_m) +
    geom_boxplot(aes(x=J, y=value, color=J)) +
    facet_wrap(~variable)
