# Author: Eric Kalosa-Kenyon
# Date: April 23, 2017
#
# Expectation maximization

x = c(125, 18, 19, 35) # Observations
C = 5
l = function(t){ # Log likelihood
    r = x[1]*log(2+t) + (x[2]+x[3])*log(1-t) + x[4]*log(t) + C
    return(r)
}
ddt = function(f, t, e=0.00000001){
    f1 = f(t-e)
    f2 = f(t)
    f3 = f(t+e)
    d1 = f2-f1
    d2 = f3-f2
    r = (d1+d2)/2
    return(r)
}
dldt = function(t){
    r = ddt(l, t)
    return(r)
}
d2ldt2 = function(t){
    r = ddt(dldt, t)
    return(r)
}

# Newton-Rhapson the 0's of the derivative of the log likelihood
tol = 0.00001
t0 = 0.5
t1 = 0.55

ts = seq(0.01,0.99,length.out=100)
# plot(ts, l(ts))

# i = 1
# while(abs(t0-t1)>tol){
#     t0 = t1
#     t1 = t1 - dldt(t1) / d2ldt2(t1)
#     print(t1)
#     i = i + 1
# }


### PR 3


## Parameterize script
K = 5000 # Number of samples
NN = (1 - exp(-2)) # Normalizing constant

fx = function(x){
    # pdf for random variable X
    if(x<0){return(0)}
    if(x>2){return(0)}

    r = exp(-x) / NN
    return(r)
}

Fx = function(x){
    # cdf for random variable X
    if(x<0){return(0)}
    if(x>2){return(0)}

    r = (1 - exp(-x)) / NN
    return(r)
}

Fu = function(u){
    # Inverse cdf for the random variable X
    if(u<0){return(NaN)}
    if(u>1){return(NaN)}

    r = -log(1 - u*(1-exp(-2)))
    return(r)
}

## Simulate
Fuv = Vectorize(Fu) # Make Fu take vector arguments cleanly
Fxv = Vectorize(Fx)
us = runif(K, 0, 1) # K samples from U(0,1)
inv_xs = Fuv(us) # Inverse sampling method

## Plot results
# Simulation results histogram with true distribution overlay
plt_xs = seq(0, 2, length.out = 100) # xs for true fx(x) pdf
fxv = Vectorize(fx)
plt_ys = fxv(plt_xs) # true fx(x) = ys

df_tru = data.frame(x_tru = plt_xs, y_tru = plt_ys)
df_sim = data.frame(x_sim = inv_xs)

plt1 = ggplot() +
    geom_histogram( # Filled in density for the simulated data
                data = df_sim,
                aes(x = x_sim, y = ..density..),
                binwidth = 0.05,
                color="royalblue"
            ) +
    geom_line(
                data = df_tru,
                aes(x = x_tru, y = y_tru),
                color="coral"
            )

### PR 4

## Parameterize
K = 500
NN = 0.6214496 # normalizing constant for fx(x) the pdf
M1 = 1/NN # Resampling constant ensuring fx(x) <= M*g1(x) for all x
M2 = pi/NN

### BEGIN: Subroutines

fx = Vectorize(function(x){
    # pdf for problem 4
    if(x<0){return(0)}
    r = exp(-x) / (1+x^2) / NN
    return(r)
})

integrate = function(ff, xs, d=0.01){
    # Reimann integrate function (ff) from 0 to each x in (xs), delta=(d)
    rs = rep(NA, length(xs))
    n = 1
    for(x in xs){
        y = 0
        r = 0
        for(i in 1:(x/d)){
            r = r + ff(y)*d
            y = y + d
        }
        rs[n] = r
        n = n + 1
    }
    return(rs) # rs is vector of doubles same length as xs
}

g1 = Vectorize(function(x){
    # First envelope density for rejection sampling
    if(x<0){return(0)}
    r = exp(-x)
    return(r)
})

g2 = Vectorize(function(x){
    # Second envelope density for rejection sampling
    if(x<0){return(0)}
    r = 2/(pi*(1+x^2))
    return(r)
})

### END: Subroutines

## Simulate
# Envelope sampling
# First g1:
# sample y from g1, u from U(0,1)
# if u < f(y)/(M_1*g_1(y), accept
N = 200 # resolution of the discretized x axis
xdisc = seq(0,5,length.out=N) # discretization of the x axis
g1i = 1 # how many loops this envelope rejection sampling takes
i = 1 # list counter
g1rs = rep(NA, K) # initialize empty result vector
while(i <= K){
    x = sample(xdisc, 1, prob=g1(xdisc)) # proposition from envelope
    u = runif(1)
    w = fx(x)/(M1*g1(x))
    if(w>u){
        # Accept sample
        g1rs[i] = x
        i = i + 1
    }
    g1i = g1i + 1
}

# Now g2
g2i = 1 # how many loops this envelope rejection sampling takes
i = 1 # list counter
g2rs = rep(NA, K) # initialize empty result vector
while(i <= K){
    x = sample(xdisc, 1, prob=g2(xdisc)) # proposition from envelope
    u = runif(1)
    w = fx(x)/(M2*g2(x))
    if(w>u){
        # Accept sample
        g2rs[i] = x
        i = i + 1
    }
    g2i = g2i + 1
}

## Plot
# simulation results
true_df = data.frame(x=xdisc,
                     g1=g1(xdisc),
                     g2=g2(xdisc),
                     f=fx(xdisc)
                     )
true_df_m = melt(true_df, id="x")
names(true_df_m) = c("x", "Distribution", "y")

gsim_df = data.frame(g1=g1rs, g2=g2rs)

plt5 = ggplot() +
    geom_line(data=true_df_m,
              aes(x=x, y=y, color=Distribution)) +
    labs(title="Target distribution and envelopes",
         x="X", y="Probability")

plt6 = ggplot() + # g1 simulation results histogram
    geom_histogram(data=gsim_df,
                   mapping=aes(x=g1,y=..density..),
                   binwidth=0.08
                   ) +
    geom_line(data=true_df_m[true_df_m$Distribution != "g2",],
              aes(x=x, y=y, color=Distribution))

# fx(x)
pltx = seq(0,5,length.out=100)
pltfx = fx(pltx)
pltFx = integrate(fx, pltx, d=0.05)
df_plt2 = melt(
                data.frame(x = pltx, fx = pltfx, Fx = pltFx),
                id = "x"
            )
plt2 = ggplot(data = df_plt2) +
    geom_line(aes(x = x, y = value, color = variable))

# g1(x)
pltg1x = g1(pltx)
pltG1x = integrate(g1, pltx, d=0.05)
df_plt3 = melt(
                data.frame(x = pltx, g1 = pltg1x, G1 = pltG1x),
                id = "x"
            )
plt3 = ggplot(data = df_plt3) +
    geom_line(aes(x = x, y = value, color = variable))

# g2(x)
pltg2x = g2(pltx)
pltG2x = integrate(g2, pltx, d=0.05)
df_plt4 = melt(
                data.frame(x = pltx, g2 = pltg2x, G2 = pltG2x),
                id = "x"
            )
plt4 = ggplot(data = df_plt4) +
    geom_line(aes(x = x, y = value, color = variable))

# ### PR 5
# ## Parameterize
# a = 1/2
# NN = 1/a # TODO: actually derive this

# ## Subroutines
# fxy = Vectorize(function(x,y){
#     if(x<0){return(0)}
#     if(y<0){return(0)}
#     if(x^2+y^2>=1){return(0)}

#     r = x^a*y/NN
#     return(r)
# })

# ## Plot density
# ns = 100
# xys = matrix(NA,ns,ns)
# rownames(xys) = ((1:ns)/ns)
# colnames(xys) = ((1:ns)/ns)
# df_plt5 = melt(xys)
# names(df_plt5) = c("x","y","z")
# df_plt5$z = fxy(df_plt5$x, df_plt5$y)
# plt5 = ggplot(df_plt5, aes(x=x,y=y)) + # Plots the fxy
#     # geom_contour(aes(z=z)) +
#     geom_raster(aes(fill=z))
