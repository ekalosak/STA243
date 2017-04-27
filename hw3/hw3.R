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
NN = 0.6214496 # normalizing constant for fx(x) the pdf
fx = Vectorize(function(x){
    # pdf for problem 4
    if(x<0){return(0)}
    r = exp(-x) / (1+x^2) / NN
    return(r)
})

Fx = Vectorize(function(x, d=0.01){
    if(x<0){return(0)}
    r = 0
    y = 0
    for(i in 1:(x/d)){
        # Reimann integration on fx(x) because I'm lazy like that
        r = r + fx(y)*d
        y = y + d
    }
    return(r)
})

M1 = 1/NN # Resampling constant ensuring fx(x) <= M*g1(x) for all x
g1 = Vectorize(function(x){
    return(NaN)
})

## Parameterize
K = 5000

## Simulate

## Plot
# fx(x)
pltx = seq(0,5,length.out=100)
pltfx = fx(pltx)
pltFx = Fx(pltx)
df_plt2 = melt(
                data.frame(x = pltx, fx = pltfx, Fx = pltFx),
                id = "x"
            )
plt2 = ggplot(data = df_plt2) +
    geom_line(aes(x = x, y = value, color = variable))
