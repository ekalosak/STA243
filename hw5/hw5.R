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

# fit_cubic_spline = function(x, y, p=3, knots, lambda){
#     # first, create the design matrix
#     # Using: source 2 (cschwarz  stat.sfu.ca)
#     # fits simple spline of degree p
#     X1 = outer(x ,1:p, "^")
#     X2 = outer(x, knots, ">") *
#         outer(x, knots, "-")^p
#     X = cbind(X1, X2)
#     model = lm(y ~ x)
# }

## Parameterize
J = 6
k = 30
n = 200
is = 1:n
xis = (is-0.5)/n
a = 0
b = 1 # bounds of inegration for lambda*int_a^b f^''(x)^2 dx penalty
knots = seq(a, b, length.out=k)

## Generate observations
fxs = f(xis)
noise_df = data.frame(x=xis, yt=fxs)
for(j in 1:J){
    yis = fxs + rnorm(n, mean=0, sd=noise_sd(j))
    colname = paste("y", j, sep="")
    noise_df[colname] = yis
}


## Select lambda for each method
# TODO
lam0 = 0.0001
lam1 = 1.5
lam2 = 2
lam3 = 2.5

# ## Fit the spline using the true data noise_df$yt
model_cv = smooth.spline(   # regular CV
    x=xis, y=noise_df$yt, nknots=k, lambda=lam0)
# model_gcv = smooth.spline(  # generalized CV
#     x=xis, y=yis, nknots=k, lambda=lam1)
# model_el = smooth.spline(   # expected loss
#     x=xis, y=yis, nknots=k, lambda=lam2)
# model_aic = smooth.spline(   # Akaike information criterion
#     x=xis, y=yis, nknots=k, lambda=lam3)

noise_df$cvys = predict(model_cv, noise_df$x)[[2]]
plot_df = melt(noise_df, id='x')
plot_df = plot_df[plot_df$variable == "yt" | plot_df$variable == "cvys",]

## Plot
plt_pdf = ggplot(data=plot_df) +
    geom_line(aes(x=x, y=value, color=variable)) #+
    # geom_point(aes(x=x, y=y1)) + #NOTE only showing one noise schedule
    # labs(x="X", y="Y", title="Raw data and true f(x)")
