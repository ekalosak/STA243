
library("ggplot2")          # fancy plotting
library("latex2exp")        # latex for plots

df = data.frame(
    x=c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10),
    y=c(47, 76, 97, 107, 123, 139, 152, 159, 191, 201, 200, 207)
)

fit = lm(y~x, data=df)
b0 = fit$coefficients[1]
b1 = fit$coefficients[2]
t1 = 1/b0
t2 = b1*t1

# Define loss functions and its derivatives
g = function(t, x, y){
    # Multivariate OLS, requires type(x & y) = list(length=n), type(t)=list(2)
    r = y - (t[1]*x)/(t[2]+x)
    return(sum(r^2))
}

dgdt1 = function(t, x, y){
    t1 = t[1]
    t2 = t[2]
    r = 2*x*(t1*x-y*(x+t2))/(x+t2)^2
    return(sum(r))
}
dgdt2 = function(t, x, y){
    t1 = t[1]
    t2 = t[2]
    r = -2*t1*x*(t1*x-y*(t2+x))/(t2+x)^3
    return(sum(r))
}
ddgdt1dt2 = function(t, x, y){
    t1 = t[1]
    t2 = t[2]
    r = -2*x*(2*t1*x-y*(x+t2))/(x+t2)^3
    return(sum(r))
}
ddgdt12 = function(t, x, y){
    t1 = t[1]
    t2 = t[2]
    r = 2*x^2/(x+t2)^2
    return(sum(r))
}
ddgdt22 = function(t, x, y){
    t1 = t[1]
    t2 = t[2]
    r = 2*t1*x*(3*t1*x-2*y*(t2+x))/(t2+x)^4
    return(sum(r))
}

gp = function(t, x, y){
    g1 = dgdt1(t, x, y)
    g2 = dgdt2(t, x, y)
    return(c(g1, g2))
}
gpp = function(t, x, y){
    g11 = ddgdt12(t, x, y)
    g12 = ddgdt1dt2(t, x, y)
    g22 = ddgdt22(t, x, y)
    m = c(g11, g12, g12, g22)
    return(matrix(m, 2, 2))
}

# Define iteration function
nr_itt = function(t, x, y){
    gp1 = gp(t, x, y)
    gp2 = gpp(t, x, y)

    t1 = tryCatch({
        gp2inv = solve(gp2)
        r = t - gp2inv %*% gp1
        return(r)

    }, error = function(e) {
        print(e)
        r = c(Inf, Inf)
        return(r)

    }, warning = function(w) {
        print(w)
        r = c(Inf, Inf)
        return(r)
    })

    return(t1)
}

# Parameterize algorithm
# eps = 0.001 # convergence criteria for multidim. grad. desc. unclear
x = df$x
y = df$y
tn = c(t1, t2) # MoM as t_0
tnp1 = nr_itt(tn, x, y)

# Gradient descent
i = 0
maxitt = 1000
track = data.frame(
                   t1=rep(0, maxitt),
                   t2=rep(0, maxitt),
                   i=rep(0, maxitt))
track[1,] = c(tn, 1)
for(i in 1:maxitt){
    # i = i + 1
    tn = tnp1
    track[i,] = c(tn, i)
    tnp1 = nr_itt(tn, x, y)
    print(tnp1)
    if(tnp1[1] == Inf){
        tnp1 = tn
        break
    }
}

# Numerical plotting the g()
ff = function(t1, t2){
    t = c(t1,t2)
    r = g(t, df$x, df$y)
}

N = 50
rngt1 = seq(-10, 500, length=N)
rngt2 = seq(-1, 1, length=N)
pltdf = data.frame(x=rep(0, N^2), y=rep(0, N^2), z=rep(0, N^2))

i = 1
for(x in rngt1){
    for(y in rngt2){
        r = c(x, y, ff(x, y))
        pltdf[i,] = r
        i = i + 1
    }
}

pltdf$loglog = log(log(pltdf$z))
plt = ggplot(pltdf, aes(x=x, y=y, color=loglog, z=z)) +
    geom_point() + # + geom_contour()
    ylab(TeX("$\\theta_2")) + xlab(TeX("$\\theta_1")) +
    scale_colour_gradientn(colours=rainbow(4)) +
    ggtitle(TeX("Loss function g($\\theta)"))
plt2 =  ggplot(track[track$i>200,],
               aes(x=t1, y=t2, color=i)) + geom_point() +
    ggtitle("Progression of Newton-Rhapson method")

inx = which(pltdf$z == min(pltdf$z))
tgrid = c(pltdf$x[inx], pltdf$y[inx])
