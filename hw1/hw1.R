library("ggplot2")
library("latex2exp")        # latex for plots

p_1 = function(x, t){
    d = pi*(1+(x-t)^2)
    return(1/d)
}
l_1 = function(x, t){
    return(log(p_1(x,t)))
}
l_n = function(xs, t){
    r = 0
    for(x in xs){
        r = r + l_1(x, t)
    }
    return(r)
}

raw = c(-13.87, -2.53, -2.44, -2.40, -1.75, -1.34, -1.05, -0.23, -0.07, 0.27,
        1.77, 2.76, 3.29, 3.47, 3.71, 3.80, 4.24, 4.53, 43.21, 56.75)
ts = seq(from=min(raw), to=max(raw), length=100)
df = data.frame(ts, l_n(xs=raw, t=ts))
names(df) = c("theta", "loglik")

g = ggplot(df, aes(x="theta", y="loglik")) +
    geom_line() +
    ggtitle("Log Likelihood") +
    ylab(TeX("l($\\theta)")) + xlab(TeX("$\\theta"))


xis = c(-11, -1, 0, 1.4, 4.1, 4.8, 7, 8, 38) # initial values for NR
eps = 0.001 # convergence criterion
delt = 0.0001 # increment parameter for derivatives

f_pdf = function(t){
    # Log likelihood given raw data as xs
    return(l_n(xs=raw, t=t))
}
fp = function(f, t, d){
    # first derivative aprox
    f0 = f(t-d/2)
    f1 = f(t+d/2)
    return((f1-f0)/d)
}
fpp = function(f, t, d){
    # second derivative aprox
    return(fp(fp(f, t, d), t, d))
}

nr_inc = function(f, x, d){
    # Increment x using the newton rhapson method
    # print(paste("Incrementing", x))
    return(x-f(x)/fp(f, x, d))
}

find_0s = function(f, x0, e, d, inc_fxn){
    # use NR to find a zero for the function
    x1 = inc_fxn(f, x0, d)
    i = 0
    while(abs(x1-x0)>e){
        # While we haven't converged, iterate the incrementing function (which
        # must have the signature [pdf, x, delta]) e.g. nr_inc
        x0 = x1
        x1 = inc_fxn(f, x0, d)
        # print(paste("i:",i,"x0:",x0,"x1:",x1))

        # Check for algorithm failures and divergences
        if(is.na(x1)){
            return(Inf)
        }
        if(abs(x1)==Inf){
            # If the algorithm diverges, return Inf
            return(Inf)
        }
        i = i+1
    }

    # If we've converged (abs(x1-x0)<e) then return the convergent value
    return(x1)
}

f_find_0 = function(t){
    return(fp(f_pdf, t, eps))
}
xmax = c() # holds results of NR algorithm, Inf when divergent
for(xi in xis){
    xmax = c(xmax, find_0s(f_find_0, xi, eps, delt, nr_inc))
}

# print("Finished NR, starting on fI")

eps = 0.001 # convergence criterion
delt = 0.001 # increment parameter for derivatives

fI = length(raw)/2  # I(\theta) = n/2 for the Cauchy dist.
fiscor_inc = function(f, x, d){
    # Increment x using the Fischer Scoring method
    return(x+f(x)/fI)
}
xmax = c() # holds results of NR algorithm, Inf when divergent
for(xi in xis){
    # print(paste("working on", xi))
    xmax = c(xmax, find_0s(f_find_0, xi, eps, delt, fiscor_inc))
}


### PR 2
raw = c(0.52, 1.96, 2.22, 2.28, 2.28, 2.46, 2.50, 2.53, 2.54, 2.99, 3.47, 3.53,
        3.70, 3.88, 3.91, 4.04, 4.06, 4.82, 4.85, 5.46)

p_1 = function(x, t){
    # Single dimensional pdf of x with parameter t

    # Domain checking
    if(x<0 || x>2*pi){
        return(0)
    }
    if(t< -pi || t>pi){
        return(0)
    }

    return((1-cos(x-t))/(2*pi))
}

l_n = function(xs, t){
    # multidimensional log likelihood for the above pdf
    r = 0
    for(x in xs){
        r = r + log(p_1(x, t))
    }
    return(r)
}

llik = function(t){
    # log likelihood above conditioned on observing the raw data above
    return(l_n(raw, t))
}

# Plot llik
N = 200
ts = seq(from=-pi, to=pi, length=N)

df = data.frame(ts, llik(ts))
names(df) = c("theta", "loglik")

rm(g)
g = qplot(df$theta, df$loglik) +
    geom_line() +
    ggtitle(paste("Log likelihood collocated at", N, "points")) +
    ylab(TeX("l($\\theta)")) + xlab(TeX("$\\theta"))
g

## Solve for mle with MoM t0
# print("Working on pr2")
eps = 0.001
delt = 0.001
t0 = asin(mean(raw) - pi)
llikp = function(t){        # d/dt l(t)
    return(fp(llik, t, delt))
}
mle = find_0s(llikp, t0, eps, delt, nr_inc)

ts = seq(from=-pi, to=pi, length.out=N)
df = data.frame(ts, llikp(ts))
names(df) = c("theta", "dldt")


N = 200
t0s = seq(from=-pi+delt, to=pi-delt, length.out=N)
mles = c()
for(t0 in t0s){
    mle = find_0s(llikp, t0, eps, delt, nr_inc)
    # print(paste("t0",t0,"mle",mle))
    mles = c(mles, mle)
}


### PROBLEM 3
print('starting problem 3')
df = data.frame(
    x=c(0.02, 0.02, 0.06, 0.06, 0.11, 0.11, 0.22, 0.22, 0.56, 0.56, 1.10, 1.10),
    y=c(47, 76, 97, 107, 123, 139, 152, 159, 191, 201, 200, 207)
)

fit = lm(y~x, data=df)
b0 = fit$coefficients[1]
b1 = fit$coefficients[2]
t1 = 1/b0
t2 = b1*t1
# print(paste(t1, t2))

# Define loss functions and its derivatives
g = function(t, x, y){
    # print("g")
    # Multivariate OLS, requires type(x & y) = list(length=n), type(t)=list(2)
    r = y - (t[1]*x)/(t[2]+x)
    s = sum(r^2)
    # print(s)
    return(s)
}

dgdt1 = function(t, x, y){
    # print("dgdt1")
    t1 = t[1]
    t2 = t[2]
    r = 2*x*(t1*x-y*(x+t2))/(x+t2)^2
    s = sum(r)
    # print(s)
    return(s)
}
dgdt2 = function(t, x, y){
    # print("dgdt2")
    t1 = t[1]
    t2 = t[2]
    r = -2*t1*x*(t1*x-y*(t2+x))/(t2+x)^3
    s = sum(r)
    # print(s)
    return(s)
}
ddgdt1dt2 = function(t, x, y){
    # print("ddgdt1dt2")
    t1 = t[1]
    t2 = t[2]
    r = -2*x*(2*t1-y*(x+t2))/(x+t2)^3
    s = sum(r)
    # print(s)
    return(s)
}
ddgdt12 = function(t, x, y){
    # print("ddgdt12")
    t1 = t[1]
    t2 = t[2]
    r = -2*x^2/(x+t2)^2
    s = sum(r)
    # print(s)
    return(s)
}
ddgdt22 = function(t, x, y){
    # print("ddgdt22")
    t1 = t[1]
    t2 = t[2]
    r = 2*t1*x*(3*t1*x-2*y*(t2+x))/(t2+x)^4
    s = sum(r)
    # print(s)
    return(s)
}

gp = function(t, x, y){
    # print("gp")
    g1 = dgdt1(t, x, y)
    g2 = dgdt2(t, x, y)
    s = c(g1, g2)
    print(s)
    return(s)
}
gpp = function(t, x, y){
    # print("gpp")
    g11 = ddgdt12(t, x, y)
    g12 = ddgdt1dt2(t, x, y)
    g22 = ddgdt22(t, x, y)
    m = c(g11, g12, g12, g22)
    s = matrix(m, 2, 2)
    # print(s)
    return(s)
}

gpn = function(t, x, y, d){
    # numeric first derivative
    t01 = c(t[1]-d, t[2])
    t11 = c(t[1]+d, t[2])
    dt1 = (g(t11, x, y)-g(t01, x, y))/(2*d)
    # print(paste("numerical dt1:", dt1))

    t02 = c(t[1], t[2]-d)
    t12 = c(t[1], t[2]+d)
    dt2 = (g(t12, x, y)-g(t02, x, y))/(2*d)
    # print(paste("numerical dt2:", dt2))

    r = c(dt1, dt2)
    return(r)
}

gppn = function(t, x, y, d){
    # numeric second derivative
    t01 = c(t[1]-d, t[2])
    t11 = c(t[1]+d, t[2])
    dt1 = (g(t11, x, y)-g(t01, x, y))/(2*d)
    # print(paste("numerical dt1:", dt1))

    t02 = c(t[1], t[2]-d)
    t12 = c(t[1], t[2]+d)
    dt2 = (g(t12, x, y)-g(t02, x, y))/(2*d)
    # print(paste("numerical dt2:", dt2))

    r = c(dt1, dt2)
    return(r)
}

# Define iteration function
nr_itt = function(t, x, y, d){

    gp1 = gp(t, x, y)
    gp2 = gpp(t, x, y)
    gp1n = gpn(t, x, y, d)
    print(paste("analytical gp1:", gp1))
    print(paste("analytical gp2:", gp2))
    print(paste("numerical gp1:", gp1n))

    t1 = t - solve(gp2) %*% gp1
    print(paste("old t:", t))
    print(paste("new t:", t1))
    return(t1)
}

# Parameterize algorithm
eps = 0.001
del = 0.0001
x = df$x
y = df$y
tn = c(t1, t2) # MoM as t_0
tnp1 = nr_itt(tn, x, y, del)

# Gradient descent
print("Entering multidimensional gradient descent")
i = 0
while(t(tn) %*% tnp1>eps){
    i = i + 1
    print(paste("working on itt", i))
    tn = tnp1
    tnp1 = nr_itt(tn, x, y, del)
    # print(tn)
    if(i>0){break}
}

# Numerical plotting the g()
ff = function(t1, t2){
    t = c(t1,t2)
    r = g(t, df$x, df$y)
}
N = 150
rng = seq(-1,1,length=N)
pltdf = data.frame(x=rep(0, N^2), y=rep(0, N^2), z=rep(0, N^2))
i = 1
for(x in rng){
    for(y in rng){
        r = c(x, y, ff(x, y))
        pltdf[i,] = r
        i = i + 1
    }
    print(i)
}

plt = ggplot(pltdf, aes(x=x, y=y, color=z, z=z)) + geom_point() + geom_contour()
df1 = data.frame(x=t1, y=t2, z=ff(t1,t2))
plt = plt + geom_point(data=df1, aes(x=x, y=y), colour="green")
