library("ggplot2")

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


xis = c(-11, -1, 0, 1.4, 4.1, 4.8, 7, 8, 38)
eps = 0.01

f_pdf = function(t){
    # Log likelihood given raw data as xs
    return(l_n(xs=raw, t=t))
}
fp = function(f, t, e){
    # first derivative aprox
    f0 = f(t-e/2)
    f1 = f(t+e/2)
    return((f1-f0)/e)
}
fpp = function(f, t, e){
    # second derivative aprox
    return(fp(fp(f, t, e), t, e))
}

nr_inc = function(f, x, e){
    # Increment x using the newton rhapson method
    return(x-f(x)/fp(f, x, e))
}

nr = function(f, x0, e){
    # use NR to find a zero for the function
    x1 = nr_inc(f, x0, e)
    i = 0
    while(abs(x1-x0)>e){
        print(i)
        print(x0)
        print(x1)
        x0 = x1
        x1 = nr_inc(f, x0, e)
        i = i + 1

        if(abs(x1)==Inf){
            # If the algorithm diverges, return False
            return(Inf)
        }
    }
    print(paste("returning after", i, "iterations"))
    return(x1)
}

xmax = c()

f_find_0 = function(t){
    return(fp(f_pdf, t, eps))
}
for(xi in xis){
    print(paste("working on", xi))
    xmax = c(xmax, nr(f_find_0, xi, eps))
}
