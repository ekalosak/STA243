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
plot(ts, l(ts))

# i = 1
# while(abs(t0-t1)>tol){
#     t0 = t1
#     t1 = t1 - dldt(t1) / d2ldt2(t1)
#     print(t1)
#     i = i + 1
# }

