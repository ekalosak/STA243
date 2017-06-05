## Supplemental code for homework 6
# Author: Eric Kalosa-Kenyon

model1 = function(x){
    t = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
    temp = 0
    for(i in 1:11) {
        temp = temp + h[i]/2 * (1 + sign(x - t[i]))
    }
    return(temp)
}

noise1 = function(x){
 return(rnorm(length(x))/3)
}
