### Author: Eric Kalosa-Kenyon
### STA 243 Hw2 - Simulated annealing

##### BEGIN

### Purpose of this algorithm
#   The goal of this program is to find an optimal solution for the traveling
#   salesman problem. The method used is simulated annealing.

### Pseudocode overview
# Generate starting sample s
# For k in k_0 : k_max
#   T = temp(k_0/k_max)
#   sample from neighborhood N(s) -> s'
#   if H(s')<H(s):
#       s <- s'
#   else:
#       if e^(-(H(s')-H(s))/T) > runif(0,1):
#           s <- s'
# return s
# Plot result

### BEGIN: Code
rm(list=ls())

## Imports
library("ggplot2")
library("dplyr")

## Define distance matrix
raw = c(
    0, 1, 2, 4, 9, 8, 3, 2, 1, 5, 7, 1, 2, 9, 3,
    1, 0, 5, 3, 7, 2, 5, 1, 3, 4, 6, 6, 6, 1, 9,
    2, 5, 0, 6, 1, 4, 7, 7, 1, 6, 5, 9, 1, 3, 4,
    4, 3, 6, 0, 5, 2, 1, 6, 5, 4, 2, 1, 2, 1, 3,
    9, 7, 1, 5, 0, 9, 1, 1, 2, 1, 3, 6, 8, 2, 5,
    8, 2, 4, 2, 9, 0, 3, 5, 4, 7, 8, 3, 1, 2, 5,
    3, 5, 7, 1, 1, 3, 0, 2, 6, 1, 7, 9, 5, 1, 4,
    2, 1, 7, 6, 1, 5, 2, 0, 9, 4, 2, 1, 1, 7, 8,
    1, 3, 1, 5, 2, 4, 6, 9, 0, 3, 3, 5, 1, 6, 4,
    5, 4, 6, 4, 1, 7, 1, 4, 3, 0, 9, 1, 8, 5, 2,
    7, 6, 5, 2, 3, 8, 7, 2, 3, 9, 0, 2, 1, 8, 1,
    1, 6, 9, 1, 6, 3, 9, 1, 5, 1, 2, 0, 5, 4, 3,
    2, 6, 1, 2, 8, 1, 5, 1, 1, 8, 1, 5, 0, 9, 6,
    9, 1, 3, 1, 2, 2, 1, 7, 6, 5, 8, 4, 9, 0, 7,
    3, 9, 4, 3, 5, 5, 4, 8, 4, 2, 1, 3, 6, 7, 0
    )
D = matrix(raw, sqrt(length(raw)), sqrt(length(raw)))

## Parameterize script
N = dim(D)[1] # number of cities
M = 10 # number of times to anneal
k_max = 1000 # steps in the annealing process

## Subroutine definitions
temp = function(r){
    ## Temperature schedule
    stopifnot(r > 0 || r <= 1)

    # when r is large, temperature should be low
    # temp should always be positive, real
    tt = -log(r)
    return(tt)
}

path_len = function(s, DM=D){
    r = 0
    N = length(s)
    for(i in 1:(N-1)){
        r = r + DM[s[i],s[i+1]]
    }
    return(r)
}

nhbd = function(s, k = 2){
    # Generate a neighborhood around path (s = e.g. [1,4,2,3] for N=4)
    #   which consists of all sigma(s) where sigma is a k-cycle rotation with k
    #   non-zero elements.

    # Check input
    N = length(s)
    stopifnot(N > 2)
    stopifnot(k >= 2) # 1 cycle neighborhood is identity
    stopifnot(k <= N) # there are no N+1 cycles in Z/Z_{N}

    # Generate k-cycle
    K = sample(1:N, k)

    # Rotate s by K
    s2 = s
    for(i in 1:(k-1)){
        s2[K[i]] = s[K[i+1]]
    }
    s2[K[k]] = s[K[1]]

    # return result
    return(s2)
}

## Tests
test_all_subroutines = function(){
    ## Unit tests

    # Test neighborhood generator
    for(i in 1:100){
        k = floor(runif(1)*18) + 2
        N = k + floor(runif(1)*10) + 1
        s = sample(1:N)
        n = nhbd(s, k)
        stopifnot(sum(n != s) == k)
    }
    print("nhbd(s,k) passed!")

}

## Anneal
results = data.frame(
    Step=integer(),
    Length=integer(),
    M=integer(),
    stringsAsFactors=F
)

# anneal M times to avoid local minima
for(m in 1:M){

    s0 = sample(1:N) # Starting path sampled uniformly at random

    for(k in 1:k_max){

        tt = temp(k/k_max) # Get temperature from temperature schedule
        s1 = nhbd(s0) # Sample uniformly at random from the neighborhood of s0
        h0 = path_len(s0) # Calculate energy of object
        h1 = path_len(s1)

        if(h1 < h0){ # If the new path is shorter
            s0 = s1 # take it
        }else{
            R = runif(1)
            P = exp(-(h1-h0)/tt)

            # if(is.nan(P)){
            #     print("error calculating P")
            #     next
            # }

            if(P>R){ # Metropolis-Hastings
                s0 = s1
            }
        }

        # record results
        results[dim(results)[1]+1,] = c(k, path_len(s0), m)

    }
}

## Plot result
plt = ggplot(results, aes(x=Step, y=Length, color=factor(M))) +
    geom_line() +
    ggtitle("Annealing the Traveling Salesman problem") +
    xlab("Step") + ylab("Path length") +
    labs(color="Annealing\nrun")
