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
#       if e^(-|H(s)-H(s')|/T) > runif(0,1):
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

## Subroutine definitions
temp = function(r){
    # Temperature schedule
    stopifnot(r >= 0 || r <= 1)
    exp(1-r)
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

## Parameterize script
N = dim(D)[1]
k_max = 1000

## Anneal
s0 = sample(1:N) # Starting path sampled uniformly at random
for(k in 1:k_max){
    tt = temp(k/k_max)
    s1 = nhbd(s0)
    h0 = path_len(s0)
    h1 = path_len(s1)

}
