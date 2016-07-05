# create random mixtures of communities

library(compositions)

#calculate possible number of permutations or combinations

perm = function(n, k) {
  return(factorial(n) / factorial(n-k))
}


comb = function(n, k) {
  return(factorial(n) / (factorial(k) * factorial(n-k)))
}



# choose(n, k) is a base function


