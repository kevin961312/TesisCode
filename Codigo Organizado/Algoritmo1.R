source("C:/Users/garfo/Downloads/Tesis/Codigo Organizado/AlgoritmoRo.R")
library(MASS)
NumVariable <- 30
Observations <- 10000
SigmaFirstDistribution <- diag(NumVariable)
Mu <- rep(0,NumVariable)
DatesNorm <- mvrnorm(Observations,Mu,SigmaFirstDistribution)
AlgorithmObjects <- AlgorithmRoMDP(DatesNorm)

