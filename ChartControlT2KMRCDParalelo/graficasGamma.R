rm(list=ls())
library(here)
library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)
library(psych)
library(expm)
library(Rfast)
library(copula)
library(lcmix)
library(plotly)


NumberVariable =2 
fun <- function(i,j) (-0.5)^(abs(i-j))
rows <- 1:NumberVariable
cols <- 1:NumberVariable
SigmaCorr = outer(rows,cols,FUN=fun)

x     <- seq(0, 5, 0.01) 
y     <- seq(0, 5, 0.01)
f <- function(x, y) dmvgamma(cbind(x, y),c(1,1),c(0.5,0.5),SigmaCorr)
fig <- plot_ly(x = x, y = y, z = outer(x, y, f), type = "surface")
fig



x     <- seq(0, 1, 0.01) 
y     <- seq(0, 1, 0.01)
f <- function(x, y) dmvgamma(cbind(x, y),c(2,2),c(5,5),SigmaCorr)
fig <- plot_ly(x = x, y = y, z = outer(x, y, f), type = "surface")
fig



x     <- seq(0, 10, 0.01) 
y     <- seq(0, 10, 0.01)
f <- function(x, y) dmvgamma(cbind(x, y),c(5,5),c(1,1),SigmaCorr)
fig <- plot_ly(x = x, y = y, z = outer(x, y, f), type = "surface")
fig

x     <- seq(0, 10, 0.01) 
y     <- seq(0, 10, 0.01)
f <- function(x, y) dmvgamma(cbind(x, y),c(5,1),c(3,0.5),SigmaCorr)
fig <- plot_ly(x = x, y = y, z = outer(x, y, f), type = "surface")
fig




