rm(list=ls())

library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(psych)

sigmaMatriz <- diag(300)
mu <- rep(0,300)
shiftMean <- rep(1,300)
NumOutliers = 10
Data <- mvrnorm(90,mu=mu,Sigma = sigmaMatriz)
DataOutlier <- mvrnorm(NumOutliers, mu = shiftMean, Sigma = sigmaMatriz) 
Data <- rbind(Data, DataOutlier)
CovMRCD <- CovMrcd(Data, alpha = 0.75)

BestSubset <- CovMRCD$best
MediaMRCD <- CovMRCD$center
SigmaMRCD <- CovMRCD$cov
DSigmaMRCD  <- diag(SigmaMRCD)
DMRCD <- diag(DSigmaMRCD)
InverseSigmaDMP <- solve(DMRCD)
CorMatrix <- sqrt(InverseSigmaDMP)%*%SigmaMRCD%*%sqrt(InverseSigmaDMP)
 
ZAlpha <- qnorm(1-0.05,mean=0,sd=1)
Phizalpha <- pnorm(ZAlpha,mean=0,sd=1)
TR2 <- tr(CorMatrix%*%CorMatrix)-300*300/100
Constantaprox <- 1+Phizalpha*sqrt(2*TR2)/(300*(1-0.05))
CP <- 1+tr(CorMatrix%*%CorMatrix)/(300^(2/3))

T2 <- mahalanobis(Data , center= MediaMRCD, cov = DMRCD)
T2Final <- Constantaprox*T2

ConstanOutlier <- 300+ZAlpha*sqrt(2*TR2*CP)
