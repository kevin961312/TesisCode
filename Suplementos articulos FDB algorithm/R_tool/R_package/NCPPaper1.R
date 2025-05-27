rm(list=ls())

library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)

set.seed(125) 

MethodUCLKernel <- function(T2, alpha)
{
  EstKernelSmooth <- density(T2,kernel="gaussian",bw="nrd")
  valuePairs <- cbind(EstKernelSmooth$x,EstKernelSmooth$y)
  valuePairs <- valuePairs[valuePairs[,1]>=0,]
  Estand <- valuePairs[,2]/sum(valuePairs[,2])
  NumberRow <- nrow(valuePairs)
  countValues <- 0
  
  for(i in 1:NumberRow)
  {
    if(countValues <  (1-alpha))
    {
      j <- i 
      countValues <- sum(Estand[1:i])
    }
  }
  
  UCL <- valuePairs[j,1]
  return(list(UCL=UCL,alpha1=(1-countValues)))
}

Observation <- 30 

NumberVariable <- 2
SigmaCorr <- diag(nrow=NumberVariable)
mu <- rep(0,NumberVariable)
numSimulation <- 5000


T2Matrix <- c()
T2Total <- c()
T2Max <- c()
for (i in 1:numSimulation)
{
  
  Data <- mvrnorm(Observation,mu=mu,Sigma = SigmaCorr)
  SigmaMRCD <- cov(as.matrix.data.frame(Data))
  MediaMRCD <- colMeans(Data)

  T2 <- mahalanobis(Data, center= MediaMRCD, cov = SigmaMRCD)
  T2Total <- c(T2Total,T2)
  T2Max <- c(T2Max,max(T2))
  
}
T2Matrix  <- matrix(T2Total, ncol = Observation)

AlphaPFA <- 0.05

UCLKernel <- MethodUCLKernel(T2Total,AlphaPFA)
UCL <- qemp(p = 1-AlphaPFA, obs = T2Total)
UCLMax <- qemp(p = 1-AlphaPFA, obs = T2Max)



Matrix <- T2Matrix
SignalCount <- 0
SignalCountMax <- 0
PFAT <- c()
PFAMax <- c()
for (i in 1: dim(Matrix)[1])
{
  RowMatrix <- Matrix[i,]
  for (j in 1:length(RowMatrix))
  {
    if (RowMatrix[j] > UCL)
    {
      SignalCount <- SignalCount + 1
    }
  }
  
  # listfilter <- Filter(function(x) x > UCL,RowMatrix)
  # SignalCount <- SignalCount + length(listfilter)
  # PFAT <- c(length(listfilter)/length(RowMatrix),PFAT)
  listfilterMax <- Filter(function(x) x > UCLMax,RowMatrix)
  SignalCountMax <- SignalCountMax + length(listfilterMax)
  PFAMax <- c(length(listfilterMax)/length(RowMatrix),PFAMax)
}

TotalValue <- dim(Matrix)[1]*dim(Matrix)[2]

SignalPro <- SignalCount/TotalValue
SignalProMax <- SignalCountMax/TotalValue





# Parte 2

percentoutliers = 5/30
shiftMean = c(1,sqrt(29))
inverse <- solve(SigmaCorr)
deltaNCP <- t(shiftMean-mu)%*%inverse%*%(shiftMean-mu)
T2Matrix2 <- c()
T2Total2 <- c()

for (i in 1:1500)
{
  
  NumOutliers = floor(percentoutliers*Observation)
  Data1 <- mvrnorm(Observation-NumOutliers,mu=mu,Sigma = SigmaCorr)
  if(NumOutliers >= 1)
  {
    DataOutlier <- mvrnorm(NumOutliers, mu = shiftMean, Sigma = SigmaCorr)
    Data1 <- rbind(Data, DataOutlier)
  }
  SigmaMRCD1 <- cov(Data1)
  MediaMRCD1 <- colMeans(Data1)
  
  T21 <- mahalanobis(DataOutlier, center= MediaMRCD1, cov = SigmaMRCD1)
  T2Total2 <- c(T2Total2,T21)
  
}
T2Matrix2  <- matrix(T2Total2, ncol = 1)



Matrix2 <- T2Matrix2
SignalCount2 <- 0
SignalCountMax2 <- 0
PFAT2 <- c()
PFAMax2 <- c()
for (i in 1: dim(Matrix2)[1])
{
  RowMatrix2 <- Matrix2[i,]
  for (j in 1:length(RowMatrix2))
  {
    if (RowMatrix2[j] > UCL)
    {
      SignalCount2 <- SignalCount2 + 1
    }
  }
  
  # listfilter <- Filter(function(x) x > UCL,RowMatrix)
  # SignalCount <- SignalCount + length(listfilter)
  # PFAT <- c(length(listfilter)/length(RowMatrix),PFAT)
  listfilterMax2 <- Filter(function(x) x > UCLMax,RowMatrix2)
  SignalCountMax2 <- SignalCountMax2 + length(listfilterMax2)
  PFAMax2 <- c(length(listfilterMax2)/length(RowMatrix2),PFAMax2)
}

TotalValue2 <- dim(Matrix2)[1]*dim(Matrix2)[2]

SignalPro2 <- SignalCount2/TotalValue2
SignalProMax2 <- SignalCountMax2/TotalValue2


