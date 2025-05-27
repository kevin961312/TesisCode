rm(list=ls())

library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)

SimulationT2MRCDChart <- function (observation, numVariables, 
                                   numSimulation, meanVector, 
                                   sigmaMatriz, shiftMean, percentoutliers = 0, alphaMRCD = 0.75 )
{

  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()
  
  for (i in 1:numSimulation)
  {
    
    NumOutliers = floor(percentoutliers*observation)
    Data <- mvrnorm(Observation-NumOutliers,mu=mu,Sigma = sigmaMatriz)
    if(NumOutliers >= 1)
    {
      DataOutlier <- mvrnorm(NumOutliers, mu = shiftMean, Sigma = sigmaMatriz)
      Data <- rbind(Data, DataOutlier)
    }
    
    
    CovMRCD <- CovMrcd(Data, alpha = alphaMRCD)
    MediaMRCD <- CovMRCD$center
    SigmaMRCD <- CovMRCD$cov
    BestSubset <- CovMRCD$best
    
    DataBestSubset <- Data[BestSubset,] 
    T2 <- mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
    T2Total <- c(T2Total,T2)
    T2Max <- c(T2Max,max(T2))
  }
  LenBestSset <- length(BestSubset)
  T2Matrix  <- matrix(T2Total, ncol = LenBestSset)
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
}

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


Observation <- 50
NumberVariable <- 100
NumSimulation <- 100
SigmaCorr <- randcorr(NumberVariable)
mu <- rep(0,NumberVariable)

# shiftmu <- sample(-5:5, size = NumberVariable, replace = TRUE)
shiftmu <- rep(0.04,NumberVariable)
# Matriz de covarianza
# DPMatriz <- genPositiveDefMat(dim = NumberVariable, covMethod = "unifcorrmat")
# SigmaCov <- DPMatriz$Sigma

inverse <- solve(SigmaCorr)
deltaNCP <- sqrt(t(shiftmu-mu)%*%inverse%*%(shiftmu-mu))

Values <- SimulationT2MRCDChart(Observation,
                                NumberVariable,
                                NumSimulation,
                                mu,
                                SigmaCorr,
                                shiftmu,
                                0.1,0.75)



#alpha = 0.005  , 200
#alpha  = 0.0027, 370
#alpha = 0.05 , 20
AlphaPFA <- 0.05

UCLKernel <- MethodUCLKernel(Values$T2Total,AlphaPFA)
UCL <- qemp(p = 1-AlphaPFA, obs = Values$T2Total) 
UCLMax <- qemp(p = 1-AlphaPFA, obs = Values$T2Max) 
UCLMeanMaxT2 <- mean(Values$T2Max)

P1 <- pemp(UCLMax,Values$T2Total)
P2 <- pemp(UCLMeanMaxT2,Values$T2Total)
T2 <- Values$T2Matrix[1,]



Matrix <- Values$T2Matrix
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

plot(ecdf(Values$T2Total),
     main = "ECDF Todos los valores T2 de la simulación",
     col = "blue")
dev.new()
plot(ecdf(Values$T2Max),
     main = "ECDF los valores maximos T2 de la simulación",
     col = "red")

dev.new()
plot(1:length(T2),
     T2,
     type='l',
     xlim=c(0,length(T2)+2),
     ylim=c(min(T2)-2, max(UCL,UCLMax,UCLMeanMaxT2)+2),
     main=expression("Carta"*~T^2),
     ylab=expression(T^2),xlab="No. Observacion",font=2)
abline(h=UCL,lty=3, col = "blue")
abline(h=UCLMax,lty=3,col = "red")
abline(h=UCLMeanMaxT2,lty=3,col = "green")




















