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
library(fftw)
library(lcmix)

source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/MethodUCLKernel.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SignalProbability.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/AlgorithmEBADI.R")) 
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/kMRCD.R"))

SimulationT2Chart <- function(observation, numVariables, numSimulation, sigma, alphaMRCD = 0.75, typeMethod = "MRCD", dfreedom = 1) {
  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()
  
  for (i in 1:numSimulation)
  {
    Data <- rmvgamma(observation,rep(1,numVariables),rep(1,numVariables),sigma)
    if (typeMethod == "KMRCD") {
      T2 <- c()
      CovKMRCD <- CovKmrcd(x = Data, kModel = "RbfKernel", alpha = alphaMRCD)
      DataBestSubset <- Data[CovKMRCD$hsubsetIndices, ]
      meanKMRCD <- colMeans(DataBestSubset)
      T2 <- mahalanobis(DataBestSubset, center = meanKMRCD, cov = CovKMRCD$cov)
      LenMatrix <- length(CovKMRCD$hsubsetIndices)
    }
    T2Total <- c(T2Total, T2)
    T2Max <- c(T2Max, max(T2))
  }
  
  T2Matrix <- matrix(T2Total, ncol = LenMatrix)
  return(list(
    "T2Matrix" = T2Matrix,
    "T2Total" = T2Total,
    "T2Max" = T2Max
  ))
}




SimulationT2ChartOutliers <- function(observation, numVariables,
                                      numSimulation, sigma, rho, rhoShift, percentoutliers = 0, alphaMRCD = 0.75, typeMethod = "MRCD", dfreedom = 1) {
  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()
  
  for (i in 1:numSimulation)
  {
    NumOutliers <- floor(percentoutliers * observation)
    if (identical(rho, rhoShift)) {
      Data <- rmvgamma(observation,rep(1,numVariables),rep(1,numVariables),sigma)
    } else {
      Data <- rmvgamma(observation - NumOutliers,rep(1,numVariables),rep(1,numVariables),sigma)
      if (NumOutliers >= 1) {
        DataOutlier <- rmvgamma(NumOutliers,rhoShift,rep(1,numVariables),sigma)
        Data <- rbind(Data, DataOutlier)
      }
    }
    if (typeMethod == "KMRCD") {
      T2 <- c()
      CovKMRCD <- CovKmrcd(x = Data, kModel = "RbfKernel",alpha = alphaMRCD)
      if (identical(rho, rhoShift)) {
        DataBestSubset <- Data[CovKMRCD$hsubsetIndices, ]
        meanKMRCD <- colMeans(DataBestSubset)
        T2 <- mahalanobis(DataBestSubset, center = meanKMRCD, cov = CovKMRCD$cov)
      } else {
        DataBestSubset <- Data[CovKMRCD$hsubsetIndices, ]
        meanKMRCD <- colMeans(DataBestSubset)
        DataBestSubsetOutliers <- rbind(DataBestSubset, DataOutlier)
        T2 <- mahalanobis(DataOutlier, center = meanKMRCD, cov = CovKMRCD$cov)
      }
      LenMatrix <- observation
    }
    T2Total <- c(T2Total, T2)
    T2Max <- c(T2Max, max(T2))
  }
  T2Matrix <- matrix(T2Total, ncol = 1)
  return(list(
    "T2Matrix" = T2Matrix,
    "T2Total" = T2Total,
    "T2Max" = T2Max
  ))
}

set.seed(123)

#Parameters Simulación
Observation = 150
NumberVariable = 200
NumSimulation = 100
fun <- function(i,j) (0.5)^(abs(i-j))

rows <- 1:NumberVariable
cols <- 1:NumberVariable

SigmaCorr = outer(rows,cols,FUN=fun)
AlphaMRCD = 0.75
AlphaPFA = 0.05
NumSimulationDelta = 100

# KMRCD
ValuesKMRCD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                SigmaCorr,
                                AlphaMRCD,
                                "KMRCD")
KernelMethodKMRCD = MethodUCLKernel(ValuesKMRCD$T2Total,AlphaPFA)
UCLKMRCD = qemp(p = 1-AlphaPFA, obs = ValuesKMRCD$T2Total) 
UCLMaxKMRCD = 0.5
UCLKernelKMRCD = KernelMethodKMRCD$UCL


SignalProbabilityT2KMRCD = SignalProbability(ValuesKMRCD$T2Matrix, UCLKMRCD,UCLMaxKMRCD, UCLKernelKMRCD)


# Parte 2 Calculo del parametro de no centralidad 

DataCopula <- rmvgamma(1000,rep(1,NumberVariable),rep(1,NumberVariable),SigmaCorr)
mu <- colMeans(DataCopula)
Sigma <- SigmaCorr
ArrayRhoShift = list(mu,
                     mu+1/100,
                     mu+5/100,
                     mu+10/100,
                     mu+15/100,
                     mu+25/100,
                     mu+50/100,
                     mu+60/100,
                     mu+75/100,
                     mu+90/100,
                     mu+1)

Inverse = solve(Sigma)
Percentoutliers = 0.1
MatrixDeltaKMRCD = matrix(,ncol = 4)
for (rhoShift in ArrayRhoShift)
{
  DeltaNCP = sqrt(t(rhoShift-rep(1,NumberVariable))%*%Inverse%*%(rhoShift-rep(1,NumberVariable)))
  
  ValuesOutliersKMRCD = SimulationT2ChartOutliers(Observation,
                                                  NumberVariable,
                                                  NumSimulationDelta,
                                                  SigmaCorr,
                                                  mu,
                                                  rhoShift,
                                                  Percentoutliers,
                                                  AlphaMRCD,
                                                  "KMRCD")
  
  SignalProbabilityT2OutliersKMRCD = SignalProbability(ValuesOutliersKMRCD$T2Matrix, UCLKMRCD,UCLMaxKMRCD, UCLMaxKMRCD)
  RowdeltaKMRCD  = c(DeltaNCP[1,1],
                     SignalProbabilityT2OutliersKMRCD$SignalPro,
                     SignalProbabilityT2OutliersKMRCD$SignalProMax,
                     SignalProbabilityT2OutliersKMRCD$SignalProKernel)
  MatrixDeltaKMRCD = rbind(MatrixDeltaKMRCD,RowdeltaKMRCD)
  
}

MatrixDeltaKMRCD = MatrixDeltaKMRCD[-1,]

p = plot(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], col='purple',pch = 15,type='b',lty = 6)

