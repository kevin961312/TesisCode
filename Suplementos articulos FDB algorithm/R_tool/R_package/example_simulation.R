library(MASS)
library(DetMCD)
library(FDB)
library(pracma)
library(rrcov)
library(randcorr)
# help(FDB)

#### example for using these functions n=400,p=40
x <- mvrnorm(39, rep(0, 40), diag(40))
pro <- FDB(x, alpha = 0.75, depth = "pro")
centerpro <- pro$center
covpro <- pro$cov
bestpro <- pro$best

# l2 <- FDB(x, alpha = 0.75, depth = "l2")
# centerl2 <- l2$center
# covl2 <- l2$cov

SimulationT2Chart <- function(observation, numVariables, numSimulation, meanVector, sigmaMatriz, alphaMRCD = 0.75, typeMethod = "MRCD") {
  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()

  for (i in 1:numSimulation)
  {
    Data <- mvrnorm(Observation, mu = mu, Sigma = sigmaMatriz)
    if (typeMethod == "MRCD") {
      CovMRCD <- CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD <- CovMRCD$center
      SigmaMRCD <- CovMRCD$cov
      BestSubset <- CovMRCD$best

      DataBestSubset <- Data[BestSubset, ]
      T2 <- mahalanobis(DataBestSubset, center = MediaMRCD, cov = SigmaMRCD)
      LenMatrix <- length(BestSubset)
    } else if (typeMethod == "pro") {
      CovPro <- FDB(Data, alpha = alphaMRCD, depth = "pro")
      MediaMRCD <- CovPro$center
      SigmaMRCD <- CovPro$cov
      BestSubset <- CovPro$best

      DataBestSubset <- Data[BestSubset, ]
      T2 <- mahalanobis(DataBestSubset, center = MediaMRCD, cov = SigmaMRCD)
      LenMatrix <- length(BestSubset)
    } else if (typeMethod == "l2") {
      Covl2 <- FDB(Data, alpha = alphaMRCD, depth = "l2")
      MediaMRCD <- Covl2$center
      SigmaMRCD <- Covl2$cov
      BestSubset <- Covl2$best

      DataBestSubset <- Data[BestSubset, ]
      T2 <- mahalanobis(DataBestSubset, center = MediaMRCD, cov = SigmaMRCD)
      LenMatrix <- length(BestSubset)
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

MethodUCLKernel <- function(T2, alpha) {
  EstKernelSmooth <- density(T2, kernel = "gaussian", bw = "nrd")
  valuePairs <- cbind(EstKernelSmooth$x, EstKernelSmooth$y)
  valuePairs <- valuePairs[valuePairs[, 1] >= 0, ]
  Estand <- valuePairs[, 2] / sum(valuePairs[, 2])
  NumberRow <- nrow(valuePairs)
  countValues <- 0

  for (i in 1:NumberRow)
  {
    if (countValues < (1 - alpha)) {
      j <- i
      countValues <- sum(Estand[1:i])
    }
  }

  UCL <- valuePairs[j, 1]
  return(list(UCL = UCL, alpha1 = (1 - countValues)))
}

# Parameters Simulación
Observation <- 100
NumberVariable <- 150
NumSimulation <- 100
SigmaCorr <- randcorr(NumberVariable)
mu <- rep(0, NumberVariable)
AlphaMRCD <- 0.75
AlphaPFA <- 0.05
NumSimulationDelta <- 100

x <- mvrnorm(Observation, mu, SigmaCorr)
pro <- FDB(x, alpha = 0.75, depth = "pro")
centerpro <- pro$center
covpro <- pro$cov
bestpro <- pro$best

# MRCD
ValuesMRCD <- SimulationT2Chart(
  Observation,
  NumberVariable,
  NumSimulation,
  mu,
  SigmaCorr,
  AlphaMRCD,
  "MRCD"
)
KernelMethodMRCD <- MethodUCLKernel(ValuesMRCD$T2Total, AlphaPFA)
UCLMRCD <- qemp(p = 1 - AlphaPFA, obs = ValuesMRCD$T2Total)
UCLMaxMRCD <- qemp(p = 1 - AlphaPFA, obs = ValuesMRCD$T2Max)
UCLKernelMRCD <- KernelMethodMRCD$UCL

# Pro
ValuesPro <- SimulationT2Chart(
  Observation,
  NumberVariable,
  NumSimulation,
  mu,
  SigmaCorr,
  AlphaMRCD,
  "pro"
)
KernelMethodPro <- MethodUCLKernel(ValuesPro$T2Total, AlphaPFA)
UCLPro <- qemp(p = 1 - AlphaPFA, obs = ValuesPro$T2Total)
UCLMaxPro <- qemp(p = 1 - AlphaPFA, obs = ValuesPro$T2Max)
UCLKernelPro <- KernelMethodPro$UCL

# L2
ValuesL2 <- SimulationT2Chart(
  Observation,
  NumberVariable,
  NumSimulation,
  mu,
  SigmaCorr,
  AlphaMRCD,
  "pro"
)
KernelMethodL2 <- MethodUCLKernel(ValuesL2$T2Total, AlphaPFA)
UCLL2 <- qemp(p = 1 - AlphaPFA, obs = ValuesL2$T2Total)
UCLMaxL2 <- qemp(p = 1 - AlphaPFA, obs = ValuesL2$T2Max)
UCLKernelL2 <- KernelMethodL2$UCL
