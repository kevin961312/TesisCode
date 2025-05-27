source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/AlgorithmEBADI.R")) 
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/kMRCD1.R"))

SimulationT2Chart <- function(observation, numVariables, numSimulation, sigma, rho, alphaMRCD = 0.75, 
typeMethod = "MRCD", beta = 1, KMRCDModel= "RbfKernel") {
  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()

  for (i in 1:numSimulation)
  {
    Data <- rmvgamma(observation,rep(rho,numVariables),rep(beta,numVariables),sigma)
    if (typeMethod == "MRCD") {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      
      DataBestSubset = Data[BestSubset,] 
      T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      LenMatrix = length(BestSubset)
    } else if (typeMethod == "T2MOD") {
      Media <- colMeans(Data)
      Sigma <- cov(Data)
      SigmaMod <- 1 / (sum(diag(Sigma)) / numVariables)

      T2 <- c()
      for (j in 1:dim(Data)[1])
      {
        Ai <- (observation / (observation - 1)) * (norm((Data[j, ] - Media) / sqrt(numVariables)^2, type = "2") / SigmaMod)
        T2 <- c(Ai, T2)
      }
      LenMatrix <- observation
    } else if (typeMethod == "EBADIUI") {
      T2 <- c()
      Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, FALSE)
      T2 <- c(Firstestimation$Ui, T2)
      LenMatrix <- observation
    } else if (typeMethod == "EBADIZI") {
      T2 <- c()
      Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, FALSE)
      T2 <- c(Firstestimation$Zi, T2)
      LenMatrix <- observation
    } else if (typeMethod == "KMRCD") {
      T2 <- c()
      CovKMRCD <- CovKmrcd(x = Data, kModel = KMRCDModel,alpha = alphaMRCD)
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
                                      numSimulation, sigma, rho, rhoShift, percentoutliers = 0, 
                                      alphaMRCD = 0.75, typeMethod = "MRCD", beta = 1, KMRCDModel= "RbfKernel") {
  T2Matrix <- c()
  T2Total <- c()
  T2Max <- c()

  for (i in 1:numSimulation)
  {
    NumOutliers <- floor(percentoutliers * observation)
    if (rhoShift == 0) {
      Data <- rmvgamma(observation,rep(rho,numVariables),rep(beta,numVariables),sigma)
    } else {
      Data <- rmvgamma(observation - NumOutliers,rep(rho,numVariables),rep(beta,numVariables),sigma)
      if (NumOutliers >= 1) {
        DataOutlier <- rmvgamma(NumOutliers,rep(rho + rhoShift,numVariables),rep(beta,numVariables),sigma)
        Data <- rbind(Data, DataOutlier)
      }
    }
    if (typeMethod == "MRCD") {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      DataBestSubset = Data[BestSubset,]
      if (rhoShift == 0) {
        T2 <- mahalanobis(DataBestSubset, center = MediaMRCD, cov = SigmaMRCD)
      } else {
        T2 <- mahalanobis(DataOutlier, center = MediaMRCD, cov = SigmaMRCD)
      }
    } else if (typeMethod == "T2MOD") {
      Media <- colMeans(Data)
      Sigma <- cov(Data)
      SigmaMod <- 1 / (sum(diag(Sigma)) / numVariables)

      T2 <- c()
      if (rhoShift == 0) {
        for (j in 1:dim(Data)[1])
        {
          Ai <- (observation / (observation - 1)) * (norm((Data[j, ] - Media) / sqrt(numVariables)^2, type = "2") / SigmaMod)
          T2 <- c(Ai, T2)
        }
      } else {
        for (j in 1:dim(DataOutlier)[1])
        {
          Ai <- (observation / (observation - 1)) * (norm((DataOutlier[j, ] - Media) / sqrt(numVariables)^2, type = "2") / SigmaMod)
          T2 <- c(Ai, T2)
        }
      }
    } else if (typeMethod == "EBADIUI") {
      T2 <- c()

      if (rhoShift == 0) {
        Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, c(), FALSE)
        T2 <- c(Firstestimation$Ui, T2)
      } else {
        Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, DataOutlier, TRUE)
        T2 <- c(Firstestimation$Ui, T2)
      }
      LenMatrix <- observation
    } else if (typeMethod == "EBADIZI") {
      T2 <- c()
      if (rhoShift == 0) {
        Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, c(), FALSE)
        T2 <- c(Firstestimation$Zi, T2)
      } else {
        Firstestimation <- AlgorithmMDPCFPart1(Data, numVariables, observation, 0.05, DataOutlier, TRUE)
        T2 <- c(Firstestimation$Zi, T2)
      }
      LenMatrix <- observation
    } else if (typeMethod == "KMRCD") {
      T2 <- c()
      CovKMRCD <- CovKmrcd(x = Data,kModel = KMRCDModel,alpha = alphaMRCD)
      if (rhoShift == 0) {
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
