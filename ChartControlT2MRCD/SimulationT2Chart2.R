source(here::here("Downloads/Tesis/TesisCode/ChartControlT2MRCD/AlgorithmEBADI.R"))

SimulationT2Chart = function (observation, numVariables, numSimulation, meanVector, sigmaMatriz, alphaMRCD = 0.75 , typeMethod = "MRCD")
{
  
  T2Matrix = c()
  T2Total = c()
  T2Max = c()
  
  for (i in 1:numSimulation)
  {
    Data = mvrnorm(Observation,mu=mu,Sigma = sigmaMatriz)
    if (typeMethod == "MRCD")
    {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      
      DataBestSubset = Data[BestSubset,] 
      T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      LenMatrix = length(BestSubset)
    }
    else if(typeMethod == "MRCD2")
    {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      
      DataBestSubset = Data[BestSubset,] 
      T2 = mahalanobis(Data , center= MediaMRCD, cov = SigmaMRCD)
      LenMatrix = length(BestSubset)
    }
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  
  T2Matrix  = matrix(T2Total, ncol = LenMatrix)
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
}




SimulationT2ChartOutliers = function (observation, numVariables, 
                                      numSimulation, meanVector, 
                                      sigmaMatriz, shiftMean, percentoutliers = 0, alphaMRCD = 0.75,  typeMethod = "MRCD")
{
  T2Matrix = c()
  T2Total = c()
  T2Max = c()
  
  for (i in 1:numSimulation)
  {
    
    NumOutliers = floor(percentoutliers*observation)
    if(identical(meanVector, shiftMean))
    {
      Data = mvrnorm(Observation,mu=mu,Sigma = sigmaMatriz)    
    }
    else 
    {
      Data = mvrnorm(Observation-NumOutliers,mu=mu,Sigma = sigmaMatriz)
      if(NumOutliers >= 1)
      {
        DataOutlier = mvrnorm(NumOutliers, mu = shiftMean, Sigma = sigmaMatriz)
        Data = rbind(Data, DataOutlier)
      }
    }
    if (typeMethod == "MRCD")
    {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      DataBestSubset = Data[BestSubset,]
      if(identical(meanVector, shiftMean))
      {
        T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      }
      else 
      {
        T2 = mahalanobis(DataOutlier, center= MediaMRCD, cov = SigmaMRCD)
      }
    }
    else if (typeMethod == "MRCD2")
    {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      DataBestSubset = Data[BestSubset,]
      if(identical(meanVector, shiftMean))
      {
        T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      }
      else 
      {
        T2 = mahalanobis(rbind(DataBestSubset, DataOutlier), center= MediaMRCD, cov = SigmaMRCD)
      }
    }else if (typeMethod == "MRCD3")
    {
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      DataBestSubset = Data[BestSubset,]
      if(identical(meanVector, shiftMean))
      {
        T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      }
      else 
      {
        T2 = mahalanobis(Data, center= MediaMRCD, cov = SigmaMRCD)
      }
    }
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  T2Matrix  = matrix(T2Total, ncol = 1)
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
}