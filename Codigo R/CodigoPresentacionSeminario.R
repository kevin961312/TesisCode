library(MASS)
library(rrcov)
library(randcorr)



Observation <-200
NumberVariable <- 400
Percentoutliers <- 0.4
AlphaMRCD <- 0.75
Mu <- rep(0,NumberVariable)
Shiftmu <- rep(6,NumberVariable)
TotalSum <- c()
for (i in 1:50)
{
  SigmaMatriz <- randcorr(NumberVariable)
  NumOutliers <- floor(Percentoutliers*Observation)
  Data <- mvrnorm(Observation-NumOutliers,mu= Mu,Sigma = SigmaMatriz)
  if(NumOutliers >= 1)
  {
    DataOutlier <- mvrnorm(NumOutliers, mu = Shiftmu, Sigma = SigmaMatriz)
    Data <- rbind(Data, DataOutlier)
  }
  
  CovMRCD <- CovMrcd(Data, alpha = AlphaMRCD)
  MediaMRCD <- CovMRCD$center
  SigmaMRCD <- CovMRCD$cov
  BestSubset <- CovMRCD$best
  
  
  DifMatCov <- SigmaMatriz-SigmaMRCD
  SumTotal <- sum(rowSums(DifMatCov%*%DifMatCov))/(NumberVariable*NumberVariable)
  TotalSum <- c(SumTotal,TotalSum)
}



