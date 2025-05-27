library(MASS)
library(expm)
library(ggplot2)
library(spatstat)
library(psych)
library(Rfast)

AlgorithmRoMDP <- function(DataSet, itertime = 100)
{
  # Dimensión del conjuto de datos
  SampleNumber <- dim(DataSet)[1]
  VariableNumber <- dim(DataSet)[2]
  MiddlePoint <- round(SampleNumber/2)+1
  TransposedDataSet <-t(DataSet)
  
  ValueInitialSubsets <- 2
  BestDet <- 0
  VecZero <- numeric(SampleNumber)
  
  #Ciclo de iteracciones por default 100
  for ( iter in 1:itertime)
  {
    Id <- sample(SampleNumber, ValueInitialSubsets, replace = FALSE)
    SubsetById <- DataSet[Id,]
    MuSubsetById <- Rfast::colmeans(SubsetById)
    VarSubsetById <- Rfast::colVars(SubsetById)
    Sama <- (TransposedDataSet-MuSubsetById)/VarSubsetById
    Distance <- Rfast::colsums(Sama)
    Criterio <- 10
    Count <- 0
    while (Criterio !=0 & Count <= 15)
    {
      Count <- Count + 1
      VectorZeroi <- numeric(SampleNumber)
      DistancePerm <- order(Distance)
      VectorZeroi[ DistancePerm[1:MiddlePoint] ] <- 1
      Criterio <- sum( abs(VectorZeroi - VecZero) )
      VecZero <- VectorZeroi
      NewDataSet <- DataSet[DistancePerm[1:MiddlePoint], ]
      MuSubsetById <- Rfast::colmeans(NewDataSet)
      VarSubsetById <- Rfast::colVars(NewDataSet)
      Sama <- (TransposedDataSet-MuSubsetById)/VarSubsetById
      Distance <- Rfast::colsums(Sama)
    }
    TempDet <- prod(VarSubsetById)
    if(BestDet == 0 | TempDet < BestDet) 
    {
      BestDet <- TempDet
      FinalVec <- VecZero
    }
  }
  SubMCD <- (1:SampleNumber)[FinalVec != 0]
  
  MuSubsetById <- Rfast::colmeans( DataSet[SubMCD, ] )
  VarSubsetById <- Rfast::colVars( DataSet[SubMCD, ] )
  Sigma <- cov(DataSet[SubMCD, ])
  return(list(MuSubsetById, VarSubsetById, Sigma, SubMCD))
}


AlgorithmMDPCFPart1 <- function(DatesNorm, NumVariable, Observations, Alpha = 0.05, OutliersFlag = FALSE)
{
  AlgorithmObjects <- AlgorithmRoMDP(DatesNorm)
  #calculo de Media y matrices de correlacion y de varianza
  MeanDMP <-unlist(AlgorithmObjects[1])    
  SigmaDMP <- matrix(diag(as.numeric(unlist(AlgorithmObjects[2]))),ncol = NumVariable)
  InverseSigmaDMP <- solve(SigmaDMP)
  SigmaOfMinimunDet <- AlgorithmObjects[3][[1]]
  CorMatrix <- sqrt(InverseSigmaDMP)%*%SigmaOfMinimunDet%*%sqrt(InverseSigmaDMP)
  
  # Estimacion objetos algoritmo Ebadi
  TraceRhoSquare <- tr(CorMatrix %^% 2) - (NumVariable**2)/MiddlePoint
  TraceRhoCubic  <- tr(CorMatrix %^% 3) - ((3*NumVariable)/MiddlePoint)*tr(CorMatrix %^% 2) + ((2*(NumVariable**3))/(MiddlePoint**2))
  
  # Estimadores Ui y Zi
  ListEstimUi <- c()
  ListEstimZi <- c()
  ListOutliers <- c()
  Zalpha <- qnorm(1-Alpha,0,1)
  ConstantMDP <- 1 + (2*NumVariable)/(Observations*sqrt(TraceRhoSquare))
  
  for( i  in 1:dim(DatesNorm)[1])
  {
    DistanceMahalanobisXi <- t((DatesNorm[i,] - MeanDMP))%*%InverseSigmaDMP%*%(DatesNorm[i,] - MeanDMP)
    Ui <- (DistanceMahalanobisXi - NumVariable)/(2*ConstantMDP*sqrt(TraceRhoSquare))
    ListEstimUi <- c(ListEstimUi,Ui)
    Zi <- Ui-(4*TraceRhoCubic*(Zalpha**2-1))/(3*(2*TraceRhoSquare)^(3/2))
    ListEstimZi <- c(ListEstimZi,Zi)
  
  }
  return(list(DatesNorm,ListEstimUi, ListEstimZi))
}

mValue <- c(20,30,40, 60, 80, 100)
PRate <- c()
for (m in mValue)
{
  # Variables locales
  NumVariable <- 30
  Observations <- 10000
  
  MiddlePoint <- floor(Observations/2)+1
  Percent <- 1
  # Matrices de Covarianza
  SigmaFirstDistribution <- diag(NumVariable)
  # Vector de Medias
  Mu <- rep(0,NumVariable)
  DatesNorm <- c()
  DatesNorm <- mvrnorm(Observations,Mu,SigmaFirstDistribution)
  Firstestimation <- AlgorithmMDPCFPart1(DatesNorm, NumVariable,Observations, 0.05, FALSE)
  
  #Funcion CDF
  # ListEstimUi <- Firstestimation[[2]]
  ListEstimZi <- Firstestimation[[3]]
  NormalStandarSeq <- seq(-4, 4, length.out= Observations)
  NormalStandarCDF <- pnorm(NormalStandarSeq, mean = 0, sd = 1)
  
  
  
  # # create sample data
  ScaleDate <- scale(ListEstimZi)
  alpha = 1-(1-0.05)^(1/m)
  countOutlier = 0 
  Zalpha <- qnorm(1-alpha,0,1)
  for ( i in ScaleDate)
  {
    if(i > Zalpha)
    {
      countOutlier = countOutlier+1
    }
  }
  print(countOutlier)
  probability <- countOutlier/dim(ScaleDate)[1]
  PRate <- c(PRate,probability)
}

CDF2 <- ecdf(ScaleDate)
# 
# # draw the cdf plot
plot(CDF2,col = "red", xlim = c(-4,4))
lines(NormalStandarSeq,NormalStandarCDF, col="black")



