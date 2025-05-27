rm(list=ls())

library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)
library(psych)
library(expm)
library(Rfast)

set.seed(123)


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


AlgorithmMDPCFPart1 <- function(DatesNorm, NumVariable, Observations, Alpha = 0.05, outliersData = c(), OutliersFlag = FALSE )
{
  MiddlePoint <- floor(Observations/2)+1
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
  
  
  if(OutliersFlag)
  {
    for( i  in 1:dim(outliersData)[1])
    {
      DistanceMahalanobisXi <- t((outliersData[i,] - MeanDMP))%*%InverseSigmaDMP%*%(outliersData[i,] - MeanDMP)
      Ui <- (DistanceMahalanobisXi - NumVariable)/(2*ConstantMDP*sqrt(TraceRhoSquare))
      ListEstimUi <- c(ListEstimUi,Ui)
      Zi <- Ui-(4*TraceRhoCubic*(Zalpha**2-1))/(3*(2*TraceRhoSquare)^(3/2))
      ListEstimZi <- c(ListEstimZi,Zi)
    }
  }
  else
  {
    for( i  in 1:dim(DatesNorm)[1])
    {
      DistanceMahalanobisXi <- t((DatesNorm[i,] - MeanDMP))%*%InverseSigmaDMP%*%(DatesNorm[i,] - MeanDMP)
      Ui <- (DistanceMahalanobisXi - NumVariable)/(2*ConstantMDP*sqrt(TraceRhoSquare))
      ListEstimUi <- c(ListEstimUi,Ui)
      Zi <- Ui-(4*TraceRhoCubic*(Zalpha**2-1))/(3*(2*TraceRhoSquare)^(3/2))
      ListEstimZi <- c(ListEstimZi,Zi)
      
    }
  }
  return(list("Ui" = ListEstimUi, "Zi"= ListEstimZi))
}

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
      tic()
      CovMRCD = CovMrcd(Data, alpha = alphaMRCD)
      exectime <- toc()
      exectime <- exectime$toc - exectime$tic
      print("MRCD tiempo")
      MediaMRCD = CovMRCD$center
      SigmaMRCD = CovMRCD$cov
      BestSubset = CovMRCD$best
      
      DataBestSubset = Data[BestSubset,] 
      T2 = mahalanobis(DataBestSubset , center= MediaMRCD, cov = SigmaMRCD)
      LenMatrix = length(BestSubset)
    }
    else if (typeMethod == "T2MOD")
    {
      tic()
      Media = colMeans(Data)
      Sigma = cov(Data)
      exectime <- toc()
      exectime <- exectime$toc - exectime$tic
      print("T2MOD tiempo")
      SigmaMod = 1/(sum(diag(Sigma))/numVariables)
      
      T2 = c()
      for (j in 1:dim(Data)[1])
      {
        Ai = (observation/(observation-1))*(norm((Data[j,]-Media)/sqrt(numVariables)^2, type = "2")/SigmaMod)
        T2 = c(Ai,T2)
      }
      LenMatrix = observation
    }
    else if (typeMethod == "EBADIUI")
    {
      T2 = c()
      tic()
      Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05, FALSE)
      exectime <- toc()
      exectime <- exectime$toc - exectime$tic
      print("EBADIUI tiempo")
      T2 = c(Firstestimation$Ui,T2)
      LenMatrix = observation
    }
    else if (typeMethod == "EBADIZI")
    {
      T2 = c()
      tic()
      Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05, FALSE)
      exectime <- toc()
      exectime <- exectime$toc - exectime$tic
      print("EBADIZI tiempo")
      T2 = c(Firstestimation$Zi,T2)
      LenMatrix = observation
    }
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  
  T2Matrix  = matrix(T2Total, ncol = LenMatrix)
  return (exectime)
}

#Parameters Simulación
Observation = 100
NumberVariable = 250
NumSimulation = 1
SigmaCorr = randcorr(NumberVariable)
mu = rep(0,NumberVariable)
AlphaMRCD = 0.75
AlphaPFA = 0.05
NumSimulationDelta = 1

ValuesMRCD = c()
for (i in 1:100) {
  # MRCD
  timeMRCD = SimulationT2Chart(Observation,
                                 NumberVariable,
                                 NumSimulation,
                                 mu,
                                 SigmaCorr,
                                 AlphaMRCD,
                                 "MRCD")
  ValuesMRCD = c(ValuesMRCD,timeMRCD)
}


ValuesT2MOD = c()
for (i in 1:100) {
#T2MOD

  timeT2MOD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                mu,
                                SigmaCorr,
                                AlphaMRCD,
                                "T2MOD")
  ValuesT2MOD = c(ValuesT2MOD,timeT2MOD)
}

ValuesEBADIUI = c()
for (i in 1:100) {
# EBADIUI
  timeEBADIUI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  mu,
                                  SigmaCorr,
                                  AlphaMRCD,
                                  "EBADIUI")
  ValuesEBADIUI = c(ValuesEBADIUI,timeEBADIUI)
}

ValuesEBADIZI = c()
for (i in 1:100) {
# EBADIZI
  timeEBADIZI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  mu,
                                  SigmaCorr,
                                  AlphaMRCD,
                                  "EBADIZI")
  ValuesEBADIZI = c(ValuesEBADIZI,timeEBADIZI)
}

meanTimeMRCD = mean(ValuesMRCD)
meanTimeT2MOD = mean(ValuesT2MOD)
meanTimeEBADIUI = mean(ValuesEBADIUI)
meanTimeEBADIZ = mean(ValuesEBADIZI)

