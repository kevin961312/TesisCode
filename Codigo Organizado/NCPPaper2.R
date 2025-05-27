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
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
}


MethodUCLKernel = function(T2, alpha)
{
  EstKernelSmooth = density(T2,kernel="gaussian",bw="nrd")
  valuePairs = cbind(EstKernelSmooth$x,EstKernelSmooth$y)
  valuePairs = valuePairs[valuePairs[,1]>=0,]
  Estand = valuePairs[,2]/sum(valuePairs[,2])
  NumberRow = nrow(valuePairs)
  countValues = 0
  
  for(i in 1:NumberRow)
  {
    if(countValues <  (1-alpha))
    {
      j = i 
      countValues = sum(Estand[1:i])
    }
  }
  
  UCL = valuePairs[j,1]
  return(list(UCL=UCL,alpha1=(1-countValues)))
}


SignalProbability = function (matrixT2, ucl, uclMax, uclKernel )
{
  
  SignalCountUCL = 0
  SignalCountMax = 0
  SignalCountKernel = 0
  
  for (i in 1: dim(matrixT2)[1])
  {
    RowMatrix = matrixT2[i,]
    
    ListFilterUCL = Filter(function(x) x > ucl,RowMatrix)
    SignalCountUCL = SignalCountUCL + length(ListFilterUCL)
    
    
    ListFilterMax = Filter(function(x) x > uclMax,RowMatrix)
    SignalCountMax = SignalCountMax + length(ListFilterMax)
    
    ListFilterKernel = Filter(function(x) x > uclKernel,RowMatrix)
    SignalCountKernel = SignalCountKernel + length(ListFilterKernel)
  }
  
  TotalValue = dim(matrixT2)[1]*dim(matrixT2)[2]
  
  SignalPro = SignalCountUCL/TotalValue
  SignalProMax = SignalCountMax/TotalValue
  SignalProKernel = SignalCountKernel/TotalValue
  
  return(list("SignalPro" = SignalPro,
              "SignalProMax" = SignalProMax,
              "SignalProKernel" = SignalProKernel))
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
    else if (typeMethod == "T2MOD")
    {
      Media = colMeans(Data)
      Sigma = cov(Data)
      SigmaMod = 1/(sum(diag(Sigma))/numVariables)
      
      T2 = c()
      if(identical(meanVector, shiftMean))
      {
        for (j in 1:dim(Data)[1])
        {
          Ai = (observation/(observation-1))*(norm((Data[j,]-Media)/sqrt(numVariables)^2, type = "2")/SigmaMod)
          T2 = c(Ai,T2)
        }
      }
      else 
      {
        for (j in 1:dim(DataOutlier)[1])
        {
          Ai = (observation/(observation-1))*(norm((DataOutlier[j,]-Media)/sqrt(numVariables)^2, type = "2")/SigmaMod)
          T2 = c(Ai,T2)
        }
      }
    }
    else if (typeMethod == "EBADIUI")
    {
      T2 = c()
      
      if(identical(meanVector, shiftMean))
      {
        Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05,c(), FALSE)
        T2 = c(Firstestimation$Ui,T2)
      }
      else 
      {
        Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05, DataOutlier, TRUE)
        T2 = c(Firstestimation$Ui,T2)
      }
      LenMatrix = observation
    }
    else if (typeMethod == "EBADIZI")
    {
      T2 = c()
      if(identical(meanVector, shiftMean))
      {
        Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05,c(), FALSE)
        T2 = c(Firstestimation$Zi,T2)
      }
      else 
      {
        Firstestimation = AlgorithmMDPCFPart1(Data, numVariables,Observation, 0.05, DataOutlier, TRUE)
        T2 = c(Firstestimation$Zi,T2)
      }
      LenMatrix = observation
    }
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  T2Matrix  = matrix(T2Total, ncol = 1)
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
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

# Parte 1 Calculo de los Limites de Probabilidad MRCD

# MRCD
ValuesMRCD = SimulationT2Chart(Observation,
                               NumberVariable,
                               NumSimulation,
                               mu,
                               SigmaCorr,
                               AlphaMRCD,
                               "MRCD")
KernelMethodMRCD = MethodUCLKernel(ValuesMRCD$T2Total,AlphaPFA)
UCLMRCD = qemp(p = 1-AlphaPFA, obs = ValuesMRCD$T2Total) 
UCLMaxMRCD = qemp(p = 1-AlphaPFA, obs = ValuesMRCD$T2Max) 
UCLKernelMRCD = KernelMethodMRCD$UCL

SignalProbabilityT2MRCD = SignalProbability(ValuesMRCD$T2Matrix, UCLMRCD,UCLMaxMRCD, UCLKernelMRCD)

#T2MOD

ValuesT2MOD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                mu,
                                SigmaCorr,
                                AlphaMRCD,
                                "T2MOD")
KernelMethodT2MOD = MethodUCLKernel(ValuesT2MOD$T2Total,AlphaPFA)
UCLT2MOD = qemp(p = 1-AlphaPFA, obs = ValuesT2MOD$T2Total) 
UCLMaxT2MOD = qemp(p = 1-AlphaPFA, obs = ValuesT2MOD$T2Max) 
UCLKernelT2MOD = KernelMethodT2MOD$UCL

SignalProbabilityT2MOD = SignalProbability(ValuesT2MOD$T2Matrix, UCLT2MOD,UCLMaxT2MOD, UCLKernelT2MOD)

# EBADIUI
ValuesEBADIUI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  mu,
                                  SigmaCorr,
                                  AlphaMRCD,
                                  "EBADIUI")
KernelMethodEBADIUI = MethodUCLKernel(ValuesEBADIUI$T2Total,AlphaPFA)
UCLEBADIUI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIUI$T2Total) 
UCLMaxEBADIUI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIUI$T2Max) 
UCLKernelEBADIUI = KernelMethodEBADIUI$UCL

SignalProbabilityT2EBADIUI = SignalProbability(ValuesEBADIUI$T2Matrix, UCLEBADIUI,UCLMaxEBADIUI, UCLKernelEBADIUI)

# EBADIZI
ValuesEBADIZI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  mu,
                                  SigmaCorr,
                                  AlphaMRCD,
                                  "EBADIZI")
KernelMethodEBADIZI = MethodUCLKernel(ValuesEBADIZI$T2Total,AlphaPFA)
UCLEBADIZI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIZI$T2Total) 
UCLMaxEBADIZI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIZI$T2Max) 
UCLKernelEBADIZI = KernelMethodEBADIZI$UCL

SignalProbabilityT2EBADIZI = SignalProbability(ValuesEBADIZI$T2Matrix, UCLEBADIZI,UCLMaxEBADIZI, UCLKernelEBADIZI)


# Parte 2 Calculo del parametro de no centralidad 
ArrayShiftmu = list(mu,
                    rep(1/100,NumberVariable),
                    rep(5/100,NumberVariable),
                    rep(10/100,NumberVariable),
                    rep(15/100,NumberVariable),
                    rep(25/100,NumberVariable),
                    rep(50/100,NumberVariable),
                    rep(60/100,NumberVariable),
                    rep(75/100,NumberVariable),
                    rep(90/100,NumberVariable),
                    rep(1,NumberVariable))

Inverse = solve(SigmaCorr)
Percentoutliers = 0.2
MatrixDeltaMRCD = matrix(,ncol = 4)
MatrixDeltaT2MOD = matrix(,ncol = 4)
MatrixDeltaEBADIUI = matrix(,ncol = 4)
MatrixDeltaEBADIZI = matrix(,ncol = 4)
for (shiftmu in ArrayShiftmu)
{
  DeltaNCP = sqrt(t(shiftmu-mu)%*%Inverse%*%(shiftmu-mu))
  
  #MRCD
  ValuesOutliersMRCD = SimulationT2ChartOutliers(Observation,
                                                 NumberVariable,
                                                 NumSimulationDelta,
                                                 mu,
                                                 SigmaCorr,
                                                 shiftmu,
                                                 Percentoutliers,
                                                 AlphaMRCD,
                                                 "MRCD")

  SignalProbabilityT2OutliersMRCD = SignalProbability(ValuesOutliersMRCD$T2Matrix, UCLMRCD,UCLMaxMRCD, UCLKernelMRCD)
  RowdeltaMRCD = c(DeltaNCP[1,1],
               SignalProbabilityT2OutliersMRCD$SignalPro,
               SignalProbabilityT2OutliersMRCD$SignalProMax,
               SignalProbabilityT2OutliersMRCD$SignalProKernel)
  MatrixDeltaMRCD = rbind(MatrixDeltaMRCD,RowdeltaMRCD)
  
  #T2MOD
  ValuesOutliersT2MOD = SimulationT2ChartOutliers(Observation,
                                                 NumberVariable,
                                                 NumSimulationDelta,
                                                 mu,
                                                 SigmaCorr,
                                                 shiftmu,
                                                 Percentoutliers,
                                                 AlphaMRCD,
                                                 "T2MOD")
  
  SignalProbabilityT2OutliersT2MOD = SignalProbability(ValuesOutliersT2MOD$T2Matrix, UCLT2MOD,UCLMaxT2MOD, UCLKernelT2MOD)
  RowdeltaT2MOD = c(DeltaNCP[1,1],
                   SignalProbabilityT2OutliersT2MOD$SignalPro,
                   SignalProbabilityT2OutliersT2MOD$SignalProMax,
                   SignalProbabilityT2OutliersT2MOD$SignalProKernel)
  MatrixDeltaT2MOD = rbind(MatrixDeltaT2MOD,RowdeltaT2MOD)
  
  #EBADIUI
  ValuesOutliersEBADIUI = SimulationT2ChartOutliers(Observation,
                                                    NumberVariable,
                                                    NumSimulationDelta,
                                                    mu,
                                                    SigmaCorr,
                                                    shiftmu,
                                                    Percentoutliers,
                                                    AlphaMRCD,
                                                    "EBADIUI")
  
  SignalProbabilityT2OutliersEBADIUI = SignalProbability(ValuesOutliersEBADIUI$T2Matrix, UCLEBADIUI,UCLMaxEBADIUI, UCLMaxEBADIUI)
  RowdeltaEBADIUI  = c(DeltaNCP[1,1],
                       SignalProbabilityT2OutliersEBADIUI$SignalPro,
                       SignalProbabilityT2OutliersEBADIUI$SignalProMax,
                       SignalProbabilityT2OutliersEBADIUI$SignalProKernel)
  MatrixDeltaEBADIUI = rbind(MatrixDeltaEBADIUI,RowdeltaEBADIUI)
  
  #EBADIZI
  ValuesOutliersEBADIZI = SimulationT2ChartOutliers(Observation,
                                                    NumberVariable,
                                                    NumSimulationDelta,
                                                    mu,
                                                    SigmaCorr,
                                                    shiftmu,
                                                    Percentoutliers,
                                                    AlphaMRCD,
                                                    "EBADIZI")
  
  SignalProbabilityT2OutliersEBADIZI = SignalProbability(ValuesOutliersEBADIZI$T2Matrix, UCLEBADIZI,UCLMaxEBADIZI, UCLMaxEBADIZI)
  RowdeltaEBADIZI  = c(DeltaNCP[1,1],
                       SignalProbabilityT2OutliersEBADIZI$SignalPro,
                       SignalProbabilityT2OutliersEBADIZI$SignalProMax,
                       SignalProbabilityT2OutliersEBADIZI$SignalProKernel)
  MatrixDeltaEBADIZI = rbind(MatrixDeltaEBADIZI,RowdeltaEBADIZI)
  
}
MatrixDeltaMRCD = MatrixDeltaMRCD[-1,]
MatrixDeltaT2MOD = MatrixDeltaT2MOD[-1,]
MatrixDeltaEBADIUI = MatrixDeltaEBADIUI[-1,]
MatrixDeltaEBADIZI = MatrixDeltaEBADIZI[-1,]

text <- c("MRCD","T2MOD","EBADIUI","EBADIZI")
plot_colors <- c("red","green","blue","orange")


library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/Imagen100x250x1000x0.05x20.svg"
svglite(urlplot, width = 8, height = 8)
par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)
lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
dev.off()


