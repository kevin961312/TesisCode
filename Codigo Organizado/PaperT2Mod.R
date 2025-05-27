rm(list=ls())

library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
library(randcorr)
library(KernSmooth)
library(psych)

set.seed(123)

SimulationT2Mod = function (observation, numVariables, numSimulation, meanVector, sigmaMatriz, alphaMRCD = 0.75 )
{
  
  T2Matrix = c()
  T2Total = c()
  T2Max = c()
  
  for (i in 1:numSimulation)
  {
    
    Data = mvrnorm(observation,mu=mu,Sigma = sigmaMatriz)
    
    
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  T2Matrix  = matrix(T2Total, ncol = observation)
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

SimulationT2ModOutliers = function (observation, numVariables, 
                                          numSimulation, meanVector, 
                                          sigmaMatriz, shiftMean, percentoutliers = 0, alphaMRCD = 0.75 )
{
  
  T2Matrix = c()
  T2Total = c()
  T2Max = c()
  
  for (i in 1:numSimulation)
  {
    
    NumOutliers = floor(percentoutliers*observation)
    Data = mvrnorm(Observation-NumOutliers,mu=mu,Sigma = sigmaMatriz)
    if(NumOutliers >= 1)
    {
      DataOutlier = mvrnorm(NumOutliers, mu = shiftMean, Sigma = sigmaMatriz)
      Data = rbind(Data, DataOutlier)
    }
    
    
    Media = colMeans(Data)
    Sigma = cov(Data)
    SigmaMod = 1/(sum(diag(Sigma))/numVariables)

    T2 = c()
    for (j in 1:dim(DataOutlier)[1])
    {
      Ai = (observation/(observation-1))*(norm((DataOutlier[j,]-Media)/sqrt(numVariables)^2, type = "2")/SigmaMod)
      T2 = c(Ai,T2)
    }
    T2Total = c(T2Total,T2)
    T2Max = c(T2Max,max(T2))
  }
  T2Matrix  = matrix(T2Total, ncol = 1)
  return(list("T2Matrix" = T2Matrix,
              "T2Total" = T2Total,
              "T2Max" = T2Max))
}




# Parte 1 Calculo de los Limites de Probabilidad
Observation = 50
NumberVariable = 100
NumSimulation = 100
SigmaCorr = randcorr(NumberVariable)
mu = rep(0,NumberVariable)

Values = SimulationT2Mod(Observation,
                               NumberVariable,
                               NumSimulation,
                               mu,
                               SigmaCorr,
                               0.75)

#alpha = 0.005  , 200
#alpha  = 0.0027, 370
#alpha = 0.05 , 20
AlphaPFA = 0.05

KernelMethod = MethodUCLKernel(Values$T2Total,AlphaPFA)
UCL = qemp(p = 1-AlphaPFA, obs = Values$T2Total) 
UCLMax = qemp(p = 1-AlphaPFA, obs = Values$T2Max) 
UCLKernel = KernelMethod$UCL

SignalProbabilityT2 = SignalProbability(Values$T2Matrix, UCL,UCLMax, UCLKernel)

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
Percentoutliers = 0.25
MatrixDelta = matrix(,ncol = 4)
for (shiftmu in ArrayShiftmu)
{
  DeltaNCP = sqrt(t(shiftmu-mu)%*%Inverse%*%(shiftmu-mu))
  
  ValuesOutliers = SimulationT2ModOutliers(Observation,
                                                 NumberVariable,
                                                 NumSimulation,
                                                 mu,
                                                 SigmaCorr,
                                                 shiftmu,
                                                 Percentoutliers,
                                                 0.75)
  
  SignalProbabilityT2Outliers = SignalProbability(ValuesOutliers$T2Matrix, UCL,UCLMax, UCLKernel)
  Rowdelta = c( DeltaNCP[1,1], 
                SignalProbabilityT2Outliers$SignalPro, 
                SignalProbabilityT2Outliers$SignalProMax, 
                SignalProbabilityT2Outliers$SignalProKernel)
  MatrixDelta = rbind(MatrixDelta,Rowdelta)
}
MatrixDelta =  MatrixDelta[-1,]

plot(MatrixDelta[,1],MatrixDelta[,2])
plot(MatrixDelta[,1],MatrixDelta[,3])
plot(MatrixDelta[,1],MatrixDelta[,4])
