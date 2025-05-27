
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


source(here::here("Downloads/Tesis/TesisCode/ChartControlT2MRCD/MethodUCLKernel.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2MRCD/SignalProbability.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2MRCD/SimulationT2Chart2.R"))
set.seed(123)

#Parameters Simulación
Observation = 150
NumberVariable = 250
NumSimulation = 10

fun <- function(i,j) (0.5)^(abs(i-j))

rows <- 1:NumberVariable
cols <- 1:NumberVariable

SigmaCorr = outer(rows,cols,FUN=fun)
# SigmaCorr = randcorr(NumberVariable)
mu = rep(0,NumberVariable)
AlphaMRCD = 0.75
AlphaPFA = 0.05
NumSimulationDelta = 10

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
Percentoutliers = 0.05
MatrixDeltaMRCD = matrix(,ncol = 4)
MatrixDeltaT2MOD = matrix(,ncol = 4)
MatrixDeltaEBADIUI = matrix(,ncol = 4)
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
                                                 "MRCD2")
  
  SignalProbabilityT2OutliersT2MOD =  SignalProbability(ValuesOutliersT2MOD$T2Matrix, UCLMRCD,UCLMaxMRCD, UCLKernelMRCD)
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
                                                    "MRCD3")
  
  SignalProbabilityT2OutliersEBADIUI =  SignalProbability(ValuesOutliersEBADIUI$T2Matrix, UCLMRCD,UCLMaxMRCD, UCLKernelMRCD)
  RowdeltaEBADIUI  = c(DeltaNCP[1,1],
                       SignalProbabilityT2OutliersEBADIUI$SignalPro,
                       SignalProbabilityT2OutliersEBADIUI$SignalProMax,
                       SignalProbabilityT2OutliersEBADIUI$SignalProKernel)
  MatrixDeltaEBADIUI = rbind(MatrixDeltaEBADIUI,RowdeltaEBADIUI)
  
}
MatrixDeltaMRCD = MatrixDeltaMRCD[-1,]
MatrixDeltaT2MOD = MatrixDeltaT2MOD[-1,]
MatrixDeltaEBADIUI = MatrixDeltaEBADIUI[-1,]

par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)

