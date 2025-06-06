
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
library(lcmix)

source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/MethodUCLKernel.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SignalProbability.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SimulationT2ChartGamma.R"))
set.seed(123)

#Parameters Simulación
Observation = 100
NumberVariable = 200
NumSimulation = 10
rho = 2
beta = 5
AlphaMRCD = 0.85
AlphaPFA = 0.05
NumSimulationDelta = 10

fun <- function(i,j) (0.5)^(abs(i-j))
rows <- 1:NumberVariable
cols <- 1:NumberVariable

SigmaCorr = outer(rows,cols,FUN=fun)

# Parte 1 Calculo de los Limites de Probabilidad MRCD

# MRCD
ValuesMRCD = SimulationT2Chart(Observation,
                               NumberVariable,
                               NumSimulation,
                               SigmaCorr,
                               rho,
                               AlphaMRCD,
                               "MRCD",beta)
KernelMethodMRCD = MethodUCLKernel(ValuesMRCD$T2Total,AlphaPFA)
UCLMRCD = qemp(p = 1-AlphaPFA, obs = ValuesMRCD$T2Total) 
UCLMaxMRCD = qemp(p = 1-AlphaPFA, obs = ValuesMRCD$T2Max) 
UCLKernelMRCD = KernelMethodMRCD$UCL

SignalProbabilityT2MRCD = SignalProbability(ValuesMRCD$T2Matrix, UCLMRCD,UCLMaxMRCD, UCLKernelMRCD)

#T2MOD

ValuesT2MOD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                SigmaCorr,
                                rho,
                                AlphaMRCD,
                                "T2MOD",beta)
KernelMethodT2MOD = MethodUCLKernel(ValuesT2MOD$T2Total,AlphaPFA)
UCLT2MOD = qemp(p = 1-AlphaPFA, obs = ValuesT2MOD$T2Total) 
UCLMaxT2MOD = qemp(p = 1-AlphaPFA, obs = ValuesT2MOD$T2Max) 
UCLKernelT2MOD = KernelMethodT2MOD$UCL

SignalProbabilityT2MOD = SignalProbability(ValuesT2MOD$T2Matrix, UCLT2MOD,UCLMaxT2MOD, UCLKernelT2MOD)

# EBADIUI
ValuesEBADIUI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  SigmaCorr,
                                  rho,
                                  AlphaMRCD,
                                  "EBADIUI",beta)
KernelMethodEBADIUI = MethodUCLKernel(ValuesEBADIUI$T2Total,AlphaPFA)
UCLEBADIUI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIUI$T2Total) 
UCLMaxEBADIUI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIUI$T2Max) 
UCLKernelEBADIUI = KernelMethodEBADIUI$UCL

SignalProbabilityT2EBADIUI = SignalProbability(ValuesEBADIUI$T2Matrix, UCLEBADIUI,UCLMaxEBADIUI, UCLKernelEBADIUI)

# EBADIZI
ValuesEBADIZI = SimulationT2Chart(Observation,
                                  NumberVariable,
                                  NumSimulation,
                                  SigmaCorr,
                                  rho,
                                  AlphaMRCD,
                                  "EBADIZI",beta)
KernelMethodEBADIZI = MethodUCLKernel(ValuesEBADIZI$T2Total,AlphaPFA)
UCLEBADIZI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIZI$T2Total) 
UCLMaxEBADIZI = qemp(p = 1-AlphaPFA, obs = ValuesEBADIZI$T2Max) 
UCLKernelEBADIZI = KernelMethodEBADIZI$UCL

SignalProbabilityT2EBADIZI = SignalProbability(ValuesEBADIZI$T2Matrix, UCLEBADIZI,UCLMaxEBADIZI, UCLKernelEBADIZI)


# KMRCD
ValuesKMRCD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                SigmaCorr,
                                rho,
                                AlphaMRCD,
                                "KMRCD",beta)
KernelMethodKMRCD = MethodUCLKernel(ValuesKMRCD$T2Total,AlphaPFA)
UCLKMRCD = qemp(p = 1-AlphaPFA, obs = ValuesKMRCD$T2Total) 
UCLMaxKMRCD = qemp(p = 1-AlphaPFA, obs = ValuesKMRCD$T2Max) 
UCLKernelKMRCD = KernelMethodKMRCD$UCL

SignalProbabilityT2KMRCD = SignalProbability(ValuesKMRCD$T2Matrix, UCLKMRCD,UCLMaxKMRCD, UCLKernelKMRCD)


# Parte 2 Calculo del parametro de no centralidad 
ArrayRhoShift = list(0,
                     1/100,
                     5/100,
                     10/100,
                     15/100,
                     25/100,
                     50/100,
                     60/100,
                     75/100,
                     90/100,
                     1)

Inverse = solve(SigmaCorr)
Percentoutliers = 0.05
MatrixDeltaMRCD = matrix(,ncol = 4)
MatrixDeltaT2MOD = matrix(,ncol = 4)
MatrixDeltaEBADIUI = matrix(,ncol = 4)
MatrixDeltaEBADIZI = matrix(,ncol = 4)
MatrixDeltaKMRCD = matrix(,ncol = 4)
for (rhoShift in ArrayRhoShift)
{
  DeltaNCP = sqrt(t(rep(rhoShift,NumberVariable))%*%Inverse%*%(rep(rhoShift,NumberVariable)))
  
  #MRCD
  ValuesOutliersMRCD = SimulationT2ChartOutliers(Observation,
                                                 NumberVariable,
                                                 NumSimulationDelta,
                                                 SigmaCorr,
                                                 rho,
                                                 rhoShift,
                                                 Percentoutliers,
                                                 AlphaMRCD,
                                                 "MRCD",beta)
  
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
                                                  SigmaCorr,
                                                  rho,
                                                  rhoShift,
                                                  Percentoutliers,
                                                  AlphaMRCD,
                                                  "T2MOD",beta)
  
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
                                                    SigmaCorr,
                                                    rho,
                                                    rhoShift,
                                                    Percentoutliers,
                                                    AlphaMRCD,
                                                    "EBADIUI",beta)
  
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
                                                    SigmaCorr,
                                                    rho,
                                                    rhoShift,
                                                    Percentoutliers,
                                                    AlphaMRCD,
                                                    "EBADIZI",beta)
  
  SignalProbabilityT2OutliersEBADIZI = SignalProbability(ValuesOutliersEBADIZI$T2Matrix, UCLEBADIZI,UCLMaxEBADIZI, UCLMaxEBADIZI)
  RowdeltaEBADIZI  = c(DeltaNCP[1,1],
                       SignalProbabilityT2OutliersEBADIZI$SignalPro,
                       SignalProbabilityT2OutliersEBADIZI$SignalProMax,
                       SignalProbabilityT2OutliersEBADIZI$SignalProKernel)
  MatrixDeltaEBADIZI = rbind(MatrixDeltaEBADIZI,RowdeltaEBADIZI)
  #KMRCD
  ValuesOutliersKMRCD = SimulationT2ChartOutliers(Observation,
                                                  NumberVariable,
                                                  NumSimulationDelta,
                                                  SigmaCorr,
                                                  rho,
                                                  rhoShift,
                                                  Percentoutliers,
                                                  AlphaMRCD,
                                                  "KMRCD",beta)
  
  SignalProbabilityT2OutliersKMRCD = SignalProbability(ValuesOutliersKMRCD$T2Matrix, UCLKMRCD,UCLMaxKMRCD, UCLMaxKMRCD)
  RowdeltaKMRCD  = c(DeltaNCP[1,1],
                     SignalProbabilityT2OutliersKMRCD$SignalPro,
                     SignalProbabilityT2OutliersKMRCD$SignalProMax,
                     SignalProbabilityT2OutliersKMRCD$SignalProKernel)
  MatrixDeltaKMRCD = rbind(MatrixDeltaKMRCD,RowdeltaKMRCD)
  
}
MatrixDeltaMRCD = MatrixDeltaMRCD[-1,]
MatrixDeltaT2MOD = MatrixDeltaT2MOD[-1,]
MatrixDeltaEBADIUI = MatrixDeltaEBADIUI[-1,]
MatrixDeltaEBADIZI = MatrixDeltaEBADIZI[-1,]
MatrixDeltaKMRCD = MatrixDeltaKMRCD[-1,]


p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)
lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
lines(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], col='purple',pch = 15,type='b',lty = 6)

