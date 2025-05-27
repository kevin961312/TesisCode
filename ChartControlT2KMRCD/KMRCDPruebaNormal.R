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

source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/MethodUCLKernel.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SignalProbability.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SimulationT2ChartNormal.R"))

set.seed(123)

#Parameters Simulación
Observation = 50
NumberVariable =250
NumSimulation = 200
fun <- function(i,j) (0.5)^(abs(i-j))

rows <- 1:NumberVariable
cols <- 1:NumberVariable

SigmaCorr = outer(rows,cols,FUN=fun)

# SigmaCorr = randcorr(NumberVariable)
mu = rep(0,NumberVariable)
AlphaMRCD = 0.75
AlphaPFA = 0.05
NumSimulationDelta = 200

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


# KMRCD
ValuesKMRCD = SimulationT2Chart(Observation,
                                NumberVariable,
                                NumSimulation,
                                mu,
                                SigmaCorr,
                                AlphaMRCD,
                                "KMRCD")
KernelMethodKMRCD = MethodUCLKernel(ValuesKMRCD$T2Total,AlphaPFA)
UCLKMRCD = qemp(p = 1-AlphaPFA, obs = ValuesKMRCD$T2Total) 
UCLMaxKMRCD = qemp(p = 1-AlphaPFA, obs = ValuesKMRCD$T2Max) 
UCLKernelKMRCD = KernelMethodKMRCD$UCL

SignalProbabilityT2KMRCD = SignalProbability(ValuesKMRCD$T2Matrix, UCLKMRCD,UCLMaxKMRCD, UCLKernelKMRCD)


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
Percentoutliers = 0.1
MatrixDeltaMRCD = matrix(,ncol = 4)
MatrixDeltaT2MOD = matrix(,ncol = 4)
MatrixDeltaEBADIUI = matrix(,ncol = 4)
MatrixDeltaEBADIZI = matrix(,ncol = 4)
MatrixDeltaKMRCD = matrix(,ncol = 4)
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
  
  SignalProbabilityT2OutliersEBADIUI = SignalProbability(ValuesOutliersEBADIUI$T2Matrix, UCLEBADIUI,UCLMaxEBADIUI, UCLKernelEBADIUI)
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
  
  SignalProbabilityT2OutliersEBADIZI = SignalProbability(ValuesOutliersEBADIZI$T2Matrix, UCLEBADIZI,UCLMaxEBADIZI, UCLKernelEBADIZI)
  RowdeltaEBADIZI  = c(DeltaNCP[1,1],
                       SignalProbabilityT2OutliersEBADIZI$SignalPro,
                       SignalProbabilityT2OutliersEBADIZI$SignalProMax,
                       SignalProbabilityT2OutliersEBADIZI$SignalProKernel)
  MatrixDeltaEBADIZI = rbind(MatrixDeltaEBADIZI,RowdeltaEBADIZI)

  #KMRCD
  ValuesOutliersKMRCD = SimulationT2ChartOutliers(Observation,
                                                  NumberVariable,
                                                  NumSimulationDelta,
                                                  mu,
                                                  SigmaCorr,
                                                  shiftmu,
                                                  Percentoutliers,
                                                  AlphaMRCD,
                                                  "KMRCD")
  
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


library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/Imagen100x250x1000x0.05x20.svg"
svglite(urlplot, width = 8, height = 8)
par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)
lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
lines(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], col='purple',pch = 15,type='b',lty = 6)

text <- c("MRCD","T2MOD","EBADIUI","EBADIZI")
plot_colors <- c("red","green","blue","orange")



lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
dev.off()