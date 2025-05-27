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
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SimulationT2Chart.R"))
set.seed(123)

#Parameters Simulación
Observation = 100
NumberVariable = 250
NumSimulation = 10
SigmaCorr = randcorr(NumberVariable)
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
Percentoutliers = 0.1
MatrixDeltaMRCD = matrix(,ncol = 4)
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
}
MatrixDeltaMRCD = MatrixDeltaMRCD[-1,]

text <- c("MRCD","T2MOD","EBADIUI","EBADIZI")
plot_colors <- c("red","green","blue","orange")


library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/Imagen100x250x1000x0.05x20.svg"
svglite(urlplot, width = 8, height = 8)
par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
dev.off()