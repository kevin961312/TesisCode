Data_VDP_A6 <- read.delim("~/Downloads/Data_VDP_A6.txt")
BoardData <- unique(Data_VDP_A6$Board)
plot(NULL, xlim=c(min(Data_VDP_A6$Depth),max(Data_VDP_A6$Depth)), ylim=c(min(Data_VDP_A6$VDP),max(Data_VDP_A6$VDP)), ylab="y label", xlab="x lablel")

MatrixValuesY <- matrix(, ncol = 314)
for (i in BoardData) {
  datafilter <- Data_VDP_A6[Data_VDP_A6$Board == i,]
  MatrixValuesY <- rbind(MatrixValuesY, datafilter$VDP)
  lines(datafilter$Depth,datafilter$VDP)
}

MatrixValuesY <- MatrixValuesY[-1,]


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
library(foreach)
library(doParallel)
library(SixSigma)
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCDParalelo/kMRCD.R"))

Data <- t(ss.data.wby)
CovKMRCD <- CovKmrcd(x = Data, kModel = "RbfKernel",alpha = 0.8)
complemento <- setdiff(1:dim(Data)[1], CovKMRCD$hsubsetIndices)
DataBestSubset <- Data[CovKMRCD$hsubsetIndices, ]
DataComplement <- Data[complemento, ]

plot(NULL, xlim=c(min(ss.data.wbx),max(ss.data.wbx)), ylim=c(min(ss.data.wby),max(ss.data.wby)), ylab="y label", xlab="x lablel")
for (i in rownames(DataBestSubset)) {
  lines(ss.data.wbx, DataBestSubset[i,])
}

for (i in rownames(DataComplement)) {
  lines(ss.data.wbx, DataComplement[i,], col='red')
}


## con MRCD
CovMRCD <- CovMrcd(Data,alpha = 0.8)

complementoMRCD <- setdiff(1:dim(Data)[1], CovMRCD$best)
DataBestSubsetMRCD <- Data[CovMRCD$best, ]
DataComplementMRCD <- Data[complementoMRCD, ]

plot(NULL, xlim=c(min(ss.data.wbx),max(ss.data.wbx)), ylim=c(min(ss.data.wby),max(ss.data.wby)), ylab="y label", xlab="x lablel")
for (i in rownames(DataBestSubsetMRCD)) {
  lines(ss.data.wbx, DataBestSubsetMRCD[i,])
}

for (i in rownames(DataComplementMRCD)) {
  lines(ss.data.wbx, DataComplementMRCD[i,], col='red')
}





