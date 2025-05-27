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
library(readr)
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCDParalelo/kMRCD.R"))
MethodUCLKernel <- function(T2, alpha) {
  EstKernelSmooth <- density(T2, kernel = "gaussian", bw = "nrd")
  valuePairs <- cbind(EstKernelSmooth$x, EstKernelSmooth$y)
  valuePairs <- valuePairs[valuePairs[, 1] >= 0, ]
  Estand <- valuePairs[, 2] / sum(valuePairs[, 2])
  NumberRow <- nrow(valuePairs)
  countValues <- 0
  
  for (i in 1:NumberRow)
  {
    if (countValues < (1 - alpha)) {
      j <- i
      countValues <- sum(Estand[1:i])
    }
  }
  
  UCL <- valuePairs[j, 1]
  return(list(UCL = UCL, alpha1 = (1 - countValues)))
}


VDP_DATA <- read_table("Downloads/Data_VDP_A6.txt")
Categori <- unique(VDP_DATA$Board)
VDP_DATA <- VDP_DATA[order(VDP_DATA$BoardNumber,VDP_DATA$Depth), ]

Data <- matrix(, nrow = 0, ncol = 314)
for ( i in Categori)
{
  VDP_DATA_Categori <- VDP_DATA[VDP_DATA$Board==i,]
  Data <- rbind(Data, VDP_DATA_Categori$VDP)
}

x <- c()
for ( i in c("A1"))
{
  VDP_DATA_Categori <- VDP_DATA[VDP_DATA$Board==i,]
  ValuesProfile <- VDP_DATA_Categori$Depth
  x <- c(x, ValuesProfile)
}

## con KMRCD
CovKMRCD <- CovKmrcd(Data, kModel = "RbfKernel",0.9)
Hsubinidice <- CovKMRCD$hsubsetIndices
complementoKMRCD <- setdiff(1:dim(Data)[1], CovKMRCD$hsubsetIndices)
DataBestSubsetKMRCD <- Data[CovKMRCD$hsubsetIndices, ]
DataComplementKMRCD <- Data[complementoKMRCD, ]
MediaKMRCD <- colmeans(DataBestSubsetKMRCD)
SigmaKMRCD <- CovKMRCD$cov

library(svglite)
urlplot ="/Users/kevin.pineda/Downloads/Tesis/TesisCode/Imagenes TG 2 KMRCD/ExampleData/VDPKMRCDPerfiles.svg"
svglite(urlplot, width = 10, height = 8)

plot(NULL, xlim=c(min(x),max(x)), ylim=c(min(Data),max(Data)), ylab="Density", xlab="Depth")
for (j in 1:dim(DataBestSubsetKMRCD)[1]) {
  lines(x, DataBestSubsetKMRCD[j,])
}

for (j in 1:dim(DataComplementKMRCD)[1]) {
  lines(x, DataComplementKMRCD[j,], col='red')
}
dev.off()

T2 <- mahalanobis(DataBestSubsetKMRCD, center = MediaKMRCD, cov = SigmaKMRCD)
T2Total<- mahalanobis(Data, center = MediaKMRCD, cov = SigmaKMRCD)
KernelMethodKMRCD <- MethodUCLKernel(T2, 0.01)
x <- c(1:dim(Data)[1])


# library(svglite)
# urlplot ="/Users/kevin.pineda/Downloads/Tesis/TesisCode/Imagenes TG 2 KMRCD/ExampleData/VDPKMRCD.svg"
# svglite(urlplot, width = 10, height = 8)
plot(x,T2Total, col="lightblue", pch=19, cex=2,)
abline(h  = KernelMethodKMRCD$UCL, col = "red") 
# dev.off()
