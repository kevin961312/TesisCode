
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/MethodUCLKernel.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SignalProbability.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/SimulationT2Chart.R"))
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/kMRCD.R"))
library(rdetools)
set.seed(123)


MatrixDeltaKMRCD = MatrixDeltaKMRCD[-1,]

par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)


DataCopula <- rCopula(60,tCopula(dim=2,0,df=1))
mu <- colMeans(DataCopula)
Sigma <- cov(DataCopula)
CovKMRCD <- CovKmrcd(x = DataCopula, kModel = "RbfKernel", alpha = 0.75)


CovKMRCD1 <- CovKmrcd(x = DataCopula, kModel = "RbfKernel", alpha = 0.75)
DataBestSubset <- DataCopula[CovKMRCD1$hsubsetIndices, ]
meanKMRCD1 <- colMeans(DataBestSubset)
T2 <- mahalanobis(DataBestSubset, center = meanKMRCD1, cov = CovKMRCD1$cov)
CovKMRCD$smd[CovKMRCD$hsubsetIndices]














