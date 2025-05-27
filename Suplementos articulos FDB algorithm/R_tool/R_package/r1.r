library(MASS)
library(DetMCD)
library(FDB)
library(pracma)
library(rrcov)
library(randcorr)
library(fds)
library(pls)



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

# Octane Data Example 485.930


data(octane)
octane <- octane[, -1]

CovMRCD <- CovMrcd(octane, alpha = 0.75)
MediaMRCD <- CovMRCD$center
SigmaMRCD <- CovMRCD$cov
BestSubset <- CovMRCD$best

DataBestSubset <- octane[BestSubset, ]
T2 <- mahalanobis(DataBestSubset, center = MediaMRCD, cov = SigmaMRCD)
T2Total<- mahalanobis(octane, center = MediaMRCD, cov = SigmaMRCD)
KernelMethodMRCD <- MethodUCLKernel(T2, 0.05)
x <- c(1:dim(octane)[1])

# library(svglite)
# urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/OctaneExample.svg"
# svglite(urlplot, width = 10, height = 8)
plot(x,T2Total, col="lightblue", pch=19, cex=2,)
abline(h  = KernelMethodMRCD$UCL, col = "red") 
text(x,T2Total, labels=rownames(octane),data=octane, cex=0.9, font=2)
# dev.off()

plot(x,T2Total,cex=2)
abline(h  = KernelMethodMRCD$UCL, col = "red") 
text(x,T2Total, labels=rownames(octane),data=octane, cex=0.9, font=2)


library(readr)
gasoline <- read_csv("/Users/kevin.pineda/Downloads/Suplementos articulos FDB algorithm/R_tool/R_package/gasolineData.csv")
gasoline[c(1,2,3,4,5,6),] <- mvrnorm(6, mu = colMeans(gasoline) + 2, Sigma = cov(gasoline))


CovMRCDG <- CovMrcd(gasoline, alpha = 0.9)
MediaMRCDG <- CovMRCDG$center
SigmaMRCDG <- CovMRCDG$cov
BestSubsetG <- CovMRCDG$best

DataBestSubsetG <- gasoline[BestSubsetG, ]

T2G <- mahalanobis(DataBestSubsetG, center = MediaMRCDG, cov = SigmaMRCDG)
T2TotalG<- mahalanobis(gasoline, center = MediaMRCDG, cov = SigmaMRCDG)
KernelMethodMRCDG <- MethodUCLKernel(T2G, 0.05)

xGasoline <- c(1:dim(gasoline)[1])

library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/GasolineExample.svg"
svglite(urlplot, width = 10, height = 8)
plot(xGasoline,T2TotalG, col="lightblue", pch=19, cex=2,)
abline(h  = KernelMethodMRCDG$UCL, col = "red") 
text(xGasoline,T2Total, labels=rownames(gasoline), cex=0.9, font=2)
dev.off()

plot()
abline(h  = KernelMethodMRCD$UCL, col = "red") 
