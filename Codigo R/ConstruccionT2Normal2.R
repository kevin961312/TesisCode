library(MASS)
library(clusterGeneration)
library(rrcov)
library(EnvStats)
Observation <- 60
NumberVariable <- 100
mu <- rep(0,NumberVariable)
DPMatriz <- genPositiveDefMat(dim = NumberVariable, covMethod = "unifcorrmat")
Sigma <- DPMatriz$Sigma
T2General <- c()
for (i in 1:120)
{
  X <- mvrnorm(Observation,mu=mu,Sigma=Sigma)
  c <- CovMrcd(X, alpha=0.75)
  Xmedia <-c$center
  Su <- c$cov
  h <- c$best
  X1 <- X[h,] 
  T2<-mahalanobis(X1,center=Xmedia,cov=Su)
  T2General <- c(T2General,max(T2))
}

plot(ecdf(T2General))
UCL <- qemp(p = 0.95, obs = T2General) 
#1/alpha=200 o 370

Observacion<-1:Observation

par(bg="cornsilk")

plot(1:length(h),
     T2,
     type='l',
     xlim=c(0,length(h)+2),
     main=expression("Carta"*~T^2),
     ylab=expression(T^2),xlab="No. Observacion",font=2)
abline(h=UCL,lty=3)

epdfPlot(T2General, epdf.col = "cyan", 
         xlab = "Value of Random Variable", 
         main = "Empirical and Theoretical PDFs",xlim = c(30,50))
