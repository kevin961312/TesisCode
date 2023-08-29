library(KernSmooth)
UCL.kernel <- function(T2,alpha)
{
  est <- density(T2,kernel="gaussian",bw="nrd")
  aux <- cbind(est$x,est$y)
  aux <- aux[aux[,1]>=0,]
  estand <- aux[,2]/sum(aux[,2])
  n1 <- nrow(aux)
  acum <- 0
  for(i in 1:n1)
  {
    if(acum<(1-alpha))
    {
      j <- i 
      acum <- sum(estand[1:i])
    }
  }
  
  yhist <- hist(T2,plot=FALSE); yhist<-yhist$density
  x.lim <- c(min(est$x),max(est$x,T2))
  y.lim<-c(0,max(est$y,yhist))
  hist(T2,freq=FALSE,xlim=x.lim,ylim=y.lim,
       main=expression(paste("Histograma de",sep=" ",T^2,sep=" ",
                             "con suavizamiento kernel")),xlab=expression(T^2),
       ylab="Densidad")
  par(new=T)
  plot(est,  main="", xlab="", ylab="",  xlim=x.lim, ylim=y.lim)
  UCL <- aux[j,1]
  
  res<-unlist(list(UCL=UCL,alpha1=(1-acum)))
}
T2 <- c(35.32785, 34.22203, 33.76553, 36.01239, 34.03818, 33.67729, 36.11540, 35.35763, 34.45803, 35.60423,
         33.36493, 35.94589, 36.21965, 36.69409, 36.74051, 35.31853, 35.24346, 35.09185, 35.63372, 34.53380,
         36.08184, 35.56808, 35.38254, 34.22530, 34.19336, 32.70254, 35.13262, 34.64065, 34.28021, 35.14733,
         35.17441, 34.46997, 34.57148, 34.75746, 35.38855, 35.02812, 32.25758, 35.71119, 33.99364, 34.32525,
         33.62640, 37.08969, 35.08300, 33.53511, 34.99878, 36.59812, 34.48782, 35.07952, 34.79494, 34.25979)

EstKernelSmooth <- density(T2,kernel="gaussian",bw="nrd")
valuePairs <- cbind(EstKernelSmooth$x,EstKernelSmooth$y)
valuePairs <- valuePairs[valuePairs[,1]>=0,]
Estand <- valuePairs[,2]/sum(valuePairs[,2])







