#Simeulación de datos:
library(MASS)
mu<-c(28.1,7.18,3.089)
Sigma<- matrix(c(140.54,49.68,1.94,49.68,72.25,3.68,1.94,3.68,0.25),ncol=3,byrow=T)
X <- mvrnorm(50,mu=mu,Sigma=Sigma)
alpha <- 0.05

n<-nrow(X)
p<-ncol(X)

UCL<-(((n-1)^2)/n)*qbeta((1-alpha),(p/2),(n-p-1)/2)

Xmedia<-apply(X,2,mean)

Su<-var(X)

T2<-mahalanobis(X,center=Xmedia,cov=Su)

Observacion<-1:n

par(bg="cornsilk")

plot(Observacion,
     T2,
     type='l',
     xlim=c(0,n+2),
     ylim=c(0,max(UCL,max(T2))+2),
     main=expression("Carta"*~T^2),
     ylab=expression(T^2),xlab="No. Observación",font=2)
abline(h=UCL,lty=3)
for(i in 1:n)
{
  temp<-ifelse((T2[i]>UCL),4,19)
  points(Observacion[i],T2[i],pch=temp)
  if(T2[i]>UCL)text(i,T2[i],labels=paste('obs=',i),pos=3,font=2,ex=0.7) 
}
text((max(Observacion)-1),
     UCL,
     paste('UCL=',round(UCL,digits=4)),
     pos=3,
     font=2,cex=0.7)
legend(locator(1),c(paste("p=",p),paste("alpha=",alpha),
                    paste("n=",n)),
       ncol=3,cex=0.7,bg='gray95')
(estimaciones<-list(medias=Xmedia,var=Su,T2=T2))


