library(MASS)
library(Rfast)
library(expm)
library(ggplot2)

AlgorithmRoMDP <- function(dataMatrix, itertime = 100) {
  dimension <- dim(dataMatrix)
  NumObservation <- dimension[1] 
  MiddlePoint <- round(NumObservation/2) + 1
  Id <- replicate( itertime, sample.int(NumObservation, 2) ) - 1
  FinalVec <- as.vector(.Call('Rfast_rmdp', PACKAGE = 'Rfast',dataMatrix,MiddlePoint,Id,itertime))
  Submcd <- seq(1, NumObservation)[FinalVec != 0]
  MDPMu <- Rfast::colmeans( dataMatrix[Submcd, ] ) 
  MDPVar <- Rfast::colVars( dataMatrix[Submcd, ] )
  Sigma <- cov(dataMatrix[Submcd, ])
  return(list(MDPMu,MDPVar,Sigma))
}

p <- 10
m <- 10000
h <- floor(m/2)+1
sigma <- diag(p)
sigma2 <- c()
for (i in 1:p){
  for (j in 1:p)
  {
    sigma2 <- c(sigma2,(0.5)^abs(i-j))
  }
}
porcent <- 0.9
sigma2 <- matrix(sigma2, ncol = p)

mu <- rep(0,p)
DatesNorm <- mvrnorm(m*porcent,mu,sigma)
DatesNorm <- rbind(DatesNorm,mvrnorm(m*(1-porcent),mu+1,sigma2))
a <-AlgorithmRoMDP(DatesNorm)

varmatrix <- as.numeric(unlist(a[2]))
muestime <-unlist(a[1])    
DMS <- matrix(diag(varmatrix),ncol = p)
DSquareinverse <- matrix(diag(sqrt(varmatrix^(-1))), ncol = p)
S <- a[3][[1]]
R <- DSquareinverse%*%S%*%DSquareinverse


trarhocuadra <- sum(diag(R %^% 2))- p^2/h
trarhoCuabico <- sum(diag(R %^% 3)) - (3*p/h)*sum(diag(R%^% 2))+2*(p^3/h^2)
zalpha <- qnorm(0.015,0,1)^2-1
cmp <- 1+(2*p)/(m*sqrt(trarhocuadra))
denominator <- sqrt(2*cmp*trarhocuadra)
listui <- c()
listzi <- c()
for( i  in 1:dim(DatesNorm)[1])
{
  DistanceMahalanobis <- t((DatesNorm[i,] - muestime))%*%DMS%*%(DatesNorm[i,] - muestime)
  
  ui <- (DistanceMahalanobis - p)/denominator
  zi <- ui - (4*trarhoCuabico*(zalpha)/(3*(2*trarhocuadra)^(3/2)))
  listui <- c(listui, ui)
  listzi <- c(listzi,zi)
}

listui <- sort(listui)
listzi <- sort(listzi)
f <- listui[listui>= 3.7]
f1 <- listzi[listzi>= 3.7]
pnorm1 <- pnorm(listui, mean = 0, sd = 1)
pnorm2 <- pnorm(listzi, mean = 0, sd = 1)
zx <- seq(min(listui[1],listzi[1]), max(listui[m],listzi[m]), length.out=m)
zpnorm <- pnorm(zx, mean = 0, sd = 1)

ggplot()+ geom_line(aes(x =listui, y = pnorm1,color="Red")) + geom_line(aes(x =listzi, y = pnorm2,color="green")) + geom_line(aes(x =zx, y = zpnorm))
ggplot()+ geom_line(aes(x =listzi, y = pnorm2,color="green"))
ggplot()+ geom_line(aes(x =zx, y = zpnorm))