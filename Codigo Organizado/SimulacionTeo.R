library(MASS)
library(IQCC)
library(moments)
library(mvtnorm)

ListTIndependent <- c()
ListNormal <- c()

NumSubgroup <- 15
NumSample <- 10
ListNormalNumbers <- c()
ListSDandVar <- c()
for(i in 1:NumSubgroup)
{
  Normalnumbers <- rnorm(n = NumSample, mean = 0, sd = 1)
  ListSDandVar <- c(sd(Normalnumbers),ListSDandVar)
  ListSDandVar <- c(var(Normalnumbers),ListSDandVar)
  ListSDandVar <- c((max(Normalnumbers)-min(Normalnumbers))/4,ListSDandVar )
  ListNormalNumbers <- c(Normalnumbers,ListNormalNumbers)
}
MatrizSamples <- matrix(ListNormalNumbers,NumSubgroup,NumSample)
MatrizSamplesSD <- matrix(ListSDandVar,3,NumSubgroup)

OverallMean <- mean(rowMeans(MatrizSamples))
#d2 y C4 constante esta en el articulo es de la tabla
d2 <-  1/d2(NumSample)
c4 <-  1/c4(NumSample)
SDStandar <- mean(MatrizSamplesSD[3,])*c4
SDRange <- mean(MatrizSamplesSD[1,])*d2
SDVariance <- sqrt(mean(MatrizSamplesSD[2,]))*c4

# rangos limites
CR <- 1.00223
KR <- 2.8792 
LR <- CR*KR*sqrt(NumSubgroup/(NumSubgroup-1))

# estandar limites
CS <- 1.00190
KS <- 2.8586 
LS <- CS*KS*sqrt(NumSubgroup/(NumSubgroup-1))

# varianza limites
KV <- 2.8575
LV <- KV*sqrt(NumSubgroup/(NumSubgroup-1))
# FAP = 0.05
#calculo probabilidad normal
PNR <- 1- (pnorm(LR, 0, 1) - pnorm(-LR, 0, 1))^NumSubgroup;PNR
PNS <- 1- (pnorm(LS, 0, 1) - pnorm(-LS, 0, 1))^NumSubgroup;PNS
PNV <- 1- (pnorm(LV, 0, 1) - pnorm(-LV, 0, 1))^NumSubgroup;PNV

# calculo probabilidad t independiente
PTIR <- 1- (pt(LR, df = 112.059, lower.tail = TRUE) - pt(-LR, df = 112.059, lower.tail = TRUE))^NumSubgroup;PTIR
PTIS <- 1- (pt(LS, df = 131.80759, lower.tail = TRUE) - pt(-LS, df = 131.80759, lower.tail = TRUE))^NumSubgroup;PTIS
PTIV <- 1- (pt(LV, df = 135, lower.tail = TRUE) - pt(-LV, df = 135, lower.tail = TRUE))^NumSubgroup;PTIV

# Calculo probabilidad t multivariada
MatrizCorrelationT <- matrix(0,NumSubgroup,NumSubgroup)
for (i in 1:NumSubgroup)
{
  for (j in 1:NumSubgroup)
  {
    if( i == j)
    {
      MatrizCorrelationT[i,j]=1
    }
    else
    {
      MatrizCorrelationT[i,j]=-1/(NumSubgroup-1)
    }
  }
}

PTMR <- 1-pmvt(lower = rep(-LR, NumSubgroup) , upper = rep(LR, NumSubgroup) , corr = MatrizCorrelationT, df = 113)[1];PTMR
PTMS <- 1-pmvt(lower = rep(-LS, NumSubgroup) , upper = rep(LS, NumSubgroup) , corr = MatrizCorrelationT, df = 132)[1];PTMS
PTMV <- 1-pmvt(lower = rep(-LV, NumSubgroup) , upper = rep(LV, NumSubgroup) , corr = MatrizCorrelationT, df = 135)[1];PTMV