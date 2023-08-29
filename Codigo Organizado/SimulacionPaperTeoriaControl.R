library(MASS)
library(IQCC)
library(mvtnorm)




simulationCSheSmall <- function(NumberSimulatio, delta = 0 ) {
  NumSubgroup <- 15
  NumObservation <- 10
  ListPFA <- c()
  for (h in 1:NumberSimulatio)
  {
    ListNormalNumbers <- c()
    ListSDandVar <- c()
    for(i in 1:NumSubgroup)
    {
      if( i>10)
      {
        Normalnumbers <- rnorm(n = NumObservation, mean = 0, sd = 1)
      }
      else
      {
        Normalnumbers <- rnorm(n = NumObservation, mean = 0+1*delta, sd = 1)
      }
      
      ListSDandVar <- c(var(Normalnumbers),ListSDandVar)
      ListNormalNumbers <- c(Normalnumbers,ListNormalNumbers)
    }
    MatrizSamples <- matrix(ListNormalNumbers,NumSubgroup,NumObservation)
    MatrizSamplesSD <- matrix(ListSDandVar,1,NumSubgroup)
    
    OverallMean <- mean(rowMeans(MatrizSamples))
    #d2 y C4 constante esta en el articulo es de la tabla
    SDVariance <- sqrt(mean(MatrizSamplesSD[1,]))/c4(NumObservation)
    
    # varianza limites
    KVC <- 2.8575*1.00190
    LV <- KVC*sqrt(NumSubgroup/(NumSubgroup-1))
    # FAP = 0.05
    # calculo probabilidad t independiente
    PTIV <- 1- (pt(LV, df = 135, lower.tail = TRUE) - pt(-LV, df = 135, lower.tail = TRUE))^NumSubgroup;
    CountP = 0 
    for (j in 1:NumSubgroup)
    {
      Date <- Filter(function(x) x > LV | -LV > x, MatrizSamples[j,])
      if(length(Date) > 0 )
      {
        CountP = CountP + 1
      }
    }
    
    PFA <- CountP/NumSubgroup
    ListPFA  <- c(ListPFA, PFA) 
  }
  listfilter <- Filter(function(x) x > 0,ListPFA)
  PFAOverall <- mean(listfilter)
  return(PFAOverall)
}


PFAWithoutOutlier <- simulationCSheSmall(10000)
PFAWithOutlier1 <- simulationCSheSmall(10000,0.1)
PFAWithOutlier2 <- simulationCSheSmall(10000,0.5)
PFAWithOutlier3 <- simulationCSheSmall(10000,1)
PFAWithOutlier4 <- simulationCSheSmall(10000,2)
PFAWithOutlier5 <- simulationCSheSmall(10000,5)