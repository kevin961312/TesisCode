library(MASS)
library(Rfast)
library(expm)
library(ggplot2)
library(psych)
library(spatstat)
# Funcion Algoritmo de Ro-2015
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


algorithmUiZi <- function(DatesNorm, NumVariable,Observations, Alpha, OutliersFlag = FALSE)
{
  # Usar Algoritmo de Ro-2015
  AlgorithmObjects <- AlgorithmRoMDP(DatesNorm)
  MeanDMP <-unlist(AlgorithmObjects[1])    
  SigmaDMP <- matrix(diag(as.numeric(unlist(AlgorithmObjects[2]))),ncol = NumVariable)
  InverseSigmaDMP <- solve(SigmaDMP)
  SigmaOfMinimunDet <- AlgorithmObjects[3][[1]]
  CorMatrix <- sqrt(InverseSigmaDMP)%*%SigmaOfMinimunDet%*%sqrt(InverseSigmaDMP)
  
  # Estimacion objetos algoritmo Ebadi
  TraceRhoSquare <- tr(CorMatrix %^% 2) - (NumVariable**2)/MiddlePoint
  TraceRhoCubic  <- tr(CorMatrix %^% 3) - ((3*NumVariable)/MiddlePoint)*tr(CorMatrix %^% 2) + ((2*(NumVariable**3))/(MiddlePoint**2))
  
  # Estimadores Ui y Zi
  ListEstimUi <- c()
  ListEstimZi <- c()
  ListOutliers <- c()
  Zalpha <- qnorm(1-Alpha,0,1)
  ZAlphaMedia <- qnorm(1-Alpha/2,0,1)
  ConstantMDP <- 1 + (2*NumVariable)/(Observations*sqrt(TraceRhoSquare))
  for( i  in 1:dim(DatesNorm)[1])
  {
    DistanceMahalanobisXi <- t((DatesNorm[i,] - MeanDMP))%*%InverseSigmaDMP%*%(DatesNorm[i,] - MeanDMP)
    Ui <- (DistanceMahalanobisXi - NumVariable)/(2*sqrt(TraceRhoSquare))
    ListEstimUi <- c(ListEstimUi,Ui)
    Zi <- Ui-(4*TraceRhoCubic*(Zalpha**2-1))/(3*(2*TraceRhoSquare)^(3/2))
    ListEstimZi <- c(ListEstimZi,Zi)
    OutlierZi <- (DistanceMahalanobisXi - NumVariable)/(2*ConstantMDP*sqrt(TraceRhoSquare))-(4*TraceRhoCubic*(ZAlphaMedia**2-1))/(3*(2*TraceRhoSquare)^(3/2))
    if(OutlierZi > ZAlphaMedia)
    {
      ListOutliers <- c(ListOutliers,i)
    }
  }
  
  if(OutliersFlag)
  {
    for (i in ListOutliers)
    {
      DatesNorm[i,] <- NA
      ListEstimUi[i] <- NA
      ListEstimZi[i] <- NA
    }
    ListEstimUi <- ListEstimUi[!is.na(ListEstimUi)]
    ListEstimZi <- ListEstimZi[!is.na(ListEstimZi)]
    DatesNorm <-  DatesNorm[!is.na(DatesNorm)]
  }
  return(list(DatesNorm,ListEstimUi,ListEstimZi))
}

algorithmPartTwo <- function(DatesNorm, NumVariable,Alpha, OutliersFlag = FALSE)
{
  
  MeanDMP <- mean(DatesNorm) 
  SigmaOfMinimunDet <- cov(DatesNorm)
  SigmaDMP <- matrix(diag(SigmaOfMinimunDet[row(SigmaOfMinimunDet)==col(SigmaOfMinimunDet)]),ncol = NumVariable)
  InverseSigmaDMP <- solve(SigmaDMP)
  SigmaOfMinimunDet <- cov(DatesNorm)
  CorMatrix <- sqrt(InverseSigmaDMP)%*%SigmaOfMinimunDet%*%sqrt(InverseSigmaDMP)
  Observations <- dim(DatesNorm)[1]
  # Estimacion objetos algoritmo Ebadi
  TraceRhoSquare <- tr(CorMatrix %^% 2) - (NumVariable**2)/Observations
  TraceRhoCubic  <- tr(CorMatrix %^% 3) - ((3*NumVariable)/Observations)*tr(CorMatrix %^% 2) + ((2*(NumVariable**3))/(Observations**2))
  
  # Estimadores Ui y Zi
  ListEstimUi <- c()
  ListEstimZi <- c()
  ListOutliers <- c()
  Zalpha <- qnorm(1-Alpha,0,1)
  ZAlphaMedia <- qnorm(1-Alpha/2,0,1)
  ConstantMDP <- 1 + (2*NumVariable)/(Observations*sqrt(TraceRhoSquare))
  for( i  in 1:dim(DatesNorm)[1])
  {
    Denominator <- 1+pnorm(ZAlphaMedia,0,1)*(NumVariable^(-1))*((1-Alpha/2)^(-1))*sqrt(2*TraceRhoSquare)
    DistanceMahalanobisXi <- (t((DatesNorm[i,] - MeanDMP))%*%InverseSigmaDMP%*%(DatesNorm[i,] - MeanDMP))/Denominator
    Ui <- (DistanceMahalanobisXi - NumVariable)/(2*sqrt(TraceRhoSquare))
    ListEstimUi <- c(ListEstimUi,Ui)
    Zi <- Ui-(4*TraceRhoCubic*(Zalpha**2-1))/(3*(2*TraceRhoSquare)^(3/2))
    ListEstimZi <- c(ListEstimZi,Zi)
    OutlierZi <- (DistanceMahalanobisXi - NumVariable)/(2*ConstantMDP*sqrt(TraceRhoSquare))-(4*TraceRhoCubic*(ZAlphaMedia**2-1))/(3*(2*TraceRhoSquare)^(3/2))
  }
  return(list(DatesNorm,ListEstimUi,ListEstimZi))
}

# Variables locales
NumVariable <- 200
Observations <- 10000
MiddlePoint <- floor(Observations/2)+1
Percent <- 1
Alpha <- 0.05
# Matrices de Covarianza
SigmaFirstDistribution <- diag(NumVariable)
SigmaSecondDistribution <- c()
for (rowindex in 1:NumVariable){
  for (columnindex in 1:NumVariable)
  {
    SigmaSecondDistribution <- c(SigmaSecondDistribution,(0.5)^abs(rowindex-columnindex))
  }
}
SigmaSecondDistribution <- matrix(SigmaSecondDistribution, ncol = NumVariable)

# Vector de Medias
Mu <- rep(0,NumVariable)
DatesNorm <- c()
# Generacion de Valores con Valores atipicos
if(Percent>0)
{
  DatesNorm <- mvrnorm(floor(Observations*Percent),Mu,SigmaFirstDistribution)
}

if(Percent < 1)
{
  DatesNorm <- rbind(DatesNorm,mvrnorm(as.integer(ceiling(Observations*(1-Percent))),Mu,SigmaSecondDistribution))
}

# Algoritmhreweight

Firstestimation <- algorithmUiZi(DatesNorm, NumVariable,Observations, Alpha, FALSE)
Secondestimation <- algorithmPartTwo(matrix(Firstestimation[[1]], ncol = NumVariable), NumVariable, Alpha)

#Funcion CDF
ListEstimUi <- Firstestimation[[2]]
ListEstimZi <- Firstestimation[[3]]
NormalStandarSeq <- seq(-4, 4, length.out= Observations)
NormalStandarCDF <- pnorm(NormalStandarSeq, mean = 0, sd = 1)



# create sample data
sample_Data <- ListEstimUi

# calculate CDF
CDF1 <- ecdf(sample_Data )
CDF2 <- ecdf(ListEstimZi)
# draw the cdf plot
plot(CDF1,col = "red", xlim = c(-4,4))
lines(CDF2,col = "blue")
lines(NormalStandarSeq,NormalStandarCDF, col="black")



