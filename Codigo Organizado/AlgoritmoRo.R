library(Rfast)

AlgorithmRoMDP <- function(DataSet, itertime = 100)
{
  # Dimensión del conjuto de datos
  SampleNumber <- dim(DataSet)[1]
  VariableNumber <- dim(DataSet)[2]
  MiddlePoint <- round(SampleNumber/2)+1
  TransposedDataSet <-t(DataSet)
  
  ValueInitialSubsets <- 2
  BestDet <- 0
  VecZero <- numeric(SampleNumber)
  
  #Ciclo de iteracciones por default 100
  for ( iter in 1:itertime)
  {
    Id <- sample(SampleNumber, ValueInitialSubsets, replace = FALSE)
    SubsetById <- DataSet[Id,]
    MuSubsetById <- Rfast::colmeans(SubsetById)
    VarSubsetById <- Rfast::colVars(SubsetById)
    Sama <- (TransposedDataSet-MuSubsetById)/VarSubsetById
    Distance <- Rfast::colsums(Sama)
    Criterio <- 10
    Count <- 0
    while (Criterio !=0 & Count <= 15)
    {
      Count <- Count + 1
      VectorZeroi <- numeric(SampleNumber)
      DistancePerm <- order(Distance)
      VectorZeroi[ DistancePerm[1:MiddlePoint] ] <- 1
      Criterio <- sum( abs(VectorZeroi - VecZero) )
      VecZero <- VectorZeroi
      NewDataSet <- DataSet[DistancePerm[1:MiddlePoint], ]
      MuSubsetById <- Rfast::colmeans(NewDataSet)
      VarSubsetById <- Rfast::colVars(NewDataSet)
      Sama <- (TransposedDataSet-MuSubsetById)/VarSubsetById
      Distance <- Rfast::colsums(Sama)
    }
    TempDet <- prod(VarSubsetById)
    if(BestDet == 0 | TempDet < BestDet) 
    {
      BestDet <- TempDet
      FinalVec <- VecZero
    }
  }
  SubMCD <- (1:SampleNumber)[FinalVec != 0]
  
  MuSubsetById <- Rfast::colmeans( DataSet[SubMCD, ] )
  VarSubsetById <- Rfast::colVars( DataSet[SubMCD, ] )
  Sigma <- cov(DataSet[SubMCD, ])
  return(list(MuSubsetById, VarSubsetById, Sigma, SubMCD))
}