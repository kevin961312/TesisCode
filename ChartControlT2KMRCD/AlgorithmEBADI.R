
source(here::here("Downloads/Tesis/TesisCode/ChartControlT2KMRCD/AlgorithmRoMDP.R"))

AlgorithmMDPCFPart1 <- function(DatesNorm, NumVariable, Observations, Alpha = 0.05, outliersData = c(), OutliersFlag = FALSE) {
  MiddlePoint <- floor(Observations / 2) + 1
  AlgorithmObjects <- AlgorithmRoMDP(DatesNorm)
  # calculo de Media y matrices de correlacion y de varianza
  MeanDMP <- unlist(AlgorithmObjects[1])
  SigmaDMP <- matrix(diag(as.numeric(unlist(AlgorithmObjects[2]))), ncol = NumVariable)
  InverseSigmaDMP <- solve(SigmaDMP)
  SigmaOfMinimunDet <- AlgorithmObjects[3][[1]]
  CorMatrix <- sqrt(InverseSigmaDMP) %*% SigmaOfMinimunDet %*% sqrt(InverseSigmaDMP)

  # Estimacion objetos algoritmo Ebadi
  TraceRhoSquare <- tr(CorMatrix %^% 2) - (NumVariable**2) / MiddlePoint
  TraceRhoCubic <- tr(CorMatrix %^% 3) - ((3 * NumVariable) / MiddlePoint) * tr(CorMatrix %^% 2) + ((2 * (NumVariable**3)) / (MiddlePoint**2))

  # Estimadores Ui y Zi
  ListEstimUi <- c()
  ListEstimZi <- c()
  ListOutliers <- c()
  Zalpha <- qnorm(1 - Alpha, 0, 1)
  ConstantMDP <- 1 + (2 * NumVariable) / (Observations * sqrt(TraceRhoSquare))


  if (OutliersFlag) {
    for (i in 1:dim(outliersData)[1])
    {
      DistanceMahalanobisXi <- t((outliersData[i, ] - MeanDMP)) %*% InverseSigmaDMP %*% (outliersData[i, ] - MeanDMP)
      Ui <- (DistanceMahalanobisXi - NumVariable) / (2 * ConstantMDP * sqrt(TraceRhoSquare))
      ListEstimUi <- c(ListEstimUi, Ui)
      Zi <- Ui - (4 * TraceRhoCubic * (Zalpha**2 - 1)) / (3 * (2 * TraceRhoSquare)^(3 / 2))
      ListEstimZi <- c(ListEstimZi, Zi)
    }
  } else {
    for (i in 1:dim(DatesNorm)[1])
    {
      DistanceMahalanobisXi <- t((DatesNorm[i, ] - MeanDMP)) %*% InverseSigmaDMP %*% (DatesNorm[i, ] - MeanDMP)
      Ui <- (DistanceMahalanobisXi - NumVariable) / (2 * ConstantMDP * sqrt(TraceRhoSquare))
      ListEstimUi <- c(ListEstimUi, Ui)
      Zi <- Ui - (4 * TraceRhoCubic * (Zalpha**2 - 1)) / (3 * (2 * TraceRhoSquare)^(3 / 2))
      ListEstimZi <- c(ListEstimZi, Zi)
    }
  }
  return(list("Ui" = ListEstimUi, "Zi" = ListEstimZi))
}
