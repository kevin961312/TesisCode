
SignalProbability <- function (matrixT2, ucl, uclMax, uclKernel )
{
  
  SignalCountUCL = 0
  SignalCountMax = 0
  SignalCountKernel = 0
  
  for (i in 1: dim(matrixT2)[1])
  {
    RowMatrix = matrixT2[i,]
    
    ListFilterUCL = Filter(function(x) x > ucl,RowMatrix)
    SignalCountUCL = SignalCountUCL + length(ListFilterUCL)
    
    
    ListFilterMax = Filter(function(x) x > uclMax,RowMatrix)
    SignalCountMax = SignalCountMax + length(ListFilterMax)
    
    ListFilterKernel = Filter(function(x) x > uclKernel,RowMatrix)
    SignalCountKernel = SignalCountKernel + length(ListFilterKernel)
  }
  
  TotalValue = dim(matrixT2)[1]*dim(matrixT2)[2]
  
  SignalPro = SignalCountUCL/TotalValue
  SignalProMax = SignalCountMax/TotalValue
  SignalProKernel = SignalCountKernel/TotalValue
  
  return(list("SignalPro" = SignalPro,
              "SignalProMax" = SignalProMax,
              "SignalProKernel" = SignalProKernel))
}

