MethodUCLKernel <- function(T2, alpha) 
{
  EstKernelSmooth = density(T2,kernel="gaussian",bw="nrd")
  valuePairs = cbind(EstKernelSmooth$x,EstKernelSmooth$y)
  valuePairs = valuePairs[valuePairs[,1]>=0,]
  Estand = valuePairs[,2]/sum(valuePairs[,2])
  NumberRow = nrow(valuePairs)
  countValues = 0
  
  for(i in 1:NumberRow)
  {
    if(countValues <  (1-alpha))
    {
      j = i 
      countValues = sum(Estand[1:i])
    }
  }
  
  UCL = valuePairs[j,1]
  return(list(UCL=UCL,alpha1=(1-countValues)))
}