rm(list=ls())
library(svglite)
urlplot ="/Users/kevin.pineda/Downloads/Tesis/TesisCode/Imagenes TG 2 KMRCD/Gamma(1-5,beta)/ImagenGamma(1-5,beta)100x200x0.2.svg"
load("/Users/kevin.pineda/Downloads/Tesis/TesisCode/KMRCD/DatosKMRCD/Parametro(1-5,beta)/SinalprobabilityGamma(1-5,beta)100x200x0.2.RData")

par(mar = c(10,5, 5, 5))
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)
lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
lines(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], col='purple',pch = 1,type='b',lty = 6)


svglite(urlplot, width = 8, height = 8)
p = plot(MatrixDeltaMRCD[,1],MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(MatrixDeltaT2MOD[,1],MatrixDeltaT2MOD[,2], col='green',pch = 17,type='b',lty = 4)
lines(MatrixDeltaEBADIUI[,1],MatrixDeltaEBADIUI[,2], col='blue',pch = 18,type='b',lty = 5)
lines(MatrixDeltaEBADIZI[,1],MatrixDeltaEBADIZI[,2], col='orange',pch = 15,type='b',lty = 6)
lines(MatrixDeltaKMRCD[,1],MatrixDeltaKMRCD[,2], col='purple',pch = 1,type='b',lty = 6)
dev.off()