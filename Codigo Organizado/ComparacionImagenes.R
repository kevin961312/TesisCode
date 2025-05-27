load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data100x200x1000x0.05x10.RData",D1x2x05x1<- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data100x200x1000x0.05x20.RData",D1x2x05x2 <- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data150x200x1000x0.05x10.RData",D15x2x05x1<- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data150x200x1000x0.05x20.RData",D15x2x05x2 <- new.env())



library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/Comparacion200x1000x0.05.svg"
svglite(urlplot, width = 8, height = 8)
par(mar = c(10,5, 5, 5))
plot(D1x2x05x1$MatrixDeltaMRCD[,1],D1x2x05x1$MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(D1x2x05x2$MatrixDeltaMRCD[,1],D1x2x05x2$MatrixDeltaMRCD[,2], col='green',pch = 17,type='b',lty = 4)
lines(D15x2x05x1$MatrixDeltaMRCD[,1],D15x2x05x1$MatrixDeltaMRCD[,2], col='blue',pch = 18,type='b',lty = 5)
lines(D15x2x05x2$MatrixDeltaMRCD[,1],D15x2x05x2$MatrixDeltaMRCD[,2], col='orange',pch = 15,type='b',lty = 6)
dev.off()



load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data100x250x1000x0.05x10.RData",D1x25x05x1<- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data100x250x1000x0.05x20.RData",D1x25x05x2 <- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data150x250x1000x0.05x10.RData",D15x25x05x1<- new.env())
load("~/Downloads/Tesis/TesisCode/Codigo Organizado/DataTG/Data150x250x1000x0.05x20.RData",D15x25x05x2 <- new.env())

library(svglite)
urlplot ="/Users/kevin.pineda/Desktop/Imagenes TG 1/Comparacion250x1000x0.05.svg"
svglite(urlplot, width = 8, height = 8)
par(mar = c(10,5, 5, 5))
plot(D1x25x05x1$MatrixDeltaMRCD[,1],D1x25x05x1$MatrixDeltaMRCD[,2], ylim = c(0,1), col='red', pch = 16,type='b', xlab ="Delta Values", ylab="Signal Probability",lty = 1)
lines(D1x25x05x2$MatrixDeltaMRCD[,1],D1x25x05x2$MatrixDeltaMRCD[,2], col='green',pch = 17,type='b',lty = 4)
lines(D15x25x05x1$MatrixDeltaMRCD[,1],D15x25x05x1$MatrixDeltaMRCD[,2], col='blue',pch = 18,type='b',lty = 5)
lines(D15x25x05x2$MatrixDeltaMRCD[,1],D15x25x05x2$MatrixDeltaMRCD[,2], col='orange',pch = 15,type='b',lty = 6)
dev.off()

