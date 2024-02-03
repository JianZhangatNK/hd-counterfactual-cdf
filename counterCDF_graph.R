### Draw graph for empirical result #####
# install packages before runing this code.
library('latex2exp')
library('plotrix')

counterCDF_graph <- function(object){
  ## Initializing...
  object1 = object$Estimate
  ylim.list = object1$evalY
  F1yhat = object1$Y1
  Fhat.dml = object1$DML
  Fhat.psdb = object1$PSDB
  StE.dml = F1yhat - Fhat.dml
  StE.psdb = F1yhat - Fhat.psdb
  ngrid = length(ylim.list)
  
  object2 = object$Bootstrap
  n = object2$n
  min_Y = object2$minY
  max_Y = object2$maxY
  F1yhat.B = object2$Y1
  Fhat.dml.B = object2$DML
  Fhat.psdb.B = object2$PSDB
  StE.dml.B = F1yhat.B - Fhat.dml.B
  StE.psdb.B = F1yhat.B - Fhat.psdb.B
  B = dim(F1yhat.B)[1]
  
  par(mfrow=c(1,1))
  ## F1(y) versus F^C(y)
  plot(ylim.list, F1yhat, type = "l", pch = 2, 
       col = "red", xlab = "y", xlim = c(min(ylim.list),max(ylim.list)), ylim = c(0, 1), main = "Treated v.s. Counterfactual")
  lines(ylim.list, Fhat.dml, type = 'l', pch = 2, col = "blue")
  lines(ylim.list, Fhat.psdb, type = 'l', pch = 2, col = "green")
  
  legend(x = "topleft", 
         legend=c("Treated", "DML", "PSDB"),
         col=c("red", "blue",  "green"), 
         lty=c(1,1), lwd = c(1,1), cex=0.8)
  
  par(mfrow=c(1,2))
  ## DML method
  sd.dml.bt = sqrt((colSums(StE.dml.B^2)-B*colMeans(StE.dml.B)^2)/(B-1))*sqrt(B/n)
  StE.dml.up = StE.dml + 1.96*sd.dml.bt
  StE.dml.lo = StE.dml - 1.96*sd.dml.bt
  plot(ylim.list, StE.dml.up, type = "l", lty = 2, pch = 2, 
       col = "blue", xlab = "y", ylab = "Stru. Effect", 
       ylim = c(min(StE.dml.lo), max(StE.dml.up)), 
       main = "Stru. Effect: DML")
  lines(ylim.list, StE.dml, type = 'l', pch = 2, col = "black")
  lines(ylim.list, StE.dml.lo, type = 'l', lty = 2,  pch = 2, col = "blue")
  abline(h=0, lty = 'dashed', col = "red")
  
  ## PSDB method
  sd.psdb.bt = sqrt((colSums(StE.psdb.B^2)-B*colMeans(StE.psdb.B)^2)/(B-1))*sqrt(B/n)
  StE.psdb.up = StE.psdb + 1.96*sd.psdb.bt
  StE.psdb.lo = StE.psdb - 1.96*sd.psdb.bt
  plot(ylim.list, StE.psdb.up, type = "l", lty = 2, pch = 2, 
       col = "blue", xlab = "y", ylab = "Stru. Effect", 
       ylim = c(min(StE.psdb.lo), max(StE.psdb.up)), 
       main = "Stru. Effect: PSDB")
  lines(ylim.list, StE.psdb, type = 'l', pch = 2, col = "black")
  lines(ylim.list, StE.psdb.lo, type = 'l', lty = 2,  pch = 2, col = "blue")
  abline(h=0, lty = 'dashed', col = "red")
  
}



