## install packages "MASS" and "glmnet" before runing this code.
library('MASS')
library('glmnet')

###### Main Functions #######
counterCDF <- function(y, X, D, y.grids, B=500, f){
  m = length(y.grids)
  n = length(y)
  X0 = X[as.logical(1-D),]
  X1 = X[D,]
  y0 = y[as.logical(1-D)]
  y1 = y[D]
  b.xi = matrix(rexp(B*n, rate = 1), B, n) #bootstrap
  datafolds.list = Rand.splitting(n, f)
  
  F0yhat = rep(NA, m)
  F1yhat = rep(NA, m)
  Fhat.dml = rep(NA, m)
  Fhat.psdb = rep(NA, m)
  F1yhat.B = matrix(NA, B, m)
  F0yhat.B = matrix(NA, B, m)
  Fhat.dml.B = matrix(NA, B, m)
  Fhat.psdb.B = matrix(NA, B, m)
  for (i in 1:m) {
    cat('eval Y =', y.grids[i], ', # =', i, '\n')
    F0yhat[i] = mean((y0<=y.grids[i]))
    F0yhat.B[,i] = (b.xi[,as.logical(1-D)]%*%(y0<=y.grids[i]))/sum(D)
    F1yhat[i] = mean((y1<=y.grids[i]))
    F1yhat.B[,i] = (b.xi[,D]%*%(y1<=y.grids[i]))/sum(D)
    
    alpha.xz = rep(NA, n)
    beta.xz = rep(NA, n)
    Fyxz = rep(NA, n)
    Myxz = rep(NA, n)
    PSxz = rep(NA, n)
    for(k in 1:f){
      cat('  fold =', k, '\n')
      foldind = datafolds.list[[k]] #datafolds[k,]
      D.nfold = D[-foldind]
      X.nfold = X[-foldind,]
      y.nfold = y[-foldind]
      
      ### Estimating alpha(x)
      ## Estimating propensity score
      d.cv.out = cv.glmnet(X.nfold, D.nfold, type.measure = 'mse',
                           family = 'binomial')
      xzcoeff = coefficients(d.cv.out, s = d.cv.out$lambda.min)
      nz_index = xzcoeff@i[-1] + 1
      e.p.ind = cbind(0,X[foldind,nz_index-1]) %*%
        c(0,xzcoeff[nz_index]) + xzcoeff[1]
      e.propensity = exp(e.p.ind)/(1+exp(e.p.ind)) 
      PSxz[foldind] = e.propensity
      
      ## alpha(x)
      alpha.xz[foldind] = e.propensity/(1-e.propensity) * 
        (1-mean(D.nfold))/(mean(D.nfold))
      
      ### Estimating conditional CDF when D=0
      
      ## propensity score
      Myxz[foldind] = (1-D[foldind])*e.propensity*(y[foldind]<=y.grids[i])/
        (1-e.propensity)
      
      ## lasso
      X0.nfold = X.nfold[as.logical(1-D.nfold),]
      X1.nfold = X.nfold[D.nfold,]
      y0.nfold = y.nfold[as.logical(1-D.nfold)]
      yvar0.nfold = as.numeric(y0.nfold<=y.grids[i])
      
      cv.out = cv.glmnet(X0.nfold, yvar0.nfold, type.measure = 'mse',
                         family = "binomial")
      y0coeff = coefficients(cv.out, s = cv.out$lambda.min)
      nz_index0 = y0coeff@i[-1] + 1
      Fyxz.index= cbind(0, X[foldind,nz_index0-1]) %*%
        c(0, y0coeff[nz_index0]) + y0coeff[1]
      Fyxz[foldind] = exp(Fyxz.index)/(1+exp(Fyxz.index))
      
      ## beta(x)
      beta.xz[foldind] = Fyxz[foldind]/
        ((1-e.propensity)*mean(D.nfold))
    }
    
    ### Estimate Counterfactual CDF (unconditional)
    Fyxz.c = Fyxz[D]
    Fyxz0 = Fyxz[as.logical(1-D)]
    Fhat.dbias = alpha.xz[as.logical(1-D)] * ((y0<=y.grids[i])-Fyxz0)
    Fhat.dml[i] = mean(Fyxz.c) + mean(Fhat.dbias)
    
    Fhat.dml.B[,i] = (b.xi[,D]%*%Fyxz.c)/sum(D) + 
      (b.xi[,as.logical(1-D)]%*%Fhat.dbias)/(n-sum(D)) ##bootstrap
    
    Fhat.dbias2 = beta.xz * (D-PSxz)
    Fhat.psdb[i] = sum(Myxz)/sum(D) + mean(Fhat.dbias2)
    
    Fhat.psdb.B[,i] = (b.xi%*%Myxz)/sum(D) + (b.xi%*%Fhat.dbias2)/n
  }
  counterCDF = list(
    evalY = y.grids,
    Y0 = F0yhat,
    Y1 = F1yhat,
    DML = Fhat.dml,
    PSDB = Fhat.psdb
  )
  bootstrap.CDF = list(
    n = n,
    minY = min(y),
    maxY = max(y),
    evalY = y.grids,
    Y0 = F0yhat.B,
    Y1 = F1yhat.B,
    DML = Fhat.dml.B,
    PSDB = Fhat.psdb.B
  )
  return(list(Estimate = counterCDF, Bootstrap = bootstrap.CDF))
}

###### Functions #######
mat.vec <- function(p, q){
  nq = dim(q)
  if (nq == 1){
    return(p*q)
  }else if (nq == 0){
    return(0)
  }else if (nq>1){
    return(p %*% q)
  }
}

## OB Decomposition
OB_decomp = function(object1){
  struc.Effect.dml = object1$Y1 - object1$DML 
  struc.Effect.psdb = object1$Y1 - object1$PSDB
  struc_object = list(
    DML = struc.Effect.dml,
    PSDB = struc.Effect.psdb
  )
  comp.Effect.dml = object1$DML - object1$Y0
  comp.Effect.psdb = object1$PSDB - object1$Y0
  comp_object = list(
    DML = comp.Effect.dml,
    PSDB = comp.Effect.psdb
  )
  return(list(StructureEffect = struc_object, 
              CompositionEffect = comp_object))
}

## basis functions
basisF <- function(x){cbind(x, x^2)}

## Random Splitting
Rand.splitting <- function(n, f){
  ndata.rd = sample(1:n, n, replace = FALSE)
  n_int = floor(n/f)*f
  datafolds = matrix(ndata.rd[1:n_int], floor(n/f), f)
  newcol = ndata.rd[(floor(n/f)*(f-1)+1):n]
  datafolds.list1 = as.list(as.data.frame(datafolds[,1:(f-1)]))
  datafolds.list2 = list(newcol)
  names(datafolds.list2) = paste0("v",f)
  datafolds.list = c(datafolds.list1, datafolds.list2)
  return(datafolds.list)
}

## A Simple Simulation Example
simulate_data = function(n, p){
  x = rnorm(n, 0, 1)
  z.var.mat = rho.mat(p, 0.6)
  z = mvrnorm(n, rep(0,p), z.var.mat) #z = replicate(p, rnorm(n, 0, 1))
  zprime = rnorm(n, 0.5, 1)
  D = (runif(n, 0, 1)>0.5)
  x = x + 0.5*D
  z[,1] = (1-D)*z[,1] + D*zprime
  z[,2] = 0.7*z[,2] + 0.7*z[,1]
  z[,3] = 0.8*z[,3] + 0.6*z[,1]
  eps = rnorm(n, 0, 1)
  y = 0.5*x + as.vector(z[,1]+0.1*z[,4]) + 0.5*D + eps
  X = cbind(x,z)
  return(list(y=y,X=X,D=D))
}

rho.mat <- function(p, rho){
  rho.mat = matrix(nrow = p, ncol = p)
  for (i in 1:p) {
    for(j in 1:p) {
      rho.mat[i,j] = rho^(abs(i-j))
    }
  }
  return(rho.mat)
}


