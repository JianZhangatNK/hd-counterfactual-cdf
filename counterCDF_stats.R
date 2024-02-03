## Calculate ASE, ATT, QTT, KS 
## Install packages "DescTools" and "dineq" before runing this code.
library('DescTools')
library('dineq')

##### Main Functions ######
counterCDF_stats <- function(object){
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
  min_Y = object2$minY
  max_Y = object2$maxY
  F1yhat.B = object2$Y1
  Fhat.dml.B = object2$DML
  Fhat.psdb.B = object2$PSDB
  StE.dml.B = F1yhat.B - Fhat.dml.B
  StE.psdb.B = F1yhat.B - Fhat.psdb.B
  B = dim(F1yhat.B)[1]
  
  ## ASE: f=approximately uniformly distributed over "ylim.list"
  ase_dml = mean(StE.dml)
  ase_psdb = mean(StE.psdb)
  
  ####
  ylim.list.expand = c(min_Y, ylim.list, max_Y)
  F1yhat.expand = c(0, F1yhat, 1)
  Fhat.dml.expand = c(0, Fhat.dml, 1)
  Fhat.psdb.expand = c(0, Fhat.psdb, 1)
  
  ## ATT
  y1.mean = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                  (F1yhat.expand[-1] - F1yhat.expand[-(ngrid+2)])/2)
  yc.dml.mean = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                      (Fhat.dml.expand[-1] - Fhat.dml.expand[-(ngrid+2)])/2)
  yc.psdb.mean = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                       (Fhat.psdb.expand[-1] - Fhat.psdb.expand[-(ngrid+2)])/2)
  att_dml = y1.mean - yc.dml.mean
  att_psdb = y1.mean - yc.psdb.mean
  
  ## QTT: quantiles=10%,25%,50%,75%,90%
  qtt_dml_10 = findQuan(ylim.list.expand, F1yhat.expand, 0.1) -
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.1)
  qtt_dml_25 = findQuan(ylim.list.expand, F1yhat.expand, 0.25) -
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.25)
  qtt_dml_50 = findQuan(ylim.list.expand, F1yhat.expand, 0.5) -
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.5)
  qtt_dml_75 = findQuan(ylim.list.expand, F1yhat.expand, 0.75) -
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.75)
  qtt_dml_90 = findQuan(ylim.list.expand, F1yhat.expand, 0.9) -
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.9)
  
  qtt_psdb_10 = findQuan(ylim.list.expand, F1yhat.expand, 0.1) -
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.1)
  qtt_psdb_25 = findQuan(ylim.list.expand, F1yhat.expand, 0.25) -
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.25)
  qtt_psdb_50 = findQuan(ylim.list.expand, F1yhat.expand, 0.5) -
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.5)
  qtt_psdb_75 = findQuan(ylim.list.expand, F1yhat.expand, 0.75) -
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.75)
  qtt_psdb_90 = findQuan(ylim.list.expand, F1yhat.expand, 0.9) -
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.9)
  
  ## Counterfactual CDF:
  ## interquantile:90-10, 90-50, 50-10
  intq_y1_90_10 = findQuan(ylim.list.expand, F1yhat.expand, 0.9) - 
    findQuan(ylim.list.expand, F1yhat.expand, 0.1)
  intq_y1_90_50 = findQuan(ylim.list.expand, F1yhat.expand, 0.9) - 
    findQuan(ylim.list.expand, F1yhat.expand, 0.5)
  intq_y1_50_10 = findQuan(ylim.list.expand, F1yhat.expand, 0.5) - 
    findQuan(ylim.list.expand, F1yhat.expand, 0.1)
  
  intq_dml_90_10 = findQuan(ylim.list.expand, Fhat.dml.expand, 0.9) - 
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.1)
  intq_dml_90_50 = findQuan(ylim.list.expand, Fhat.dml.expand, 0.9) - 
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.5)
  intq_dml_50_10 = findQuan(ylim.list.expand, Fhat.dml.expand, 0.5) - 
    findQuan(ylim.list.expand, Fhat.dml.expand, 0.1)
  
  intq_psdb_90_10 = findQuan(ylim.list.expand, Fhat.psdb.expand, 0.9) - 
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.1)
  intq_psdb_90_50 = findQuan(ylim.list.expand, Fhat.psdb.expand, 0.9) - 
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.5)
  intq_psdb_50_10 = findQuan(ylim.list.expand, Fhat.psdb.expand, 0.5) - 
    findQuan(ylim.list.expand, Fhat.psdb.expand, 0.1)
  
  ## Standard deviation
  #ylim.list.expand1 = c(min(y[D]), ylim.list, max(y[D]))
  std_y1 = findSTD2(ylim.list.expand, F1yhat.expand)
  std_dml = findSTD2(ylim.list.expand, Fhat.dml.expand)
  std_psdb = findSTD2(ylim.list.expand, Fhat.psdb.expand)
  
  ## Gini index
  Gini_y1 = findGINI(ylim.list.expand, F1yhat.expand)
  Gini_dml = findGINI(ylim.list.expand, Fhat.dml.expand)
  Gini_psdb = findGINI(ylim.list.expand, Fhat.psdb.expand)
  
  ## Theil index
  Theil_y1 = findTheil(ylim.list.expand, F1yhat.expand)
  Theil_dml = findTheil(ylim.list.expand, Fhat.dml.expand)
  Theil_psdb = findTheil(ylim.list.expand, Fhat.psdb.expand)
  
  ### Kolmogorov-Smirnov Test
  KShat.dml = max(abs(StE.dml))
  KShat.psdb = max(abs(StE.psdb))
  
  
  ## Standard deviation by bootstrapping:
  # structure effects
  ase_dml.B = rep(NA, B)
  ase_psdb.B = rep(NA, B)
  att_dml.B = rep(NA, B)
  att_psdb.B = rep(NA, B)
  qtt_dml_10.B = rep(NA, B)
  qtt_dml_25.B = rep(NA, B)
  qtt_dml_50.B = rep(NA, B)
  qtt_dml_75.B = rep(NA, B)
  qtt_dml_90.B = rep(NA, B)
  qtt_psdb_10.B = rep(NA, B)
  qtt_psdb_25.B = rep(NA, B)
  qtt_psdb_50.B = rep(NA, B)
  qtt_psdb_75.B = rep(NA, B)
  qtt_psdb_90.B = rep(NA, B)
  # Counterfactuals
  intq_y1_50_10.B = rep(NA, B)
  intq_y1_90_10.B = rep(NA, B)
  intq_y1_90_50.B = rep(NA, B)
  intq_dml_50_10.B = rep(NA, B)
  intq_dml_90_10.B = rep(NA, B)
  intq_dml_90_50.B = rep(NA, B)
  intq_psdb_50_10.B = rep(NA, B)
  intq_psdb_90_10.B = rep(NA, B)
  intq_psdb_90_50.B = rep(NA, B)
  std_y1.B = rep(NA, B)
  std_dml.B = rep(NA, B)
  std_psdb.B = rep(NA, B)
  Gini_y1.B = rep(NA, B)
  Gini_dml.B = rep(NA, B)
  Gini_psdb.B = rep(NA, B)
  Theil_y1.B = rep(NA, B)
  Theil_dml.B = rep(NA, B)
  Theil_psdb.B = rep(NA, B)
  # K-S test
  KShat.dml.B = rep(NA, B)
  KShat.psdb.B = rep(NA, B)
  for (b in 1:B) {
    ## ASE: Bootstrap version
    ase_dml.B[b] = mean(as.vector(StE.dml.B[b,]))
    ase_psdb.B[b] = mean(as.vector(StE.psdb.B[b,]))
    
    ####
    F1yhat.B.expand = c(0, as.vector(F1yhat.B[b,]), 1)
    Fhat.dml.B.expand = c(0, as.vector(Fhat.dml.B[b,]), 1)
    Fhat.psdb.B.expand = c(0, as.vector(Fhat.psdb.B[b,]), 1)
    
    ## ATT: Bootstrap version
    y1.mean.B = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                      (F1yhat.B.expand[-1] - F1yhat.B.expand[-(ngrid+2)])/2)
    yc.dml.mean.B = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                          (Fhat.dml.B.expand[-1] - 
                             Fhat.dml.B.expand[-(ngrid+2)])/2)
    yc.psdb.mean.B = sum((ylim.list.expand[-1]+ylim.list.expand[-(ngrid+2)])*
                           (Fhat.psdb.B.expand[-1] - 
                              Fhat.psdb.B.expand[-(ngrid+2)])/2)
    att_dml.B[b] = y1.mean.B - yc.dml.mean.B
    att_psdb.B[b] = y1.mean.B - yc.psdb.mean.B
    
    ## QTT: quantiles=10%,25%,50%,75%,90%
    qtt_dml_10.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.1) -
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.1)
    qtt_dml_25.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.25) -
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.25)
    qtt_dml_50.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.5) -
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.5)
    qtt_dml_75.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.75) -
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.75)
    qtt_dml_90.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.9) -
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.9)
    
    qtt_psdb_10.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.1) -
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.1)
    qtt_psdb_25.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.25) -
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.25)
    qtt_psdb_50.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.5) -
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.5)
    qtt_psdb_75.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.75) -
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.75)
    qtt_psdb_90.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.9) -
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.9)
    
    ## Counterfactual CDF:
    ## interquantile:90-10, 90-50, 50-10
    intq_y1_90_10.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.9) - 
      findQuan(ylim.list.expand, F1yhat.B.expand, 0.1)
    intq_y1_90_50.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.9) - 
      findQuan(ylim.list.expand, F1yhat.B.expand, 0.5)
    intq_y1_50_10.B[b] = findQuan(ylim.list.expand, F1yhat.B.expand, 0.5) - 
      findQuan(ylim.list.expand, F1yhat.B.expand, 0.1)
    
    intq_dml_90_10.B[b] = findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.9) - 
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.1)
    intq_dml_90_50.B[b] = findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.9) - 
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.5)
    intq_dml_50_10.B[b] = findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.5) - 
      findQuan(ylim.list.expand, Fhat.dml.B.expand, 0.1)
    
    intq_psdb_90_10.B[b] = findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.9) - 
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.1)
    intq_psdb_90_50.B[b] = findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.9) - 
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.5)
    intq_psdb_50_10.B[b] = findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.5) - 
      findQuan(ylim.list.expand, Fhat.psdb.B.expand, 0.1)
    
    ## Standard deviation
    std_y1.B[b] = findSTD2(ylim.list.expand, F1yhat.B.expand)
    std_dml.B[b] = findSTD2(ylim.list.expand, Fhat.dml.B.expand)
    std_psdb.B[b] = findSTD2(ylim.list.expand, Fhat.psdb.B.expand)
    
    ## Gini index
    Gini_y1.B[b] = findGINI(ylim.list.expand, F1yhat.B.expand)
    Gini_dml.B[b] = findGINI(ylim.list.expand, Fhat.dml.B.expand)
    Gini_psdb.B[b] = findGINI(ylim.list.expand, Fhat.psdb.B.expand)
    
    ## Theil index
    Theil_y1.B[b] = findTheil(ylim.list.expand, F1yhat.B.expand)
    Theil_dml.B[b] = findTheil(ylim.list.expand, Fhat.dml.B.expand)
    Theil_psdb.B[b] = findTheil(ylim.list.expand, Fhat.psdb.B.expand)
    
    ### K-S TEST
    KShat.dml.B[b] = max(abs(StE.dml.B[b,] - StE.dml))
    KShat.psdb.B[b] = max(abs(StE.psdb.B[b,] - StE.psdb))
  }
  ### K-S test
  pvalue.KS.dml = mean((KShat.dml<=KShat.dml.B))
  pvalue.KS.psdb = mean((KShat.psdb<=KShat.psdb.B))
  
  ### OUTPUT
  ASE = list(
    DML = ase_dml,
    PSDB = ase_psdb,
    sd_DML = sd(ase_dml.B)*sqrt(B/n),
    sd_PSDB = sd(ase_psdb.B)*sqrt(B/n)
  )
  ATT = list(
    DML = att_dml,
    PSDB = att_psdb,
    sd_DML = sd(att_dml.B)*sqrt(B/n),
    sd_PSDB = sd(att_psdb.B)*sqrt(B/n)
  )
  QTT = list(
    DML10 = qtt_dml_10,
    DML25 = qtt_dml_25,
    DML50 = qtt_dml_50,
    DML75 = qtt_dml_75,
    DML90 = qtt_dml_90,
    sd_DML10 = sd(qtt_dml_10.B)*sqrt(B/n),
    sd_DML25 = sd(qtt_dml_25.B)*sqrt(B/n),
    sd_DML50 = sd(qtt_dml_50.B)*sqrt(B/n),
    sd_DML75 = sd(qtt_dml_75.B)*sqrt(B/n),
    sd_DML90 = sd(qtt_dml_90.B)*sqrt(B/n),
    PSDB10 = qtt_psdb_10,
    PSDB25 = qtt_psdb_25,
    PSDB50 = qtt_psdb_50,
    PSDB75 = qtt_psdb_75,
    PSDB90 = qtt_psdb_90,
    sd_PSDB10 = sd(qtt_psdb_10.B)*sqrt(B/n),
    sd_PSDB25 = sd(qtt_psdb_25.B)*sqrt(B/n),
    sd_PSDB50 = sd(qtt_psdb_50.B)*sqrt(B/n),
    sd_PSDB75 = sd(qtt_psdb_75.B)*sqrt(B/n),
    sd_PSDB90 = sd(qtt_psdb_90.B)*sqrt(B/n)
  )
  KS_pvalue = list(
    DML = pvalue.KS.dml,
    PSDB = pvalue.KS.psdb
  )
  Interquantile = list(
    DML90_10 = intq_dml_90_10,
    DML90_50 = intq_dml_90_50,
    DML50_10 = intq_dml_50_10,
    sd_DML90_10 = sd(intq_dml_90_10.B)*sqrt(B/n),
    sd_DML90_50 = sd(intq_dml_90_50.B)*sqrt(B/n),
    sd_DML50_10 = sd(intq_dml_50_10.B)*sqrt(B/n),
    PSDB90_10 = intq_psdb_90_10,
    PSDB90_50 = intq_psdb_90_50,
    PSDB50_10 = intq_psdb_50_10,
    sd_PSDB90_10 = sd(intq_psdb_90_10.B)*sqrt(B/n),
    sd_PSDB90_50 = sd(intq_psdb_90_50.B)*sqrt(B/n),
    sd_PSDB50_10 = sd(intq_psdb_50_10.B)*sqrt(B/n)
  )
  Std.err = list(
    Y1 = std_y1,
    DML = std_dml,
    PSDB = std_psdb,
    sd_Y1 = sd(std_y1.B)*sqrt(B/n),
    sd_DML = sd(std_dml.B)*sqrt(B/n),
    sd_PSDB = sd(std_psdb.B)*sqrt(B/n)
  )
  Gini = list(
    Y1 = Gini_y1,
    DML = Gini_dml,
    PSDB = Gini_psdb,
    sd_Y1 = sd(Gini_y1.B)*sqrt(B/n),
    sd_DML = sd(Gini_dml.B)*sqrt(B/n),
    sd_PSDB = sd(Gini_psdb.B)*sqrt(B/n)
  )
  Theil = list(
    Y1 = Theil_y1,
    DML = Theil_dml,
    PSDB = Theil_psdb,
    sd_Y1 = sd(Theil_y1.B)*sqrt(B/n),
    sd_DML = sd(Theil_dml.B)*sqrt(B/n),
    sd_PSDB = sd(Theil_psdb.B)*sqrt(B/n)
  )
  
  ## OUTPUT TO TXT FILE
  outfile <- "CZZ_OUTPUT.txt"
  cat("************* Jun Cai, Jian Zhang, Yahong Zhou(2024) OUTPUT RESULTS ***********", file=outfile, append=TRUE, "\n")
  cat("**********ASE:*********",file=outfile, append=TRUE, "\n")
  cat("ASE_DML = ", sprintf("%0.3f", ase_dml), "(", sprintf("%0.3f", sd(ase_dml.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("ASE_PSDB = ", sprintf("%0.3f", ase_psdb), "(", sprintf("%0.3f", sd(ase_psdb.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("******ATT: ************",file=outfile, append=TRUE, "\n")
  cat("ATT_DML = ", sprintf("%0.3f", att_dml), "(", sprintf("%0.3f", sd(att_dml.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("ATT_PSDB = ", sprintf("%0.3f", att_psdb), "(", sprintf("%0.3f", sd(att_psdb.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("*****QTT: quantiles=10%,25%,50%,75%,90%*****", file=outfile, append=TRUE, "\n")
  cat("DML 10%:", sprintf("%0.3f", qtt_dml_10), "(", sprintf("%0.3f", sd(qtt_dml_10.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    25%:", sprintf("%0.3f", qtt_dml_25), "(", sprintf("%0.3f", sd(qtt_dml_25.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    50%:", sprintf("%0.3f", qtt_dml_50), "(", sprintf("%0.3f", sd(qtt_dml_50.B)*sqrt(B/n)), ")", file=outfile, append=TRUE,"\n")
  cat("    75%:", sprintf("%0.3f", qtt_dml_75), "(", sprintf("%0.3f", sd(qtt_dml_75.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    90%:", sprintf("%0.3f", qtt_dml_90), "(", sprintf("%0.3f", sd(qtt_dml_90.B)*sqrt(B/n)), ")", file=outfile, append=TRUE,"\n")
  cat("PSDB10%:", sprintf("%0.3f", qtt_psdb_10), "(", sprintf("%0.3f", sd(qtt_psdb_10.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    25%:", sprintf("%0.3f", qtt_psdb_25), "(", sprintf("%0.3f", sd(qtt_psdb_25.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    50%:", sprintf("%0.3f", qtt_psdb_50), "(", sprintf("%0.3f", sd(qtt_psdb_50.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    75%:", sprintf("%0.3f", qtt_psdb_75), "(", sprintf("%0.3f", sd(qtt_psdb_75.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("    90%:", sprintf("%0.3f", qtt_psdb_90), "(", sprintf("%0.3f", sd(qtt_psdb_90.B)*sqrt(B/n)), ")", file=outfile, append=TRUE, "\n")
  cat("***********************Counterfactuals**********************", file=outfile, append=TRUE, "\n")
  cat("*********Treated: Y|D=1*********Counterfactual(DML)***********Counterfactual(PSDB)*****************", file=outfile, append=TRUE, "\n")
  cat("        90-10:", sprintf("%0.3f", intq_y1_90_10), "(", sprintf("%0.3f", sd(intq_y1_90_10.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_dml_90_10), "(", sprintf("%0.3f", sd(intq_dml_90_10.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_psdb_90_10), "(", sprintf("%0.3f", sd(intq_psdb_90_10.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("        90-50:", sprintf("%0.3f", intq_y1_90_50), "(", sprintf("%0.3f", sd(intq_y1_90_50.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_dml_90_50), "(", sprintf("%0.3f", sd(intq_dml_90_50.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_psdb_90_50), "(", sprintf("%0.3f", sd(intq_psdb_90_50.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("        50-10:", sprintf("%0.3f", intq_y1_50_10), "(", sprintf("%0.3f", sd(intq_y1_50_10.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_dml_50_10), "(", sprintf("%0.3f", sd(intq_dml_50_10.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", intq_psdb_50_10), "(", sprintf("%0.3f", sd(intq_psdb_50_10.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("Std(log wage):", sprintf("%0.3f", std_y1), "(", sprintf("%0.3f", sd(std_y1.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", std_dml), "(", sprintf("%0.3f", sd(std_dml.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", std_psdb), "(", sprintf("%0.3f", sd(std_psdb.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("         Gini:", sprintf("%0.3f", Gini_y1), "(", sprintf("%0.3f", sd(Gini_y1.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", Gini_dml), "(", sprintf("%0.3f", sd(Gini_dml.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", Gini_psdb), "(", sprintf("%0.3f", sd(Gini_psdb.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("        Theil:", sprintf("%0.3f", Theil_y1), "(", sprintf("%0.3f", sd(Theil_y1.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", Theil_dml), "(", sprintf("%0.3f", sd(Theil_dml.B)*sqrt(B/n)), ")", " && ",
      sprintf("%0.3f", Theil_psdb), "(", sprintf("%0.3f", sd(Theil_psdb.B)*sqrt(B/n)), ")",
      file=outfile, append=TRUE, "\n")
  cat("**********Kolmogorov-Smirnov Test*********",file=outfile, append=TRUE, "\n")
  cat("  P-value(DML): ", sprintf("%0.3f", pvalue.KS.dml), file=outfile, append=TRUE, "\n")
  cat("P-value(PSDB): ", sprintf("%0.3f", pvalue.KS.psdb), file=outfile, append=TRUE, "\n")
  cat(" ", file=outfile, append=TRUE, "\n")
  cat(" ", file=outfile, append=TRUE, "\n")
  
  return(list(
    ASE = ASE,
    ATT = ATT,
    QTT = QTT,
    KS_pvalue = KS_pvalue,
    Interquantile = Interquantile,
    Std_err = Std.err,
    Gini = Gini,
    Theil = Theil
  ))
}

##### Functions ######
# 1. Find quantile of a discrete distribution
findQuan = function(ylist, dist, tau){
  l_y = length(ylist) # dim(ylist) = dim(dist)
  dist_min = dist[1]
  dist_max = dist[l_y]
  if(tau>=dist_max){
    y.quantile = ylist[l_y]
  }
  if(tau<dist_min){
    y.quantile = ylist[1]
  }
  if((tau>=dist_min)*(tau<dist_max)){
    dist0 = max(dist[dist<=tau])
    dist1 = min(dist[dist>tau])
    y0 = max(ylist[dist<=tau])
    y1 = min(ylist[dist>tau])
    y.quantile = y0 + ((tau - dist0)/(dist1 - dist0))*(y1-y0)
  }
  return(y.quantile)
}

# 2. Calculate Standard deviation from distribution function:
findSTD = function(ylist, ydist){
  ylen = length(ylist)
  mean.y = sum((ylist[-1] + ylist[-ylen])*(ydist[-1] - ydist[-ylen])/2)
  mean.y2 = sum((ylist[-1]^2 + ylist[-ylen]^2 + ylist[-1]*ylist[-ylen])*
                  (ydist[-1] - ydist[-ylen])/3)
  std.y = sqrt(mean.y2  - (mean.y)^2)
  return(std.y)
}

findSTD2 = function(ylist, ydist){
  y.intq = findQuan(ylist, ydist, 0.75) - findQuan(ylist, ydist, 0.25)
  std.y = y.intq/1.349
  return(std.y)
}

# 3. Calculate Gini coefficients
findGINI = function(ylist, ydist){
  # Lorentz curve: L(p)
  # pi:
  pi = ydist
  
  # Li
  ylen = length(ylist)
  mean.y = sum((ylist[-1] + ylist[-ylen])*(ydist[-1] - ydist[-ylen])/2)
  Li = rep(NA, ylen)
  Li[1] = 0
  Li[ylen] = 1
  for (jj in 2:(ylen-1)) {
    Li[jj] = sum((ylist[2:jj] + ylist[1:(jj-1)])*
                   (ydist[2:jj] - ydist[1:(jj-1)])/2)/mean.y
  }
  
  # Gini
  xi = abs((Li[-1] - Li[-ylen])/(pi[-1]-pi[-ylen]))
  wi = abs(pi[-1] - pi[-ylen])
  y.Gini = Gini(xi, weights = wi)
  return(y.Gini)
}

# 4. Calculate Gini coefficients
findTheil = function(ylist, ydist){
  # Lorentz curve: L(p)
  # pi
  pi = ydist
  
  # Li
  ylen = length(ylist)
  mean.y = sum((ylist[-1] + ylist[-ylen])*(ydist[-1] - ydist[-ylen])/2)
  Li = rep(NA, ylen)
  Li[1] = 0
  Li[ylen] = 1
  for (jj in 2:(ylen-1)) {
    Li[jj] = sum((ylist[2:jj] + ylist[1:(jj-1)])*
                   (ydist[2:jj] - ydist[1:(jj-1)])/2)/mean.y
  }
  
  # Theil
  xi = abs((Li[-1] - Li[-ylen])/(pi[-1]-pi[-ylen]))
  wi = abs(pi[-1] - pi[-ylen])
  y.theil = theil.wtd(xi, weights = wi)
  return(y.theil)
}




