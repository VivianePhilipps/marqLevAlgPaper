### code of section 6 : comparison with other optimization algorithms ###


## prothro data from JM package
library("JM")

prothro$t0 <- as.numeric(prothro$time == 0)


## initial models
lmeFit <- lme(pro ~ treat * (time + t0), random = ~ time | id, data = prothro)
survFit <- coxph(Surv(Time, death) ~ treat, data = prothros, x = TRUE)

##derivForm
dFpro <- list(fixed = ~ 1 + treat, indFixed = c(3, 5), random = ~ 1, indRandom = 2)


## 0.1 scaling
prothro$pro01 <- prothro$pro/10
lmeFit01 <- lme(pro01 ~ treat * (time + t0), random = ~ time | id, data = prothro)

## 10 scaling
prothro$pro10 <- prothro$pro*10
lmeFit10 <- lme(pro10 ~ treat * (time + t0), random = ~ time | id, data = prothro)



## estimation by BFGS alogorithm from optim function, scaling 1 :

jm_pro_val_GH15_BFGS <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_slo_GH15_BFGS <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_both_GH15_BFGS <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by BFGS alogorithm from optim function, scaling 0.1 :

jm_pro01_val_GH15_BFGS <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_slo_GH15_BFGS <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_both_GH15_BFGS <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by BFGS alogorithm from optim function, scaling 10 :

jm_pro10_val_GH15_BFGS <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro10_slo_GH15_BFGS <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro10_both_GH15_BFGS <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",method="BFGS",GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by optimParallel function, scaling 1 :

jm_pro_val_GH15_LBFGSB <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_slo_GH15_LBFGSB <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_both_GH15_LBFGSB <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by optimParallel function, scaling 0.1 :

jm_pro01_val_GH15_LBFGSB <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_slo_GH15_LBFGSB <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_both_GH15_LBFGSB <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by optimParallel function, scaling 10 :

jm_pro10_val_GH15_LBFGSB <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)
if(!exists("jm_pro10_val_GH15_LBFGSB"))
{
    jm_pro10_val_GH15_LBFGSB <- vector("list",length(jm_pro_val_GH15_LBFGSB))
    names(jm_pro10_val_GH15_LBFGSB) <- names(jm_pro_val_GH15_LBFGSB)
    jm_pro10_val_GH15_LBFGSB$convergence <- 1
    jm_pro10_val_GH15_LBFGSB$iters <- 0
    jm_pro10_val_GH15_LBFGSB$logLik <- NA
    jm_pro10_val_GH15_LBFGSB$coefficients <- list(Dalpha=NA, alpha=NA)
    jm_pro10_val_GH15_LBFGSB$qNtime <- NA
}

jm_pro10_slo_GH15_LBFGSB <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro10_both_GH15_LBFGSB <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="optim",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)



## estimation by EM, scaling 1 :

jm_pro_val_GH15_EM <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro_slo_GH15_EM <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro_both_GH15_EM <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)


## estimation by EM, scaling 0.1 :

jm_pro01_val_GH15_EM <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro01_slo_GH15_EM <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro01_both_GH15_EM <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)


## estimation by EM, scaling 10 :

jm_pro10_val_GH15_EM <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro10_slo_GH15_EM <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)

jm_pro10_both_GH15_EM <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(GHk = 15,lng.in.kn = 1,iter.EM=1000,iter.qN=0),verbose = TRUE)


## estimation by marqLevAlg, scaling 1 :

jm_pro_val_GH15_marq <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_slo_GH15_marq <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro_both_GH15_marq <- jointModel(lmeFit, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by marqLevAlg, scaling 0.1 :

jm_pro01_val_GH15_marq <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_slo_GH15_marq <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro01_both_GH15_marq <- jointModel(lmeFit01, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## estimation by marqLevAlg, scaling 10 :

jm_pro10_val_GH15_marq <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro10_slo_GH15_marq <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="slope",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)

jm_pro10_both_GH15_marq <- jointModel(lmeFit10, survFit, timeVar = "time",method = "spline-PH-aGH",parameterization="both",derivForm=dFpro,control=list(optimizer="marq",nproc=3,GHk = 15,lng.in.kn = 1,iter.EM=0,iter.qN=1000),verbose = TRUE)


## create Table 4 of the manuscript
a <- expand.grid(c("", "01", "10"), c("val", "slo", "both"), "GH15", c("BFGS", "LBFGSB", "EM", "marq"), stringsAsFactors=FALSE)
a <- a[order(a[,2], decreasing=TRUE),]
mod <- apply(a, 1, function(x) paste("jm_pro", paste(x, collapse = "_"), sep = ""))


res <- matrix(NA, 27+9, 7)
table4 <- data.frame(dep=rep(c("value","slope","both"), each=12),
                   algorithm = a[,4], scaling = "1", loglik=0, value = 0, slope = 0,
                   iterations=0, time=0)

k <- 0
for(m in mod)
{
  k <- k + 1
  jm <- get(m)
  
  sc <- 1
  if(length(grep("pro10", m))) sc <- 10
  if(length(grep("pro01", m))) sc <- 0.1

  res[k, 1] <- jm$convergence
  res[k, 2] <- jm$iters
  res[k, 3] <- jm$logLik
  res[k, 4] <- jm$logLik + nrow(prothro) * log(sc)
  if(!length(jm$coefficients$Dalpha)) res[k, 7] <- NA else res[k, 7] <- as.numeric(jm$coefficients$Dalpha)
  if(!length(jm$coefficients$alpha)) res[k, 6] <- NA else  res[k, 6] <- as.numeric(jm$coefficients$alpha)
  if(!is.null(jm$qNtime)) res[k, 5] <- jm$qNtime[3]
  if(length(grep("EM", m))) res[k, 5] <- jm$EMtime

  table4$scaling[k] <- sc
  table4$loglik[k] <- res[k,4]
  table4$value[k] <- res[k,6]*sc
  table4$slope[k] <- res[k,7]*sc
  table4$iterations[k] <- res[k,2]
  table4$time[k] <- res[k,5]
}

colnames(res) <- c("convergence", "iterations", "LogLik", "scaledLogLik", "time", "alpha", "Dalpha")

table4$value[1:12] <- 100 * (table4$value[1:12] - table4$value[10]) / table4$value[10]
table4$slope[13:24] <- 100 * (table4$slope[13:24] - table4$slope[22]) / table4$slope[22]

table4$value[25:36] <- 100 * (table4$value[25:36] - table4$value[34]) / table4$value[34]
table4$slope[25:36] <- 100 * (table4$slope[25:36] - table4$slope[34]) / table4$slope[34]

table4


