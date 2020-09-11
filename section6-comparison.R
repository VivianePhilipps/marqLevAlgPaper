### code of section 6 : comparison with other optimization algorithms ###


## prothro data from JM package
library(JM)

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


## create Table 4 of the manuscript
a <- expand.grid(c("", "01", "10"), c("val", "slo", "both"), "GH15", c("BFGS", "marq", "EM"), stringsAsFactors=FALSE)
mod <- apply(a, 1, function(x) paste("jm_pro", paste(x, collapse = "_"), sep = ""))

abis <- expand.grid(c("", "01", "10"), c("val"), "GH15", c("BFGS", "marq", "EM"), stringsAsFactors = FALSE)
modbis <- apply(abis, 1, function(x) paste("jm_pro", paste(x, collapse = "_"), sep = ""))
ater <- expand.grid(c("", "01", "10"), c("slo"), "GH15", c("BFGS", "marq", "EM"), stringsAsFactors = FALSE)
modter <- apply(ater, 1, function(x) paste("jm_pro", paste(x, collapse = "_"), sep = ""))

res <- matrix(NA, 27, 7)
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
  if(m %in% modbis) res[k, 7] <- NA
  else res[k, 7] <- as.numeric(jm$coefficients$Dalpha)
  if(m %in% modter) res[k, 6] <- NA
  else  res[k, 6] <- as.numeric(jm$coefficients$alpha)
  if(!is.null(jm$qNtime)) res[k, 5] <- jm$qNtime[3]
  if(!is.null(jm$EMtime)) res[k, 5] <- jm$EMtime
}

colnames(res) <- c("convergence", "iterations", "LogLik", "scaledLogLik", "time", "alpha", "Dalpha")

res2 <- data.frame(res, algorithm = 0, scaling = 0, association = 0, value = 0, slope = 0)

res2$algorithm[c(1:9)] <- "BFGS"
res2$algorithm[c(10:18)] <- "MLA"
res2$algorithm[c(19:27)] <- "EM"
res2$scaling[c(1, 4, 7, 10, 13, 16, 19, 22, 25)] <- "1"
res2$scaling[c(2, 5, 8, 11, 14, 17, 20, 23, 26)] <- "0.1"
res2$scaling[c(3, 6, 9, 12, 15, 18, 21, 24, 27)] <- "10"
res2$association[c(1:3, 10:12, 19:21)] <- "value"
res2$association[c(4:6, 13:15, 22:24)] <- "slope"
res2$association[c(7:9, 16:18, 25:27)] <- "both"

res2$alpha[c(2, 5, 8, 20, 23, 26, 11, 14, 17)] <- res2$alpha[c(2, 5, 8, 20, 23, 26, 11, 14, 17)] / 10
res2$alpha[c(3, 6, 9, 21, 24, 27, 12, 15, 18)] <- res2$alpha[c(3, 6, 9, 21, 24, 27, 12, 15, 18)] * 10
res2$Dalpha[c(2, 5, 8, 20, 23, 26, 11, 14, 17)] <- res2$Dalpha[c(2, 5, 8, 20, 23, 26, 11, 14, 17)] / 10
res2$Dalpha[c(3, 6, 9, 21, 24, 27, 12, 15, 18)] <- res2$Dalpha[c(3, 6, 9, 21, 24, 27, 12, 15, 18)] * 10


res2$value[c(1, 2, 3, 19, 20, 21, 10, 11, 12)] <- 100 * (res2$alpha[c(1, 2, 3, 19, 20, 21, 10, 11, 12)] - res2$alpha[10]) / res2$alpha[10]
res2$slope[c(1, 2, 3, 19, 20, 21, 10, 11, 12)] <- 100 * (res2$Dalpha[c(1, 2, 3, 19, 20, 21, 10, 11, 12)] - res2$Dalpha[10]) / res2$Dalpha[10]

res2$value[c(4, 5, 6, 22, 23, 24, 13, 14, 15)] <- 100 * (res2$alpha[c(4, 5, 6, 22, 23, 24, 13, 14, 15)] - res2$alpha[13]) / res2$alpha[13]
res2$slope[c(4, 5, 6, 22, 23, 24, 13, 14, 15)] <- 100 * (res2$Dalpha[c(4, 5, 6, 22, 23, 24, 13, 14, 15)] - res2$Dalpha[13]) / res2$Dalpha[13]

res2$value[c(7, 8, 9, 25, 26, 27, 16, 17, 18)] <- 100 * (res2$alpha[c(7, 8, 9, 25, 26, 27, 16, 17, 18)] - res2$alpha[16]) / res2$alpha[16]
res2$slope[c(7, 8, 9, 25, 26, 27, 16, 17, 18)] <- 100 * (res2$Dalpha[c(7, 8, 9, 25, 26, 27, 16, 17, 18)] - res2$Dalpha[16]) / res2$Dalpha[16]

table4 <- res2[c(1, 2, 3, 19, 20, 21, 10, 11, 12, 4, 5, 6, 22, 23, 24, 13, 14, 15, 7, 8, 9, 25, 26, 27, 16, 17, 18),
              c("association", "algorithm", "scaling", "scaledLogLik", "value", "slope", "iterations", "time")]

table4
