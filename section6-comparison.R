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

