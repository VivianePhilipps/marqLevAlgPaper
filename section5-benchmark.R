### code of section 5 : benckmark ###

## This section of the manuscript presents a benchmark on 7 models with 100 replications.
## We give here the code for one replication.

#########################################
############### JM + mla ################
#########################################

library(JM)



#### load data and initial models on 5000 subjects ####

load("data_mla_5000.RData")
load("lmeFit.RData")
load("coxFit.RData")
JM <- NULL
sortie <- NULL

core <- c(1,2,3,4,6,8,10,15,20,25,30)

### numerical derivatives

sortie <- NULL

for(k in core) {

 JM <- jointModel(lmeObject = lmeFit,
               survObject = coxFit,
               timeVar = "t",
               parameterization = "value",
               method = "spline-PH-aGH",
               control = list(optimizer="marq",GHk = 3, lng.in.kn = 1,
                iter.EM=0,iter.qN=50,nproc=k, numeriDeriv=TRUE),
               verbose = TRUE)

sortie <- rbind(sortie,c(k,JM$CPUtime,JM$logLik,JM$iters,JM$convergence))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="JM_numeriDeriv.RData")
}



rm(sortie,JM,k)



### anaytical derivatives

sortie <- NULL

for(k in core) {

 JM <- jointModel(lmeObject = lmeFit,
               survObject = coxFit,
               timeVar = "t",
               parameterization = "value",
               method = "spline-PH-aGH",
               control = list(optimizer="marq",GHk = 3, lng.in.kn = 1,
                iter.EM=0,iter.qN=50,nproc=k, numeriDeriv=FALSE),
               verbose = TRUE)

sortie <- rbind(sortie,c(k,JM$CPUtime,JM$logLik,JM$iters,JM$convergence))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="JM_analyDeriv.RData")
}


rm(list=ls())




#########################################
############## hlme + mla ###############
#########################################

library(lcmmMLA)


## load simulated data set
load("data_mla.RData")

core <- c(1,2,3,4,6,8,10,15,20,25,30)


### one latent class

sortie <- NULL

for(k in core) {
    mG1 <- hlme(Y~1+I(t/10)+I((t/10)^2), random=~1+I(t/10)+I((t/10)^2),
                subject="i",ng=1,data=data_mla,verbose=FALSE,maxiter=30,nproc=k)

sortie <- rbind(sortie,c(k,mG1$time,mG1$loglik,mG1$ni,mG1$istop))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="hlme_G1.RData")
}


### 2 latent classes

b2 <- c(0, 0, 0, 0, 54.2713165484272, 51.9692401385805, 2.4708605906538, 
5.53807280570749, -8.82863376679862, -8.18035913550552, 329.06447101401, 
-108.035593294684, 153.357722977144, -11.2339964792991, -166.00077515535, 
244.634691319755, 1, 4.34225696290442)

sortie <- NULL

for(k in core) {
    mG2 <- hlme(Y~1+I(t/10)+I((t/10)^2), random=~1+I(t/10)+I((t/10)^2),
                mixture=~1+I(t/10)+I((t/10)^2),nwg=TRUE,classmb=~X1+X2+X3,
                subject="i",ng=2,data=data_mla,verbose=FALSE,maxiter=30,nproc=k,B=b2)

sortie <- rbind(sortie,c(k,mG2$time,mG2$loglik,mG2$ni,mG2$istop))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="hlme_G2.RData")
}



### 3 latent classes


b3 <- c(0, 0, 0, 0, 0, 0, 0, 0, 53.0945627380509, 53.2399420396394, 
54.4269074001362, 5.26465623011422, 3.79143012314735, 6.23839975066632, 
-9.44143356198573, -10.4950552552757, -10.1529978271412, 292.751768720462, 
-80.2567239174774, 132.206422427329, -9.92168321839977, -152.906018466692, 
220.23371919073, 1, 1, 4.40372813987142)

sortie <- NULL

for(k in core) {
    mG3 <- hlme(Y~1+I(t/10)+I((t/10)^2), random=~1+I(t/10)+I((t/10)^2),
                mixture=~1+I(t/10)+I((t/10)^2),nwg=TRUE,classmb=~X1+X2+X3,
                subject="i",ng=3,data=data_mla,verbose=FALSE,maxiter=30,nproc=k,B=b3)

sortie <- rbind(sortie,c(k,mG3$time,mG3$loglik,mG3$ni,mG3$istop))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="hlme_G3.RData")
}




### 4 latent classes

b4 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52.8490592741525, 55.6394754968402, 
54.1545598758033, 52.9982705556955, 3.35731963419959, 6.97751925567454, 
4.9458820485929, 5.14394239823614, -6.64671457489409, -10.4877076311514, 
-8.03560028321582, -9.41373051406899, 248.988298845135, -81.2369505771433, 
145.484846690813, -3.58765944283958, -172.145995642558, 252.738512339079, 
1, 1, 1, 4.5231699453073)

sortie <- NULL

for(k in core) {
    mG4 <- hlme(Y~1+I(t/10)+I((t/10)^2), random=~1+I(t/10)+I((t/10)^2),
                mixture=~1+I(t/10)+I((t/10)^2),nwg=TRUE,classmb=~X1+X2+X3,
                subject="i",ng=4,data=data_mla,verbose=FALSE,maxiter=30,nproc=k,B=b4)

sortie <- rbind(sortie,c(k,mG4$time,mG4$loglik,mG4$ni,mG4$istop))
colnames(sortie) <- c("core","CPUtime","loglik","iters","convergence")
#save(sortie,file="hlme_G4.RData")
}


rm(list=ls())






#########################################
############# CInLPN + mla ##############
#########################################


library(CInLPN)

core <- c(1,2,3,4,6,8,10,15,20,25,30)


epsa <- 0.0001
epsb <- 0.0001
epsd <- 0.0001
paras.ini <- NULL
paras.ini <- c(0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,rep(0,9),rep(1,9))

indexparaFixeUser <- c(1,2,3,6+c(1, 2,3,5,6, 7, 8,9,11, 12, 13,14,17,18,20))
paraFixeUser <- c(0,0,0,1,rep(0,4),1,rep(0,3),1,rep(0,5))

DeltaT <- 1


sortie <- NULL

for(k in core) {
    mod <- CInLPN(structural.model = list(fixed.LP0 = ~1|1|1,fixed.DeltaLP = L1 | L2 | L3 ~ 1|1|1,random.DeltaLP = ~1|1|1,trans.matrix = ~1,delta.time = DeltaT),
                   measurement.model = list(link.functions = list(links = c(NULL, NULL, NULL),knots = list(NULL, NULL, NULL))),
                   parameters = list(paras.ini = paras.ini, Fixed.para.index = indexparaFixeUser, Fixed.para.values = paraFixeUser),
                   option = list(parallel = TRUE, nproc = k, print.info = FALSE, epsa = epsa, epsb = epsb, epsd = epsd,maxiter=50,univarmaxiter = 7, makepred = FALSE),
                  Time = "time", subject = "id", data = data)

sortie <- rbind(sortie,c(k,mod$MLAtime,mod$loglik,mod$niter,mod$conv))
colnames(sortie) <- c("core","CPUtime","logLik","iters","convergence")
#save(sortie,file="CInLPN.RData")
}
    
}



q("no")
