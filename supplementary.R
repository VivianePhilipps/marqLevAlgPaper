#######################################################
############### supplementary materials ###############
#######################################################

library("marqLevAlg")
library("funconstrain")
library("parallel")
library("optimParallel")
library("minpack.lm")
library("nlmrt")
library("microbenchmark")
library("rgenoud")
library("DEoptim")
library("hydroPSO")
library("GA")


############### code of section 1 : examples from the litterature ###############


## 35 problems implemented in funconstrain package
pb <- list(rosen(), freud_roth(), powell_bs(), brown_bs(), beale(), jenn_samp(m=10), helical(), bard(), gauss(), meyer(), gulf(m=99), box_3d(m=20), powell_s(), wood(), kow_osb(), brown_den(m=20), osborne_1(), biggs_exp6(m=13), osborne_2(), watson(), ex_rosen(), ex_powell(), penalty_1(), penalty_2(), var_dim(), trigon(), brown_al(), disc_bv(), disc_ie(), broyden_tri(),broyden_band(), linfun_fr(m=100), linfun_r1(m=100), linfun_r1z(m=100), chebyquad())

## initial values
init <- lapply(pb,"[[", name="x0")

vi <- lapply(init, function(x) {if(class(x)=="function") do.call(x,list()) else x})
vi[[20]] <- rep(0, 6) # watson n=6
vi[[23]] <- 1:4 #penalty_1 n=4
vi[[24]] <- rep(0.5, 4) #penalty_2 n=4
n <- 35
vi[[28]] <- (1:n * (1/(n+1))) * ( (1:n * (1/(n+1))) -1)
vi[[29]] <- (1:n * (1/(n+1))) * ( (1:n * (1/(n+1))) -1)
vi[[35]] <- 1:7/8

## objective function's value at optimum
solutionf <- c(0, 48.9842, 0, 0, 0, 124.362, 0, 8.21487e-3, 1.12793e-8, 87.9458, 0, 0, 0, 0, 3.07505e-4, 85822.2, 5.46489e-5, 5.65565e-3, 4.01377e-2, 2.28767e-3,0, 0, 2.24997e-5, 9.37629e-6, 0, 0, 0, 0,0,0,0, 55, (100*99)/(2*(200+1)),(100^2+300-6)/(2*(200-3)), 0)

## set default cluster (needed for optimParallel)
cl <- parallel::makeCluster(1)
setDefaultCluster(cl)

## initialize the matrix of results
res35 <- matrix(NA, nrow=35, ncol=8)
rownames(res35) <- c("rosen", "freud_roth", "powell_bs", "brown_bs","beale", "jenn_samp", "helical", "bard", "gauss", "meyer", "gulf", "box_3d", "powell_s", "wood","kow_osb", "brown_den", "osborne_1", "biggs_exp6","osborne_2", "watson", "ex_rosen", "ex_powell", "penalty_1","penalty_2", "var_dim", "trigon", "brown_al", "disc_bv","disc_ie", "broyden_tri", "broyden_band", "linfun_fr", "linfun_r1", "linfun_r1z", "chebyquad")
colnames(res35) <- c("solution", "marqLevAlg", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "optimParallel", "nlminb")
res35[,1] <- solutionf

## run the optimizations
for(k in 1:35)
{
    resmla <- mla(b=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, maxiter=500)
    if(resmla$istop!=1) res35[k,2] <- NA else res35[k,2] <- resmla$fn.value

    resoptim <- optim(par=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, method="Nelder-Mead",control=list(maxit=500), hessian=TRUE)
    if(resoptim$convergence!=0) res35[k,3] <- NA else res35[k,3] <- resoptim$value
    rm(resoptim)
    
    resoptim <- optim(par=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, method="BFGS",control=list(maxit=500), hessian=TRUE)
    if(resoptim$convergence!=0) res35[k,4] <- NA else res35[k,4] <- resoptim$value
    rm(resoptim)
    
    resoptim <- optim(par=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, method="CG",control=list(maxit=500), hessian=TRUE)
    if(resoptim$convergence!=0) res35[k,5] <- NA else res35[k,5] <- resoptim$value
    rm(resoptim)
    
    try(resoptim <- optim(par=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, method="L-BFGS-B",control=list(maxit=500), hessian=TRUE), silent=TRUE)
    res35[k,6] <- NaN
    if(exists("resoptim"))
    {
        if(resoptim$convergence!=0) res35[k,6] <- NA else res35[k,6] <- resoptim$value
        rm(resoptim)
    }
    
    try(resoptp <- optimParallel(par=vi[[k]], fn=pb[[k]]$fn, gr=pb[[k]]$gr, control=list(maxit=500), hessian=TRUE, parallel=cl), silent=TRUE)
    res35[k,7] <- NaN
    if(exists("resoptp"))
    {
        if(resoptp$convergence!=0) res35[k,7] <- NA else res35[k,7] <- resoptp$value
        rm(resoptp)
    }
        
    resnlm <- nlminb(vi[[k]], pb[[k]]$fn, pb[[k]]$gr, control=list(iter.max=500))
    if(resnlm$convergence!=0) res35[k,8] <- NA else res35[k,8] <- resnlm$objective
    rm(resnlm)
}


## results : absolute biases
round(sweep(res35[,-c(1)], 1, res35[,1]),5) # in table (Table 1)
boxplot(log(1+round(sweep(res35[,-c(1)], 1, res35[,1]),5)), ylab="log(1 + bias)") # in boxplot (Figure 1)



############### code of section 2 : nonlinear least squares ###############

## examples from nlmrt package ##

traceval <- FALSE

## Problem in 1 parameter to ensure methods work in trivial case
set.seed(1)
nobs <- 8
tt <- seq(1,nobs)
dd <- 1.23*tt + 4*runif(nobs)
df <- data.frame(tt, dd)

runtime1 <- microbenchmark(
    n1 <- nlxb(dd ~ a*tt, start=c(a=1), data=df),
    m1 <- mla(b=1, fn=function(b,dd,tt){sum((b[1]*tt-dd)^2)}, dd=df$dd, tt=df$tt))


## Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972) # for testing
y  <-  ydat  # for testing
tdat  <-  seq_along(ydat) # for testing
weeddata1  <-  data.frame(y=ydat, tt=tdat)

eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
escal  <-   y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))

suneasy  <-  c(b1=200, b2=50, b3=0.3)
ssceasy  <-  c(b1=2, b2=5, b3=3)
st1scal  <-  c(b1=100, b2=10, b3=0.1)
start1  <-  c(b1=1, b2=1, b3=1)

shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
                        # This variant uses looping
     if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
     y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
              38.558, 50.156, 62.948, 75.995, 91.972)
     tt  <-  1:12
     res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}
      
shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
     jj  <-  matrix(0.0, 12, 3)
     tt  <-  1:12
     yy  <-  exp(-0.1*x[3]*tt)
     zz  <-  100.0/(1+10.*x[2]*yy)
     jj[tt,1]   <-   zz
     jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
     jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
     return(jj)
}

grhobbs <- function(b,y,tt) {
    .expr8 <- 1 + 10 * b[2] * exp(-0.1 * b[3] * tt)
    g1 <- 100/.expr8

    .expr1 <- 100 * b[1]
    .expr6 <- exp(-0.1 * b[3] * tt)
    .expr8 <- 1 + 10 * b[2] * .expr6
    g2 <- -(.expr1 * (10 * .expr6)/.expr8^2)

    .expr1 <- 100 * b[1]
    .expr2 <- 10 * b[2]
    .expr6 <- exp(-0.1 * b[3] * tt)
    .expr8 <- 1 + .expr2 * .expr6
    g3 <- .expr1 * (.expr2 * (.expr6 * (0.1 * tt)))/.expr8^2

    c(sum(-2 * g1 * (y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))),
      sum(-2 * g2 * (y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))),
      sum(-2 * g3 * (y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))))
}


runtime2 <- microbenchmark(
    n2 <- nlxb(eunsc, start=start1, trace=traceval, data=weeddata1),
    m2 <- mla(b=start1, fn=function(b,y,tt){sum((y - b[1]/(1+b[2]*exp(-b[3]*tt)))^2)}, y=weeddata1$y, tt=weeddata1$tt),
    n3 <- nlxb(eunsc, start=suneasy, trace=traceval, data=weeddata1),
    m3 <- mla(b=suneasy, fn=function(b,y,tt){sum((y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))^2)}, y=weeddata1$y, tt=weeddata1$tt),
    n4 <- nlxb(escal, start=start1, trace=traceval, data=weeddata1),
    m4 <- mla(b=start1, fn=function(b,y,tt){sum((y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))^2)}, y=weeddata1$y, tt=weeddata1$tt),
    n5 <- nlxb(escal, start=ssceasy, trace=traceval, data=weeddata1),
    m5 <- mla(b=ssceasy, fn=function(b,y,tt){sum((y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))^2)}, y=weeddata1$y, tt=weeddata1$tt),
    n6 <- nlxb(escal, start=st1scal, trace=traceval, data=weeddata1),
    m6 <- mla(b=st1scal, fn=function(b,y,tt){sum((y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))^2)}, y=weeddata1$y, tt=weeddata1$tt),
    n7 <- nlfb(start1, shobbs.res, shobbs.jac, trace=traceval),
    m7 <- mla(b=start1, fn=function(b,y,tt){sum((y - 100*b[1]/(1+10*b[2]*exp(-0.1*b[3]*tt)))^2)}, gr=grhobbs, y=weeddata1$y, tt=weeddata1$tt))


## Gabor Grothendieck problem
DF <- data.frame(x = c(5, 4, 3, 2, 1), y = c(1, 2, 3, 4, 5))

runtime3 <- microbenchmark(
    n8 <- nlxb(y ~ A * x + B, data = DF, start = c(A = 1, B = 6), trace=traceval),
    m8 <- mla(b=c(1,6), fn=function(b, x, y){sum((y - b[1]*x - b[2])^2)}, x=DF$x, y=DF$y))

meantime1 <- summary(runtime1, unit="us")[,4]
meantime2 <- summary(runtime2, unit="us")[,4]
meantime3 <- summary(runtime3, unit="us")[,4]

rssn <- sapply(paste("n", 1:8, sep=""), function(x){get(x)$ssquares})

rssm <- sapply(paste("m", 1:8, sep=""), function(x){
    m <- get(x)
    ifelse(m$istop==1, m$fn.value, NA)})

res <- matrix(cbind(rssn, c(meantime1[1], meantime2[seq(1,12,2)], meantime3[1]), rssm, c(meantime1[2], meantime2[seq(2,12,2)], meantime3[2])), 8, 4)
colnames(res) <- rep(c("objective function", "runtime"), 2)
rownames(res) <- c("One parameter problem", paste("Hobbs problem", c("unscaled - start1", "unscaled - easy", "scaled - start1", "scaled - easy", "scaled - hard", "scaled - start1 - gradient")), "Gabor Grothendieck problem")

res # Table 2



## examples from minpack.lm package ##

## exemple 1

set.seed(1)     
## values over which to simulate data 
x <- seq(0,5,length=100)
     
## model based on a list of parameters 
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c 
     
## parameter values used to simulate data
pp <- list(a=9,b=-1, c=6) 
     
## simulated data, with noise  
simDNoisy <- getPred(pp,x) + rnorm(length(x),sd=.1)
      
## residual function 
residFun <- function(p, observed, xx) observed - getPred(p,xx)
     
## starting values for parameters  
parStart <- list(a=3,b=-.001, c=1)
     
## for mla
residFunSum2 <- function(p, observed, xx) sum( (observed - getPredVec(p,xx) )^2)
getPredVec <- function(parS, xx) parS[1] * exp(xx * parS[2]) + parS[3] 
parStartVec <- c(3, -0.001, 1)


runtime1 <- microbenchmark(
 nls1 <- nls.lm(par=parStart, fn = residFun, observed = simDNoisy,
                xx = x, control = nls.lm.control(maxiter=500)),
 mla1 <- marqLevAlg(b = parStartVec, fn = residFunSum2,
                observed = simDNoisy, xx = x, maxiter=500))


## exemple 2

## function to simulate data 
f <- function(TT, tau, N0, a, f0) {
         expr <- expression(N0*exp(-TT/tau)*(1 + a*cos(f0*TT)))
         eval(expr)
}
     
## helper function for an analytical gradient 
j <- function(TT, tau, N0, a, f0) {
expr <- expression(N0*exp(-TT/tau)*(1 + a*cos(f0*TT)))
c(eval(D(expr, "tau")), eval(D(expr, "N0" )),
  eval(D(expr, "a"  )), eval(D(expr, "f0" )))
}
     
## values over which to simulate data 
TT <- seq(0, 8, length=501)
     
## parameter values underlying simulated data  
p <- c(tau = 2.2, N0 = 1000, a = 0.25, f0 = 8)
     
## get data 
Ndet <- do.call("f", c(list(TT = TT), as.list(p)))
## with noise
set.seed(1)
N <- Ndet +  rnorm(length(Ndet), mean=Ndet, sd=.01*max(Ndet))

## define a residual function 
fcn <- function(p, TT, N, fcall, jcall){
         (N - do.call("fcall", c(list(TT = TT), as.list(p))))
}     
 
## define analytical expression for the gradient 
fcn.jac <- function(p, TT, N, fcall, jcall){
         -do.call("jcall", c(list(TT = TT), as.list(p)))
}     

## starting values 
guess <- c(tau = 2.2, N0 = 1500, a = 0.25, f0 = 10)
     
gr <- function(b, TT, N){
    .expr3 <- exp(-TT/b[1])
    .expr8 <- 1 + b[3] * cos(b[4] * TT)
    g1 <- b[2] * (.expr3 * (TT/b[1]^2)) * .expr8

    .expr3 <- exp(-TT/b[1])
    .expr8 <- 1 + b[3] * cos(b[4] * TT)
    g2 <- .expr3 * .expr8

    .expr4 <- b[2] * exp(-TT/b[1])
    .expr6 <- cos(b[4] * TT)
    g3 <- .expr4 * .expr6
    
    .expr4 <- b[2] * exp(-TT/b[1])
    .expr5 <- b[4] * TT
    g4 <- -(.expr4 * (b[3] * (sin(.expr5) * TT)))

    c(sum(-2 * g1 * (N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))),
      sum(-2 * g2 * (N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))),
      sum(-2 * g3 * (N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))),
      sum(-2 * g4 * (N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))))    
}


runtime2 <- microbenchmark(
    nls2 <- nls.lm(par = guess, fn = fcn, jac = fcn.jac,fcall = f, jcall = j,TT = TT, N = N),
    mla2 <- mla(b=guess, fn=function(b, TT, N){sum((N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))^2)}, N=N, TT=TT, gr=gr),
    nls2gr <- nls.lm(par = guess, fn = fcn,fcall = f,TT = TT, N = N),
    mla2gr <- mla(b=guess, fn=function(b, TT, N){sum((N - b[2]*exp(-TT/b[1])*(1 + b[3]*cos(b[4]*TT)))^2)}, N=N, TT=TT))


meantime1 <- summary(runtime1, unit="us")[,4]
meantime2 <- summary(runtime2, unit="us")[,4]


res <- cbind(c(nls1$deviance, nls2$deviance, nls2gr$deviance), 
             c(meantime1[1], meantime2[c(1,3)]), 
             c(mla1$fn.value, mla2$fn.value, mla2gr$fn.value),
             c(meantime1[2], meantime2[c(2,4)]))
colnames(res) <- rep(c("objective function", "runtime"), 2)
rownames(res) <- c("Example1", "Example2", "Example2 - gradient")

res # Table 3


############### code of section 3 : other parallelized algorithms ###############

## sequential estimations
cl1 <- parallel::makeCluster(1)
clusterExport(cl1, list("fLMM","loglikLMM","X", "Y", "ni"))
setDefaultCluster(cl1)
set.seed(123)

time1proc <- microbenchmark(
    mla1 <- marqLevAlg(b = binit, fn = loglikLMM, minimize = FALSE, X = X, Y = Y, ni = ni),
    genoud1   <- genoud(loglikLMM, nvars=length(binit), max=TRUE, X = X, Y = Y, ni = ni),
    DEoptim1 <- DEoptim(fn=fLMM, lower=low, upper=up, Z=X, Y=Y, ni=ni),
    pso1 <- hydroPSO(binit, fn=fLMM, Y=Y, Z=X, ni=ni, lower=low, upper=up),
    ga1 <- ga(type="real-valued", fitness=loglikLMM, lower=low, upper=up, Y=Y, X=X, ni=ni),
    opt1 <- optimParallel(par=binit, fn=fLMM, control=list(maxit=500), hessian=TRUE, parallel=cl1, Y=Y, Z=X, ni=ni),
    times=10)


## in parallel mode
cl <- parallel::makeCluster(2)
clusterExport(cl, list("fLMM","loglikLMM","X", "Y", "ni"))
setDefaultCluster(cl)
set.seed(123)

time2proc <- microbenchmark(
    mla2 <- marqLevAlg(b = binit, fn = loglikLMM, minimize = FALSE, nproc = 2, X = X, Y = Y, ni = ni),
    genoud2   <- genoud(loglikLMM, nvars=length(binit), max=TRUE, X = X, Y = Y, ni = ni, cluster=cl),
    DEoptim2 <- DEoptim(fn=fLMM, lower=low, upper=up, Z=X, Y=Y, ni=ni, control=DEoptim.control(cluster=cl)),
    pso2 <- hydroPSO(binit, fn=fLMM, Y=Y, Z=X, ni=ni, lower=low, upper=up, control=list(parallel="parallel", par.nnodes=2)),
    ga2 <- ga(type="real-valued", fitness=loglikLMM, lower=low, upper=up, Y=Y, X=X, ni=ni, parallel=2),
    opt2 <- optimParallel(par=binit, fn=fLMM, control=list(maxit=500), hessian=TRUE, parallel=cl, Y=Y, Z=X, ni=ni),
    times=10)


## results
t1 <- summary(time1proc)
t2 <- summary(time2proc)
round(cbind(t1[,4], t2[,4], t1[,4] / t2[,4]), 2) # Table 4



############### code of section 4 : sensitivity to initial values ###############


## comparison with minpack.lm ##


## example 1 in the help of nls.lm function
set.seed(123)

## values over which to simulate data 
x <- seq(0,5,length=100)

## model based on a list of parameters 
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c 

## parameter values used to simulate data
pp <- list(a=9,b=-1, c=6) 
     
   
## residual function 
residFun <- function(p, observed, xx) observed - getPred(p,xx)

## sum squares residuals function
residFunSum2 <- function(p, observed, xx) sum( (observed - getPredVec(p,xx) )^2)

## model based on a list of parameters 
getPredVec <- function(parS, xx) parS[1] * exp(xx * parS[2]) + parS[3] 


## doone : data simulation and estimation with nls.lm and marqLevAlg with 100 different starting points
doone <- function(seed)
{
    set.seed(seed)
    
    ## simulated data, plus noise  
    simD <- getPred(pp,x) + rnorm(length(x),sd=.05)

    res <- matrix(NA,4,100)
    for(i in 1:100) {
        
        ## set starting values at random in interval -10 -- 10 
        param <- c(runif(1, -10, 10), runif(1, -10, 10), runif(1, -10, 10))
        
        ## starting values for parameters  
        parStart <- list(a=param[1],b=param[2], c=param[3])
        parStartVec <- c(a=param[1],b=param[2], c=param[3])
        
        ## perform fit 
        nls.out <- nls.lm(par = parStart, fn = residFun, observed = simD,
                          xx = x, control = nls.lm.control(nprint=0, maxiter=500))
        ml.out <- marqLevAlg(b = parStartVec, fn = residFunSum2,
                             observed = simD, xx = x, maxiter=500)
        
        res[1,i] <- nls.out$info
        res[2,i] <- nls.out$deviance
        res[3,i] <- ml.out$istop
        res[4,i] <- ml.out$fn.value
    }

    return(res)
}

## 100 different seeds
set.seed(123)
seeds <- round(runif(100, 0, 1000))

## run 100 replicates
res100 <- sapply(seeds, doone)


## results
convnls <- res100[seq(1,397,4),]
convmla <- res100[seq(3,399,4),]

table(convnls, convmla) # Table 5

fnls <- res100[seq(2,398,4),]
fmla <- res100[seq(4,400,4),]

ok1 <- length(which(fnls[which(convnls==1)]<1)) / length(which(convnls==1))
ok2 <- length(which(fnls[which(convnls==2)]<1)) / length(which(convnls==2))
ok3 <- length(which(fnls[which(convnls==3)]<1)) / length(which(convnls==3))
okm <- length(which(fmla[which(convmla==1)]<1)) / length(which(convmla==1))

par(mfrow=c(2,2), mgp=c(2,0.5,0), mar=c(5,4,2,2), tcl=-0.5, xpd=NA)
hist(log(1+fmla[which(convmla==1)]), ylim=c(0,1), freq=FALSE, breaks=seq(0,45,1), xlab="log(1+f)", main="marqLevAlg - istop=1 (N=5155)")
text(0.7, okm+0.05, paste(round(okm*100),"%", sep=""))
hist(log(1+fnls[which(convnls==1)]), ylim=c(0,1), breaks=seq(0,45,1), freq=FALSE, xlab="log(1+f)", main="nls.lm - info=1 (N=4730)")
text(0.7, ok1+0.05, paste(round(ok1*100,1),"%", sep=""))
text(6.6, 1-ok1+0.05, paste(round((1-ok1)*100,1),"%", sep=""))
hist(log(1+fnls[which(convnls==2)]), ylim=c(0,1), freq=FALSE, breaks=seq(0,45,1), xlab="log(1+f)", main="nls.lm - info=2 (N=1274)")
text(0.7, ok2+0.05, paste(round(ok2*100,1),"%", sep=""))
hist(log(1+fnls[which(convnls==3)]), ylim=c(0, 1), freq=FALSE, breaks=seq(0,45,1), xlab="log(1+f)", main="nls.lm - info=3 (N=594)")
text(0.7, ok3+0.05, paste(round(ok3*100,1),"%", sep=""))
text(6.5, 1-ok3+0.05, paste(round((1-ok3)*100,1),"%", sep="")) # Figure 2



## global optimization ##

## the Wild function
fw <- function (x){ 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x+80}
plot(fw, -50, 50, n = 1000, ylab="fw")
points(-15.81515, fw(-15.81515), col=2, lwd=3, pch=3) # Figure 4

## optimization using SANN
resSANN <- optim(50, fw, method = "SANN",
                 control = list(maxit = 20000, temp = 20, parscale = 20))

## optimization using marqLevAlg with grid search
res200 <- sapply(seq(-50,50, length.out=200), function(x){
    z <- mla(b=x, fn=fw)
    res <- c(NA, NA)
    if(z$istop==1) res <- c(z$fn.value, z$b)
    return(res)}) 

## results
res <- cbind(c(resSANN$value, resSANN$par), res200[,which.min(res200[1,])] )
colnames(res) <- c("SANN", "grid search MLA")
rownames(res) <- c("minimum", "param")

res # Table 6

