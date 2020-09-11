### code of section 4 : example ###


library(marqLevAlg)


## We first define the quantities to include as argument in loglikLMM function: 

Y <- dataEx$Y
X <- as.matrix(cbind(1,dataEx[,c("t","X1","X3")],dataEx$t*dataEx$X1))
ni <- as.numeric(table(dataEx$i))

binit <- c(0,0,0,0,0,1,1)


## maximum likelihood estimation of the linear mixed model in sequential mode :
estim <- marqLevAlg(b=binit, fn=loglikLMM, minimize=FALSE, X=X, Y=Y, ni=ni)
estim


## exact same model estimated in parallel mode using FORK implementation of parallelism : 
estim2 <- marqLevAlg(b=binit, fn=loglikLMM, minimize=FALSE, 
                     nproc = 2, clustertype = "FORK", X=X, Y=Y, ni=ni)


## estimation by using analytical gradients :
estim3 <- marqLevAlg(b=binit, fn=loglikLMM, gr=gradLMM, minimize=FALSE,
                     X=X, Y=Y, ni=ni)



## results :
res <- function(x){
    res <- data.frame(loglik=x$fn,
                      iterations=x$ni,
                      criterion1=x$ca,
                      criterion2= x$cb,
                      criterion3=x$rdm)
    rownames(res) <- deparse(substitute(x))
    return(t(res))
}

cbind(res(estim),res(estim2),res(estim3))


coef <- function(x){
    coef <- cbind(x$b,sqrt(x$v[c(1,3,6,10,15,21,28)]))
    colnames(coef) <- c(paste(deparse(substitute(x)),":b",sep=""),
                        paste(deparse(substitute(x)),":se",sep=""))
    rownames(coef) <- paste("coef",1:7)  
    return(round(coef,digits = 4))
}

## Table 1 in the manuscript
cbind(res(estim, 1, "no"), res(estim2, 2, "no"), res(estim3, 1, "yes"))

## Table 2 in the manuscript
cbind(coef(estim), coef(estim2), coef(estim3))
