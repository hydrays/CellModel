convertThetaIter  <- function(theta, thetabar)
{
    dtheta  <- theta - thetabar
    theta[dtheta > .5*pi]  <-  theta[dtheta > .5*pi] - pi
    theta[dtheta < -.5*pi]  <- pi + theta[dtheta < -.5*pi]
    if ( any(abs(dtheta) > .5*pi) )
    {
        theta <- convertThetaNew(theta, thetabar)
    }
    return(theta)
}

convertThetaNew  <- function(theta, thetabar)
{
    dtheta  <- theta - thetabar
    if ( any(abs(dtheta) > pi ) )
    {
        cat("something wrong with theta!\n")
        exit()
    }
    theta[dtheta > .5*pi]  <-  theta[dtheta > .5*pi] - pi
    theta[dtheta < -.5*pi]  <- pi + theta[dtheta < -.5*pi]
    if ( any(abs(dtheta) > pi) )
    {
        cat("something wrong with theta 2!\n")
        exit()
    }
    return(theta)
}

convertTheta  <- function(theta)
{
    if ( any(abs(theta) > pi ) )
    {
        cat("something wrong with theta!\n")
        exit()
    }
    theta[theta > .5*pi]  <-  theta[theta > .5*pi] - pi
    theta[theta < -.5*pi]  <- pi + theta[theta < -.5*pi]
    if ( any(abs(theta) > pi) )
    {
        cat("something wrong with theta 2!\n")
        exit()
    }
    return(theta)
}

simu  <- function(lambda2, k1, k2, alpha, ec, A, index)
{
    dt <- 0.01
    a0 <- 0.01
    e0 <- 1.0
    lambda1 <- 1
    Ncell <- 50
    Nstep <- 40000
    v <- a0*A
    DispWindow  <- seq(28000, 40000)
    ##Selection  <- 28000 + 2400*seq(5)
    Selection  <- 28000 + 100*seq(120)

    nt <- Ncell
    lt <- 1
    output1 <- NULL

    e <- seq(1.2, 1.2, length.out=Ncell)
    theta <- runif(Ncell, -pi/2, pi/2)
    thetabar <- 0
    S <- 0
    
    ## ############### model I
    for ( i in seq(Nstep) )
    {
        if ( i %% 20 == 0 )
        {
            n <- cbind(cos(theta), sin(theta))
            Q <- cbind(n[,1]*n[,1]-0.5, n[,1]*n[,2], n[,2]*n[,1], n[,2]*n[,2]-0.5)
            Q <- colSums((e-1)*Q)/sum((e-1))
            M <- matrix(Q, nrow=2, ncol=2, byrow=TRUE)
            ##M <- matrix(colMeans(Q), nrow=2, ncol=2, byrow=TRUE)
            eigenvec <- eigen(M)$vectors[,1]
            S <- 2.0*eigen(M)$values[1]
            if(eigenvec[1] < 0)
            {
                eigenvec  <- -eigenvec
            }
            thetabar <- sign(eigenvec[2])*atan(abs(eigenvec[2])/pmax(0.000001, abs(eigenvec[1])))
        }
        Lt <- exp(v*dt*i)
        Force <- alpha*(Lt-lt)
        e <- e + dt*(Force*abs(cos(theta)) - lambda1*(e-1)) + k1*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        e <- pmax(e, 1)
        theta <- theta - dt*lambda2*S*(e-1)*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        ##theta <- theta - dt*lambda2*S*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        theta <- convertThetaIter(theta, thetabar)
        ##theta <- atan2(sin(theta), cos(theta))
        
        theta1 <- atan(abs(sin(theta))/pmax(0.000001, abs(cos(theta))))
        ndiv=rbinom(1, nt, dt*a0)
        divtheta <- sample(theta1, ndiv)
        lt<-lt*(1+sum(cos(divtheta))/nt)
        nt<-nt+ndiv
        theta<-c(theta, runif(ndiv, -pi/2, pi/2))
        e<-c(e, rep(e0, ndiv))
        A1 <- mean(cos(theta1))
        A2 <- mean(sin(theta1)) - mean(cos(theta1))

        output1 <- rbind(output1, c(Lt, lt, A1, A2, Force, nt, thetabar, S))
    }

    eModelI <- e
    thetaModelI <- theta
    
################################ model II
    output2 <- NULL
    nt=Ncell
    lt=1
    e <- seq(1.2, 1.2, length.out=Ncell)
    theta <- runif(Ncell, -pi/2, pi/2)

    for ( i in seq(Nstep) )
    {
        if ( i %% 20 == 0 )
        {
            n <- cbind(cos(theta), sin(theta))
            Q <- cbind(n[,1]*n[,1]-0.5, n[,1]*n[,2], n[,2]*n[,1], n[,2]*n[,2]-0.5)
            Q <- colSums((e-1)*Q)/sum((e-1))
            M <- matrix(Q, nrow=2, ncol=2, byrow=TRUE)
            ##M <- matrix(colMeans(Q), nrow=2, ncol=2, byrow=TRUE)
            eigenvec <- eigen(M)$vectors[,1]
            S <- 2.0*eigen(M)$values[1]
            if(eigenvec[1] < 0)
            {
                eigenvec  <- -eigenvec
            }
            thetabar <- sign(eigenvec[2])*atan(abs(eigenvec[2])/pmax(0.000001, abs(eigenvec[1])))
        }
        Lt <- exp(v*dt*i)
        Force <- alpha*(Lt-lt)
        e <- e + dt*(Force*abs(cos(theta)) - lambda1*(e-1)) + k1*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        e <- pmax(e, 1)
        theta <- theta - dt*lambda2*S*(e-1)*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        ##theta <- theta - dt*lambda2*S*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        theta <- convertThetaIter(theta, thetabar)
        ##theta <- atan2(sin(theta), cos(theta))
        
        theta1 <- atan(abs(sin(theta))/pmax(0.000001, abs(cos(theta))))
        ecut <- e < ec
        theta2 <- theta1 * (1 - ecut) + ecut * runif(nt, 0, pi/2)
        
        ndiv=rbinom(1, nt, dt*a0)
        
        divtheta <- sample(theta2, ndiv)
        lt<-lt*(1+sum(cos(divtheta))/nt)
        nt<-nt+ndiv
        theta<-c(theta, runif(ndiv, -pi/2, pi/2))
        e<-c(e, rep(e0, ndiv))
        A1 <- mean(cos(theta1))
        A2 <- mean(sin(theta1)) - mean(cos(theta1))
        
        output2<-rbind(output2, c(Lt, lt, A1, A2, Force, nt, thetabar, S))
    }


    ##if ( index <= 3 )
    {
        png(paste("plotlt", "_A_", A, "_alpha_", alpha,
                  "_lambda2_", lambda2,
                  "_k1_", k1, "_k2_", k2,
                  "_", index, ".png", sep=''),
            width=1200, height=600)
        par(mfrow=c(2,5))
        plot(eModelI, thetaModelI*180/pi, ylim=c(-90,90))
        abline(h=c(-30,30))
        plot(output1[DispWindow, 1]-output1[DispWindow, 2],
             col = 'blue', ylab="L(t)-l(t)",
             ylim=c(-1,1),
             xlab="timesteps", cex=0.5
             )##, ylim=c(-1,1), cex=0.5)
        abline(h=0, lwd = 2, col='grey')
        plot(output1[DispWindow,5], cex=0.25)
        plot(output1[DispWindow, 3], col='blue', ylim=c(-1, 1))
        lines(output1[DispWindow, 4], col='red')
        lines(output2[DispWindow, 7], col='green', pch=5, type='p')        
        abline(h=c(-0.3))
        plot(output1[DispWindow, 8], col='blue', ylim=c(0, 1))        

        plot(e, theta*180/pi, ylim=c(-90,90))
        abline(h=c(-30,30))
        plot(output2[DispWindow, 1]-output2[DispWindow, 2],
             col = 'blue', ylab="L(t)-l(t)",
             ylim=c(-1,1),
             xlab="timesteps", cex=0.5
             )
        abline(h=0, lwd = 2, col='grey')
        plot(output2[DispWindow,5], cex=0.25)
        plot(output2[DispWindow, 3], col='blue', ylim=c(-1, 1))
        lines(output2[DispWindow, 4], col='red')
        lines(output2[DispWindow, 7], col='green', pch=5, type='p')
        abline(h=c(-0.3))
        plot(output2[DispWindow, 8], col='blue', ylim=c(0, 1))        

        dev.off()
    }

    simu <- NULL
    simu$model1 <- output1[Selection, 1] - output1[Selection, 2]
    simu$model2 <- output2[Selection, 1] - output2[Selection, 2]
    return(simu)
}

baseModels  <- function(lambda2, k1, k2, ec, index, Force)
{
    dt <- 0.01
    a0 <- 0.01
    e0 <- 1.0
    lambda1 <- 1
    Ncell <- 2000
    Nstep <- 10000

    nt <- Ncell
    e <- seq(1.2, 1.2, length.out=Ncell)
    theta <- 0.3+runif(Ncell, -pi/2, pi/2)
    thetabar <- 0.3
    S  <- 0
    
    ## ############### model I
    for ( i in seq(Nstep) )
    {
        if ( i %% 20 == 0 )
        {
            n <- cbind(cos(theta), sin(theta))
            Q <- cbind(n[,1]*n[,1]-0.5, n[,1]*n[,2], n[,2]*n[,1], n[,2]*n[,2]-0.5)
            M <- matrix(colMeans(Q), nrow=2, ncol=2, byrow=TRUE)
            eigenvec <- eigen(M)$vectors[,1]
            S <- 2.0*eigen(M)$values[1]
            if(eigenvec[1] < 0)
            {
                eigenvec  <- -eigenvec
            }
            thetabar <- sign(eigenvec[2])*atan(abs(eigenvec[2])/pmax(0.000001, abs(eigenvec[1])))
            ## cat(S, thetabar, '\n')
        }
        e <- e + dt*(Force*cos(theta) - lambda1*(e-1)) + k1*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        e <- pmax(e, 1)
        ##theta <- theta - dt*lambda2*S*(e-1)*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        ##theta <- theta - dt*lambda2*S*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        theta <- theta - dt*lambda2*(e-1)*S*(theta-thetabar) + k2*sqrt(dt)*rnorm(nt, mean=0, sd=1)
        theta <- convertThetaIter(theta, thetabar)
        ##theta <- atan2(sin(theta), cos(theta))
    }
    
    theta1 <- atan(abs(sin(theta))/pmax(0.000001, abs(cos(theta))))
    ecut <- e < ec
    theta2 <- theta1 * (1 - ecut) + ecut * runif(nt, 0, pi/2)
    A1 <- mean(cos(theta1))
    A2 <- mean(sin(theta1)) - mean(cos(theta1))
    B1 <- mean(cos(theta2))
    B2 <- mean(sin(theta2)) - mean(cos(theta2))

    output <- NULL
    output$values <- c(A1, A2, B1, B2, S, thetabar)
    output$e <- e
    output$theta <- theta
        
    baseModels  <- output
    return(baseModels)
}

library(foreach)
library(doParallel)

lambda2 <- 0.1
k1 <- 0.2
k2 <- 0.1
alpha <- 2
ec <- 1.6
A  <- 0.8

Alist <- seq(0.6, 1, 0.05)
lambda2list <- seq(0.1, 3, 0.2)
eclist <- seq(1.1, 2, 0.1)
alphalist <- c(0.1, 0.2, 0.3, 0.5, 1, 1.5, 2, 5, 10, 20, 40, 60, 100)
k1list <- seq(0.1, 2, 0.2)
k2list <- seq(0.1, 2, 0.2)

singleRun <- function(Nsample)
{
    resAll <- NULL
    for ( i in seq(Nsample))
    {
        cat(i, '\n')
        res <- simu(lambda2, k1, k2, alpha, ec, A, i)
        resAll$model1 <- c(resAll$model1, res$model1)
        resAll$model2 <- c(resAll$model2, res$model2)
    }
    singleRun <- c(sd(resAll$model1), sd(resAll$model2))
    return(singleRun)
}

paraRun <- function(Nsample)
{
    paraRun <- foreach ( i=seq(Nsample), .export=c("convertTheta", "simu", "lambda2", "k1", "k2", "alpha", "ec", "A"), .combine = list, .multicombine = TRUE ) %dopar%
        simu(lambda2, k1, k2, alpha, ec, A, i)
    return(paraRun)
}

computeSDpara  <- function(nthreads, Nsample)
{
    zz <- file("sd.txt", "w")
    close(zz)
    
    cl<-makeCluster(nthreads)
    registerDoParallel(cl)
    ##for ( k1 in k1list )
    {
        ##for ( k2 in k2list )
        {
            res <- paraRun(Nsample)
            sd1 <- sd(unlist(lapply(res, "[[", 1)))
            sd2 <- sd(unlist(lapply(res, "[[", 2)))
            cat(lambda2, k1, k2, alpha, ec, A, sd1, sd2, '\n',
                file="sd.txt", sep=',', append=TRUE)
        }
    }
    stopCluster(cl)
}

figure2 <- function(Nsample)
{
    mycol = c("#99CCFF","darkgrey","#FFCCCC","#FF99FF","#99FFCC")
    ##mycol = c("#6666FF","darkgrey","#FFCC99","#FF99FF","#99FFCC")
    mylty = c(1, 1, 1, 1, 1)
    mypch = c(0,1,5,16,2)
    
    pdf("figure2A.pdf", height=7, width=9)
    par(mar=c(4,6,4,5))
    plot(0, 0, type='n', 
         col = 'blue', ylab="L(t)-l(t)",
         ylim=c(-.75,.75),
         xlim=c(0, 120),
         xlab="t", cex.lab=2, cex.axis=1.5
         )
    for ( i in seq(Nsample))
    {
        cat(i, '\n')
        res <- simu(lambda2, k1, k2, alpha, ec, A, i)
        lines(res$model1, pch=mypch[i], lty=mylty[i], col=mycol[i], type='o', lwd=1.5, cex=0.6)
    }
    dev.off()

    pdf("figure2B.pdf", height=7, width=9)
    par(mar=c(4,6,4,5))
    plot(0, 0, type='n', 
         col = 'blue', ylab="L(t)-l(t)",
         ylim=c(-.75,.75),
         xlim=c(0, 120),
         xlab="t", cex.lab=2, cex.axis=1.5
         )
    for ( i in seq(Nsample))
    {
        cat(i, '\n')
        res <- simu(lambda2, k1, k2, alpha, ec, A, i)
        lines(res$model2, pch=mypch[i], lty=mylty[i], col=mycol[i], type='o', lwd=1.5, cex=0.6)
    }
    dev.off()
}

## Force <- 0.5
## Slist  <- NULL
## ##for (lambda2 in lambda2list )
## {
##     res <- baseModels(lambda2, k1, k2, ec, i, Force)
##     Slist  <- c(Slist, res$values[5])
##     cat(res$values, '\n')
## }

figure2(5)

##computeSDpara(2, 2)

##system("cat run.r > sd.txt")

