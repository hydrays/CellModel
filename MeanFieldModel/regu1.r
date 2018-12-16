dt <- 0.01
Nstep <- 10000
k1 <- 0.1
lambda1 <- 1.0
lambda2 <- 0.5
k2 <- 0.2
Ncell <- 600
F0 <- 0.1
p0 <- 0.7
## F0 <- 0.8
## p0 <- 0.8
ec <- 1.6
sig1 <- 0
sig2 <- 0
alpha <- 0.5
sig <- 0.0

thetabarList <- NULL
for (nsample in seq(1))
{
e <- seq(1.2, 1.2, length.out=Ncell)
theta <- runif(Ncell, -pi/2, pi/2)

output <- NULL

thetabar <- 0
for ( i in seq(Nstep) )
{
    if ( i %% 20 == 0 )
    {
        n <- cbind(cos(theta), sin(theta))
        Q <- cbind(n[,1]*n[,1]-0.5, n[,1]*n[,2], n[,2]*n[,1], n[,2]*n[,2]-0.5)
        M <- matrix(colMeans(Q), nrow=2, ncol=2, byrow=TRUE)
        eigenvec <- eigen(M)$vectors[,1]
        thetabar <- atan(abs(eigenvec[2])/pmax(0.000001, abs(eigenvec[1])))
    }
    F <- max(0, F0 - alpha * sig/(1+abs(sig)))
    S <- mean(cos(2.0*theta))
    e <- e + dt*(F*cos(theta) - lambda1*(e-1)) + k1*sqrt(dt)*rnorm(Ncell, mean=0, sd=1)
    e <- pmax(e, 1)
    theta <- theta - dt*lambda2*(e-1)*(theta - thetabar) + k2*sqrt(dt)*rnorm(Ncell, mean=0, sd=1)
    ##theta <- theta - dt*S*lambda2*(theta-thetabar) + k2*sqrt(dt)*rnorm(Ncell, mean=0, sd=1)
    ##theta <- atan(abs(sin(theta))/pmax(0.000001, abs(cos(theta))))
    #cat(i, S, '\n')

    if ( i%%10 == 0 )
    {
        theta1 <- theta
        theta1 <- atan(abs(sin(theta1))/pmax(0.000001, abs(cos(theta1))))
        ecut <- e < ec
        theta2 <- theta1 * (1 - ecut) + ecut * runif(Ncell, 0, pi/2)
        A1 <- mean(sin(theta1) - cos(theta1))
        A2 <- mean(sin(theta2) - cos(theta2))
        ES <- mean((1-p0)*cos(theta1) - p0*sin(theta1))
        #cat(A1, A2, S, ES, F, sig, thetabar*180/pi, '\n')

        divtheta <- sample(theta1, 1)
        sig1 <- sig1 + cos(divtheta)
        sig2 <- sig2 + sin(divtheta)
        sig <- ((1-p0)*sig1 - p0*sig2)

        output <- rbind(output, c(sig, F, A1, A2, ES, sig1, sig2, divtheta))
    }
}

cat(A1, A2, S, ES, F, sig, thetabar*180/pi, '\n')
thetabarList <- c(thetabarList, thetabar*180/pi)
}
##png("output1.png")
pdf("output1.pdf")
par(mfrow=c(2, 1))
par(mar=c(4,6,4,4))
plot(output[200:800,1], type='o', col='blue', ylim=c(-6, 6),
     xlab="Number of divisions",cex.lab=1.5, cex.axis=1.5,
     ylab="Signal")
abline(h=0, lty=2)
##lines(output[,2], col='green', lty=1)
##lines(output[,3], col='blue')
plot(output[200:800,5], col='red')
dev.off()

sig0 <- mean(output[,1])
F0rec = F0 - alpha * sig0/(1+abs(sig0))
## theta1 <- theta
## theta1 <- atan(abs(sin(theta1))/pmax(0.000001, abs(cos(theta1))))
## ecut <- e < ec
## theta2 <- theta1 * (1 - ecut) + ecut * runif(Ncell, 0, pi/2)
## A1 <- mean(sin(theta1) - cos(theta1))
## A2 <- mean(sin(theta2) - cos(theta2))
## A1List <- c(A1List, A1)
## A2List <- c(A2List, A2)
## cat(A1, A2, S, thetabar*180/pi, '\n')

