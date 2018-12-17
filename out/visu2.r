require(lattice)
cellIDset = c(56, 70, 30, 20)
##cellIDset = NULL
edgesOld <- read.csv(file=paste("cell_edges", 0, ".txt", sep=''), header=FALSE)
ellipsesOld <- read.csv(file=paste("ellipses", 0, ".txt", sep=''), header=FALSE)
r0 <- 1.53

ratioAll <- NULL

##for (file_index in seq(170, 172) )
for (file_index in seq(0,10100,by=10) )
{
    cat(file_index, "\n")
    nodes <- read.csv(file=paste("cell_nodes", file_index, ".txt", sep=''), header=FALSE)
    edges <- read.csv(file=paste("cell_edges", file_index, ".txt", sep=''), header=FALSE)
    ##force <- read.csv(file=paste("force", file_index, ".txt", sep=''), header=FALSE)
    ##force2 <- read.csv(file=paste("force2", file_index, ".txt", sep=''), header=FALSE)
    ellipses <- read.csv(file=paste("ellipses", file_index, ".txt", sep=''), header=FALSE)
    ##sumforce <- read.csv(file=paste("sum_force", file_index, ".txt", sep=''), header=FALSE)
    ##eigenshape <- read.csv(file=paste("eigen", file_index, ".txt", sep=''), header=FALSE)
    ##subset <- force[,1] %in% cellIDset
    ##subset2 <- sumforce[,1] %in% cellIDset
    ##eigensubset <- eigenshape[,1]==1
    ##eigensubset <- eigenshape[,1]>0

    type <- ellipses[,12]
    ratio <- ellipses[,8]/ellipses[,9]
    ratioAll <- c(ratioAll, sum(ratio>r0)/length(ratio))
    
    ## png(paste("ar", file_index, ".png", sep=''), height=600, width=800)
    ## ##pdf(paste("ar", file_index, ".pdf", sep=''), height=7, width=7)
    ## par(mar=c(5,6,5,5))
    ## ## hist(ellipses[type<10,6]/ellipses[type<10,7], breaks=15, freq=FALSE,
    ## ##      xlim=c(1,2.6),
    ## ##      xlab="aspect ratio", ylab="frequency", cex.lab=2, cex.axis=1.5,
    ## ##      ##main=paste("mean", mean(ratio), "var", var(ratio)))
    ## ##      main="")
    ## ar <- ellipses[type<10,6]/ellipses[type<10,7]
    ## ##ar <- ar + rnorm(length(ar), 0, 0.1)
    ## F <- ecdf(ar)
    ## plot(F, xlim=c(1,2.6))
    ## abline(v=1.53, col='red', lwd=4, lty=2)
    ## abline(h=sum(ar<1.53)/length(ar), col='blue', lwd=4, lty=2)
    ## ##cat(file="meanar.txt", mean(ellipses[,6]/ellipses[,7]), '\n', append=T)
    ## dev.off()
    
    png(paste("plot", file_index, ".png", sep=''),
        height=600, width=1200)
    par(mfrow=c(1,2))
    ##pdf(paste("plot", file_index, ".pdf", sep=''),
    ##    height=6, width=8)
    plot(0, 0, type='n', xlim=c(-30, 30), ylim=c(-30, 30), asp=1, axes=FALSE,
         ##main=paste("p1 =", format(round(p0, 2), nsmall = 2))
         ##main="", cex.lab=1.5, cex.axis=1.5, xlab="longitudinal direction", ylab="circumferential direction")
         main="", cex.lab=1.5, cex.axis=1.5, xlab="", ylab="") 
    ##points(nodes[,2], nodes[,3])
    ##segments(edgesOld[, 4], edgesOld[, 5], edgesOld[, 6], edgesOld[, 7], col='gray')
    segments(edges[, 4], edges[, 5], edges[, 6], edges[, 7], col='black')

    ## arrows(force[subset,2], force[subset,3], force[subset,4], force[subset,5],
    ##        length=0.1, col='blue')
    ## arrows(force2[subset,2], force2[subset,3], force2[subset,4], force2[subset,5],
    ##        length=0.1, col='red')
    ## arrows(sumforce[subset2,2], sumforce[subset2,3], sumforce[subset2,4], sumforce[subset2,5], length=0.1, col='red')
    
    ##eigensubset <- eigenshape[,1]==1
    ##eigensubset <- eigenshape[,1]>0
    ##points(ellipses[ratio>r0, 2:3], col='gold', pch=16, cex=2)
    ##points(ellipses[ratio<=r0, 2:3], col='blue', pch=16, cex=2)
    points(ellipses[(type<10) & (ratio>r0), 2:3], col='darkgray', pch=1, cex=1)
    points(ellipses[(type<10) & (ratio<=r0), 2:3], col='black', pch=16, cex=0.75)
    points(ellipses[type>=10, 2:3], col='red', pch=16, cex=0.75)
    ## arrows(eigenshape[eigensubset,2],
    ##        eigenshape[eigensubset,3],
    ##        eigenshape[eigensubset,6],
    ##        eigenshape[eigensubset,7], length=0.1, col='blue')
    ## arrows(eigenshape[eigensubset,2],
    ##        eigenshape[eigensubset,3],
    ##        eigenshape[eigensubset,8],
    ##        eigenshape[eigensubset,9], length=0.1, col='red')
    arrows(ellipses[ratio>r0,2],
           ellipses[ratio>r0,3],
           ellipses[ratio>r0,10],
           ellipses[ratio>r0,11], length=0.1, col='gold')
    arrows(ellipses[ratio<=r0,2],
           ellipses[ratio<=r0,3],
           ellipses[ratio<=r0,10],
           ellipses[ratio<=r0,11], length=0.1, col='blue')
    ##abline(h = c(Hmin, Hmax))
    ##abline(v = c(Lmin, Lmax))
    par(mar=c(14, 6, 12, 8))
    plot(ratioAll, col='black', xlim=c(0, 2000), ylim=c(0,1),
         xlab="time", ylab="percentage of elongated cells", cex.axis=1.5, cex.lab=1.5, type='l', lwd=2)
    abline(h=c(0.2, 0.4, 0.6, 0.8), lty=4, col='gray')
    points(file_index+1, ratioAll[file_index+1], col='red', cex=1.5, pch=16)
    dev.off()
}
