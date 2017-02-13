require(lattice)
#cellIDset = c(56)##, 70, 30, 20)
cellIDset = NULL
edgesOld <- read.csv(file=paste("cell_edges", 100, ".txt", sep=''), header=FALSE)
ellipsesOld <- read.csv(file=paste("ellipses", 100, ".txt", sep=''), header=FALSE)
r0 <- 1.6

for (file_index in seq(100, 1100) )
{
    cat(file_index, "\n")
    nodes <- read.csv(file=paste("cell_nodes", file_index, ".txt", sep=''), header=FALSE)
    edges <- read.csv(file=paste("cell_edges", file_index, ".txt", sep=''), header=FALSE)
    force <- read.csv(file=paste("force", file_index, ".txt", sep=''), header=FALSE)
    force2 <- read.csv(file=paste("force2", file_index, ".txt", sep=''), header=FALSE)
    ellipses <- read.csv(file=paste("ellipses", file_index, ".txt", sep=''), header=FALSE)
    sumforce <- read.csv(file=paste("sum_force", file_index, ".txt", sep=''), header=FALSE)
    eigenshape <- read.csv(file=paste("eigen", file_index, ".txt", sep=''), header=FALSE)
    subset <- force[,1] %in% cellIDset
    subset2 <- sumforce[,1] %in% cellIDset
    eigensubset <- eigenshape[,1]==1
    eigensubset <- eigenshape[,1]>0
    ratio <- ellipses[,6]/ellipses[,7]
    p0 <- sum(ratio>r0)/length(ratio)
    
    png(paste("plot", file_index, ".png", sep=''),
        height=600, width=800)
    plot(0, 0, type='n', xlim=c(-30, 30), ylim=c(-20, 20), asp=1,
         main=paste("p0 =", format(round(p0, 2), nsmall = 2))) 
    ##points(nodes[,2], nodes[,3])
    segments(edgesOld[, 4], edgesOld[, 5], edgesOld[, 6], edgesOld[, 7], col='gray')
    segments(edges[, 4], edges[, 5], edges[, 6], edges[, 7], col='black')
    arrows(force[subset,2], force[subset,3], force[subset,4], force[subset,5],
           length=0.1, col='blue')
    arrows(force2[subset,2], force2[subset,3], force2[subset,4], force2[subset,5],
           length=0.1, col='red')
    arrows(sumforce[subset2,2], sumforce[subset2,3], sumforce[subset2,4], sumforce[subset2,5], length=0.1, col='red')
    eigensubset <- eigenshape[,1]==1
    eigensubset <- eigenshape[,1]>0
    points(ellipses[ratio>r0, 2:3], col='gold', pch=16, cex=2)
    points(ellipses[ratio<=r0, 2:3], col='blue', pch=16, cex=2)
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
    abline(h = c(-10, 10))
    abline(v = c(-15, 15))
    dev.off()
}
