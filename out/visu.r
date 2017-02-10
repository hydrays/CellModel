require(lattice)
nodes <- read.csv(file="cell_nodes.txt", header=FALSE)
edges <- read.csv(file="cell_edges.txt", header=FALSE)
ellipses <- read.csv(file="ellipses100.txt", header=FALSE)
##new_ellipses <- read.csv(file="ellipses105.txt", header=FALSE)
old_ellipses <- read.csv(file="cellPos.txt", header=FALSE)
eigenshape <- read.csv(file="eigen.txt", header=FALSE)
partition_table <- as.matrix(read.csv(file="partition_table.txt", header=FALSE))
## plot(0, 0, type='n', xlim=c(-20, 20), ylim=c(-15, 15), asp=1)
## #plot(0, 0, type='n', xlim=c(-14, -8), ylim=c(-1, 5), asp=1)
## points(nodes[, 2:3])
## segments(edges[, 4], edges[, 5], edges[, 6], edges[, 7])
## points(new_ellipses[, 2:3], col='green', pch=15)
## points(ellipses[, 2:3], col='red', pch=15)

## ellipses <- read.csv(file="ellipses100.txt", header=FALSE)
## cat(dim(ellipses), '\n')
## plot(ellipses[,2], ellipses[,3], xlim=c(-30, 30), ylim=c(-20, 20), asp=1)
## ellipses <- read.csv(file="ellipses118.txt", header=FALSE)
## cat(dim(ellipses), '\n')
## points(ellipses[,2], ellipses[,3], col='red', pch=3)
## abline(h = c(-10, 10))
## abline(v = c(-15, 15))

plot(ellipses[,2], ellipses[,3], type='p', xlim=c(-30, 30), ylim=c(-20, 20), asp=1)
##plot(ellipses[,2], ellipses[,3], type='p', xlim=c(-15, -5), ylim=c(5, 15), asp=1)
points(old_ellipses[,1], old_ellipses[,2], col='red', pch=3)
segments(edges[, 4], edges[, 5], edges[, 6], edges[, 7])
#points(nodes[, 2:3])
force <- read.csv(file="force189.txt", header=FALSE)
sumforce <- read.csv(file="sum_force120.txt", header=FALSE)
subset <- force[,1]==103
subset2 <- sumforce[,1]==103
eigensubset <- eigenshape[,1]==1
eigensubset <- eigenshape[,1]>0
arrows(force[subset,2], force[subset,3], force[subset,4], force[subset,5],
       length=0.1, col='blue')
arrows(sumforce[subset2,2], sumforce[subset2,3], sumforce[subset2,4], sumforce[subset2,5], length=0.1, col='red')
arrows(eigenshape[eigensubset,2],
       eigenshape[eigensubset,3],
       eigenshape[eigensubset,6],
       eigenshape[eigensubset,7], length=0.1, col='green')
arrows(eigenshape[eigensubset,2],
       eigenshape[eigensubset,3],
       eigenshape[eigensubset,8],
       eigenshape[eigensubset,9], length=0.1, col='green')
abline(h = c(-10, 10))
abline(v = c(-15, 15))
