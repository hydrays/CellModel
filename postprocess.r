require(lattice)
output_file_prefix <- "basemodel2_a135"

A = NULL
rep <- 0
force <- c(1, 3, 7, 10, 14, 16)
for ( i in seq(6) )
{
    rep = rep + 1
    cat("\n ***********run********* ", rep, "\n")
    FolderName <- paste("case_", output_file_prefix, "_", as.character(rep), "/out/", sep='')
    cmd = paste("ls ", FolderName, "stat* > ", FolderName, "filelist.txt", sep='')
    system(cmd)
    fileList <- readLines(paste(FolderName, "filelist.txt", sep=''))
    for ( fileid in fileList )
    {
        FileName <- fileid
        stat <- read.csv(file=FileName, header=FALSE)
        L <- length(stat[,1])
        if ( L > 200 )
        {
            A <- rbind(A, c(force[rep], stat[L,10]))
        }
        cat(L, " ")
    }
    cat("\n")
}
plot(A[,1], A[,2])

## colList <- c("red", "blue", "black", "green", "cyan")
## colList <- rep(colList, 5)
## png(paste("output_", output_file_prefix, ".png", sep=''), height=600, width=1200)
## par(mar=c(5,6,4,4))
## plot(0, 0, type='n', xlim=c(0, 2000),
##      ##ylim=c(-10, 10),
##      ylim=c(-1, 0),
##      xlab="number of divisions", ylab="deviation X(t)",
##      cex.axis=1.5, cex.lab=1.5)
## abline(h=c(0, -5, 5),col='blue')
## dev.off()


## hist(atan(abs(sin(totalDividingAngleList)/cos(totalDividingAngleList))))

## png('AspectRatio.png', height=600, width=600)
## ##par(mfrow=c(1,2))
## ##angle <- angle + rnorm(length(angle), 0, pi*15.37/180)
## angle=atan(abs(sin(angle)/cos(angle)))*180/pi
## t = hist(angle, breaks=c(0, 30, 60, 90), plot=F)
## barplot(t$counts/sum(t$counts), ylim=c(0, 1))
## ##t = hist(angle[braf==1], breaks=c(0, 30, 60, 90))
## #tt = hist(angle[braf==0], breaks=c(0, 30, 60, 90))
## #d1 = t$density/sum(t$density)
## #d0 = tt$density/sum(tt$density)
## dev.off()

## write.table(t$counts/sum(t$counts), file="as.csv", sep=',')

## cat('Fcal: ', var(Fcal, na.rm=T), '\n')
## cat('Gcal: ', var(Gcal, na.rm=T), '\n')

## pdf("p1.pdf")
## ##png("p1.png")
## par(mar=c(5,6,4,4))
## plot(0, 0, type='n', xlim=c(0, 2000), ylim=c(0, 1),
##      xlab="time", ylab="percentage of elongated cells",
##      cex.axis=1.5, cex.lab=2)
## abline(h=c(0.37), col='black', lwd=4, lty=2)
## rep = rep0
## n0 = 600
## n <- 0
## for ( i in seq(5) )
## {
##     rep = rep + 1
##     cat("\n ***********run********* ", rep, "\n")
##     FolderName <- paste("case",as.character(rep), "/out/", sep='')
##     cmd = paste("ls ", FolderName, "stat* > ", FolderName, "filelist.txt", sep='')
##     system(cmd)
##     fileList <- readLines(paste(FolderName, "filelist.txt", sep=''))
##     for ( fileid in fileList )
##     {
##         n <- n + 1
##         FileName <- fileid
##         A <- read.csv(file=FileName, header=FALSE)
##         if ( length(A[,10])/10 >= Nsample )
##         {
##             ##lines(A[,10], col=colList[i])
##             lines(A[,10], col='red')
##         }
##     }
## }

## rep = 20
## n0 = 600
## n <- 0
## for ( i in seq(10) )
## {
##     rep = rep + 1
##     cat("\n ***********run********* ", rep, "\n")
##     FolderName <- paste("../360_InitialStates_CellModel_MitotiSpindle/case",as.character(rep), "/out/", sep='')
##     cmd = paste("ls ", FolderName, "stat* > ", FolderName, "filelist.txt", sep='')
##     system(cmd)
##     fileList <- readLines(paste(FolderName, "filelist.txt", sep=''))
##     for ( fileid in fileList )
##     {
##         n <- n + 1
##         FileName <- fileid
##         A <- read.csv(file=FileName, header=FALSE)
##         if ( length(A[,10])/10 >= Nsample )
##         {
##             ##lines(A[,10], col=colList[i])
##             lines(A[,10], col='blue')
##         }
##     }
## }
## dev.off()

## plot(0, 0, type='n', xlim=c(0, 1000), ylim=c(0, 1),
##      xlab="number of divisions", ylab="p",
##      cex.axis=1.5, cex.lab=1.5)
## abline(h=c(0, 0.37),col='blue')

## rep = rep0
## n0 = 1
## Fcal <- NULL
## n <- 0
## for ( i in seq(10) )
## ##for ( i in 10 )
## {
##     rep = rep + 1
##     cat("\n ***********run********* ", rep, "\n")
##     FolderName <- paste("case",as.character(rep), "/out/", sep='')
##     cmd = paste("ls ", FolderName, "stat* > ", FolderName, "filelist.txt", sep='')
##     system(cmd)
##     fileList <- readLines(paste(FolderName, "filelist.txt", sep=''))
##     for ( fileid in fileList )
##     {
##         n <- n + 1
##         FileName <- fileid
##         stat <- read.csv(file=FileName, header=FALSE)
##         z <- 100*(stat[n0:dim(stat)[1],9] - 0.37)
##         ##z <- z - 0.1*seq(length(z))
##         ##lines(stat[,8], type='o', col=colList[i])
##         lines(z, col='green')
##     }
## }


## rep = rep0
## plot(0, 0, type='n', xlim=c(0, 200), ylim=c(-15, 15))
## abline(h=c(0, -5, 5),col='blue')
## n0 = 1
## Fcal <- NULL
## for ( i in seq(40) ) {
##   rep = rep + 1
##   cat("\n ***********run********* ", rep, "\n")
##   FolderName <- gsub("(\ )", "", paste("case",as.character(rep)))
##   FileName <- paste(FolderName, "/out/statistics.txt", sep='')
##   stat <- read.csv(file=FileName, header=FALSE)
##   z <- stat[n0:dim(stat)[1],7] - stat[n0,7]
##   z <- z - 0.0*seq(length(z))
##   ##lines(stat[,8], type='o', col=colList[i])
##   lines(z, col=colList[i])
##   if ( length(z) >= Nsample )
##   {
##       Fcal[i] <- z[Nsample]
##   }
##   else
##     {
##         Fcal[i] <- NA
##     }
##   ##lines(15*(z[,8]-0.37), col=colList[i])
## }

##legend("bottomleft", c("random", "model"), lwd=1, col=c('green', 'blue'), cex=1.5)

##dev.off()

## plot(0, 0, type='n', xlim=c(0, 200), ylim=c(0, 1))
## abline(h=0.37,col='blue')
## rep = rep0
## n0 = 1
## for ( i in seq(10) ) {
##   rep = rep + 1
##   cat("\n ***********run********* ", rep, "\n")
##   FolderName <- gsub("(\ )", "", paste("case",as.character(rep)))
##   FileName <- paste(FolderName, "/out/statistics.txt", sep='')
##   stat <- read.csv(file=FileName, header=FALSE)
##   ## z <- stat[n0:dim(stat)[1],9] - stat[n0,7]
##   ## z <- z + 0.02*seq(length(z))
##   lines(stat[,9], col=colList[i])
##   ##lines(15*(z[,8]-0.37), col=colList[i])
## }

