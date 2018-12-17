require(lattice)
OP <- function(angles)
{
    OP <- mean(cos(2*angles))
    return(OP)
}

output_file_prefix <- "testrun_2model"

rep0 = 0

#rep0 = 110
p0 <- 0.75
pbraf <- 0.0
Nsample <- 200

colList <- c("red", "blue", "black", "green", "cyan")
colList <- rep(colList, 5)
##par(mfrow=c(2,1))

png(paste("output_", output_file_prefix, ".png", sep=''), height=600, width=1200)
par(mar=c(5,6,4,4))
##par(mfrow=c(2,2))

plot(0, 0, type='n', xlim=c(0, 2000),
     ##ylim=c(-10, 10),
     ylim=c(-1, 0),
     xlab="number of divisions", ylab="deviation X(t)",
     cex.axis=1.5, cex.lab=1.5)
abline(h=c(0, -5, 5),col='blue')

angle = NULL
braf = NULL
A = NULL
rep = rep0
n0 = 600
Fcal <- NULL
n <- 0
totalDividingAngleList <- NULL
OrderParameterList <- NULL
for ( i in seq(16) )
##for ( i in 10 )
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
        z <- stat[n0:dim(stat)[1],9] - stat[n0,9]
        ##z <- z - 0.1*seq(length(z))
	##lines(stat[,8], type='o', col=colList[i])
	##lines(stat[,11], col=colList[floor((i-1)/1)+1])
        ##lines((1-p0)*((stat[,5]-stat[,4])-40) - p0*2.0*(stat[,6]-17), col=colList[floor((i-1)/1)+1])
        lines(stat[,11], col=colList[floor((i-1)/1)+1])
        ##lines(stat[,10], col=colList[floor((i-1)/5) + 1])
        ##lines(cumsum(stat[stat[,3]>=0,9])-0.15*seq_along(stat[stat[,3]>=0,9]), col=colList[floor((i-1)/5) + 1])
        if ( length(stat[,11]) >= 450 )
        {
            stableStat <- stat[-seq(600),]
            dividingAngleList <- stableStat[stableStat[,3]>0,12]
            totalDividingAngleList <- c(totalDividingAngleList, dividingAngleList)
            OrderParameterList <- rbind(OrderParameterList, c(i, stat[450,11]))
            cat('order parameter: ', OP(dividingAngleList), '\n')
        }
        if ( length(stat[,11]) == 600 )
        {
            n <- n + 1
            ##A[rep] <- stat[length(stat[,11]),11]
            A <- rbind(A, c(rep, stat[560,11]))
            Fcal[n] <- z[10*Nsample]
            ##lines(z[1:(10*Nsample)], col="blue")
            angle = c(angle, stat[10*seq(n0/10, length.out=Nsample),11])
            ##braf = c(braf, stat[10*seq(20, Nsample),11])
        }
    }
}
dev.off()

#hist(atan(abs(sin(totalDividingAngleList)/cos(totalDividingAngleList))))

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

