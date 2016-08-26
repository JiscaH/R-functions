# create a stacked histogram out of a variable V with K levels given in vector L

StackHist <- function(V,L=NULL, br=seq(.99,1.15,by=.005)-.0025, beside=F,...) 
{
  if (is.vector(V)) {
    L <- as.factor(L)
    HIST <- matrix(,length(levels(L)),length(br)-1)
    for (k in 1:length(levels(L))) {
      HIST[k,] <- hist(as.numeric(V[L==levels(L)[k]]), breaks=br, plot=F)$counts
    }
  } else if (is.matrix(V)) {
    if (is.null(L)) L <- 1:ncol(V)  # order of columns
    HIST <- matrix(, ncol(V), length(br)-1)
    for (k in 1:ncol(V)) {
      Vtmp <- as.numeric(V[, L[k]])
      Vtmp[Vtmp <= min(br)] <- NA
      Vtmp[Vtmp >= max(br)] <- NA
      HIST[k,] <- hist(Vtmp, breaks=br, plot=F)$counts
    }
  }

  if (beside==FALSE) {
    # get pretty labels for x-axis
    hist(Vtmp, breaks=br, xlab="", ylab="", yaxt="n", main="", border=0)
    ticks <- pretty(axTicks(1))
    par(new=TRUE)
    barplot(HIST, beside=F, space=0, las=1, width=(br[2]-br[1]), ...)
    axis(1, at=ticks-min(br), labels=ticks)
  } else {
    b <- barplot(HIST, beside=TRUE)#, ...) 
    lbls.all <- (br[-length(br)]+br[-1])/2
    lbls.pretty <- pretty((br[-length(br)]+br[-1])/2)
    axis(1, at=colMeans(b)[lbls.all %in% lbls.pretty], labels=lbls.pretty)
  }

}