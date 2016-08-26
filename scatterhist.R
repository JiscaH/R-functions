scatterhist = function(x, y, br=25, 
                       xlab, ylab, LINE=NULL,
                       colX=colours()[81], colY=colours()[26],
                       layW = c(2/3, 1/3), layH = c(1/3, 2/3),
                       marB=0.6, marL=0.7,  # margin bottom, left
                       ...){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=layW, heights=layH)
  xhist = hist(x, plot=FALSE, breaks=br)
  yhist = hist(y, plot=FALSE, breaks=br)
  top = c(x=max(xhist$counts), y=max(yhist$counts))
  par(mai=c(marB,marL,0,0))
  plot(x,y, xlim=range(xhist$breaks), ylim=range(yhist$breaks), xlab="", ylab="",
       cex.axis=1.5, las=1, ...)
  #  abline(h=0,v=0, col=2)
  if (!is.null(LINE)) abline(LINE, lty=2,lwd=2)
  mtext(xlab, side=1, line=3, outer=FALSE, adj=.5,cex=1.5) 
  mtext(ylab, side=2, line=3, outer=FALSE, adj=.5,cex=1.5)
  par(mai=c(0,marL,.1,0))
  bp.width <- xhist$breaks[2]-xhist$breaks[1]
  barplot(xhist$counts, axes=FALSE, space=0, width=bp.width, ylim=c(0, top["x"]), col=colX)
  par(mai=c(marB,0,0,.1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top["y"]), space=0, width=bp.width, 
          horiz=TRUE, col=colY)
}