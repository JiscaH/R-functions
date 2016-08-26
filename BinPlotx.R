Table <- function(...) table(..., useNA="ifany")
l <- function(x=1, a=.5) adjustcolor(x, alpha.f=a)
LB <- adjustcolor(1, alpha.f=0.5)

# for x continuous variable:

BinPlot <- function(x,y,
                    Pr=seq(0,1,by=.1), 
                    Bins=NULL,  # give either quantile prob. or bins
                    printN=FALSE,   # return numbers underlying each dot
                    col=1,   		# colour of outline 
                    PCH=21,
          					bg=NULL,
          					Offs=0,  		# horizontal offset (to decrease overlap)
                    hln = NA,       # horizontal line
                    xloc = "avg",   # avg = average, eq = equally spaced 
                    atX = NULL,
                    Nsize = TRUE,  # scale point size by no. datapoints underlying it
                    Scale = NULL,   # scaling factor
                    DifSmall = NULL,   # use diff. symbols for else too-small points (n<=DifSmall)
					          add = FALSE,    # add to existing plot?
          					ErB = TRUE,     # add error bars?
                    ...)            # ohter graphical parameters
  {
  x2 <- x[!is.na(x) & !is.na(y)]
  y2 <- y[!is.na(x) & !is.na(y)]
  x <- x2
  y <- y2
  if(is.null(Bins)) {
    Bins <- quantile(x,probs=Pr,na.rm=T)
    if(any(duplicated(Bins))) {
      cat("warning: breaks are not unique, collapsing. \n")
      Bins <- unique(Bins)
    }
  }
  xb <- cut(x, breaks=Bins, include.lowest=T,right=F)
  if (xloc=="avg") {    
    xb.m <- tapply(x, xb, mean, na.rm=T)
    xb.se <- tapply(x, xb, sd, na.rm=T)/sqrt(table(xb))
  } else if (xloc=="eq") {
    if (is.null(atX)) {
      xb.m <- 1:(length(Bins)-1)
    } else {
      xb.m <- atX
    }
    xb.se <- rep(.001, length(xb.m))  # else arrows gives warnings
  }
  
  yb.m <- tapply(y, xb, mean, na.rm=T)
  yb.se <- tapply(y, xb, sd, na.rm=T)/sqrt(table(xb))
  
  if (!add) plot(xb.m, yb.m, pch=21,las=1, cex.axis=1.3,cex.lab=1.5, type="n",  ...)
  if (!is.na(hln)) abline(h=hln, col=l())
  if (is.null(bg)) bg=l(col)
  
  if(ErB) {
    arrows(xb.m+Offs, y0=yb.m-yb.se, y1=yb.m+yb.se, angle=90, code=3, length=0, col=l(col))  
    arrows(x0=xb.m-xb.se, x1=xb.m+xb.se, y0=yb.m, angle=90, code=3, length=0, col=l(col)) 
  }
  if (!Nsize) {
    points(xb.m+Offs, yb.m, pch=PCH, col=col, bg=bg) #, ...)
  } else {
    n <- sqrt(table(xb))  # scale surface proportionally to N (rather than diameter)
    if(is.null(Scale)) Scale <- max(n)
    if(is.null(DifSmall)) {  # | max(n)/min(n[n>0]) < 10 |
   #   PCH <- 21     
     } else {
       PCH <- ifelse(n>DifSmall, 21, 24)
#       CEX <- ifelse(n^3>DifSmall, n/Scale, n/Scale*2)
     }
    CEX <- n/Scale
    points(xb.m+Offs, yb.m, pch=PCH, cex=CEX, col=col, bg=bg) #, ...)
  }
  
  if(printN) print(Table(xb))
}
  
# for backwards compatibility:
BinPoints <- function(...) BinPlot(..., add=TRUE)


################################################
# for x categorical variable:

MeanPlot <- function(x,y, ...) {
  X <- x[!is.na(x) & !is.na(y)]
  Y <- y[!is.na(x) & !is.na(y)]
  x <- as.factor(X)
  y <- as.numeric(Y)
  y.m <- tapply(y, x, mean, na.rm=T)
  y.se <- tapply(y, x, sd, na.rm=T)/sqrt(table(x))
  y.se[is.na(y.se)] <- 0
  xx <- as.numeric(levels(x))
#  YLIM <- c(0.98*min(y.m-y.se), 1.02*max(y.m+y.se))
  plot(xx, y.m, pch=16,las=1, cex.axis=1.3,cex.lab=1.5,...)  # ylim=YLIM
  arrows(xx, y0=y.m-y.se, y1=y.m+y.se, angle=90, code=3, length=0)  
}

MeanPoints <- function(x,y, off=0.1, ...) {
  X <- x[!is.na(x) & !is.na(y)]
  Y <- y[!is.na(x) & !is.na(y)]
  x <- as.factor(X)
  y <- as.numeric(Y)
  y.m <- tapply(y, x, mean, na.rm=T)
  y.se <- tapply(y, x, sd, na.rm=T)/sqrt(table(x))
  y.se[is.na(y.se)] <- 0
  xx <- as.numeric(levels(x))
  points(xx+off, y.m, pch=16, ...)  
  arrows(xx+off, y0=y.m-y.se, y1=y.m+y.se, angle=90, code=3, length=0, ...)  
}