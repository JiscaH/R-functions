
# functions to summarise ASREML results

##############

Table <- function(...) table(..., useNA="ifany")

PropV <- function(ASR, ndec=5) {
  VC <- summary(ASR)$varcomp
  PrV <-  setNames(VC[,2] / sum(VC[,2]), rownames(VC))
  n <- nrow(VC)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  vcov <- matrix(NA,n,n)
  vcov[cbind(i,j)] <- ASR$ai
  vcov[upper.tri(vcov)] <- t(vcov)[upper.tri(vcov)]
  vcovX <- vcov[-n,-n]   # exclude residual variance; irrelevant for proportions
  PrV.SE <- numeric()
  if (n>2) {
    for (x in 1:(n-1)) {
      PrV.SE[x] <- sqrt(1/sum(VC$gamma)^4 * 
                          (sum(VC$gamma[-x])^2 * diag(vcovX)[x] + 
                             VC$gamma[x]^2 * sum(vcovX[-x,-x]) -    # diagonal & off-diagonal 
                             2*VC$gamma[x] * sum(VC$gamma[-x]) * sum(vcovX[x,-x])))
    }
  } else {
    PrV.SE[1] <- sqrt(1/sum(VC$gamma)^4 * vcovX)  
  }
  PrV.SE[n] <- sqrt(1/sum(VC$gamma)^4 * sum(vcovX))
  round(cbind(PrV, PrV.SE), ndec)
}


Prand <- function(ASR, ndec=5) {
  RanCall <- as.list(summary(ASR)$call)$random
  RanNames <- attr(terms.formula(RanCall), "term.labels")
  RanNames <- lapply(RanNames, function(x) as.formula(paste("~. - ", x)))
  PR <- numeric(length(RanNames))
  for (i in 1:length(RanNames)) {
    ASR.tmp <- update.asreml(ASR, random = RanNames[[i]], trace=F)
    lldif <- summary(ASR)$loglik - summary(ASR.tmp)$loglik 
    PR[i] <- 1-pchisq(2*lldif,1)
  }
  signif(PR, 5) 
}


ASReml.EstEffects <- function(model){
  x <- summary(model)$varcomp
  totvar   <- sum(x[,1])
  x$Effect <- x$gamma/totvar
  x$SE <- NA
  
  object <- model
  pframe <- as.list(object$gammas) 
  
  denominator <- "1 "
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    denominator <- paste(denominator, "+ `", effname, "` ", sep = "")
  }
  
  denominator <- paste("(", denominator, ")", sep = "")
  
  for(effname in names(pframe)[1:length(pframe)-1]){
    transform <- eval(parse(text = paste("`placeholder` ~ `", effname, "`/", denominator, sep = "")))
    
    tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe) 
    X <- as.vector(attr(tvalue, "gradient")) 
    tname <- if (length(transform) == 3) {
      transform[[2]] 
    } else "" 
    V <- object$ai 
    n <- length(pframe) 
    i <- rep(1:n, 1:n) 
    if (length(V) != length(i)) 
      stop("vcov matrix incompatible with\nnumber of variance components") 
    j <- sequence(1:n) 
    k <- 1 + (i > j) 
    se <- sqrt(sum(V * X[i] * X[j] * k))
    x[effname,"SE"] <- se
  }
  
  x
}

PropSS <- function(ASR,ndec=5) {
  WA <- wald.asreml(ASR)
  SoS <- WA[,"Sum of Sq"]
  names(SoS) <- rownames(WA)
  signif(SoS/sum(SoS),ndec)
}

FixSE <- function(ASR, ndec=5) {
  round(sqrt(ASR$vcoeff$fixed * ASR$sigma2),ndec)
}

Summary.ASR <- function(ASR.name, ndec=5, Pfix=TRUE, Pran=FALSE, ToFile=TRUE) {
  
  ASR <- eval(parse(text=ASR.name))
  if(!ASR$converge) stop("did not converge!")
  
  fixef <- cbind(coef(ASR)$fixed, FixSE(ASR,ndec))
  if (nrow(fixef)>1) {
    fixef <- fixef[nrow(fixef):1,]
  }
  
  nfix <- which(names(ASR$nolev)=="(Intercept)")
  fix.nlev <- ASR$nolev[nfix:1]
  ranef <- summary(ASR)$varcomp
  
  names.f <- rownames(fixef)
  names.r <- sapply(strsplit(rownames(ranef),"!"), function(x) x[1])
  
  est.f <- round(fixef,ndec)
  est.r <- round(ranef$component,ndec)
  
#  Props.f <- PropSS(ASR, ndec)
 # Props.f <- rep(Props.f[1:length(fix.nlev)],fix.nlev)  # last=residual
  Props.r <- PropV(ASR, ndec)
  
  if (Pfix) {
    Wald.cond <- wald.asreml(ASR, ssType="conditional", denDF="numeric", trace=F) 
    Pr.f <- setNames(signif(Wald.cond$Wald$Pr, ndec), rownames(Wald.cond$Wald))
    Pr.f <- Pr.f[rep(names(fix.nlev),fix.nlev)]
  }else {
    Pr.f <- NA
  }

  se.r <- round(ranef$std.error,ndec)
  if (Pran) {
    Pr.r <- c(Prand(ASR, ndec=5),NA)
  } else {
    Pr.r <- NA
  }
  
  OUT <- list(fixed=data.frame(Est=est.f[,1],SE=est.f[,2], Pr=Pr.f, row.names=names.f),  # Prop=Props.f,
              random=data.frame(Est=est.r, StErr=se.r, Prop=Props.r[,1], Prop.SE=Props.r[,2], 
                                Pr=Pr.r, row.names=names.r),
              loglik=ASR$loglik)  
  assign(paste0(ASR.name,".summary"), OUT)
  
  if (ToFile) { 
#    save(list= ASR.name, file= paste0("ASR_", ASR.name, ".RData"))
    save(list=paste0(ASR.name,".summary"), file=paste0("ASR_", ASR.name, "_sum.RData"))
  } else {
    return(OUT)
  }
}


loadsum <- function(ASR.name) {
  load(paste0("ASR_", ASR.name, "_sum.RData"))
  eval(parse(text=paste0(ASR.name, ".summary")))
}
################################
