#=======================================================
#' @title Simulated genotypes
#'
#' @description Simulate SNP genotype data from a pedigree, with optional
#'   missingess and errors.
#'
#' @details There are lots of details.
#'
#' @param Pedigree  Dataframe, pedigree with the first three columns being id -
#'   dam - sire. Column names are ignored, as are additional columns, with the
#'   exception of a 'Sex' column when Inherit is not 'autosomal'.
#' @param nSnp  number of SNPs to simulate.
#' @param ParMis  Single number or vector length two with proportion of parents
#'   with fully missing genotype. Ignored if CallRate is a named vector.
#' @param MAF  minimum minor allele frequency, and allele frequencies will be
#'   sampled uniformly between this minimum and 0.5, OR a vector with minor
#'   allele frequency at each locus. In both cases, this is the MAF among
#'   pedigree founders, the MAF in the sample will deviate due to drift.
#' @param CallRate Either a single number for the mean call rate (genotyping
#'   success), OR a vector with the call rate at each SNP, OR a named vector
#'   with the call rate for each individual. In the third case, ParMis is
#'   ignored, and individuals in the pedigree (as id or parent) not included in
#'   this vector are presumed non-genotyped.
#' @param SnpError  mean per-locus genotyping error rate across SNPs, and a
#'   beta-distribution will be used to simulate the number of missing cases per
#'   SNP, OR a vector with the genotyping error for each SNP.
#' @param ErrorFM  function taking the error rate (scalar) as argument and
#'   returning a 3x3 matrix with probabilities that actual genotype i (rows) is
#'   observed as genotype j (columns). Inbuilt ones are as used in sequoia
#'   'version2.0', 'version1.3', or 'version1.1'. See details.
#' @param ReturnStats in addition to the genotype matrix, return the input
#'   parameters and mean & quantiles of MAF, error rate and call rates.
#' @param OutFile  file name for simulated genotypes. If NA (default), return
#'   results within R.
#' @param Inherit  inheritance pattern, scalar or vector of length nSnp,
#'   Defaults to 'autosomal'. An excel file included in the package has
#'   inheritance patterns for the X and Y chromosome and mtDNA, and allows
#'   custom inheritance patterns. Note that these are NOT currently supported by
#'   the pedigree reconstruction
#' @param InheritFile  file name for excel file with inheritance patterns,
#'   requires library xlsx.
#' @param quiet suppress messages.
#' @param PropLQ  [deprecated] proportion of low-quality samples.
#' @param MisHQ  [deprecated] average missingness for high-quality samples,
#'   assuming a beta-distribution with alpha = 1.
#' @param MisLQ  [deprecated] average missingness in low-quality samples.
#' @param ErHQ  [deprecated] error rate in high quality samples (defaults to
#'   0.005).
#' @param ErLQ  [deprecated] error rate in low quality samples.
#'
#' @return if \code{ReturnStats=FALSE} (the default), a matrix with genotype data in
#'   sequoia's input format, encoded as 0/1/2/-9.
#'
#'   If \code{ReturnStats=TRUE}, a named list with three elements: list
#'   'ParamsIN', matrix 'SGeno', and list 'StatsOUT':
#'   \item{AF}{Frequency in 'observed' genotypes of '1' allele}
#'   \item{AF.act}{Allele frequency in 'actual' (without genotyping errors &
#'     missingness)}
#'   \item{SnpError}{Error rate per SNP (actual /= observed AND observed /=
#'     missing)}
#'   \item{SnpCallRate}{Non-missing per SNP}
#'   \item{IndivError}{Error rate per individual}
#'   \item{IndivCallRate}{Non-missing per individual}
#'
#' @seealso other stuff.
#'
#' @author Jisca Huisman, \email{jisca.huisman@gmail.com}
#'
#' @section Disclaimer: This simulation is highly simplistic and assumes that
#'   all SNPs segregate completely independently, that the SNPs are in
#'   Hardy-Weinberg equilibrium in the pedigree founders. It assumes that
#'   genotyping errors are not due to heritable mutations of the SNPs, and that
#'   missingness is random and not e.g. due to heritable mutations of SNP
#'   flanking regions. Results based on this simulated data will provide an
#'   minimum estimate of the number of SNPs required, and an optimistic estimate
#'   of pedigree reconstruction performance.
#'
#' @importFrom stats rbinom runif rbeta setNames
#' @importFrom plyr laply
#'
#' @export

SimGeno <- function(Pedigree,
                    nSnp = 400,
                    ParMis = 0.4,
                    MAF = 0.3,
                    CallRate = 0.99,
                    SnpError = 5e-4,
                    ErrorFM = "version2.0",
					          ReturnStats = FALSE,
					          OutFile = NA,
					          Inherit = "autosomal",
					          InheritFile = NA,
					          quiet = FALSE,
                    PropLQ, # 0,    # backwards compatability
                    MisHQ, # 0.005,
                    MisLQ, #0.30,
                    ErHQ, #5e-4,
                    ErLQ) #5e-3,
{
  if (is.null(OutFile)) stop("'OutFile' must be filename or NA")

  if (missing(Pedigree)) stop("please provide a pedigree to simulate from")
  if (!class(Pedigree) %in% c("data.frame", "matrix"))  stop("Pedigree should be a dataframe with (at least) 3 columns")
  if (ncol(Pedigree) < 3)  stop("Pedigree should be a dataframe with at least 3 columns (id - dam - sire)")

  if ( !missing(PropLQ)) if (PropLQ!=0) stop("PropLQ, MisLQ, and ErLQ are deprecated, use CallRate")
  if (!missing(MisHQ)) {
     CallRate <- 1 - MisHQ
     warning("NOTE: 'MisHQ' will be deprecated, please use 'CallRate'", immediate.=TRUE)
  }
  if (!missing(ErHQ)) {
    SnpError <- ErHQ
    warning("NOTE: 'ErHQ' will be deprecated, please use 'SnpError'", immediate.=TRUE)
  }

  # unavoidable partial matching when 'Err' is specified instead of 'SnpError'
  if (is.numeric(ErrorFM)) {
    SnpError <- ErrorFM
    ErrorFM <- "version2.0"   # problem doesn't occur when 'ErrorFM' specified; restore default
  }

  #================================
  # check input

  if(!(is.numeric(nSnp)) || nSnp<=0)  stop("nSnp should be a number greater than 0")

  if(length(ParMis)==1) ParMis <- rep(ParMis, 2)

  params <- list(ParMis1 = ParMis[1], ParMis2 = ParMis[2], MAF = MAF,
                 SnpError = SnpError, CallRate = CallRate)
  for (p in seq_along(params)) {
    if (p==1 & length(params[[p]]) != 1) {
      stop("Length ", names(params)[p], " should be 1")
    } else if (p<4 & !length(params[[p]]) %in% c(1, nSnp)) {
      stop("Length ", names(params)[p], " should be 1 or nSnp")
    } else if (is.null(params[[p]]) | all(is.na(params[[p]]))) {
      stop("Please provide ", names(params)[p])
    } else if ( ! (is.numeric(params[[p]]) & all(params[[p]]>=0) & all(params[[p]]<=1)) ) {
      stop(names(params)[p], " must be a number between 0 and 1")
    }
  }

  if (length(CallRate) > 1) {
    if (is.null(names(CallRate))) {
      if (length(CallRate) != nSnp)  stop("CallRate should be length 1 or nSnp, or a named vector")
    } else {
      if (length(intersect(names(CallRate), Pedigree[,1]))==0) stop("names of CallRate vector do not match pedigree")
    }
  }

  # Error matrix: rows = actual, columns = observed
	ErFunc <- ErrToM(Err=0.1, flavour=ErrorFM, Return="function")  # 0.1 is test value only

  if(!is.na(OutFile) & ReturnStats)  stop("Cannot write Return Stats to OutFile")

  if (interactive() & !quiet & !is.na(OutFile)) {
    if (file.exists(OutFile)) {
      ANS <- readline(prompt = paste("WARNING: ", OutFile,
                                     "will be overwritten.",
                                     "Press <N> to abort, or any other key to continue."))
    } else {
      ANS <- readline(prompt = paste("Genotypes will be written to ", OutFile,
                                     ". Press <N> to abort, or any other key to continue."))
    }
    if (substr(ANS, 1, 1) %in% c("N", "n")) stop()
  }

  ParamsIN <- as.list(environment())

  #================================
  # minor allele frequencies (among founders)
  if (length(MAF)==1) {
    Q <- round(runif(nSnp, min=MAF, max=0.5),3)
  } else {
    Q <- as.numeric(MAF)
  }

  #================================
  # check & prep
#  Ped <- PedPolish(Pedigree, ZeroToNA=TRUE)
  Ped <- Pedigree
  for (x in 1:3) Ped[,x] <- as.character(Ped[,x])
  nInd <- nrow(Ped)
  if (any(round(Q*nInd) %in% c(0,1)))  warning("some simulated SNPs have fixed alleles")


  #================================
  # simulate genotypes
  # founders: random draw of alleles under HWE
  # non-founders: following Mendelian inheritance

  # rownumber of dam & sire
  Ped$damIDx <- sapply(Ped[,2], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))
  Ped$sireIDx <- sapply(Ped[,3], function(x) ifelse(is.na(x), 0, which(Ped[,1]==x)))

  #~~~~~~~~~~~~~~
  # divide pedigree into `generations': the parents of an individual must come
  # from earlier cohorts than itself, or from the founder population (gen 0)
  Ped$gen <- getGenerations(Ped)[, "gen"]
  nGen <- max(Ped$gen)


  #~~~~~~~~~~~~~~
  if (all(Inherit == "autosomal")) {

    # founders
    SGeno <- matrix(NA, nInd, nSnp,
                    dimnames=list(Ped[,1], NULL))
    for (i in which(Ped$gen==0)) {
      SGeno[i, ] <- rbinom(nSnp, 2, prob=Q)
    }

    # non-founders
    getHaplo <- function(j, G, Q) {
      if (j >0) {
        rbinom(ncol(G), 1, prob=G[j,]/2)
      } else {
        rbinom(length(Q), 1, prob=Q)
      }
    }

    for (g in 1:nGen) {
      for (i in which(Ped$gen==g)) {
        SGeno[i, ] <- rowSums( cbind( getHaplo(Ped$damIDx[i], SGeno, Q),
                                      getHaplo(Ped$sireIDx[i], SGeno, Q) ))
      }
    }

  #~~~~~~~~~~~~~~
  } else {

    if (!"Sex" %in% names(Ped))  stop("Inherit other than autosomal requires 'Sex' column in Ped")
    Ped$Sex[is.na(Ped$Sex)] <- 3
    if (any(Ped$Sex == 4)) stop("GenoSim not yet implemented for hermaphrodites (Sex=4)")
    Ped$Sex[!Ped$Sex %in% 1:2] <- 3

    # inheritance patterns
    if (is.na(InheritFile) | is.null(InheritFile)) {
      utils::data(Inherit)
      INHA <- Inherit
      rm(Inherit)
#    } else {
#      INHA <- ReadSpecialInherit(InheritFile)  # array: inherit - off sex - geno off - dam - sire
    }

    # founders
    founderProp <- apply(INHA, 1:3, sum)
    SGeno.4 <- matrix(NA, nInd, nSnp, dimnames=list(Ped[,1], NULL))
    for (i in which(Ped$gen==0)) {
      SGeno.4[i, ] <- sapply(1:nSnp, function(l) sample.int(4, size=1,
                                prob=c(Q[l]^2, Q[l]*(1-Q[l]), Q[l]*(1-Q[l]), (1-Q[l])^2) *
                founderProp[Inherit[l], Ped$Sex[i], ]) )   # CHECK
    }

    # non-founders
    Gprob <- matrix(NA, nSnp, 4,
                    dimnames=list(1:nSnp, c("aa", "aA", "Aa", "AA")))
    for (g in 1:nGen) {
      for (i in which(Ped$gen==g)) {
        if (Ped$damIDx[i] != 0 & Ped$sireIDx[i] != 0) {
          for (x in 1:4) {  #
            Gprob[,x] <- INHA[,Ped$Sex[i],x,,][cbind(Inherit, SGeno[Ped$damIDx[i],], SGeno[Ped$sireIDx[i],]) ]
          }
        } else {
          stop("Non-autosomal single parents not implemented yet!!")
        }
        SGeno.4[i,] <- apply(Gprob, 1, function(p) sample.int(4, size=1, prob=p))
      }
    }

    SGeno <- apply(SGeno.4, 1, function(v) c(0,1,1,2)[v])   # CHECK
  }


  #================================
  # genotyping errors & missing values:

  SGeno.actual <- SGeno
  SGeno <- MkGenoErrors(SGeno, CallRate, SnpError, ErFunc)

  #================================
  # Non-genotyped individuals

  NotSampled <- which(apply(SGeno, 1, function(x) all(x==-9)))

  if (any(ParMis>0) & is.null(names((CallRate)))==1) {
    for (p in 1:2) {
      if (ParMis[p]>0) {
        IsParent <- which(Ped[,1] %in% Ped[,p+1])
      }
      if (round(length(IsParent)*ParMis[p]) > 0) {
        NotSampled <- c(NotSampled,
                        sample(IsParent, round(length(IsParent)*ParMis[p]),
                               replace=FALSE) )
      }
    }
  }
  if (length(NotSampled)>0) {
    SGeno <- SGeno[-NotSampled, ]
    SGeno.actual <- SGeno.actual[-NotSampled, ]
    nInd.g <- nInd -length(NotSampled)
  } else {
    nInd.g <- nInd
  }

  #================================
  # output

  if (!is.na(OutFile)) {
    utils::write.table(SGeno, OutFile, quote=FALSE, col.names=FALSE)

  } else if (ReturnStats) {

    Params <- c(ParamsIN[c("nSnp", "ParMis", "MAF", "CallRate", "SnpError",
                           "Inherit", "ErrorFM", "OutFile", "InheritFile")],
                list(nInd = nrow(SGeno)))

    StatsOUT <- list(AF = apply(SGeno, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9))),
                     AF.act = apply(SGeno.actual, 2, function(x) sum(x[x!=-9])/(2*sum(x!=-9))),
                     SnpError = sapply(1:nSnp, function(l) {
                       sum(SGeno[,l] != SGeno.actual[,l] & SGeno[,l]!=-9)/nInd.g }),
                     SnpCallRate = apply(SGeno, 2, function(x) sum(x!=-9))/nInd.g,
                     IndivError = sapply(1:nInd.g, function(i) {
                       sum(SGeno[i,] != SGeno.actual[i,] & SGeno[i,]!=-9)/nSnp }),
                     IndivCallRate = apply(SGeno, 1, function(x) sum(x!=-9))/nSnp )

    return( list(ParamsIN = ParamsIN, StatsOUT = StatsOUT, SGeno = SGeno) )  # Ped = Ped,

  } else {
    return( SGeno )
  }
}


#=============================================================================
#' @title Simulate genotyping errors
#'
#' @description Generate errors and missing values in a (simulated) genotype
#'   matrix
#'
#' @param SGeno  Matrix with genotype data in Sequoia's format: 1 row per
#'   individual, 1 column per SNP, and genotypes coded as 0/1/2.
#' @param CallRate Either a single number for the mean call rate (genotyping
#'   success), OR a vector with the call rate at each SNP, OR a named vector
#'   with the call rate for each individual. In the third case, ParMis is
#'   ignored, and individuals in the pedigree (as id or parent) not included in
#'   this vector are presumed non-genotyped.
#' @param SnpError  mean per-locus genotyping error rate across SNPs, and a
#'   beta-distribution will be used to simulate the number of missing cases per
#'   SNP, OR a vector with the genotyping error for each SNP.
#' @param ErrorFM  function taking the error rate (scalar) as argument and
#'   returning a 4x4 or 3x3 matrix with probabilities that actual genotype i
#'   (rows) is observed as genotype j (columns).
#' @param Error.shape first shape parameter (alpha) of beta-distribution of
#'   per-SNP error rates. A higher value results in a flatter distribution.
#' @param CallRate.shape as Error.shape, for per-SNP call rates.
#'
#' @return  The input genotype matrix, with some genotypes replaced, and some
#'   set to missing (-9)
#'
#'
#' @export

MkGenoErrors <- function(SGeno,
                         CallRate = 0.99,
                         SnpError = 5e-4,
                         ErrorFM = function(E) {
                           matrix(c(1-E-(E/2)^2, E, (E/2)^2,
                                    E/2, 1-E, E/2,
                                    (E/2)^2, E, 1-E-(E/2)^2),
                                  3,3, byrow=TRUE) },
                         Error.shape=0.5,
                         CallRate.shape=1)
{
  nSnp <- ncol(SGeno)
  nInd <- nrow(SGeno)

#  if (! all(SGeno %in% c(0,1,2)) ) stop("SGeno may only contain 0, 1, or 2")

  #~~~~~~~~~
  if (any(SnpError >0)) {
    if (length(SnpError)==1) {
      El <- rbeta(nSnp, shape1=Error.shape, shape2=Error.shape*(1/SnpError -1))
    } else if (length(SnpError) == ncol(SGeno)) {
      El <- SnpError
    } else {
      stop("length of SnpError should equal 1 or number of SNPs")
    }

    shrinkET <- function(M) {
      N <- matrix(NA, 3,3)
      N[c(1,3), c(1,3)] <- M[c(1,4), c(1,4)]
      N[2, c(1,3)] <- M[2, c(1,4)]+M[3, c(1,4)]
      N[c(1,3), 2] <- M[c(1,4), 2] + M[c(1,4), 2]
      N[2, 2] <- sum(M[2:3, 2:3])
      return( N )
    }

    if (all(dim(ErrorFM(0.1))==3)) {
      RealToObs <- laply(El, ErrorFM)
    } else if (all(dim(ErrorFM(0.1))==4)) {
      RealToObs <- laply(El, function(e) shrinkET(ErrorFM(e)))
    } else {
      stop("ErrorFM(E) should return a 4x4 or 3x3 matrix")
    }

#    for (l in 1:nSnp) {
#      SGeno[,l] <- sapply(SGeno[,l], function(x) sample.int(3, 1, prob=RealToObs[l,x+1,]) -1 )
#    } # rather slow; implemented in Fortran instead:
    SGeno <- DoErrors(SGeno, RealToObs)
  }

  #~~~~~~~~~
  if (any(CallRate <1)) {
    CRtype <- ifelse(length(CallRate)==1, "mean",
                ifelse(!is.null(names(CallRate)), "Indiv", "SNP"))
    MisX <- matrix(FALSE, nInd, nSnp)

    if (CRtype == "Indiv") {
      IndivCallRate <- setNames(CallRate[rownames(SGeno)], rownames(SGeno))
      IndivCallRate[is.na(IndivCallRate)] <- 0
      lmis <- round((1-IndivCallRate) *nSnp) # no. missing SNPs per indiv
      for (i in 1:nInd) {
        MisX[i, sample.int(nSnp, lmis[i])] <- TRUE   # ERROR
      }
    } else {
      if (CRtype == "mean") {
        imis <- round(rbeta(nSnp, CallRate.shape, CallRate.shape*(1/(1-CallRate) -1)) *nInd)
        # no. missing indiv per SNP
      } else if (CRtype == "SNP") {
        imis <- round((1-CallRate) *nInd)
      }
      for (l in 1:nSnp) {
        MisX[sample.int(nInd, imis[l]), l] <- TRUE
      }
    }
    SGeno[MisX] <- -9
  }

  #~~~~~~~~~
  return( SGeno )
}


#=============================================================================
#' @title Simulate genotyping errors
#'
#' @description Wrapper for Fortran function to simulate genotyping errors.
#'
#' @param SGeno matrix with genotype data, size nInd x nSnp
#' @param RealToObs array with conditional probability of observing genotype i
#'   conditional on actual genotype j, size nSnp x 3 x 3
#'
#' @return SGeno with errors
#'
#' @useDynLib SimGE, .registration = TRUE
#'
#' @keywords internal

DoErrors <- function(SGeno, RealToObs) {
   dnames <- dimnames(SGeno)
   TMP <- .Fortran(mkerrors,
                   nind = as.integer(nrow(SGeno)),
                   nsnp = as.integer(ncol(SGeno)),
                   genofr = as.integer(SGeno),
                   eprobfr = as.double(RealToObs))
   return( matrix(TMP$genofr, nrow(SGeno), ncol(SGeno), dimnames=dnames) )
}


#=============================================================================

#' @title Count generations
#'
#' @description For each individual in a pedigree, count the number of
#'   generations since the most distant pedigree founder
#'
#' @param Ped  Dataframe, pedigree with the first three columns being ID - dam -
#'   sire. Column names are ignored, as are additional columns.
#'
#' @return  The input pedigree with three columns added: gen, gen.dam, and
#'   gen.sire
#'
#' @keywords internal

getGenerations <- function(Ped) {
  for (p in 2:3) {
    Ped[which(Ped[,p]==0), p] <- NA
  }

  Ped$gen <- NA  # individual's generation
  Ped$gen.dam <- NA
  Ped$gen.sire <- NA
  Ped$gen[is.na(Ped[,2]) & is.na(Ped[,3])] <- 0
  for (x in 0:1000) {
    Ped$gen.dam[is.na(Ped$gen.dam) & Ped[,2] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen.sire[is.na(Ped$gen.sire) & Ped[,3] %in% Ped[which(Ped$gen<=x), 1]] <- x
    Ped$gen[which(is.na(Ped$gen) &
                   (Ped$gen.dam<=x | is.na(Ped[,2])) &
                   (Ped$gen.sire<=x | is.na(Ped[,3])))] <- x+1
    if (!any(is.na(Ped$gen)))  break
  }
  return( Ped[, c("gen", "gen.dam", "gen.sire")] )
}


#=============================================================================
