#' @title Generate error matrix
#'
#' @description Generate a matrix with the probabilities of observed genotypes
#'   (columns) conditional on actual genotypes (rows), or return a function to
#'   generate such matrices (using a single value Err as input to that function)
#'
#' @details By default (\code{flavour} = "SNPchip"), \code{Err} is interpreted
#'   as a locus-level error rate (rather than allele-level), and equals the
#'   probability that an actual heterozygote is observed as either homozygote
#'   (i.e., the probability that it is observed as AA = probability that
#'   observed as aa = \code{Err}/2). The probability that one homozygote is
#'   observed as the other is (\code{Err}/2\eqn{)^2}.
#'
#' The inbuilt 'flavours' correspond to the presumed and simulated error
#' structures, which have changed with sequoia versions. The most appropriate
#' error structure will depend on the genotyping platform; 'version0.9' and
#' 'version1.1' were inspired by SNP array genotyping while 'version1.3' and
#' 'version2.0' are intended to be more general.
#'
#' \emph{version2.0:}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{(1-E/2)^2} \tab \eqn{E(1-E/2)} \tab \eqn{(E/2)^2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{(E/2)^2}   \tab \eqn{E(1-E/2)} \tab \eqn{(1-E/2)^2} \cr
#' }
#'
#' \emph{version1.3}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E-(E/2)^2} \tab \eqn{E} \tab \eqn{(E/2)^2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{(E/2)^2}   \tab \eqn{E} \tab \eqn{1-E-(E/2)^2} \cr
#' }
#'
#' \emph{version1.1}
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E} \tab \eqn{E/2} \tab \eqn{E/2} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{E/2}   \tab \eqn{E/2} \tab \eqn{1-E} \cr
#' }
#'
#' \emph{version0.9} (not recommended)
#' \tabular{lccc}{
#'     \tab \strong{0} \tab \strong{1} \tab \strong{2} \cr
#'  \strong{0}  \tab \eqn{1-E} \tab \eqn{E} \tab \eqn{0} \cr
#'  \strong{1}  \tab \eqn{E/2}       \tab \eqn{1-E}      \tab \eqn{E/2}    \cr
#'  \strong{2}  \tab \eqn{0}   \tab \eqn{E} \tab \eqn{1-E} \cr
#' }
#'
#' @param Err estimated genotyping error rate, as a single number or 3x3 or 4x4
#'   matrix. If a single number, an error model is used that aims to deal with
#'   scoring errors typical for SNP arrays. If a matrix, this should be the
#'   probability of observed genotype (columns) conditional on actual genotype
#'   (rows). Each row must therefore sum to 1. If \code{Return='function'}, this
#'   may be \code{NA}.
#' @param flavour matrix-generating function, or one of 'version2.0',
#'   'version1.3' (='SNPchip'), 'version1.1' (='version111'), referring to the
#'   sequoia version in which it was used as default. Ignored if \code{Err} is a
#'   matrix and \code{Return='matrix'} (in which case the matrix will only be
#'   checked for validity).
#' @param Return output, 'matrix' (always 3x3) or 'function'.
#'
#' @return either a 3x3 matrix, or a function generating a 3x3 matrix.
#'
#' @export

ErrToM <- function(Err = NA,
                   flavour = "version2.0",
                   Return = "matrix")
{
  if (length(Err)==1 && is.na(Err) && Return == "function")  Err <- 0.1   # only used for testing
  if (!is.atomic(Err) || !length(Err) %in% c(1,9))  stop("'Err' must be a single number or 3x3 matrix")
	if (any(Err<0 | Err>1) || !is.double(Err)) stop("'Err' must be a number between 0 and 1")
  # ErrM: observed (columns) conditional on actual (rows)
  ErrM <- NULL
  ErFunc <- NULL
  if (is.matrix(Err)) {
    ErrM <- Err

  } else if (is.function(flavour)) {
    ErrM <- flavour(Err)
    if (!is.matrix(ErrM))  stop("ErFunc(E) should return a 3x3 or 4x4 matrix")
    if (!(all(dim(ErrM)==4) | all(dim(ErrM)==3)) )  stop("ErFunc(E) should return a 4x4 or 3x3 matrix")
    ErrM.B <- flavour(Err+0.1)
    if (all(ErrM.B == ErrM))  stop("ErFunc(E) is not a function of error rate E")
    ErFunc <- flavour

  } else if (flavour == "version2.0") {
    ErFunc <- function(E) {
      matrix(c((1-E/2)^2, E*(1-E/2), (E/2)^2,
              E/2, 1-E, E/2,
              (E/2)^2, E*(1-E/2), (1-E/2)^2),
             3,3, byrow=TRUE)
    }
  } else if (flavour %in% c("version1.3", "SNPchip")) {
    ErFunc <- function(E) {
      matrix(c(1-E-(E/2)^2, E, (E/2)^2,
               E/2, 1-E, E/2,
               (E/2)^2, E, 1-E-(E/2)^2),
             3,3, byrow=TRUE)
    }
  } else if (flavour %in% c("version111", "version1.1")) {
    ErFunc <- function(E) {
      matrix(c(1-E, E/2, E/2,
              E/2, 1-E, E/2,
              E/2, E/2, 1-E),
             3,3, byrow=TRUE)
    }
  } else if (flavour %in% c("version0.9", "version0.7")) {
    ErFunc <- function(E) {
      matrix(c(1-E, E, 0,
              E/2, 1-E, E/2,
              0, E, 1-E),
             3,3, byrow=TRUE)
    }
  } else {
    stop("Unknown error matrix flavour, choose 'version2.0', 'version1.3', 'version1.1', \n",
         "or specify matrix(-generating function)")
  }

  if (is.null(ErrM))  ErrM <- ErFunc(Err)

  if (!is.double(ErrM) || any(ErrM<0 | ErrM>1)) {
    stop("Error matrix values must be between 0 and 1")
  }
  if (nrow(ErrM)==4 & ncol(ErrM)==4) {
    ErrM <- shrinkEM(ErrM)
  } else if (nrow(ErrM)!=3 | ncol(ErrM)!=3){
    stop("Error matrix should be a 3x3 or 4x4 matrix")
  }
  if (!all(abs(rowSums(ErrM) - 1) < sqrt(.Machine$double.eps))) {
    stop("Error matrix rows must sum to 1")
  }

  if (Return == "matrix") {
    return( ErrM )
  } else if (Return == "function") {
    if (!is.null(ErFunc)) {
      return( ErFunc )
    } else {
      stop("Don't know how to make error function from error matrix")
    }
  } else {
    stop("Unknown Return format")
  }
}


#===============================================================================
#===============================================================================


shrinkEM <- function(EM4) {
  EM3 <- matrix(NA, 3,3)
  EM3[c(1,3), c(1,3)] <- EM4[c(1,4), c(1,4)]
  EM3[2, c(1,3)] <- EM4[2, c(1,4)]+EM4[3, c(1,4)]
  EM3[c(1,3), 2] <- EM4[c(1,4), 2] + EM4[c(1,4), 2]
  EM3[2, 2] <- sum(EM4[2:3, 2:3])
  return( EM3 )
}
