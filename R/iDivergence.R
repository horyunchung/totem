#' iDivergence
#'
#' Computes the I-divergence (a.k.a. as Kullback-Leibler divergence) of distribution
#' \code{p} from \code{q}
#'
#' @param p a numeric vector for the candidate distribution
#' @param q a numeric vector for the reference distribution
#' (if it is \code{NULL} the entropy is calculated)
#'
#' @returns the I-divergence
#' @export
iDivergence <- function(p, q = NULL){
  ## check input
  if (! is.vector(p)){
    stop("p is not a vector")
  } else{
    if (class(p) != "numeric"){
      stop("p is not numeric")
    }
  }
  if (! is.null(q)){
    if (! is.vector(q)){
      stop("q is not a vector")
    } else{
      if (class(q) != "numeric"){
        stop("q is not numeric")
      }
    }
  }

  sel <- which(p > 0)
  ## entropy
  if (is.null(q)){
    -sum(p[sel] * log(p[sel]))
  } else{
    sum(p[sel] * (log(p[sel]) - log(q[sel])))
  }
}
