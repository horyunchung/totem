#' iProjector
#'
#' Compute the I-projection given constraints in C and expected values in targets
#'
#' @param C a matrix with D constraints in the rows for |E| entities in the columns, must be of full row rank
#' @param targets a matrix/vector of with D rows/entries with the expected values for the constraints in C
#' @param v the reference distribution, defaults to the uniform distribution
#' @param eps the tolerance
#' @param maxIter the maximal number of iterations, defaults to 10 000
#'
#' @returns a list with following entries
#'
#' @export
iProjector <- function(C, targets, v = NULL, eps = .Machine$double.eps, maxIter = 10000L){

  ## check whether M is a matrix
  if(!is.matrix(C)){
    stop("M is not a matrix")
  }
  ## the number of entities
  numberEntities <- ncol(C)

  ## coerce the targets into a column matrix
  targets <- matrix(targets, ncol = 1)

  ## check whether the number of targets is equal to the number of rows in M
  if (nrow(targets) != nrow(C)){
    stop("The number of targets does not match the number of rows in C: # rows provided by C = ", nrow(C), " and # targets = ", nrow(targets))
  }

  ## reference distribution
  ### if not given set it to the uniform
  if (is.null(v)){
    iProjection <- matrix(
      data = rep(1 / numberEntities, numberEntities),
      ncol = 1
    )
    ### otherwise take the provided one
  } else{
    #### check dimensions
    if (length(v) != numberEntities){
      stop("The number of entities do not match: # provided by C = ", numberEntities, " and # provided by reference distribution = ", length(v))
    }
    iProjection <- matrix(
      data = v,
      ncol = 1
    )
  }

  ## precompute transpose of C
  CT <- t(C)

  ## set the status
  converged = FALSE
  status = NULL
  error = NULL
  for (iter in seq_len(maxIter)){
    update <- tryCatch(
      {
        inverseJacobian <- solve(C %*% (CT * iProjection[,1]))
        CT %*% inverseJacobian %*% ((C %*% iProjection) - targets)

      },
      error = function(e){
        e
      }
    )
    if(!is.matrix(update)){
      error = update
      break
    }
    updatedProjection<- iProjection * exp(-update)
    cat(isTRUE(all.equal(updatedProjection[,1], iProjection[,1], tolerance = eps)), "\n")
    if (isTRUE(all.equal(updatedProjection[,1], iProjection[,1], tolerance = eps))){
      updatedtargets <- C %*% updatedProjection
      if (isTRUE(all.equal(updatedtargets[,1], targets[,1]))){
        status <- "converged"
        converged <- TRUE
        iProjection <- updatedProjection
      } else{
        status <- "targets not fulfilled"
        iProjection <- updatedProjection
      }
      break
    }
    iProjection <- updatedProjection
    cat(iDivergence(iProjection, tab$f), "\n")
  }

  list(
    p = iProjection[,1],
    converged = converged,
    info = list(
      status = status,
      iter = iter,
      error = error
    )
  )
}
