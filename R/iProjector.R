#' iProjector
#'
#'
#'
#' @export
iProjector <- function(M, targets, v = NULL, eps = .Machine$double.eps, maxIter = 10000L){

  ## check whether M is a matrix
  if(!is.matrix(M)){
    stop("M is not a matrix")
  }
  ## the number of entities
  numberEntities <- ncol(M)

  ## coerce the targets into a column matrix
  targets <- matrix(targets, ncol = 1)

  ## check whether the number of targets is equal to the number of rows in M
  if (nrow(targets) != nrow(M)){
    stop("The number of targets does not match the number of rows in M: # rows provided by M = ", nrow(M), " and # targets = ", nrow(targets))
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
      stop("The number of entities do not match: # provided by M = ", numberEntities, " and # provided by reference distribution = ", length(v))
    }
    iProjection <- matrix(
      data = v,
      ncol = 1
    )
  }

  ## precompute transpose of M
  MT <- t(M)

  ## set the status
  converged = FALSE
  status = NULL
  error = NULL
  for (iter in seq_len(maxIter)){
    #cat(iter, "\n")
    ## calculate the Jacobian, note since R stores matrix in the
    ## column major format, we multiply iProjection with
    ## M transpose
    update <- tryCatch(
      {
        inverseJacobian <- solve(M %*% (MT * iProjection[,1]))
        MT %*% inverseJacobian %*% ((M %*% iProjection) - targets)

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
    if (isTRUE(all.equal(updatedProjection[,1], iProjection[,1], tolerance = eps))){
      updatedtargets <- M %*% updatedProjection
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
