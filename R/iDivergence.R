#' iDivergence
#' @export
iDivergence <- function(p, q = NULL){
  sel <- which(p > 0)
  ## entropy
  if (is.null(q)){
    -sum(p[sel] * log(p[sel]))
  } else{
    sum(p[sel] * (log(p[sel]) - log(q[sel])))
  }
}
