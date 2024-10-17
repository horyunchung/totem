#' toEntities
#'
#' Converts a data table to the entity representation
#'
#' @import dplyr
#'
#' @param data a data.frame to be converted
#' @param select a vector of column names (attributes). Default is to take all columns.
#' @param cartesian a logical indicating whether the entities should comprise
#'  only the observed tupples of manifestations (\code{FALSE}, default) or
#'  the cartesian product (all possible combinations) of the attribute
#'  manifestations (\code{TRUE})
#'
#' @return A data.frame with unique rows of the data table with selected columns with
#' a new column \code{empirical} denoting the empirical relative frequency distribution.
#' The attribute \code{attr(*, "N")} corresponds to the number of data points (rows).
#' @export
#'
toEntities <- function(data, select = NULL, cartesian = FALSE){
  ## check input
  if (!is.data.frame(data)){
    stop("data is not a data.frame")
  }

  entities <- if (cartesian == TRUE){
    res <- if (is.null(select)){
      data.frame(table(data))
    } else{
      if (all(select %in% colnames(data))){
        data.frame(table(data[, select]))
      } else{
        stop("selected colnames not found:", select[!select %in% colnames(data)])
      }
    }
    for (col in colnames(res)[-NCOL(res)]){
      if (class(data[, col]) == "numeric"){
        res[, col] = as.numeric(as.character(res[, col]))
      }
    }
    colnames(res)[which(colnames(res) == "Freq")] <- "empirical"
    res
  } else{
    if (is.null(select)){
      data %>% count(!!!syms(colnames(data)), name = "empirical")
    } else{
      if (all(select %in% colnames(data))){
        data %>% count(!!!syms(select), name = "empirical")
      } else{
      stop("selected colnames not found:", select[!select %in% colnames(data)])
      }
    }
  }
  attr(entities, "N") = nrow(data)
  entities$empirical <- entities$empirical / attr(entities, "N")

  entities
}
