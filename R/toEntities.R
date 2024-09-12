#' toEntities
#'
#' Converts a data table to the entity representation
#'
#' @import dplyr
#'
#' @param data A data.frame to be converted
#' @param select A vector of column names (attributes). Default is to take all columns.
#' @return A data.frame with unique rows of the data table with selected columns with
#' a new column \code{empirical} denoting the empirical relative frequency distribution.
#' The attribute \code{attr(*, "N")} corresponds to the number of data points (rows).
#' @export
#'
toEntities <- function(data, select = NULL){
  entities <- if (is.null(select)){
    data %>% count(!!!syms(colnames(data)), name = "empirical")
  } else{
    if (all(select %in% colnames(data))){
      data %>% count(!!!syms(select), name = "empirical")
    } else{
      stop("selected colnames not found:", select[!select %in% colnames(data)])
    }
  }
  attr(entities, "N") = nrow(data)
  entities$empirical <- entities$empirical / attr(entities, "N")

  entities
}
