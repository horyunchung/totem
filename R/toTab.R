#' toTab
#'
#'
#' @export
#' @import dplyr
toTab <- function(df, select = NULL){
  tab <- if (is.null(select)){
    data %>% count(!!!syms(colnames(data)), name = "Freq")
  } else{
    data %>% count(!!!syms(select), name = "Freq")
  }
  tab$f <- tab$Freq / sum(tab$Freq)
  tab
}
