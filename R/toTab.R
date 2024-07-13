#' toTab
#'
#'
#' @export
#' @import dplyr
toTab <- function(df, select = NULL){
  tab <- if (is.null(select)){
    df %>% count(!!!syms(colnames(df)), name = "Freq")
  } else{
    if (all(select %in% colnames(df))){
      df %>% count(!!!syms(select), name = "Freq")
    } else{
      stop("selected colnames not found:", select[!select %in% colnames(df)])
    }
  }
  tab$f <- tab$Freq / sum(tab$Freq)
  tab
}
