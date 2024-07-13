#' @export
i.test <- function(x, ...) UseMethod("i.test")

#' @export
i.test.character <- function(x, y = NULL, ...){
  x <- as.factor(x)
  i.test(x, y, ...)
}
#' @export
i.test.numeric <-
  function(x, y = NULL, alternative = c("two.sided", "less", "greater"),
           mu = 0, paired = FALSE, conf.level = 0.95)
{
  alternative <- match.arg(alternative)
  ## check the conf.level: should be a single number between 0 and 1
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  ## check y
  if (is.null(y)){
    ## one sample test
    dname <- deparse(substitute(x))
    if (paired) stop("'y' is missing for paired test")
    xok <- !is.na(x)
    yok <- NULL
    df <- data.frame(
      x = x[xok],
      group = factor(rep("X", sum(xok)))
    )
  } else{
    if (class(y) == "character"){
      y <- as.factor(y)
    }
    if (class(y) == "factor"){
      ## we have a group attribute
      if (paired) stop("paired design not possible for groups")
      if (nlevels(y) < 2) stop("only", nlevels(y), "groups. Need at least 2")

      dname = paste(levels(y), collapse = ",")
      xok <- yok <- complete.cases(x,y)
      df <- data.frame(
        x = x[xok],
        group = y[yok]
      )
    } else if (class(y) == "numeric"){
      ## two sample test
      dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
      if (paired){
        xok <- yok <- complete.cases(x,y)
        df <- data.frame(
          x = x[xok] - y[yok],
          group = factor(rep("X", sum(xok)))
        )
      } else{
        yok <- !is.na(y)
        xok <- !is.na(x)
        df <- data.frame(
          x = c(x[xok], y[yok]),
          group = factor(c(rep("X", sum(xok)), rep("Y", sum(yok))))
        )
      }

    }
  }
  tab <- toTab(df)
  nGroups <- nlevels(tab)
  if (nGroups < NROW(tab))
    stop("not enough observations")
  if (nGroups == 1){
    ## one sample test
    C <- rbind(
      tab$X,
      rep(1, NROW(tab))
    )
    targets <- c(mu, 1)
  } else
    groupF <- tapply(tab$f, tab$group, sum)
    if (nlevels(tab$group) == 2){
      C <- rbind(
        ifelse(tab$group == levels(tab$group)[1],
               tab$X, -tab$X) / groupF[tab$group],
        t(sapply(
          levels(tab$group),
          function(group)
            ifelse(tab$group == group, 1, 0)
        ))
      )
      targets <- c(mu, groupF)

    } else{
        if (alternative != "two.sided"){
          stop("for more than 2 groups only `two.sided` alternative is allowed")
        }

  }
}
#' @export
i.test.factor <- function(x, y = NULL){

}

#' @export
i.test.formula <- function(formula, data, ...){

}
## want to implement the generic i-test
## 1. x = numeric, y = NULL, mu, one sample t-test
## 2. x = numeric, y = factor, two+ sample t-test/anova
## 3. x = factor, y = NULL, mu, one sample binomial test
## 4. x = factor, y = factor two+ sample barnard test
## 5. x = factro, y = factor two+ sample fisher's exact test
## 6. x = numeric, y = numeric, model test

