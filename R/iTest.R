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
  if(!missing(mu) && any(is.na(mu)))
    stop("'mu' must be a (vector) of numbers")
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
  cat(nlevels(tab$group), "\n")
  nGroups <- nlevels(tab$group)
  if (nGroups > NROW(tab))
    stop("not enough observations")
  if (nGroups == 1){
    ## one sample test
    C <- rbind(
      tab$X,
      rep(1, NROW(tab))
    )
    targets <- c(mu, 1)
    dof <- 1
  } else {
    if (nlevels(tab$group) > 2 && alternative != "two.sided"){
      stop("for more than 2 groups only `two.sided` alternative is allowed")
    }
    groupF <- tapply(tab$f, tab$group, sum)
    C <- rbind(
      t(sapply(
        levels(tab$group)[-1],
        function(group){
          ifelse(
            tab$group == levels(tab$group)[1],
            tab$x,
            ifelse(tab$group == group,
                   -tab$x, 0
            )
          )
        }
      )),
      t(sapply(
        levels(tab$group),
        function(group)
          ifelse(tab$group == group, 1, 0)
      ))
    )
    dof <- nlevels(tab$group) - 1
    if (length(mu) == 1){
      mu <- rep(mu, dof)
    }
    targets <- c(mu, groupF)

  }
  p0 <- iProjector(C, targets, v = tab$f)
  statistic = if (p0$converged){
    p1 <- iProjector(C, C %*% tab$f, v = p0$p)
    if (p1$converged){
      iDivergence(p1$p, p0$p)
    } else{
      -2
    }
  } else{
    -1
  }
  p.value <- if (statistic >= 0){
    if (alternative == "two.sided"){
      pchisq(2 * sum(tab$Freq) * statistic, df = dof, lower.tail = FALSE)
    } else{
      if (alternative == "greater"){
        #hmm
        if (C[1,] %*% tab$f > 0){
          pchisq(2 * sum(tab$Freq) * statistic, df = dof, lower.tail = FALSE) / 2
        } else{
          0.5 + pchisq(2 * sum(tab$Freq) * statistic, df = dof, lower.tail = TRUE) / 2
        }
      } else{
        if (C[1,] %*% tab$f < 0){
          pchisq(2 * sum(tab$Freq) * statistic, df = dof, lower.tail = FALSE) / 2
        } else{
          0.5 + pchisq(2 * sum(tab$Freq) * statistic, df = dof, lower.tail = TRUE) / 2
        }
      }

    }
  } else{
    0
  }
  p.value
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

