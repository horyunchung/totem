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

  nGroups <- nlevels(tab$group)
  if (nGroups > NROW(tab))
    stop("not enough observations")
  if (nGroups == 1){
    ## one sample test
    method <- if (paired)  "Paired i-test" else "One Sample i-test"
    mx <- sum(tab$f * tab$x)
    estimate <- setNames(mx, if (paired) "mean of difference" else "mean of X")
    C <- rbind(
      tab$x,
      rep(1, NROW(tab))
    )
    targets <- c(mu, 1)
    dof <- 1
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
  } else {
    if (nlevels(tab$group) > 2 && alternative != "two.sided"){
      stop("for more than 2 groups only `two.sided` alternative is allowed")
    }
    method <- if (nGroups  == 2) "Two sample i-test" else paste(nGroups, "sample i-test")
    estimate <- tapply(tab$x * tab$f, tab$group, sum)
    names(estimate) = paste("mean of", names(estimate))
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
      mus <- rep(mu, dof)
    } else{
      if (length(mu) == dof){
        mus <- mu
      } else{
        stop("number of provided mean differences do not match the number of differences")
      }
    }


    targets <- c(mus, groupF)

  }
  p0 <- iProjector(C, targets, v = tab$f)
  idiv = if (p0$converged){
    p1 <- iProjector(C, C %*% tab$f, v = p0$p)
    if (p1$converged){
      iDivergence(p1$p, p0$p)
    } else{
      -2
    }
  } else{
    -1
  }
  p.value <- if (idiv >= 0){
    statistic <- 2 * sum(tab$Freq) * idiv
    if (alternative == "two.sided"){
      pchisq(statistic, df = dof, lower.tail = FALSE)
    } else{
      if (alternative == "greater"){
        #hmm
        if (C[1,] %*% tab$f > 0){
          pchisq(statistic, df = dof, lower.tail = FALSE) / 2
        } else{
          0.5 + pchisq(statistic, df = dof, lower.tail = TRUE) / 2
        }
      } else{
        if (C[1,] %*% tab$f < 0){
          pchisq(statistic, df = dof, lower.tail = FALSE) / 2
        } else{
          0.5 + pchisq(statistic, df = dof, lower.tail = TRUE) / 2
        }
      }

    }
  } else{
    statistic <- Inf
    0
  }
  method <- paste(method, " of means")
  names(mu) <- if (nGroups < 3 || length(mu) == 1)
    if(paired || !is.null(y)) "difference in means" else "mean"
  else
    paste("difference in means", levels(tab$group)[1], "and", levels(tab$group)[-1])
  names(statistic) <- "i"

  rval <- list(
    statistic = statistic,
    parameter = c(N = sum(tab$Freq), iDivergence = idiv, df = dof),
    p.value = p.value,
    estimate = estimate,
    null.value = mu,
    alternative = alternative,
    method = method,
    data.name = dname
  )
  class(rval) <- "htest"
  return(rval)
}
#' @export
i.test.factor <-
  function(x, group, alternative = c("two.sided", "less", "greater"),
           p = 0.5, fix.margins = c("group", "both", "none"), conf.level = 0.95){
    alternative <- match.arg(alternative)
    fix.margins <- match.arg(fix.margins)
    ## check the conf.level: should be a single number between 0 and 1
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    if(!missing(p) && any(is.na(p)))
      stop("'p' must be a (vector) of numbers")
    if (class(group) == "character"){
      group <- as.factor(group)
    }

    if (length(x) != length(group))
      stop("length of 'x' should match length of 'group'")
    if (nlevels(x) != 2)
      stop("there should be only two levels for the outcome 'x'")
    if (nlevels(group) < 2)
      stop("there should be at least two groups")

    xok <- groupok <- complete.cases(x,group)

    data <- data.frame(
      x = x[xok],
      group= group[groupok]
    )

    tab <- toTab(data)

    ## fix nothing,
    ## Pearson's chisquared test

    cat(fix.margins, "\n")

    ## fix the prevalence of the groups
    ## Barnard's test
    dof <- nlevels(group) - 1
    if (fix.margins == "group"){
      groupF <- tapply(tab$f, tab$group, sum)
      C <- rbind(
        t(
          sapply(
            levels(tab$group)[-1],
            function(group){
              ifelse(
                tab$group == levels(tab$group)[1],
                ifelse(
                  tab$x == levels(tab$x)[1],
                  1,
                  0
                ),
                ifelse(
                  tab$group == group,
                  ifelse(
                    tab$x == levels(tab$x)[1],
                    -1,
                    0
                  ),
                  0
                )
              ) / groupF[tab$group]
            } ## function(group)
          ) ## sapply

        ), ## t
        t(
          sapply(
            levels(tab$group),
            function(group){
              ifelse(
                tab$group == group,
                1,
                0
              )
            }
          )
        )
      )
      targets <- c(rep(0, nlevels(group) - 1), groupF)
    } else if(fix.margins == "both"){
      ## fix the prevalence of the groups
      ## and the x
      ## Fisher's exact test
      groupF <- tapply(tab$f, tab$group, sum)
      xF <- tapply(tab$f, tab$x, sum)
      C <- rbind(
        t(
          sapply(
            levels(tab$group)[-1],
            function(group){
              ifelse(
                tab$group == levels(tab$group)[1],
                ifelse(
                  tab$x == levels(tab$x)[1],
                  1,
                  0
                ),
                ifelse(
                  tab$group == group,
                  ifelse(
                    tab$x == levels(tab$x)[1],
                    -1,
                    0
                  ),
                  0
                )
              ) / groupF[tab$group]
            } ## function(group)
          ) ## sapply

        ), ## t
        t(
          sapply(
            levels(tab$group),
            function(group){
              ifelse(
                tab$group == group,
                1,
                0
              )
            }
          )
        ),
        ifelse(tab$x == levels(tab$x)[1], 1, 0)
      )
      targets <- c(rep(0, nlevels(group) - 1), groupF, xF[1])

    }

    p0 <- iProjector(C, targets, v = tab$f)
    cat(C %*% p0$p, "\n")
    cat(C %*% tab$f, "\n")
    idiv = if (p0$converged){
      p1 <- iProjector(C, C %*% tab$f, v = p0$p)
      if (p1$converged){
        iDivergence(p1$p, p0$p)
      } else{
        -2
      }
    } else{
      -1
    }
    p.value <- if (idiv >= 0){
      statistic <- 2 * sum(tab$Freq) * idiv
      if (alternative == "two.sided"){
        pchisq(statistic, df = dof, lower.tail = FALSE)
      } else{
        if (alternative == "greater"){
          #hmm
          if (C[1,] %*% tab$f > 0){
            pchisq(statistic, df = dof, lower.tail = FALSE) / 2
          } else{
            0.5 + pchisq(statistic, df = dof, lower.tail = TRUE) / 2
          }
        } else{
          if (C[1,] %*% tab$f < 0){
            pchisq(statistic, df = dof, lower.tail = FALSE) / 2
          } else{
            0.5 + pchisq(statistic, df = dof, lower.tail = TRUE) / 2
          }
        }

      }
    } else{
      statistic <- Inf
      0
    }
    p.value



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

