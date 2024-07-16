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
           mu = 0, paired = FALSE, conf.level = 0.95, tolerance = .Machine$double.eps)
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
  idiv <- computeStatistic(C, targets, tab, tolerance = tolerance)
  N = sum(tab$Freq)
  p.value <- computePvalue(idiv, N, dof, alternative)
  method <- paste(method, " of means")
  names(mu) <- if (nGroups < 3 || length(mu) == 1)
    if(paired || !is.null(y)) "difference in means" else "mean"
  else
    paste("difference in means", levels(tab$group)[1], "and", levels(tab$group)[-1])


  rval <- list(
    statistic = p.value$statistic,
    parameter = c(N = sum(tab$Freq), iDivergence = if (idiv >= 0) idiv else NA, df = dof),
    p.value = p.value$p.value,
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
  function(x, y, alternative = c("two.sided", "less", "greater"),
           mu = 0, fix.margins = c("group", "both", "none"), conf.level = 0.95,
           tolerance = .Machine$double.eps){
    alternative <- match.arg(alternative)
    fix.margins <- match.arg(fix.margins)
    ## check the conf.level: should be a single number between 0 and 1
    if(!missing(conf.level) &&
       (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    if(!missing(mu) && any(is.na(mu)))
      stop("'mu' must be a (vector) of numbers")
    dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    if (class(y) == "character"){
      y <- factor(y)
    }

    if (length(x) != length(y))
      stop("length of 'x' should match length of 'group'")
    if (nlevels(x) < 2)
      stop("there should be at least two levels for the outcome 'x'")
    if (nlevels(y) < 2)
      stop("there should be at least two groups")


    xok <- yok <- complete.cases(x,y)

    data <- data.frame(
      x = x[xok],
      group= y[yok]
    )

    tab <- toTab(data)

    ## fix the prevalence of the groups
    ## Barnard's test
    dof <- (nlevels(tab$group) - 1) * (nlevels(tab$x) - 1)
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
      targets <- c(rep(0, nlevels(tab$group) - 1), groupF)
      p0 <- NULL
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
      targets <- c(rep(0, nlevels(tab$group) - 1), groupF, xF[1])

      f1 <- tapply(tab$f, tab$x, sum)
      f2 <- tapply(tab$f, tab$group, sum)

      p0 <- list(
        p = f1[tab$x] * f2[tab$group],
        converged = TRUE
      )
    }

    idiv <- computeStatistic(C, targets, tab, p0, tolerance = tolerance)
    N = sum(tab$Freq)
    p.value <- computePvalue(idiv, N, dof, alternative)
    estimate <- with(subset(tab, subset = tab$x == levels(tab$x)[1]), tapply(f, group, sum)) / tapply(tab$f, tab$group, sum)
    names(mu) = paste0("difference in p(X='", levels(tab$x)[1], "'|group)")
    method = "i-test"

    rval <- list(
      statistic = p.value$statistic,
      parameter = c(N = N, iDivergence = idiv, df = dof),
      p.value = p.value$p.value,
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
i.test.formula <- function(formula, data, ...){

}
## want to implement the generic i-test
## 1. x = numeric, y = NULL, mu, one sample t-test
## 2. x = numeric, y = factor, two+ sample t-test/anova
## 3. x = factor, y = NULL, mu, one sample binomial test
## 4. x = factor, y = factor two+ sample barnard test
## 5. x = factro, y = factor two+ sample fisher's exact test
## 6. x = numeric, y = numeric, model test



computeStatistic <- function(C, targets, tab, p0 = NULL, tolerance = .Machine$double.eps){
  if (is.null(p0))
    p0 <- iProjector(C, targets, v = tab$f)

  if (p0$converged){
    mask <- which(p0$p >= tolerance)
    if (length(mask) != ncol(C)){
      warning("the probability of some entities was set to zero")
    }
    p1 <- iProjector(C[, mask], C %*% tab$f, v = p0$p[mask])
    if (p1$converged){
      iDivergence(p1$p, p0$p)
    } else{
      warning("could not project back to the data")
      -2
    }
  } else{
    warning("the hypothesis is inconsistent with the data")
    -1
  }
}

computePvalue <- function(idiv, N, dof, alternative){
  p.value <- if (idiv >= 0){
    statistic <- 2 * N * idiv
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
  } else if (idiv == -2){
    statistic <- Inf
    0
  } else if (idiv == -1){
    statistic <- NA
    NA
  }
  names(statistic) = "i"
  list(
    p.value = p.value, statistic = statistic
  )
}
