#' i.test
#'
#' Performs a two way $I$-test.
#'
#' @export
i.test <- function(formula, data, alternative = c("two.sided", "less", "greater"),
                   mu = 0, fix = FALSE
                   )
{
  alternative <- match.arg(alternative)

  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) > 1L))
    stop("'formula' missing or incorrect")

  ## extract the attribute manifestations
  dname <- paste(attr(attr(terms(formula), "factors"), "dimnames")[[1]], collapse = " by ")
  yx <- eval(attr(terms(formula), "variables"), envir = data)
  yx <- as.data.frame(yx, col.names = if (length(yx) == 2) c("y", "x") else "y")
  ## retain only complete rows
  yx = subset(yx, subset = complete.cases(yx))

  ## convert to entities
  entities <- toEntities(yx)

  ## determine the attribute types
  classY <- NULL
  classX <- NULL

  gY <- NULL
  gX <- NULL

  if (is.numeric(yx$y)){
    classY <- "numeric"
    gY <- entities$y
  } else{
    ## convert characters to factors
    if (is.character(yx$y))
      yx$y <- factor(yx$y)

    if(is.factor(yx$y)){
      if (nlevels(yx$y) != 2){
        stop("categorical 'response' should have exactly two levels")
      }
      classY <- "factor"
      gY <- ifelse(entities$y == levels(entities$y)[2], 1, 0)
    } else{
      stop("'response' is neither numeric nor categorical")
    }
  }
  class(yx$y) <- classY

  oneSample <- NCOL(yx) == 1

  if (oneSample == FALSE){
    if(is.numeric(yx$x)){
      classX <- "numeric"
      gX <- entities$x
    } else{
      if (is.character(yx$x))
        yx$x <- factor(yx$x)

      if(is.factor(yx$x)){
        if (nlevels(yx$x) != 2)
          stop("catergorical 'predictors' should have exactly two levels")
        classX <- "factor"
        gX <- ifelse(entities$x == levels(entities$x)[2], 1, 0)
      } else{
        stop("'predictor' is neither numeric nor categorical")
      }

    }
    class(yx$x) = classX
    if(is.logical(fix)){
      if (fix == TRUE){
        if (classX == "factor"){
          fix <- sum(entities$f * gX)
        } else{
          fix <- sum(entities$f * gY / gX)
        }
      }
    } else if (!is.numeric(fix) || length(fix) != 1L){
      stop("invalid value for fix")
    }
  }


  if(classY == "factor" && classX == "numeric")
    stop("categorical 'response' with numeric 'predictor' is not supported")

  ## so what do we do?
  idiv <- NULL
  estimate <- NULL
  if (fix == FALSE && !oneSample){
    suppressWarnings(
      optH <- optimize(
        function(muX){
          i.test_internal(
            entities = entities,
            v = entities$empirical,
            mu = mu,
            classX = classX,
            gY = gY, gX = gX,
            oneSample = oneSample,
            fix = muX
          )$idiv
        },
        interval = if (classX == "factor") c(0, 1) else range(gY / gX),
        tol = .Machine$double.eps
      )
    )
    pH <- i.test_internal(
      entities = entities,
      v = entities$empirical,
      mu = mu,
      classX = classX,
      gY = gY, gX = gX,
      oneSample = FALSE,
      fix = optH$minimum
    )
    estimate <- pH$estimate

    if (is.finite(optH$objective)){
      suppressWarnings(
        optA <- optimize(
          function(muX){
            i.test_internal(
              entities = entities,
              v = pH$p$p,
              mu = mu,
              classX = classX,
              gY = gY, gX = gX,
              oneSample = FALSE,
              fix = muX
            )$idiv
          },
          interval = if (classX == "factor") c(0, 1) else range(gY / gX),
          tol = .Machine$double.eps
        ))
      pA <- i.test_internal(
        entities = entities,
        v = pH$p$p,
        mu = estimate,
        classX = classX,
        gY = gY, gX = gX,
        oneSample = FALSE,
        fix = optA$minimum
      )
      idiv <- pA$idiv
    } else{
      idiv <- NA
    }

  } else{
    ## hypothesis distribution
    pH <- i.test_internal(
      entities = entities,
      v = entities$empirical,
      mu = mu,
      classX = classX,
      gY = gY, gX = gX,
      oneSample = oneSample,
      fix = fix
    )
    estimate <- pH$estimate
    if (is.finite(pH$idiv)){
      pA <- i.test_internal(
        entities = entities,
        v = pH$p$p,
        mu = estimate,
        classX = classX,
        gY = gY, gX = gX,
        oneSample = oneSample,
        fix = fix
      )
      idiv <- pA$idiv
    } else{
      idiv <- NA
    }

  }
  N = attr(entities, "N")
  statistic <- 2 * N * idiv
  names(statistic) = "i"
  names(mu) = if (oneSample){
    "mean"
  } else{
    if (classY == "factor"){
      "difference of rates"
    } else{
      if (classX == "factor"){
        "difference of means"
      } else{
        "bias"
      }
    }
  }
  p.value <- if(alternative == "two.sided"){
    pchisq(statistic, df = 1, lower.tail = FALSE)
  } else if (alternative == "less"){
    if (estimate < mu){
      pchisq(statistic, df = 1, lower.tail = FALSE) / 2
    } else{
      0.5 + pchisq(statistic, df = 1, lower.tail = TRUE) / 2
    }
  } else{
    if (estimate < mu){
      pchisq(statistic, df = 1, lower.tail = FALSE) / 2
    } else{
      0.5 + pchisq(statistic, df = 1, lower.tail = TRUE) / 2
    }
  }

  method <- if (oneSample){
    if(classY == "factor"){
      "One Sample I-test for a rate"
    } else{
      "One Sample I-test for a mean"
    }
  } else{
    if (classY == "numeric" && classX == "numeric"){
      "I-test for a two-factor relationship"
    } else{
      paste("Two Sample I-test for",
          if (classY == "factor"){
            "rates"
          } else{
            "means"
          }

      )
    }
  }



  rval <- list(
    statistic = statistic,
    parameter = c(iDivergence = idiv),
    p.value = p.value,
    estimate = estimate,
    null.value = mu,
    alternative = alternative,
    method = method,
    data.name = dname,
    sample.size = attr(entities, "N")
  )
  class(rval) <- "htest"
  return(rval)

}


#' @export
i.test_internal <- function(entities, v, mu, classX, gY, gX, oneSample, fix){

  C <- NULL
  targets <- NULL
  estimate <- NULL
  ## one sample test
  if (oneSample){
    C <- rbind(
      gY,
      1
    )
    estimate = sum(entities$empirical * gY)
    names(estimate) = "mean Y"
    targets <- c(mu, 1)
  ## two samples
  } else{
    if (classX == "factor"){
      C <- rbind(
        ifelse(entities$x == levels(entities$x)[2], gY / fix, -gY / (1 - fix)),
        gX,
        1
      )
      marginal <- sum(entities$empirical * gX)
      estimate <- sum(entities$empirical * ifelse(entities$x == levels(entities$x)[2], gY / marginal, -gY / (1 - marginal)))
      names(estimate) = paste("difference", if (class(gY) == "numeric") "means" else "rates")
      targets <- c(mu, fix, 1)
    } else{
      C <- rbind(
        gY - fix * gX,
        gY / gX,
        1
      )
      marginal <- sum(entities$empirical * gY / gX)
      estimate <- sum(entities$empirical * (gY - marginal * gX))
      names(estimate) = "bias"
      targets <- c(mu, fix, 1)
    }

  }

  p <- iProjector(C, targets, v)
  idiv <- if(p$converged){
    iDivergence(p$p, v)
  } else{
    Inf
  }
  list(
    p = p,
    idiv = idiv,
    estimate = estimate
  )
}

#' @export
i.test.character <- function(y, x = NULL, ...){
  y <- as.factor(y)
  i.test(y, x, ...)
}
#' @export
i.test.numeric <-
  function(y, x = NULL, alternative = c("two.sided", "less", "greater"), fix = c("y", "both", "none"),
           mu = 0, paired = FALSE, conf.level = NULL, tolerance = .Machine$double.eps)
{
  alternative <- match.arg(alternative)
  fix <- match.arg(fix)

  ## check the conf.level: should be a single number between 0 and 1
  if(!missing(conf.level) &&
     (length(conf.level) != 1 || !is.finite(conf.level) ||
      conf.level < 0 || conf.level > 1))
    stop("'conf.level' must be a single number between 0 and 1")
  ## check mu: should be a single number
  if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
    stop("'mu' must be a single number")

  ## check the independent variable
  ## one sample I-test
  if (is.null(x)){
    dname <- deparse(substitute(y))
    if (paired) stop("'x' is missing for paired test")
    yok <- !is.na(y)
    xok <- NULL
    data <- data.frame(
      y = y[xok]
    )
  ## two sample I-test
  } else{
    if (class(x) == "character"){
      x <- as.factor(x)
    }
    if (class(x) == "factor"){
      ## we have a group attribute
      if (paired) stop("paired design not possible for groups")
      if (nlevels(x) != 2) stop("Currently, only two groups are supported.")

      dname = paste(levels(x), collapse = ",")
      yok <- xok <- complete.cases(y,x)
      data <- data.frame(
        y = y[yok],
        group = x[xok]
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

    groupF <- tapply(tab$f, tab$group, sum)
    estimate <- tapply(tab$x * tab$f, tab$group, sum) / groupF
    names(estimate) = paste("mean of", names(estimate))
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
  p.value <- if(alternative == "two.sided"){
    pchisq(2 * N * idiv, df = dof, lower.tail = FALSE)
  } else if (alternative == "less"){
    if (paired){
      if (estimate < 0){
        pchisq(2 * N * idiv, df = dof, lower.tail = FALSE) / 2
      } else{
        0.5 + pchisq(2 * N * idiv, df = dof, lower.tail = TRUE) / 2
      }
    } else{
      if (estimate < mu){
        pchisq(2 * N * idiv, df = dof, lower.tail = FALSE) / 2
      } else{
        0.5 + pchisq(2 * N * idiv, df = dof, lower.tail = TRUE) / 2
      }
    }
  } else{
    if (paired){
      if (estimate > 0){
        pchisq(2 * N * idiv, df = dof, lower.tail = FALSE) / 2
      } else{
        0.5 + pchisq(2 * N * idiv, df = dof, lower.tail = TRUE) / 2
      }
    } else{
      if (estimate > mu){
        pchisq(2 * N * idiv, df = dof, lower.tail = FALSE) / 2
      } else{
        0.5 + pchisq(2 * N * idiv, df = dof, lower.tail = TRUE) / 2
      }
    }
  }
  method <- paste(method, " of means")
  names(mu) <- if (nGroups < 3 || length(mu) == 1)
    if(paired || !is.null(y)) "difference in means" else "mean"
  else
    paste("difference in means", levels(tab$group)[1], "and", levels(tab$group)[-1])


  rval <- list(
    statistic = 2 * N * idiv,
    parameter = c(N = sum(tab$Freq), iDivergence = if (idiv >= 0) idiv else NA, df = dof),
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
  ## canonicalize the arguments
  formula <- as.formula(formula)

  if(!is.list(data) && !is.environment(data))
    stop("'data' must be a list or an environment")

  mf <- cl <- match.call()		# for creating the model frame
  varNames <- all.vars(formula)
  ## left hand side is the response
  ## right hand side is the model
  ## we should have no parameters
  form2 <- formula; form2[[2L]] <- 0
  varNamesRHS <- all.vars(form2)

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()

  ## all names in formula should represent actual
  ## variables
  ## This aux.function needs to be as complicated because
  ## exists(var, data) does not work (with lists or dataframes):
  lenVar <- function(var) tryCatch(length(eval(as.name(var), data, env)),
                                   error = function(e) -1L)
  if (length(varNames)){
    n <- vapply(varNames, lenVar, 0)
    if(any(not.there <- n == -1L)) {
      stop("")
    }
  }




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

computePvalue <- function(idiv, N, dof, alternative, direction){
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

test <- function(formula, data){
  mf <- cl <- match.call()
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf)
  x <- model.matrix(mf)
  print(str(x))
}

