#' i.test
#'
#' Performs a two way $I$-test.
#'
#' @param formula an object of class "formula"
#' @param data a data.frame containg the variables in the formula
#' @param alternative a character string specifying the alternative hypothesis,
#' must be one of "two.sided" (default), "greater" or "less". You can specify
#' just the initial letter.
#' @param mu a number indicating the true value of the rate/mean difference or
#' the bias
#' @param fix a logical indicating whether to fix the predictor marginal
#' to the sample marginal (fix = TRUE) or not (fix = FALSE). Additionally also
#' the value of the predictor marginal can be given

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
  if (is.character(yx$x))
    yx$x <- factor(yx$x)

  ## determine the attribute types
  classY <- NULL
  classX <- NULL

  if (is.numeric(yx$y)){
    classY <- "numeric"
  } else{
    ## convert characters to factors
    if (is.character(yx$y))
      yx$y <- factor(yx$y)

    if(is.factor(yx$y)){
      if (nlevels(yx$y) != 2){
        stop("categorical 'response' should have exactly two levels")
      }
      classY <- "factor"
    } else{
      stop("'response' is neither numeric nor categorical")
    }
  }
  class(yx$y) <- classY


  oneSample <- NCOL(yx) == 1

  if (oneSample == FALSE){
    if(is.numeric(yx$x)){
      classX <- "numeric"
    } else{
      if(is.factor(yx$x)){
        if (nlevels(yx$x) != 2)
          stop("catergorical 'predictors' should have exactly two levels")
        classX <- "factor"
      } else{
        stop("'predictor' is neither numeric nor categorical")
      }

    }
    class(yx$x) = classX
  }
  entities <- toEntities(yx)

  gY <- if (classY == "numeric"){
    entities$y
  } else{
    ifelse(entities$y == levels(entities$y)[2], 1, 0)
  }
  gX <- NULL
  if (! oneSample){
    gX <- if (classX == "numeric"){
      entities$x
    } else{
      ifelse(entities$x == levels(entities$x)[2], 1, 0)
    }

    if(is.logical(fix)){
      if (fix == TRUE){
        fix <- sum(entities$empirical * gX)
        cat("asd", fix, "\n")
      }
    } else if (!is.numeric(fix) || length(fix) != 1L){
      stop("invalid value for fix")
    }
  }

  #print("YY: ", str(gY))
  #print("XX: ", str(gX))
  if(classY == "factor" && classX == "numeric")
    stop("categorical 'response' with numeric 'predictor' is not supported")

  ## so what do we do?
  idiv <- NULL
  estimate <- NULL

  if (fix == FALSE && !oneSample){
    message("estimating marginal")
    pH <- estimateMarginal(entities, entities$empirical, mu, classX, gY, gX, oneSample)

    cat("H: estimated marginal:", pH$estimatedMarginal, "\n")
    estimate <- pH$estimate

    if (pH$p$converged == TRUE){
      pA <- estimateMarginal(entities, pH$p$p, estimate, classX, gY, gX, oneSample)
      cat("A: estimated marginal:", pA$estimatedMarginal, "\n")
      idiv <- pA$idiv
    } else{
      idiv <- NA
    }

  } else{
    if (!oneSample) message("fixing marginals")
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
    if (pH$p$converged == TRUE){
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
    if (estimate > mu){
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
      paste("I-test for a two-factor relationship", if (fix == FALSE) "for observational study design" else "for experimental study design")
    } else{
      paste("Two Sample I-test for",
          if (classY == "factor"){
            "rates"
          } else{
            "means"
          },
          if (fix == FALSE) "for observational study design" else "for experimental study design"

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
        gY - gY/gX * fix,
        gX,
        1
      )
      marginal <- sum(entities$empirical * gX)
      estimate <- sum(entities$empirical * (gY - gY/gX * marginal))
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

estimateMarginal <- function(entities, v, mu, classX, gY, gX, oneSample){



  f <- function(x){
    i.test_internal(
      entities = entities,
      v = v,
      mu = mu,
      classX = classX,
      gY = gY, gX = gX,
      oneSample = oneSample,
      fix = x
    )$idiv
  }
  lower <- upper <- NULL
  if (classX == "numeric"){
    xx <- sort(unique(gX))
    i <- 1
    while(!is.finite((fLower <- f(xx[i]))) && i <= length(xx)){

      i <- i + 1
      lower <- xx[i]
    }


    ## what do we do if fLower is infinite?
    ## there may be still a solution
    if (!is.finite(fLower)){
      res <-  i.test_internal(
        entities = entities,
        v = v,
        mu = mu,
        classX = classX,
        gY = gY, gX = gX,
        oneSample = oneSample,
        fix = lower
      )
      res$estimatedMarginal <- NA
      return(res)
    }

    i <- length(xx)

    while(!is.finite((fUpper <- f(xx[i]))) && i > 0){
      i <- i - 1
      upper <- xx[i]
    }
  } else{
    lower <- 0
    upper <- 1
  }

  opt <- optimize(
    f,
    interval = c(lower, upper),
    tol = .Machine$double.eps
  )
  cat(opt$objective, lower, upper, "\n")
  res <-  i.test_internal(
    entities = entities,
    v = v,
    mu = mu,
    classX = classX,
    gY = gY, gX = gX,
    oneSample = oneSample,
    fix = opt$minimum
  )
  res$estimatedMarginal <- opt$minimum
  return(res)

}
