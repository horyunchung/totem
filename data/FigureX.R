library(readxl)
library(caper)
tree <- read.nexus("data/mammalia.nexus")
data <- read_xlsx("data/brv12350-sup-0003-tables2.xlsx", skip = 2)

data <- subset(data, subset = SELECT == 1, select = c("Species Nexus", "BMR", "body mass"))

colnames(data) <- c("Species nexus", "BMR", "Mass")

data$BMR <- as.numeric(data$BMR)
data$Mass <- as.numeric(data$Mass)


tab <- toTab(data)

## test pure allowmetric laws with fixed exponents
testAllometric <- function(b){
  ## <x^b>
  meanX <- sum(tab$f * tab$Mass^b)

  ## constraints
  C <- rbind(
    ## hypothesis constraint
    tab$BMR - tab$BMR / tab$Mass^b * meanX,
    ## structural constraint
    ### <x^b>_f = <x^b>_p
    tab$Mass^b,
    ### normalization
    rep(1, nrow(tab))
  )
  ## expected values if the theory is
  ## correct
  targets <- c(0, meanX, 1)
  ## compute hypothesis distribution
  p0 <- iProjector(C, targets, v = tab$f)
  ## compute I-divergence
  idiv <- totem:::computeStatistic(C, targets, tab)
  ## compute p-value
  p.value <- totem:::computePvalue(idiv, sum(tab$Freq), 1, alternative = "two.sided")
  ## compute expected A
  A0 <- sum(p0$p * tab$BMR / tab$Mass^b)
  ## return result
  list(iDivergence = idiv, p.value, A0 = A0 , p0 = p0$p)
}
## test b = 2/3
allometricTwoThird = testAllometric(2/3)

## test b = 3/4
allometricThreeFourth = testAllometric(3/4)

## test the thermodynamic model
testThermodynamic <- function(kPrime = 0.1, f = 0.21){
  ## <x>
  meanX <- sum(tab$f * tab$Mass)
  ## constraints
  C <- rbind(
    ## hypothesis constraint
    tab$BMR - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3) - (tab$BMR  / tab$Mass - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3 - 1)) * meanX,
    ## structural constraint
    ### <x>_f = <x>_p
    tab$Mass,
    ### normalization
    rep(1, nrow(tab))
  )
  targets <- c(0, meanX, 1)
  p0 = iProjector(C, targets, v = tab$f)
  idiv <- totem:::computeStatistic(C, targets, tab)
  p.value <- totem:::computePvalue(idiv, sum(tab$Freq), 1, alternative = "two.sided")
  list(iDivergence = idiv, p.value, K = sum(p0$p * (tab$BMR - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3)) / (f * tab$Mass)), p0 = p0$p)
}



## calculate 95% CI
optF <- function(b){
  res <- testAllometric(b)
  res$iDivergence
}

rootF <- function(b){
  res <- testAllometric(b)
  res$iDivergence - qchisq(0.95, 1, lower.tail = FALSE)
}
opt = optimize(optF, interval = c(0.74, 0.8), tol = .Machine$double.eps)
lower = uniroot(rootF, interval = c(0.74, opt$minimum), tol = .Machine$double.eps)
upper = uniroot(rootF, interval = c(opt$minimum, 0.8), tol = .Machine$double.eps)


## todo: nls, ols, pgls analysis

plot(BMR ~ Mass, data = data, log = 'xy')

myF = function(x, kPrime = 0.1, f = 0.21){
  f * 0.4397837 * x + (1 - f) * kPrime / 0.0201 * x^(2/3)
}
curve(myF, add = TRUE)

myG = function(x){
  3.089168 * x^0.75
}
myH = function(x){
  4.63997 * x^(2/3)
}
curve(myG, add = TRUE, col = 2)
curve(myH, add = TRUE, col = 3)
