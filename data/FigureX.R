library(readxl)
data <- read_xlsx("data/brv12350-sup-0003-tables2.xlsx", skip = 2)

data <- subset(data, subset = SELECT == 1, select = c("BMR", "body mass"))

colnames(data) <- c("BMR", "Mass")

data$BMR <- as.numeric(data$BMR)
data$Mass <- as.numeric(data$Mass)


tab <- toTab(data)

## test allowmetric laws with fixed exponents

testAllometric <- function(b){
  ## <x^b>
  meanX <- sum(tab$f * tab$Mass^b)

  ## constraints
  C <- rbind(
    ## hypothesis constraint
    tab$BMR - tab$BMR / tab$Mass^b * meanX,
    ## structural constraint
    tab$Mass^b,
    rep(1, nrow(tab))
  )

  targets <- c(0, meanX, 1)
  p0 = iProjector(C, targets, v = tab$f)
  idiv <- totem:::computeStatistic(C, targets, tab)

  p.value <- totem:::computePvalue(idiv, sum(tab$Freq), 1, alternative = "two.sided")

  list(iDivergence = idiv, p.value, A = sum(p0$p * tab$BMR / tab$Mass^b), p0 = p0$p)
}

testThermodynamic <- function(kPrime = 0.1, f = 0.21){
  meanX <- sum(tab$f * tab$Mass)
  ## constraints
  C <- rbind(
    ## hypothesis constraint
    tab$BMR - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3) - (tab$BMR  / tab$Mass - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3 - 1)) * meanX,
    tab$Mass,
    rep(1, nrow(tab))
  )
  targets <- c(0, meanX, 1)
  p0 = iProjector(C, targets, v = tab$f)
  idiv <- totem:::computeStatistic(C, targets, tab)
  p.value <- totem:::computePvalue(idiv, sum(tab$Freq), 1, alternative = "two.sided")
  list(iDivergence = idiv, p.value, K = sum(p0$p * (tab$BMR - (1 - f) * kPrime / 0.0201 * tab$Mass^(2/3)) / (f * tab$Mass)), p0 = p0$p)
}

## pure allometric models
allometricTwoThird = testAllometric(2/3)
allometricThreeFourth = testAllometric(3/4)

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
