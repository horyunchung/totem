library(readxl)
library(caper)
library(totem)
# Prepare data
## load the nexus tree
tree <- read.nexus("data/mammalia.nexus")
## load the data
data <- read_xlsx("data/brv12350-sup-0003-tables2.xlsx", skip = 2, na = "NA")
## only the SELECT subset
data <- subset(data, subset = SELECT == 1, select = c("Species Nexus", "BMR", "body mass"))

## set the column names
colnames(data) <- c("Species", "BMR", "Mass")

## reformat the species names to match the species names in the nexus file
## (substitute ` ` with `_` between genus and species)
data$Species <- data$Species
data$Species <- gsub("\\s", "\\_", data$Species)

# TOTEM analysis of theories

## empirical entity distribution for the attributes
## BMR and Mass
tab <- toTab(subset(data, select = c("Species", "BMR", "Mass")))

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
### test b = 2/3
allometricTwoThird <- testAllometric(2/3)

### test b = 3/4
allometricThreeFourth <- testAllometric(3/4)

### calculate 95% CI
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

allometricOpt <- testAllometric(opt$minimum)



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

thermodynamic <- testThermodynamic()

### sensitivity analysis for kPrime and f?

kPrimes <- seq(0.05, 0.15, by = 0.0005)
fs <- seq(0.01, 0.42, by = 0.001)

sensitivityAnalysis <- sapply(
  fs,
  function(f){
    sapply(
      kPrimes,
      function(kPrime){
        testThermodynamic(kPrime, f)[[2]]
      }
    )
  }
)
image(kPrimes, fs, sensitivityAnalysis, useRaster = TRUE)
contour(kPrimes, fs, sensitivityAnalysis, levels = 0, add =TRUE)
points(0.1, 0.21, pch = 3, col = "black")
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


nls1 <- nls(BMR ~ A * Mass^b, data = data, start = c(A = 4, b = 0.7))
nls0.75 <- nls(BMR ~ A * Mass^0.75, data = data, start = c(A = 4))
nls0.66<- nls(BMR ~ A * Mass^0.66, data = data, start = c(A = 4))

nls0.75.anova <- anova(nls0.75, nls1)
nls0.66.anova <- anova(nls0.66, nls1)
nls1.shapiro <- shapiro.test(residuals(nls1))
nls0.75.shapiro <- shapiro.test(residuals(nls0.75))
nls0.66.shapiro <- shapiro.test(residuals(nls0.66))

nls1.ks <- ks.test(residuals(nls1), pnorm)
nls0.75.ks <- ks.test(residuals(nls0.75), pnorm)
nls0.66.ks <- ks.test(residuals(nls0.66), pnorm)

linm1 <- lm(log10(BMR) ~ log10(Mass), data = data)
linm0.75 = lm(log10(BMR) ~ offset(3/4 * log10(Mass)), data = data)
linm0.66 = lm(log10(BMR) ~ offset(2/3 * log10(Mass)), data = data)

linm0.75.anova <- anova(linm0.75, linm1)
linm0.66.anova <- anova(linm0.66, linm1)
linm1.shapiro <- shapiro.test(residuals(linm1))
linm0.75.shapiro <- shapiro.test(residuals(linm0.75))
linm0.66.shapiro <- shapiro.test(residuals(linm0.66))

linm1.ks <- ks.test(residuals(linm1), pnorm)
linm0.75.ks <- ks.test(residuals(linm0.75), pnorm)
linm0.66.ks <- ks.test(residuals(linm0.66), pnorm)

## pgls analysis
dataSpeciesComperativeData <- comparative.data(
  data = data.frame(
    data,
    log10BMR = log10(data$BMR),
    log10Mass = log10(data$Mass),
    log10BMR0.75 = log10(data$BMR) - 3/4 * log10(data$Mass),
    log10BMR0.66 = log10(data$BMR) - 2/3 * log10(data$Mass)
  ),
  phy = tree, names.col = "Species"
)
## the pgls function does not allow to fit an intercept only model
## with a specified power law exponent, so we do it by hand
pglsAnalysis1 <- pgls(log10BMR ~ log10Mass, data = dataSpeciesComperativeData, lambda = "ML")
pglsAnalysis0.75 <- pgls(log10BMR0.75 ~ 1, data = dataSpeciesComperativeData, lambda = pglsAnalysis1$param["lambda"])
pglsAnalysis0.66 <- pgls(log10BMR0.66 ~ 1, data = dataSpeciesComperativeData, lambda = pglsAnalysis1$param["lambda"])




## moreover, we cannot use the anova function, because we fitted different models
## so, we calculate the F-statistic by hand
residuals1 = residuals(pglsAnalysis1, phylo = TRUE)
residuals0.75 = residuals(pglsAnalysis0.75, phylo = TRUE)
residuals0.66 = residuals(pglsAnalysis0.66, phylo = TRUE)

pgls1.shapiro <- shapiro.test(residuals1)
pgls0.75.shapiro <- shapiro.test(residuals0.75)
pgls0.66.shapiro <- shapiro.test(residuals0.66)


rss1 <- sum(residuals1^2)
rss0.75 <- sum(residuals0.75^2)
rss0.66 <- sum(residuals0.66^2)

F0.75 <- (rss0.75 - rss1) / (rss1 / 547)
F0.66 <- (rss0.66 - rss1) / (rss1 / 547)


results <- data.frame(
  analysis = c(
    "NLS b = 0.66",
    "NLS b = 0.75",
    "OLS b = 0.66",
    "OLS b = 0.75",
    "pgls b = 0.66",
    "pgls b = 0.75",
    "I-test b = 0.66",
    "I-test b = 0.75"
  ),
  H0.A = c(
    coefficients(nls0.66)["A"],
    coefficients(nls0.75)["A"],
    10^coefficients(linm0.66)[1],
    10^coefficients(linm0.75)[1],
    10^coefficients(pglsAnalysis0.66)[1],
    10^coefficients(pglsAnalysis0.75)[1],
    allometricTwoThird$A,
    allometricThreeFourth$A
  ),
  H0.logA = c(
    log10(coefficients(nls0.66)["A"]),
    log10(coefficients(nls0.75)["A"]),
    coefficients(linm0.66)[1],
    coefficients(linm0.75)[1],
    coefficients(pglsAnalysis0.66)[1],
    coefficients(pglsAnalysis0.75)[1],
    log10(allometricTwoThird$A),
    log10(allometricThreeFourth$A)
  ),
  H0.b = rep(c("2/3", "3/4"), 4),
  H0.bias = c(
    mean(residuals(nls0.66)),
    mean(residuals(nls0.75)),
    mean(data$BMR - 10^(fitted(linm0.66))),
    mean(data$BMR - 10^(fitted(linm0.75))),
    mean(data$BMR - 10^(fitted(pglsAnalysis0.66))),
    mean(data$BMR - 10^(fitted(pglsAnalysis0.75))),
    sum(allometricThreeFourth$p0 * (tab$BMR - tab$BMR / tab$Mass^0.75 * sum(allometricThreeFourth$p0 * tab$Mass^0.75))),
    sum(allometricTwoThird$p0 * (tab$BMR - tab$BMR / tab$Mass^(2/3) * sum(allometricTwoThird$p0 * tab$Mass^(2/3))))
  ),
  H0.normality = c(
    nls0.66.shapiro$p.value,
    nls0.75.shapiro$p.value,
    linm0.66.shapiro$p.value,
    linm0.75.shapiro$p.value,
    pgls0.66.shapiro$p.value,
    pgls0.75.shapiro$p.value,
    NA,
    NA
  ),
  H1.A = rep(
    c(coefficients(nls1)["A"], 10^coefficients(linm1)[1], 10^coefficients(pglsAnalysis1)[1], NA), each = 2
  ),
  H1.logA = rep(
    c(NA, coefficients(linm1)[1], coefficients(pglsAnalysis1)[1], NA), each = 2
  ),
  H1.b = rep(
    c(coefficients(nls1)["b"], coefficients(linm1)[2], coefficients(pglsAnalysis1)[2], NA), each = 2
    ),
  H1.bias = rep(
    c(
      mean(residuals(nls1)),
      mean(data$BMR - 10^fitted(linm1)),
      mean(data$BMR - 10^fitted(pglsAnalysis1)),
      NA),
    each = 2
  ),
  H1.normality = rep(
    c(nls1.shapiro$p.value, linm1.shapiro$p.value, pgls1.shapiro$p.value, NA), each = 2
  ),
  statistic = c(
    nls0.66.anova$`F value`[2],
    nls0.75.anova$`F value`[2],
    linm0.66.anova$F[2],
    linm0.75.anova$F[2],
    F0.66,
    F0.75,
    allometricTwoThird[[2]]$statistic,
    allometricThreeFourth[[2]]$statistic
  ),
  p.value = c(
    nls0.66.anova$`Pr(>F)`[2],
    nls0.75.anova$`Pr(>F)`[2],
    linm0.66.anova$`Pr(>F)`[2],
    linm0.75.anova$`Pr(>F)`[2],
    pf(F0.66, 1, 547, lower.tail = FALSE),
    pf(F0.75, 1, 547, lower.tail = FALSE),
    allometricTwoThird[[2]]$p.value,
    allometricThreeFourth[[2]]$p.value
  )

)



