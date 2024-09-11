#' @export
i.test <- function(x, ...) UseMethod("i.test")

## want to implement the generic i-test
## 1. x = numeric, y = NULL, mu, one sample t-test
## 2. x = numeric, y = factor, two+ sample t-test/anova
## 3. x = factor, y = NULL, mu, one sample binomial test
## 4. x = factor, y = factor two+ sample barnard test
## 5. x = factro, y = factor two+ sample fisher's exact test
## 6. x = numeric, y = numeric, model test



