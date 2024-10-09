#' Basic metabolic rate/body mass data
#'
#' Data from Genoud et al. (2018) for the basic metabolic rate and body mass of
#' 549 mammals.
#'
#' @docType data
#'
#' @usage data(bmrMass)
#'
#' @format A data.frame with 549 observations of 2 variables (BMR, and Mass)
#'
#' @keywords dataset
#'
#' @references Genoud et al. (2018) Comparative analyses of basal rate of metabolism in
#' mammals: data selection does matter. Biological Reviews 93 (1), 404–438.
#'
#' @examples
#' data(bmrMass)
#' ## testing Kleiber's law
#' i.test(BMR ~ I(Mass^(3/4)), data = bmrMass)
"bmrMass"
