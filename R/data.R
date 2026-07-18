#' National Longitudinal Survey of Youth 1997 (27 items)
#'
#' Responses of 3,458 individuals to 27 mixed-type items, with 1.12% missing
#' data. Carried over from the LAWBL package (which `vbpm` supersedes) for
#' empirical illustration of the estimators.
#'
#' @format A list with three components:
#' \describe{
#'   \item{dat}{A `3458 x 27` data frame of item responses (17 polytomous
#'     items, 10 continuous), containing `NA`s.}
#'   \item{Q}{A `27 x 3` initial design matrix with two to three specified
#'     loadings per factor, in the `-1/0/1` coding used by [vbfa()].}
#'   \item{cati}{Indices of the 17 polytomous items.}
#' }
#'
#' @source Bureau of Labor Statistics, U.S. Department of Labor. National
#'   Longitudinal Survey of Youth 1997 cohort. Previously distributed in the
#'   LAWBL package (Chen, 2022).
#'
#' @examples
#' data(nlsy27)
#' dim(nlsy27$dat)
#' nlsy27$Q
"nlsy27"
