#' Virtual dataset Case 29
#'
#' A dataset containing Ages (\code{Tnew}), depth values (\code{D}), stable
#' oxygen isotope values (\eqn{\delta^(18)}{δ18}O) and clumped isotope values
#' \eqn{\Delta_{47}}{Δ47} of a simulated carbonate record based on environmental
#' parameters following Case 29 and emplying a sampling resolution of
#' \code{0.1 mm}, \code{0.2 mm}, \code{0.45 mm}, \code{0.75 mm}, \code{1.55 mm}
#' and \code{3.25 mm}.
#' 
#' Case 29 describes an ideal temperature sinusoid without distortion by changes
#' in growth rate, but with seasonal variability in \eqn{\delta^(18)O_w}{δ18Ow}
#' with a minimum in spring season and a tiny (1 days SD) error on age model
#' 
#' Generated using the code in "Generate_Case29.r" in \code{data–raw}
#'
#' @format A data frame with 1200 rows and 30 variables:
#' \describe{
#'   \item{SR_ 0.1}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.1 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_ 0.2}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.2 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_ 0.45}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.45 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_ 0.75}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.75 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_ 1.55}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 1.55 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_ 3.25}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 3.25 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   ...
#' }
#' @source See code to generate data in \code{data–raw}
#' Details on how these example cases are defined is provided in:
#' 
#' de Winter, N. J., Agterhuis, T., Ziegler, M., Optimizing sampling strategies
#' in high–resolution paleoclimate records, _Climate of the Past Discussions_
#' **2020**, 1–52.
#'     \url{https://doi.org/fpc4}
"Case29"