#' Virtual dataset Case 16
#'
#' A dataset containing Ages (/code{Tnew}), depth values (/code{D}), stable
#' oxygen isotope values (/eqn{\delta^(18)}{δ18}O) and clumped isotope values
#' \eqn{\Delta_{47}}{Δ47} of a simulated carbonate record based on environmental
#' parameters following Case 16 and emplying a sampling resolution of
#' /code{0.1 mm}, /code{0.2 mm}, /code{0.45 mm}, /code{0.75 mm}, /code{1.55 mm}
#' and /code{3.25 mm}.
#' 
#' Case simulating a lagoonal ecosystem - Seasonal change in growth rate with
#' slower growth in summer and linear growth decrease. Pulse of increased
#' /eqn{\delta^(18)O_w}{δ18Ow} in summer (evaporation).
#' 
#' Generated using the code in "Generate_Case16.r" in \code{data-raw}
#'
#' @format A data frame with 1200 rows and 30 variables:
#' \describe{
#'   \item{SR_0.1}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.1 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_0.2}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.2 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_0.45}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.45 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_0.75}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 0.75 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_1.55}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 1.55 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   \item{SR_3.25}{Empty column denoting the start of the record sampled at a
#'                  sampling resolution of 3.25 mm}
#'   \item{Tnew}{Age, in years relative to the start of the record}
#'   \item{D}{Depth, in mm along the virtual record}
#'   \item{d18Oc}{stable oxygen isotope value, in permille VPDB}
#'   \item{D47}{clumped isotope value, in permille}
#'   ...
#' }
#' @source See code to generate data in \code{data-raw}
#' Details on how these example cases are defined is provided in:
#' 
#' de Winter, N.J., Agterhuis, T.A., Ziegler, M., Optimizing sampling strategies
#' in high-resolution paleoclimate records, _Climate of the Past Discussions_ 
#' **2020**, _cp-2020-118_.
#' \url{http://dx.doi.org/10.5194/cp-2020-118}
"Case16"