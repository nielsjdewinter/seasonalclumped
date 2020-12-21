#' Function for oxygen isotope based seasonality reconstructions.
#' 
#' Uses records of stable oxygen isotope ratios (d18Oc) through subannually
#' resolved carbonate archives (e.g. mollusk shells or corals) to reconstruct
#' monthly variability in temperature and salinity (assuming user provided
#' variability in d18O composition of the precipitation fluid).
#'
#' @param d18Oc Vector containing subannually resolved d18Oc data 
#' @param ages Vector containing ages for of all samples in years relative to
#' the shell chronology
#' @param SD_d18Oc Error on the d18Oc measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.1 permille).
#' @param d18Ow Vector containing d18O values (in permille VSMOW) of the
#' precipitation fluid used to calculate temperatures. If only a single value is
#' provided, the d18O of the fluid is presumed constant at this value.
#' Default = 0 permille VSMOW.
#' @param d18O_fun String containing the name of the transfer function used to
#' convert temperature and d18Ow to d18Oc data (for example: \code{"KimONeil97"}
#' or \code{"GrossmanKu86"}). Defaults to Kim and O'Neil (1997).
#' @param export Export table summary of result (CSV format)? \code{TRUE/FALSE}
#' @return A data frame containing monthly reconstructions of temperature,
#' d18O of the precipitation fluid and d18Oc.
#' @references Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
#' aragonite: temperature effects, _Chemical Geology_ **1986**, _59.1_, 59-74.
#'     \url{http://dx.doi.org/10.1016/0168-9622(86)90057-6}
#' Kim, S., O'Niel, J.R., Equilibrium and nonequilibrium oxygen
#' isotope effects in synthetic carbonates, _Geochimica et Cosmochimica Acta_
#' **1997**, _61.16_, 3461-3475.
#'     \url{http://dx.doi.org/10.1016/S0016-7037(97)00169-5}
#' Dettman, D.L., Reische, A.K., Lohmann, K.C., Controls on the stable isotope
#' composition of seasonal growth bands in aragonitic fresh-water bivalves
#' (Unionidae), _Geochimica et Cosmochimica Acta_ **1999**, _63.7-8_, 1049-1057.
#'     \url{http://dx.doi.org/10.1016/S0016-7037(99)00020-4}
#' Brand, W.A., Coplen, T.B., Vogl, J., Rosner, M., Prohaska, T., Assessment of
#' international reference materials for isotope-ratio analysis (IUPAC Technical
#' Report), _Pure and Applied Chemistry_ **2014**, _86.3_, 425-467.
#'     \url{http://dx.doi.org/10.1515/pac-2013-1023}
#' de Winter, N. J., Agterhuis, T., Ziegler, M., Optimizing sampling strategies
#' in high-resolution paleoclimate records, _Climate of the Past Discussions_
#' **2020**, 1-52.
#'     \url{http://dx.doi.org/10.5194/cp-2020-118}
#' @examples
#' \donttest{
#'     # find attached dummy data
#'     Case1 <- seasonalclumped::Case1
#'     d18Oc <- Case1[, 29]
#'     ages <- Case1[, 27]
#'     # Run function
#'     monthly <- optimization_seasonality(d18Oc,
#'     ages,
#'     0.1,
#'     "KimONeil97",
#'     FALSE)
#'     }
#' @export
oxygen_isotope_seasonality <- function(d18Oc, # Sub-annually resolved d18Oc data 
    ages, # Vector containing ages for of all samples in years relative to the shell chronology
    SD_d18Oc = 0.1, # Error (1 SD) on d18Oc data 
    d18Ow = 0, # Vector containing d18O values of the precipitation fluid.
    d18O_fun = "KimONeil97",
    export = FALSE # Should the result be exported? 
    ){
    
    # Prepare data
    # Check if data has equal length
    if(length(d18Oc) != length(ages)){
        return("ERROR: Vectors 'd18Oc' and 'ages' should have equal length")
    }
    if(length(SD_d18Oc) == 1){
        SD_d18Oc <- rep(SD_d18Oc, length(d18Oc)) # Duplicate SD of d18Oc error through entire record length if only a single value is given (constant uncertainty)
    }
    if(length(d18Ow) == 1){
        d18Ow <- rep(d18Ow, length(d18Oc)) # Duplicate d18Ow value through entire record length if only a single value is given (constant d18O of precipitation fluid)
    }

    # Group data into monthly bins
    resultmat <- data.frame(d18Oc = d18Oc, # Group d18Oc and age data
        month = ceiling((ages %% 1) * 12) # Use age data to group results into monthly bins
    )

    # Calculate monthly statistics of all d18Oc values
    cat("Grouping d18Oc data into monthly bins: ", "\r")
    d18Oc_monthly <- data.frame(d18Oc_mean = vapply(1:12, function(x) mean(resultmat$d18Oc[which(resultmat$month == x)]), 1),
        d18Oc_median = vapply(1:12, function(x) stats::median(resultmat$d18Oc[which(resultmat$month == x)]), 1),
        d18Oc_SD = vapply(1:12, function(x) stats::sd(resultmat$d18Oc[which(resultmat$month == x)]), 1),
        d18Oc_SDint = sqrt(vapply(1:12, function(x) sum(SD_d18Oc[which(resultmat$month == x)] ^ 2) / length(SD_d18Oc[which(resultmat$month == x)]), 1))
    )
    d18Oc_monthly$d18Oc_SDtot = sqrt(d18Oc_monthly$d18Oc_SDint ^ 2 + d18Oc_monthly$d18Oc_SD ^ 2)
    d18Oc_monthly$d18Oc_SE <- d18Oc_monthly$d18Oc_SDtot / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))

    # Calculate monthly statistics of all d18Ow values
    cat("Grouping d18Oc data into monthly bins: ", "\r")
    d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(d18Ow[which(resultmat$month == x)]), 1),
        d18Ow_median = vapply(1:12, function(x) stats::median(d18Ow[which(resultmat$month == x)]), 1),
        d18Ow_SD = vapply(1:12, function(x) stats::sd(d18Ow[which(resultmat$month == x)]), 1)
    )
    d18Ow_monthly$d18Oc_SE <- d18Ow_monthly$d18Ow_SD / sqrt(vapply(1:12, function(x) length(d18Ow[which(resultmat$month == x)]), 1))

    # Calculate monthly statistics of all temperature reconstructions
    cat("Grouping d18Oc data into monthly bins: ", "\r")
    if(d18O_fun == "KimONeil97"){ # Use transfer function by Kim and O'Neil (1997)
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(18.03 * 10 ^ 3 / (log((resultmat$d18Oc[which(resultmat$month == x)] - (0.97002 * d18Ow[which(resultmat$month == x)] - 29.98)) / 1000 + 1) * 1000 + 32.42) - 273.15), 1),
            T_median = vapply(1:12, function(x) stats::median(18.03 * 10 ^ 3 / (log((resultmat$d18Oc[which(resultmat$month == x)] - (0.97002 * d18Ow[which(resultmat$month == x)] - 29.98)) / 1000 + 1) * 1000 + 32.42) - 273.15), 1),
            T_SD = vapply(1:12, function(x) stats::sd(18.03 * 10 ^ 3 / (log((resultmat$d18Oc[which(resultmat$month == x)] - (0.97002 * d18Ow[which(resultmat$month == x)] - 29.98)) / 1000 + 1) * 1000 + 32.42) - 273.15), 1)
        )
    }else if(d18O_fun == "GrossmanKu86"){ # Use transfer function by Grossman and Ku (1986) adapted by Dettman et al., (1999)
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(20.6 - 4.34 * (resultmat$d18Oc[which(resultmat$month == x)] - d18Ow[which(resultmat$month == x)] - 0.2)), 1),
            T_median = vapply(1:12, function(x) stats::median(20.6 - 4.34 * (resultmat$d18Oc[which(resultmat$month == x)] - d18Ow[which(resultmat$month == x)] - 0.2)), 1),
            T_SD = vapply(1:12, function(x) stats::sd(20.6 - 4.34 * (resultmat$d18Oc[which(resultmat$month == x)] - d18Ow[which(resultmat$month == x)] - 0.2)), 1)
        )
    }
    T_monthly$T_SE = T_monthly$T_SD / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))

    monthly<-cbind(d18Oc_monthly,
        d18Ow_monthly,
        T_monthly)

    # Export results of monthly grouped data
    if(export == TRUE){
        utils::write.csv(monthly, paste("Monthly_results.csv"))
    }

    return(monthly)
    }