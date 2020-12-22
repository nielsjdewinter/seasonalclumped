#' Function for monthly binning based clumped isotope seasonality
#' reconstruction.
#' 
#' Combines records of stable oxygen isotope ratios (d18Oc) and
#' clumped isotope ratios (D47) through sub annually resolved carbonate archives
#' (e.g. mollusk shells or corals) to reconstruct monthly variability in
#' temperature and salinity (through the d18O composition of the precipitation
#' fluid), using the monthly binning method detailed in de Winter et
#' al., 2020 (Climate of the Past).
#'
#' @param d18Oc Vector containing subannually resolved d18Oc data 
#' @param D47 Vector containing subannually resolved D47 data
#' @param ages Vector containing ages for of all samples in years relative to
#' the shell chronology
#' @param SD_d18Oc Error on the d18Oc measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.1 permille).
#' @param SD_D47 Error on the D47 measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.04 permille).
#' @param N Number of datapoints for Monte Carlo simulation (defailts to 1000)
#' @param binsize Size of the bins in which records are subdivided. \code{month}
#' and \code{season} (period of three months) is currently supported.
#' @param d18O_fun String containing the name of the transfer function used to
#' convert temperature and d18Ow to d18Oc data (for example: \code{"KimONeil97"}
#' or \code{"GrossmanKu86"}). Defaults to Kim and O'Neil (1997).
#' @param D47_fun String containing the name of the transfer function used to
#' convert temperature to D47 data (for example: \code{"Bernasconi18"} or
#' \code{"Jautzy20"}). Defaults to Bernasconi et al., 2018).
#' @param export Export table summary of result (CSV format)? \code{TRUE/FALSE}
#' @return A data frame containing monthly reconstructions of D47, temperature,
#' d18O of the precipitation fluid and d18Oc.
#' @references Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
#' aragonite: temperature effects, _Chemical Geology_ **1986**, _59.1_, 59–74.
#'     \url{http://dx.doi.org/10.1016/0168-9622(86)90057-6}
#'
#' Kim, S., O'Niel, J.R., Equilibrium and nonequilibrium oxygen
#' isotope effects in synthetic carbonates, _Geochimica et Cosmochimica Acta_
#' **1997**, _61.16_, 3461–3475.
#'     \url{http://dx.doi.org/10.1016/S0016-7037(97)00169-5}
#'
#' Dettman, D.L., Reische, A.K., Lohmann, K.C., Controls on the stable isotope
#' composition of seasonal growth bands in aragonitic fresh–water bivalves
#' (Unionidae), _Geochimica et Cosmochimica Acta_ **1999**, _63.7–8_, 1049–1057.
#'     \url{http://dx.doi.org/10.1016/S0016-7037(99)00020-4}
#'
#' Brand, W.A., Coplen, T.B., Vogl, J., Rosner, M., Prohaska, T., Assessment of
#' international reference materials for isotope–ratio analysis (IUPAC Technical
#' Report), _Pure and Applied Chemistry_ **2014**, _86.3_, 425–467.
#'     \url{http://dx.doi.org/10.1515/pac-2013-1023}
#'
#' Kele, S., Breitenbach, S. F., Capezzuoli, E., Meckler, A. N., Ziegler, M.,
#' Millan, I. M., Kluge, T., Deák, J., Hanselmann, K. and John, C. M.,
#' Temperature dependence of oxygen– and clumped isotope fractionation in
#' carbonates: a study of travertines and tufas in the 6–95 C temperature range,
#' _Geochimica et Cosmochimica Acta_ **2015**, 168, 172–192.
#'     \url{http://dx.doi.org/10.1016/j.gca.2015.06.032}
#'
#' Bernasconi, S.M., Müller, I.A., Bergmann, K.D., Breitenbach, S.F., Fernandez,
#' A., Hodell, D.A., Jaggi, M., Meckler, A.N., Millan, I. and Ziegler, M.,
#' Reducing uncertainties in carbonate–clumped isotope analysis through
#' consistent carbonate based standardization. _Geochemistry, Geophysics,
#' Geosystems_ **2018**, 19–9, 2895–2914.
#'     \url{http://dx.doi.org/10.1029/2017GC007385}
#'
#' Petersen, S. V., Defliese, W. F., Saenger, C., Daëron, M., Huntington, K. W.,
#' John, C. M., Kelson, J. R., Bernasconi, S. M., Colman, A. S., Kluge, T.,
#' Olack, G. A., Schauer, A. J., Bajnai, D., Bonifacie, M., Breitenbach, S. F.
#' M., Fiebig, J., Fernandez, A. B., Henkes, G. A., Hodell, D., Katz, A., Kele,
#' S., Lohmann, K. C., Passey, B. H., Peral, M. Y., Petrizzo, D. A., Rosenheim,
#' B. E., Tripati, A., Venturelli, R., Young, E. D. and Winkelstern, I. Z.,
#' Effects of Improved 17O Correction on Interlaboratory Agreement in Clumped
#' Isotope Calibrations, Estimates of Mineral–Specific Offsets, and Temperature
#' Dependence of Acid Digestion Fractionation, _Geochemistry, Geophysics,
#' Geosystems_ **2019*, 20–7, 3495–3519.
#'     \url{http://dx.doi.org/10.1029/2018GC008127}
#'
#' Jautzy, J. J., Savard, M. M., Dhillon, R. S., Bernasconi, S. M. and Smirnoff,
#' A., Clumped isotope temperature calibration for calcite: Bridging theory and
#' experimentation, _Geochemical Perspectives Letters_ **2020**, 14, 36–41.
#'     \url{http://dx.doi.org/10.7185/geochemlet.2021}
#'
#' de Winter, N. J., Agterhuis, T., Ziegler, M., Optimizing sampling strategies
#' in high–resolution paleoclimate records, _Climate of the Past Discussions_
#' **2020**, 1–52.
#'     \url{http://dx.doi.org/10.5194/cp-2020-118}
#' @examples
#' \donttest{
#'     # find attached dummy data
#'     Case1 <- seasonalclumped::Case1
#'     d18Oc <- Case1[, 29]
#'     d18Oc <- d18Oc[-which(is.na(d18Oc))]
#'     D47 <- Case1[, 30]
#'     D47 <- D47[-which(is.na(D47))]
#'     ages <- Case1[, 27]
#'     ages <- ages[-which(is.na(ages))]
#'     # Run function
#'     binned <- binning_seasonality(d18Oc,
#'     D47,
#'     ages,
#'     0.1,
#'     0.04,
#'     1000,
#'     "month",
#'     "KimONeil97",
#'     "Bernasconi18",
#'     FALSE)
#'     }
#' @export
binning_seasonality <- function(d18Oc, # Subannually resolved d18Oc data 
    D47, # Subannually resolved D47 data
    ages, # Vector containing ages for of all samples in years relative to the shell chronology
    SD_d18Oc = 0.1, # Error (1 SD) on d18Oc data 
    SD_D47 = 0.04, # Error (1 SD) on D47 data 
    N = 1000, # Number of Monte Carlo simulations for binning
    binsize = "month",
    d18O_fun = "KimONeil97",
    D47_fun = "Bernasconi18",
    export = FALSE # Should the result be exported? 
    ){
    
    # Prepare data
    # Check if data has equal length
    if(length(unique(c(length(d18Oc), length(D47), length(ages)))) > 1){
        return("ERROR: Vectors 'd18Oc', 'D47' and 'ages' should have equal length")
    }
    if(length(SD_d18Oc) == 1){
        SD_d18Oc <- rep(SD_d18Oc, length(d18Oc)) # Duplicate SD of d18Oc error through entire record length if only a single value is given (constant uncertainty)
    }
    if(length(SD_D47) == 1){
        SD_D47 <- rep(SD_D47, length(D47)) # Duplicate SD of D47 error through entire record length if only a single value is given (constant uncertainty)
    }

    d18Omat <- as.data.frame(matrix(stats::rnorm(N * length(d18Oc), d18Oc, SD_d18Oc), ncol = N)) # Randomly resample d18O data using measurement uncertainty
    colnames(d18Omat) <- paste("Sim", seq(1, N, 1), sep = "")
    D47mat <- as.data.frame(matrix(stats::rnorm(N * length(D47), D47, SD_D47), ncol = N)) # Randomly resample D47 data using measurement uncertainty
    colnames(D47mat) <- paste("Sim", seq(1, N, 1), sep = "")

    # Use age model to bin D47, T and d18Osw
    if(binsize == "month"){
        resultmat <- data.frame(d18Oc = d18Oc, # Group d18Oc, D47 and age data
            D47 = D47,
            bin = ceiling((ages %% 1) * 12) # Use age data to group results into binned bins
        )
        bins = 1:12
    }else if(binsize == "season"){
            resultmat <- data.frame(d18Oc = d18Oc, # Group d18Oc, D47 and age data
            D47 = D47,
            bin = ceiling((ages %% 1) * 4) # Use age data to group results into binned bins
        )
        bins = 1:4
    }

    # Calculate binned statistics of all d18Oc values
    cat("Grouping d18Oc data into binned bins: ", "\r")
    d18Oc_binned <- data.frame(d18Oc_mean = vapply(bins, function(x) mean(as.matrix(d18Omat[which(resultmat$bin == x), ])), 1),
        d18Oc_median = vapply(bins, function(x) stats::median(as.matrix(d18Omat[which(resultmat$bin == x), ])), 1),
        d18Oc_SD = vapply(bins, function(x) stats::sd(as.matrix(d18Omat[which(resultmat$bin == x), ])), 1)
    )
    d18Oc_binned$d18Oc_SE <- d18Oc_binned$d18Oc_SD / sqrt(vapply(bins, function(x) length(resultmat$d18Oc[which(resultmat$bin == x)]), 1))

    # Calculate binned statistics of all D47 values using the d18Oc measurements and the D47-d18Oc slopes of all successful simulations
    cat("Grouping D47 data into binned bins: ", "\r")
    D47_binned <- data.frame(D47_mean = vapply(bins, function(x) mean(as.matrix(D47mat[which(resultmat$bin == x), ])), 1),
        D47_median = vapply(bins, function(x) stats::median(as.matrix(D47mat[which(resultmat$bin == x), ])), 1),
        D47_SD = vapply(bins, function(x) stats::sd(as.matrix(D47mat[which(resultmat$bin == x), ])), 1)
    )
    D47_binned$D47_SE <- D47_binned$D47_SD / sqrt(vapply(bins, function(x) length(resultmat$D47[which(resultmat$bin == x)]), 1))
    
    # Repeat for binned temperature reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping Temperature data into binned bins: ", "\r")
    if(D47_fun == "Bernasconi18"){
        T_binned <- data.frame(T_mean = vapply(bins, function(x) mean(sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15), 1),
            T_median = vapply(bins, function(x) stats::median(sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15), 1),
            T_SD = vapply(bins, function(x) stats::sd(sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15), 1)
        )
    }else if(D47_fun == "Jautzy20"){
        T_binned <- data.frame(T_mean = vapply(bins, function(x) mean(sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15), 1),
            T_median = vapply(bins, function(x) stats::median(sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15), 1),
            T_SD = vapply(bins, function(x) stats::sd(sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15), 1)
        )
    }
    T_binned$T_SE = T_binned$T_SD / sqrt(vapply(bins, function(x) length(resultmat$d18Oc[which(resultmat$bin == x)]), 1))

    # Repeat also for the case where D47 is first binned and then converted to temperatures (not the correct solution)
    if(D47_fun == "Bernasconi18"){
        Trev_binned <- data.frame(T_mean = sqrt((0.0449 * 10 ^ 6) / (D47_binned$D47_mean - 0.167)) - 273.15,
            T_median = sqrt((0.0449 * 10 ^ 6) / (D47_binned$D47_median - 0.167)) - 273.15,
            T_SD = sqrt((0.0449 * 10 ^ 6) / (D47_binned$D47_mean - 0.167)) - 273.15 - (sqrt((0.0449 * 10 ^ 6) / (D47_binned$D47_mean + D47_binned$D47_SD - 0.167)) - 273.15)
        )
    }else if(D47_fun == "Jautzy20"){
        Trev_binned <- data.frame(T_mean = sqrt((0.0433 * 10 ^ 6) / (D47_binned$D47_mean - 0.119 - 0.066)) - 273.15,
            T_median = sqrt((0.0433 * 10 ^ 6) / (D47_binned$D47_median - 0.119 - 0.066)) - 273.15,
            T_SD = sqrt((0.0433 * 10 ^ 6) / (D47_binned$D47_mean - 0.119 - 0.066)) - 273.15 - (sqrt((0.0433 * 10 ^ 6) / (D47_binned$D47_mean + D47_binned$D47_SD - 0.119 - 0.066)) - 273.15)
        )
    }
    Trev_binned$T_SE = Trev_binned$T_SD / sqrt(vapply(bins, function(x) length(resultmat$d18Oc[which(resultmat$bin == x)]), 1))

    # Repeat for binned d18Ow reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping d18Ow data into binned bins: ", "\r")
    if(d18O_fun == "KimONeil97"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_binned <- data.frame(d18Ow_mean = vapply(bins, function(x) mean(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1),
                d18Ow_median = vapply(bins, function(x) stats::median(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1),
                d18Ow_SD = vapply(bins, function(x) stats::sd(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1)
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_binned <- data.frame(d18Ow_mean = vapply(bins, function(x) mean(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1),
                d18Ow_median = vapply(bins, function(x) stats::median(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1),
                d18Ow_SD = vapply(bins, function(x) stats::sd(((as.matrix(d18Omat[which(resultmat$bin == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92), 1)
            )
        }
    }else if(d18O_fun == "GrossmanKu86"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_binned <- data.frame(d18Ow_mean = vapply(bins, function(x) mean(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_median = vapply(bins, function(x) stats::median(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_SD = vapply(bins, function(x) stats::sd(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1)
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_binned <- data.frame(d18Ow_mean = vapply(bins, function(x) mean(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_median = vapply(bins, function(x) stats::median(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_SD = vapply(bins, function(x) stats::sd(as.matrix(d18Omat[which(resultmat$bin == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(D47mat[which(resultmat$bin == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1)
            )
        }
    }
    d18Ow_binned$d18Ow_SE = d18Ow_binned$d18Ow_SD / sqrt(vapply(bins, function(x) length(resultmat$d18Oc[which(resultmat$bin == x)]), 1))

    bin_stats <- cbind(bins,
        d18Oc_binned,
        D47_binned,
        T_binned,
        Trev_binned,
        d18Ow_binned)
    colnames(bin_stats)[1] <- "bin #"
    
    # Export results of binned data (OPTIONAL)
    if(export == TRUE){
        utils::write.csv(bin_stats, paste("binned_results.csv"))
    }

    return(bin_stats)
}