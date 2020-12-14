#' Function for sample size optimization-based clumped isotope seasonality
#' reconstruction.
#' 
#' Combines records of stable oxygen isotope ratios (d18Oc) and
#' clumped isotope ratios (D47) through sub-annually resolved carbonate archives
#' (e.g. mollusk shells or corals) to reconstruct monthly variability in
#' temperature and salinity (through the d18O composition of the precipitation
#' fluid), using the sample size optimization method detailed in de Winter et
#' al., 2020 (Climate of the Past).
#'
#' @param d18Oc Vector containing sub-annually resolved d18Oc data 
#' @param D47 Vector containing sub-annually resolved D47 data
#' @param ages Vector containing ages for of all samples in years relative to
#' the shell chronology
#' @param SD_d18OC Error on the d18Oc measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.1 permille).
#' @param SD_D47 Error on the D47 measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.04 permille).
#' @param N Number of datapoints for Monte Carlo simulation (defailts to 1000)
#' @param p Threshold value for the p-value of separating summer from winter
#' reconstructions. Defaults to 0.05 (95% confidence level)
#' @param d18O_fun String containing the name of the transfer function used to
#' convert temperature and d18Ow to d18Oc data (for example: \code{"KimONeil97"}
#' or \code{"GrossmanKu86"}). Defaults to Kim and O'Neil (1997).
#' @param D47_fun String containing the name of the transfer function used to
#' convert temperature to D47 data (for example: \code{"Bernasconi18"} or
#' \code{"Jautzy20"}). Defaults to Bernasconi et al., 2018).
#' @param export_raw Export tables containing all raw model
#' results before being merged into tidy tables? \code{TRUE/FALSE}
#' @return A data frame containing monthly reconstructions of D47, temperature,
#' d18O of the precipitation fluid and d18Oc.
#' @references package dependencies: tidyverse, TTR
#' Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
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
#' Kele, S., Breitenbach, S. F., Capezzuoli, E., Meckler, A. N., Ziegler, M.,
#' Millan, I. M., Kluge, T., Deák, J., Hanselmann, K. and John, C. M.,
#' Temperature dependence of oxygen-and clumped isotope fractionation in
#' carbonates: a study of travertines and tufas in the 6–95 C temperature range,
#' _Geochimica et Cosmochimica Acta_ **2015**, 168, 172–192.
#'     \url{http://dx.doi.org/10.1016/j.gca.2015.06.032}
#' Bernasconi, S.M., Müller, I.A., Bergmann, K.D., Breitenbach, S.F., Fernandez,
#' A., Hodell, D.A., Jaggi, M., Meckler, A.N., Millan, I. and Ziegler, M.,
#' Reducing uncertainties in carbonate clumped isotope analysis through
#' consistent carbonate‐based standardization. _Geochemistry, Geophysics,
#' Geosystems_ **2018**, 19-9, 2895-2914.
#'     \url{http://dx.doi.org/10.1029/2017GC007385}
#' Petersen, S. V., Defliese, W. F., Saenger, C., Daëron, M., Huntington, K. W.,
#' John, C. M., Kelson, J. R., Bernasconi, S. M., Colman, A. S., Kluge, T.,
#' Olack, G. A., Schauer, A. J., Bajnai, D., Bonifacie, M., Breitenbach, S. F.
#' M., Fiebig, J., Fernandez, A. B., Henkes, G. A., Hodell, D., Katz, A., Kele,
#' S., Lohmann, K. C., Passey, B. H., Peral, M. Y., Petrizzo, D. A., Rosenheim,
#' B. E., Tripati, A., Venturelli, R., Young, E. D. and Winkelstern, I. Z.,
#' Effects of Improved 17O Correction on Interlaboratory Agreement in Clumped
#' Isotope Calibrations, Estimates of Mineral-Specific Offsets, and Temperature
#' Dependence of Acid Digestion Fractionation, _Geochemistry, Geophysics,
#' Geosystems_ **2019*, 20-7, 3495–3519.
#'     \url{http://dx.doi.org/10.1029/2018GC008127}
#' Jautzy, J. J., Savard, M. M., Dhillon, R. S., Bernasconi, S. M. and Smirnoff,
#' A., Clumped isotope temperature calibration for calcite: Bridging theory and
#' experimentation, _Geochemical Perspectives Letters_ **2020**, 14, 36–41.
#'     \url{http://dx.doi.org/10.7185/geochemlet.2021}
#' de Winter, N. J., Agterhuis, T., Ziegler, M., Optimizing sampling strategies
#' in high-resolution paleoclimate records, _Climate of the Past Discussions_
#' **2020**, 1-52.
#'     \url{http://dx.doi.org/10.5194/cp-2020-118}
#' @examples
#' \donttest{
#'     # find attached dummy data
#'     Case1 <- seasonalclumped:::Case1
#'     d18Oc <- Case1[, 29]
#'     D47 <- Case1[, 30]
#'     ages <- Case1[, 27]
#'     # Run function
#'     monthly <- optimization_seasonality(d18Oc,
#'     D47,
#'     ages,
#'     0.1,
#'     0.04,
#'     1000,
#'     0.05,
#'     "KimONeil97",
#'     "Bernasconi18",
#'     FALSE}
#' @export
optimization_seasonality <- function(d18Oc, # Sub-annually resolved d18Oc data 
    D47, # Sub-annually resolved D47 data
    ages, # Vector containing ages for of all samples in years relative to the shell chronology
    SD_d18Oc = 0.1, # Error (1 SD) on d18Oc data 
    SD_D47 = 0.04, # Error (1 SD) on D47 data 
    N = 1000, # Number of Monte Carlo simulations for optimization
    p = 0.05, # p-value threshold for considering successful separation of seasons
    d18O_fun = "KimONeil97",
    D47_fun = "Bernasconi18",
    export_raw = FALSE # Should raw data of successful individual simulations be exported (WARNING: Files can get large!)
    ){
    
    # Prepare data
    # Check if data has equal length
    if(length(unique(c(length(d18Oc), length(D47), length(ages)))) > 1){
        return("ERROR: Vectors 'd18Oc', 'D47' and 'ages' should have equal length")
    }
    Popt <- matrix(NA, ncol = 5, nrow = N * length(d18Oc)) # Create matrix with maximum length
    colnames(Popt) <- c("Samplesize",
        "dOsum",
        "dOwin",
        "Dsum",
        "Dwin") # Create template for storing simulation results
    if(length(SD_d18Oc) == 1){
        SD_d18Oc <- rep(SD_d18Oc, length(d18Oc)) # Duplicate SD of d18Oc error through entire record length if only a single value is given (constant uncertainty)
    }
    if(length(SD_D47) == 1){
        SD_D47 <- rep(SD_D47, length(D47)) # Duplicate SD of D47 error through entire record length if only a single value is given (constant uncertainty)
    }

    d18Omat <- as.data.frame(matrix(rep(rnorm(length(d18Oc), d18Oc, SD_d18Oc), N), ncol = N)) # Randomly resample d18O data using measurement uncertainty
    colnames(d18Omat) <- paste("Sim", seq(1, N, 1), sep = "")
    D47mat <- as.data.frame(matrix(rep(rnorm(length(D47), D47, SD_D47), N), ncol = N)) # Randomly resample D47 data using measurement uncertainty
    colnames(D47mat) <- paste("Sim", seq(1, N, 1), sep = "")

    win <- seq(1, length(d18Oc), 1) # Create vector of sample size windows
    row = 1
    
    # MONTE CARLO SIMULATION
    for(i in 1:N){ # Loop through all Monte Carlo simulations
        # Progress
        cat(paste("Sample size optimization Monte Carlo Iteration: ", i),"\r")
        flush.console()
        
        # Expanding window seasonality and T-test
        # Keep record of summer and winter values for successful sample windows
        # Isolate simulated d18O and D47 data
        X <- cbind(d18Omat[, i], D47mat[, i])
        # Sort by d18O (summer side of the record)
        X <- X[order(X[, 1]), ]
        dsum <- TTR::runMean(X[, 1] ,  1, cumulative = TRUE) # Calculate average summer d18O for progressively large sample size windows
        Dsum <- TTR::runMean(X[, 2], 1, cumulative = TRUE) # Calculate average summer D47 for progressively large sample size windows
        Dsumsd <- TTR::runSD(X[, 2], 1, cumulative = TRUE) # Calculate standard deviation within summer D47 values for progressively large sample size windows
        # Inverse sort by d18O (winter side of the record)
        X <- X[order(X[, 1], decreasing = TRUE), ]
        dwin <- TTR::runMean(X[, 1], 1, cumulative = TRUE) # Calculate average winter d18O for progressively large sample size windows
        Dwin <- TTR::runMean(X[, 2], 1, cumulative = TRUE) # Calculate average winter D47 for progressively large sample size windows
        Dwinsd <- TTR::runSD(X[, 2], 1, cumulative = TRUE) # Calculate standard deviation within winter D47 values for progressively large sample size windows
        # Calculate statistics of summer and winter separation for each sampling window
        SDpool <- sqrt((Dsumsd ^ 2 + Dwinsd ^ 2) / 2) # Calculate pooled standard deviation for each window
        T <- (Dsum - Dwin) / (SDpool * sqrt(2 / win)) # Calculate two-sample T-value for each window (equal sample size, equal variance)
        Pval <- pt(T, win - 1) # Calculate p-value for each window
        res<-cbind(win[which(Pval < p)],
            dsum[which(Pval < p)],
            dwin[which(Pval < p)],
            Dsum[which(Pval < p)],
            Dwin[which(Pval < p)]) # Combine results
        Popt[row:(row + length(which(Pval < p)) - 1), ] <- res # Add results to running matrix of optimized simulations
        row <- row + length(which(Pval < p)) - 1 # Increment row number for efficient data storage
    }

    # POST PROCESSING
    Popt <- as.data.frame(Popt[complete.cases(Popt), ]) # Remove rows with NAs
    Popt <- Popt[-which(Popt$Samplesize > (length(d18Oc) / 2)), ] # Remove simulations taking more than half of the samples (summer and winter samples overlap)

    # Add temperature calculations of optimal runs
    if(D47_fun == "Bernasconi18"){
        Popt$Tsum <- sqrt((0.0449 * 10 ^ 6) / (Popt$Dsum - 0.167)) - 273.15 # Calculate summer and winter temperatures for each successful simulation according to Kele et al., 2015 modified by Bernasconi et al., 2018
        Popt$Twin <- sqrt((0.0449 * 10 ^ 6) / (Popt$Dwin - 0.167)) - 273.15
    }else if(D47_fun == "Jautzy20"){
        Popt$Tsum <- sqrt((0.0433 * 10 ^ 6) / (Popt$Dsum - 0.119 + 0.066)) - 273.15 # Calculate summer and inter temperatures for each successful simulation according to Jautzy et al., 2020 brought into 25 degrees CDES reference frame using 70-25 acid fractionation factor by Petersen et al., 2019
        Popt$Twin <- sqrt((0.0433 * 10 ^ 6) / (Popt$Dwin - 0.119 + 0.066)) - 273.15
    }else{
        return("ERROR: Supplied D47 transfer function is not recognized")
    }
    # Add seawater d18O calculations of optimal runs
    if(d18O_fun == "KimONeil97"){
        Popt$dOwsum <- (Popt$dOsum - (exp(((18.03 * 10 ^ 3) / (Popt$Tsum + 273.15) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92 # Calculate d18O of the precipitation fluid (dOw) for summer and winter simulations using Kim and O'Neil, 1997 with conversion to PDB (following Brand et al., 2014)
        Popt$dOwwin <- (Popt$dOwin - (exp(((18.03 * 10 ^ 3) / (Popt$Twin + 273.15) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92
    }else if(d18O_fun == "GrossmanKu86"){
        Popt$dOwsum <- (Popt$Tsum - 20.6) / 4.34 + Popt$dOsum - 0.2 # Calculate d18O of the precipitation fluid (dOw) for summer and winter simulations using Grossmann and Ku (1986) modified by Dettmann et al. (1999)
        Popt$dOwwin <- (Popt$Twin - 20.6) / 4.34 + Popt$dOwin - 0.2
    }else{
        return("ERROR: Supplied d18Oc transfer function is not recognized")
    }
    # Add slopes and intercepts for D47-d18O conversion
    Popt$D_dO_slope <- (Popt$Dsum - Popt$Dwin) / (Popt$dOsum - Popt$dOwin)
    Popt$D_dO_int <- ((Popt$Dsum + Popt$Dwin) / 2) - Popt$D_dO_slope * ((Popt$dOsum + Popt$dOwin) / 2)

    # OPTIONAL: Export results of optimized sample sizes
    if(export_raw == TRUE){
        write.csv(Popt, "Optimized_simulations.csv")
    }

    # Use Popt matrix to model monthly D47, T and d18Osw from d18O data
    resultmat <- data.frame(d18Oc = d18Oc, # Group d18Oc and age data
        month = ceiling((ages %% 1) * 12) # Use age data to group results into monthly bins
    )

    # Calculate monthly statistics of all d18Oc values
    cat("Grouping d18Oc data into monthly bins: ", "\r")
    d18Oc_monthly <- data.frame(d18Oc_mean = vapply(1:12, function(x) mean(resultmat$d18Oc[which(resultmat$month == x)]), 1),
        d18Oc_median = vapply(1:12, function(x) median(resultmat$d18Oc[which(resultmat$month == x)]), 1),
        d18Oc_SD = vapply(1:12, function(x) sd(resultmat$d18Oc[which(resultmat$month == x)]), 1)
    )
    d18Oc_monthly$d18Oc_SE <- d18Oc_monthly$d18Oc_SD / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))

    # Calculate monthly statistics of all D47 values using the d18Oc measurements and the D47-d18Oc slopes of all successful simulations
    cat("Grouping D47 data into monthly bins: ", "\r")
    D47_monthly <- data.frame(D47_mean = vapply(1:12, function(x) mean(outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int), 1),
        D47_median = vapply(1:12, function(x) median(outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int), 1),
        D47_SD = vapply(1:12, function(x) sd(outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int), 1)
    )
    D47_monthly$D47_SE = D47_monthly$D47_SD / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))
    # Repeat for monthly temperature reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping Temperature data into monthly bins: ", "\r")
    if(D47_fun == "Bernasconi18"){
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15), 1),
            T_median = vapply(1:12, function(x) median(sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15), 1),
            T_SD = vapply(1:12, function(x) sd(sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15), 1)
        )
    }else if(D47_fun == "Jautzy20"){
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15), 1),
            T_median = vapply(1:12, function(x) median(sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15), 1),
            T_SD = vapply(1:12, function(x) sd(sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15), 1)
        )
    }
    T_monthly$T_SE = T_monthly$T_SD / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))
    # Repeat for monthly d18Ow reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping d18Ow data into monthly bins: ", "\r")
    if(d18O_fun == "KimONeil97"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1),
                d18Ow_median = vapply(1:12, function(x) median((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1),
                d18Ow_SD = vapply(1:12, function(x) sd((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1)
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1),
                d18Ow_median = vapply(1:12, function(x) median((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1),
                d18Ow_SD = vapply(1:12, function(x) sd((resultmat$d18Oc[which(resultmat$month == x)] - (exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 32.42) / 1000) - 1) * 1000) * 1.03092 + 30.92), 1)
            )
        }
    }else if(d18O_fun == "GrossmanKu86"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_median = vapply(1:12, function(x) median(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_SD = vapply(1:12, function(x) sd(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0449 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2), 1)
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_median = vapply(1:12, function(x) median(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1),
                d18Ow_SD = vapply(1:12, function(x) sd(resultmat$d18Oc[which(resultmat$month == x)] + ((sqrt((0.0433 * 10 ^ 6) / (outer(resultmat$d18Oc[which(resultmat$month == x)], Popt$D_dO_slope) + Popt$D_dO_int - 0.119 + 0.066)) - 273.15) - 20.6) / 4.34 - 0.2), 1)
            )
        }
    }
    d18Ow_monthly$d18Ow_SE = d18Ow_monthly$d18Ow_SD / sqrt(vapply(1:12, function(x) length(resultmat$d18Oc[which(resultmat$month == x)]), 1))

    monthly<-cbind(d18Oc_monthly,
        D47_monthly,
        T_monthly,
        d18Ow_monthly)
    # Export results of monthly grouped data
    write.csv(monthly, paste("Monthly_results.csv"))

    return(monthly)
    }