#' Function for sample size optimization based clumped isotope seasonality
#' reconstruction.
#' 
#' Combines records of stable oxygen isotope ratios (\eqn{\delta^{18}O_{w}}{δ18Ow})
#' and clumped isotope ratios (D47) through subannually resolved carbonate
#' archives (e.g. mollusk shells or corals) to reconstruct monthly variability
#' in temperature and salinity (through the \eqn{\delta^{18}O}{δ18O} composition
#' of the precipitation fluid), using the moving average method detailed in de
#' Winter et al., 2020 (Climate of the Past).
#'
#' @param d18Oc Vector containing subannually resolved
#' \eqn{\delta^{18}O_{w}}{δ18Ow} data 
#' @param D47 Vector containing subannually resolved D47 data
#' @param ages Vector containing ages for of all samples in years relative to
#' the shell chronology
#' @param SD_d18Oc Error on the \eqn{\delta^{18}O_{w}}{δ18Ow} measurements.
#' Either a single value (constant uncertainty) or a vector of length equal to
#' the period in SST data (365 days by default) containing information about the
#' error of each datapoint (1 standard deviation; default = 0.1 permille).
#' @param SD_D47 Error on the D47 measurements. Either a single value
#' (constant uncertainty) or a vector of length equal to the period in SST data 
#' (365 days by default) containing information about the error of each 
#' datapoint (1 standard deviation; default = 0.04 permille).
#' @param window Either supply the size of the window used for moving average
#' calculation (integer with values between 2 and the length of the record), or
#' enter the term "optimize" to let the function find the optimum window size
#' for the record through a Monte Carlo approach.
#' @param N Number of datapoints for Monte Carlo simulation (defaults to 1000)
#' @param p Threshold value for the p value of separating summer from winter
#' reconstructions. Defaults to 0.05 (95% confidence level)
#' @param d18O_fun String containing the name of the transfer function used to
#' convert temperature and \eqn{\delta^{18}O_{w}}{δ18Ow} to
#' \eqn{\delta^{18}O_{w}}{δ18Ow} data (for example: \code{"KimONeil97"} or
#' \code{"GrossmanKu86"}). Defaults to Kim and O'Neil (1997).
#' @param D47_fun String containing the name of the transfer function used to
#' convert temperature to D47 data (for example: \code{"Bernasconi18"} or
#' \code{"Jautzy20"}). Defaults to Bernasconi et al., 2018).
#' @param export Export table summary of result (CSV format)? \code{TRUE/FALSE}
#' @param export_raw Export tables containing all raw model
#' results before being merged into tidy tables? \code{TRUE/FALSE}
#' @return A data frame containing monthly reconstructions of D47, temperature,
#' \eqn{\delta^{18}O}{δ18O} of the precipitation fluid and
#' \eqn{\delta^{18}O_{w}}{δ18Ow}.
#' @references package dependencies: TTR
#' Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
#' aragonite: temperature effects, _Chemical Geology_ **1986**, _59.1_, 59–74.
#'     \url{https://doi.org/bvpzws}
#'
#' Kim, S., O'Niel, J.R., Equilibrium and nonequilibrium oxygen
#' isotope effects in synthetic carbonates, _Geochimica et Cosmochimica Acta_
#' **1997**, _61.16_, 3461–3475.
#'     \url{https://doi.org/c7bwbp}
#'
#' Dettman, D.L., Reische, A.K., Lohmann, K.C., Controls on the stable isotope
#' composition of seasonal growth bands in aragonitic fresh–water bivalves
#' (Unionidae), _Geochimica et Cosmochimica Acta_ **1999**, _63.7–8_, 1049–1057.
#'     \url{https://doi.org/cbb7zc}
#'
#' Brand, W.A., Coplen, T.B., Vogl, J., Rosner, M., Prohaska, T., Assessment of
#' international reference materials for isotope–ratio analysis (IUPAC Technical
#' Report), _Pure and Applied Chemistry_ **2014**, _86.3_, 425–467.
#'     \url{https://doi.org/fpc2}
#'
#' Kele, S., Breitenbach, S. F., Capezzuoli, E., Meckler, A. N., Ziegler, M.,
#' Millan, I. M., Kluge, T., Deák, J., Hanselmann, K. and John, C. M.,
#' Temperature dependence of oxygen– and clumped isotope fractionation in
#' carbonates: a study of travertines and tufas in the 6–95 C temperature range,
#' _Geochimica et Cosmochimica Acta_ **2015**, 168, 172–192.
#'     \url{https://doi.org/f7sgn6}
#'
#' Bernasconi, S.M., Müller, I.A., Bergmann, K.D., Breitenbach, S.F., Fernandez,
#' A., Hodell, D.A., Jaggi, M., Meckler, A.N., Millan, I. and Ziegler, M.,
#' Reducing uncertainties in carbonate–clumped isotope analysis through
#' consistent carbonate based standardization. _Geochemistry, Geophysics,
#' Geosystems_ **2018**, 19–9, 2895–2914.
#'     \url{https://doi.org/gfmjrw}
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
#'     \url{https://doi.org/ggrc39}
#'
#' Jautzy, J. J., Savard, M. M., Dhillon, R. S., Bernasconi, S. M. and Smirnoff,
#' A., Clumped isotope temperature calibration for calcite: Bridging theory and
#' experimentation, _Geochemical Perspectives Letters_ **2020**, 14, 36–41.
#'     \url{https://doi.org/fpc3}
#'
#' de Winter, N. J., Agterhuis, T., Ziegler, M., Optimizing sampling strategies
#' in high–resolution paleoclimate records, _Climate of the Past Discussions_
#' **2020**, 1–52.
#'     \url{https://doi.org/fpc4}
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
#'     monthly <- smoothing_seasonality(d18Oc,
#'     D47,
#'     ages,
#'     0.1,
#'     0.04,
#'     "optimize",
#'     1000,
#'     0.05,
#'     "KimONeil97",
#'     "Bernasconi18",
#'     FALSE,
#'     FALSE)
#'     }
#' @export
smoothing_seasonality <- function(d18Oc, # Sub–annually resolved d18Oc data 
    D47, # Sub–annually resolved D47 data
    ages, # Vector containing ages for of all samples in years relative to the shell chronology
    SD_d18Oc = 0.1, # Error (1 SD) on d18Oc data 
    SD_D47 = 0.04, # Error (1 SD) on D47 data 
    window = "optimize", # Either the size of the moving average window, or the word "optimize" (by default) to indicate that the optimal sampling window should be selected using a Monte Carlo approach
    N = 1000, # Number of Monte Carlo simulations for optimization
    p = 0.05, # p–value threshold for considering successful separation of seasons
    d18O_fun = "KimONeil97",
    D47_fun = "Bernasconi18",
    export = FALSE, # Should the result be exported? 
    export_raw = FALSE # Should raw data of successful individual simulations be exported (WARNING: Files can get large!)
    ){
    
    # Prepare data
    # Check if data has equal length
    if(length(unique(c(length(d18Oc), length(D47), length(ages)))) > 1){
        return("ERROR: Vectors 'd18Oc', 'D47' and 'ages' should have equal length")
    }
    Popt <- matrix(NA, ncol = 7, nrow = N) # Create matrix with maximum length
    colnames(Popt) <- c("optimal sample size",
        "number of successful simulations",
        "Optimum p-value",
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

    d18Omat <- as.data.frame(matrix(stats::rnorm(N * length(d18Oc), d18Oc, SD_d18Oc), ncol = N)) # Randomly resample d18O data using measurement uncertainty
    colnames(d18Omat) <- paste("Sim", seq(1, N, 1), sep = "")
    D47mat <- as.data.frame(matrix(stats::rnorm(N * length(D47), D47, SD_D47), ncol = N)) # Randomly resample D47 data using measurement uncertainty
    colnames(D47mat) <- paste("Sim", seq(1, N, 1), sep = "")

    if(window == "optimize"){
        # Expanding moving average window
        if(length(d18Oc) / ceiling(diff(range(ages))) < 2){ # If there are less than 2 datapoints per year, the window is allowed to become as long as the entire record
            win <- seq(2, length(d18Oc), 1)
        }else{
            win <- seq(2, length(d18Oc) / ceiling(diff(range(ages))), 1) # Create vector of sample size windows (max sample size is 1 year)
        }
        
        recond <- as.data.frame(matrix(NA, ncol = N, nrow = length(ages))) # Create matrix for d18O results
        reconD <- as.data.frame(matrix(NA, ncol = N, nrow = length(ages))) # Create matrix for D47 results
        colnames(recond) <- colnames(reconD) <- paste("Sim", seq(1, N, 1), sep = "")

        for(i in 1:N){
            # Progress
            cat(paste("Optimizing moving average iteration: ",i),"\r")
            utils::flush.console()

            Dwin <- vapply(win, function(x) max(TTR::runMean(D47mat[, i], x), na.rm = TRUE), 1) # Find highest D47 values (winter)
            Dsum <- vapply(win, function(x) min(TTR::runMean(D47mat[, i], x), na.rm = TRUE), 1) # Find lowest D47 values (summer)
            dwin <- vapply(win, function(x) max(TTR::runMean(d18Omat[, i], x), na.rm = TRUE), 1) # Find highest d18O values (winter)
            dsum <- vapply(win, function(x) min(TTR::runMean(d18Omat[, i], x), na.rm = TRUE), 1) # Find lowest d18O values (summer)
            Dwinsd <- vapply(win, function(x) TTR::runSD(D47mat[, i], x)[which.max(TTR::runMean(D47mat[, i], x))], 1) # Find standard deviation of winter D47 values
            Dsumsd <- vapply(win, function(x) TTR::runSD(D47mat[, i], x)[which.min(TTR::runMean(D47mat[, i], x))], 1) # Find standard deviation of summer D47 values
            dwinsd <- vapply(win, function(x) TTR::runSD(d18Omat[, i], x)[which.max(TTR::runMean(d18Omat[, i], x))], 1) # Find standard deviation of winter d18O values
            dsumsd <- vapply(win, function(x) TTR::runSD(d18Omat[, i], x)[which.min(TTR::runMean(d18Omat[, i], x))], 1) # Find standard deviation of summer d18O values

            SDpool <- sqrt((Dsumsd ^ 2 + Dwinsd ^ 2) / 2) # Calculate pooled standard deviation for each moving average window
            T <- (Dsum - Dwin) / (SDpool * sqrt(2 / win)) # Calculate two–sample T–value for each window (equal sample size, equal variance)
            Pval <- stats::pt(T, win - 1) # Calculate p–value for each window
            Popt[i, ] <- c(win[which(Pval == min(Pval, na.rm = TRUE))],
                length(which(Pval < 0.05)),
                min(Pval, na.rm = TRUE),
                dsum[which(Pval == min(Pval, na.rm = TRUE))],
                dwin[which(Pval == min(Pval, na.rm = TRUE))],
                Dsum[which(Pval == min(Pval, na.rm = TRUE))],
                Dwin[which(Pval == min(Pval, na.rm = TRUE))])
                
            # Shift reconstructed time series to correct for offset due to right–ordered moving average
            shift <- round(win[which(Pval == min(Pval, na.rm = TRUE))] / 2, 0)
            recond[, i] <- c(TTR::runMean(d18Omat[, i], win[which(Pval == min(Pval, na.rm=TRUE))])[(shift + 1):length(d18Omat[, i])], rep(NA, shift)) # Collect best moving window d18Oc data into matrix for export
            reconD[, i] <- c(TTR::runMean(D47mat[, i], win[which(Pval == min(Pval, na.rm=TRUE))])[(shift + 1):length(D47mat[, i])], rep(NA, shift)) # Collect best moving window D47 data into matrix for export
        }

        # OPTIONAL: Export results of optimized sample sizes
        if(export_raw == TRUE){
            utils::write.csv(Popt, "Optimized_moveing average_simulations.csv")
        }
    }else if(is.numeric(window)){ # Use fixed moving window to construct moving average
        # Shift reconstructed time series to correct for offset due to right–ordered moving average
        shift <- round(window / 2, 0)
        recond <- apply(d18Omat, 2, TTR::runMean, n = window)
        recond <- rbind(recond[(shift + 1):nrow(recond), ], matrix(NA, ncol = ncol(recond), nrow = shift))
        reconD <- apply(D47mat, 2, TTR::runMean, n = window)
        reconD <- rbind(reconD[(shift + 1):nrow(reconD), ], matrix(NA, ncol = ncol(reconD), nrow = shift))
    }else{
        print("ERROR: Window input not supported")
        return("ERROR: Window input not supported")
    }

    # Use age data to group results into monthly bins
    month <- ceiling((ages %% 1) * 12)

    # Calculate monthly statistics of all d18Oc values
    cat("Grouping d18Oc data into monthly bins: ", "\r")
    utils::flush.console()
    d18Oc_monthly <- data.frame(d18Oc_mean = vapply(1:12, function(x) mean(as.matrix(recond[which(month == x), ]), na.rm = TRUE), 1),
        d18Oc_median = vapply(1:12, function(x) stats::median(as.matrix(recond[which(month == x), ]), na.rm = TRUE), 1),
        d18Oc_SD = vapply(1:12, function(x) stats::sd(as.matrix(recond[which(month == x), ]), na.rm = TRUE), 1)
    )
    d18Oc_monthly$d18Oc_SE <- d18Oc_monthly$d18Oc_SD / sqrt(vapply(1:12, function(x) length(as.matrix(recond[which(month == x), 1])), 1))

    # Calculate monthly statistics of all D47 values using the d18Oc measurements and the D47–d18Oc slopes of all successful simulations
    cat("Grouping D47 data into monthly bins: ", "\r")
    utils::flush.console()
    D47_monthly <- data.frame(D47_mean = vapply(1:12, function(x) mean(as.matrix(reconD[which(month == x), ]), na.rm = TRUE), 1),
        D47_median = vapply(1:12, function(x) stats::median(as.matrix(reconD[which(month == x), ]), na.rm = TRUE), 1),
        D47_SD = vapply(1:12, function(x) stats::sd(as.matrix(reconD[which(month == x), ]), na.rm = TRUE), 1)
    )
    D47_monthly$D47_SE <- D47_monthly$D47_SD / sqrt(vapply(1:12, function(x) length(as.matrix(reconD[which(month == x), 1])), 1))

    # Repeat for monthly temperature reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping Temperature data into monthly bins: ", "\r")
    utils::flush.console()
    if(D47_fun == "Bernasconi18"){
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15, na.rm = TRUE), 1),
            T_median = vapply(1:12, function(x) stats::median(sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15, na.rm = TRUE), 1),
            T_SD = vapply(1:12, function(x) stats::sd(sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15, na.rm = TRUE), 1)
        )
    }else if(D47_fun == "Jautzy20"){
        T_monthly <- data.frame(T_mean = vapply(1:12, function(x) mean(sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15, na.rm = TRUE), 1),
            T_median = vapply(1:12, function(x) stats::median(sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15, na.rm = TRUE), 1),
            T_SD = vapply(1:12, function(x) stats::sd(sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15, na.rm = TRUE), 1)
        )
    }
    T_monthly$T_SE = T_monthly$T_SD / sqrt(vapply(1:12, function(x) length(as.matrix(reconD[which(month == x), 1])), 1))

    # Repeat for monthly d18Ow reconstructions by calculating temperatures for each combination before averaging
    cat("Grouping d18Ow data into monthly bins: ", "\r")
    utils::flush.console()
    if(d18O_fun == "KimONeil97"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1),
                d18Ow_median = vapply(1:12, function(x) stats::median(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1),
                d18Ow_SD = sqrt((vapply(1:12, function(x) stats::sd(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1)) ^ 2 + (1.03092 * d18Oc_monthly$d18Oc_SD) ^ 2) # Include MC simulated error on d18Oc in the analysis
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1),
                d18Ow_median = vapply(1:12, function(x) stats::median(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1),
                d18Ow_SD = sqrt((vapply(1:12, function(x) stats::sd(((as.matrix(recond[which(month == x), ]) / 1000 + 1) / exp(((18.03 * 10 ^ 3) / sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 32.42) / 1000) - 1) * 1000 * 1.03092 + 30.92, na.rm = TRUE), 1)) ^ 2 + (1.03092 * d18Oc_monthly$d18Oc_SD) ^ 2) # Include MC simulated error on d18Oc in the analysis
            )
        }
    }else if(d18O_fun == "GrossmanKu86"){
        if(D47_fun == "Bernasconi18"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1),
                d18Ow_median = vapply(1:12, function(x) stats::median(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1),
                d18Ow_SD = sqrt((vapply(1:12, function(x) stats::sd(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0449 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.167)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1)) ^ 2 + d18Oc_monthly$d18Oc_SD ^ 2) # Include MC simulated error on d18Oc in the analysis
            )
        }else if(D47_fun == "Jautzy20"){
            d18Ow_monthly <- data.frame(d18Ow_mean = vapply(1:12, function(x) mean(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1),
                d18Ow_median = vapply(1:12, function(x) stats::median(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1),
                d18Ow_SD = sqrt((vapply(1:12, function(x) stats::sd(as.matrix(recond[which(month == x), ]) + ((sqrt((0.0433 * 10 ^ 6) / (as.matrix(reconD[which(month == x), ]) - 0.119 - 0.066)) - 273.15) - 20.6) / 4.34 - 0.2, na.rm = TRUE), 1)) ^ 2 + d18Oc_monthly$d18Oc_SD ^ 2) # Include MC simulated error on d18Oc in the analysis
            )
        }
    }
    d18Ow_monthly$d18Ow_SE = d18Ow_monthly$d18Ow_SD / sqrt(vapply(1:12, function(x) length(as.matrix(recond[which(month == x), 1])), 1))

    monthly<-cbind(1:12,
        d18Oc_monthly,
        D47_monthly,
        T_monthly,
        d18Ow_monthly)

    colnames(monthly)[1] <- "month"

    # Export results of monthly grouped data
    if(export == TRUE){
        utils::write.csv(monthly, paste("Monthly_moving average_results.csv"))
    }

    return(monthly)
}