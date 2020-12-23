#' Function that produces d18O and D47 records
#' 
#' Takes vectors of time, temperature, growth rate and d18O of the fluid and
#' converts them into a d18O and D47 record. The d18O and D47 values are
#' calculated for every depth value provided in the "D" vector. By default, the
#' empirical transfer function by Kim and O'Neil (1997) is used to produce the
#' d18O record, but other transfer functions (e.g. Grossman and Ku, 1986) are
#' also supported. The default transfer function for converting temperature data
#' to D47 data is based on Bernasconi et al. (2018), but other transfer
#' functions (e.g. Jautzy et al., 2020) are also supported.
#' 
#' @param time Time vector (values in years)
#' @param SST A vector containing temperature data (values in degrees C; length
#' must be equal to that of \code{time})
#' @param GR Growth rate vector (values in same time unit as \code{time} (years);
#' length must be equal to that of \code{time})
#' @param d18Ow A vector containing data on the d18O value of the precipitation
#' fluid (values in permille VSMOW; length must be equal to that of \code{time})
#' @param D Depth vector (values in same depth unit as \code{GR})
#' @param d18O_fun String containing the name of the transfer function used to
#' convert temperature and d18Ow to d18Oc data (for example: \code{"KimONeil97"}
#' or \code{"GrossmanKu86"}). Defaults to Kim and O'Neil (1997).
#' @param D47_fun String containing the name of the transfer function used to
#' convert temperature to D47 data (for example: \code{"Bernasconi18"} or
#' \code{"Jautzy20"}). Defaults to Bernasconi et al., 2018).
#' @param AV Should the subsampling take into account the mean value within the
#' sample interval? \code{TRUE/FALSE} If \code{FALSE}, the interpolated value
#' corresponding to the exact position is used instead of the mean of the interval
#' @param plot Should the result be plotted? \code{TRUE/FALSE}
#' @return A matrix containing subsampled time, depth, d18Oc and D47 values:
#' \code{"Tnew"}): New time vector after subsampling
#' \code{"D"}): New depth vector after subsampling
#' \code{"d18Oc"}): Vector listing d18Oc values for each sample
#' \code{"D47"}): Vector listing D47 values for each sample
#' @references function dependencies: subsample, subsample_mean
#' Grossman, E.L., Ku, T., Oxygen and carbon isotope fractionation in biogenic
#' aragonite: temperature effects, _Chemical Geology_ **1986**, _59.1_, 59–74.
#'     \url{http://doi.org/bvpzws}
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
#' @examples
#' # Create test data (= ideal case)
#' # Set boundary conditions
#' Td <- seq(1, 12 * 365, 1) # Create timeline of 12 years in days
#' Ty <- Td / 365 # Convert to years
#' MAT <- 20 # Set mean annual temperature
#' Amp <- 10 # Set seasonal amplitude
#' Sext <- 2 * Amp # Calculate extent of seasonal variability
#' TSD <- 1.5 # Set the degree of random non–seasonal noise on the SST curve
#' # ("weather")
#' SST <- rnorm(length(Ty), MAT + Amp * sin(2 * pi * Ty), TSD) # Create virtual
#' # daily SST data
#' GR <- rep(10 / 365, length(Ty)) # Set growth rate to 10 mm/yr and create daily
#' # GR vector
#' DSD <- 0.6 # Set the degree of random non–seasonal noise on the d18Osw curve
#' # ("salinity fluctuations")
#' d18Osw<-rnorm(length(Ty), rep(0, length(Ty)), DSD) # Set d18Osw to 0 permille
#' # VSMOW, create daily d18Osw vector
#' SR <- 0.75 # Set sampling resolution to 0.75 mm
#' # Create vector for all samples along entire shell length by applying constant
#' # sampling resolution
#' D <- seq(SR, sum(GR), SR)
#' # Calculate virtual data
#' newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV = TRUE)
#' @export
carbmodel<-function(time,
    SST,
    GR,
    d18Ow,
    D,
    d18O_fun = "KimONeil97",
    D47_fun = "Bernasconi18",
    AV=FALSE,
    plot=FALSE
    ){
    D_cum <- cumsum(GR) # Create vector linking days to depth values
    SSTnew <- subsample(SST, D_cum, D, AV = AV) # Subsample SST along the new sample set, using mean values if AV = TRUE
    d18Ownew <- subsample(d18Ow, D_cum, D, AV = AV) # Subsample d18Ow along the new sample set, using mean values if AV = TRUE
    Tnew <- subsample(time, D_cum, D, AV = AV) # Subsample time (yr) along the new sample set, using mean values if AV = TRUE
    if(d18O_fun == "KimONeil97"){
        alpha = exp((18.03 * 1000 / (SSTnew + 273.15) - 32.42) / 1000) # Calculate alpha of calcite fractionation
        d18Ow_PDB = (0.97002 * d18Ownew - 29.98) # Convert d18Ow to PDB (following Brand et al., 2014)
        d18Oc = ((alpha * (d18Ow_PDB / 1000 + 1)) - 1) * 1000 # Calculate d18O of calcite for each sample according to Kim and O'Neil, 1997
    }else if(d18O_fun == "GrossmanKu86"){
        d18Oc <- (20.6 - SSTnew) / 4.34 + d18Ow + 0.2 # Use Grossmann and Ku (1986) modified by Dettmann et al. (1999)
    }else{
        return("ERROR: Supplied d18Oc transfer function is not recognized")
    }
    if(D47_fun == "Bernasconi18"){
        D47 <- (0.0449 * 10 ^ 6) / (SSTnew + 273.15) ^ 2 + 0.167 # Calculate D47 of calcite for each sample according to Kele et al., 2015 modified by Bernasconi et al., 2018
    }else if(D47_fun == "Jautzy20"){
        D47 <- (0.0433 * 10 ^ 6) / (SSTnew + 273.15) ^ 2 + 0.119 + 0.066 # Calculate D47 of calcite for each sample according to Jautzy et al., 2020 brought into 25 degrees CDES reference frame using 70–25 acid fractionation factor by Petersen et al., 2019
    }else{
        return("ERROR: Supplied D47 transfer function is not recognized")
    }
    if(plot == TRUE){ # Create plots of new data if requested
        plot(D, d18Oc, col = "blue")
        graphics::lines(D, d18Oc, col = "blue")
        graphics::par(new = TRUE)
        plot(D, D47, axes = FALSE, bty = "n", xlab = "", ylab = "", col = "red")
        graphics::lines(D, D47, bty = "n", xlab = "", ylab = "", col = "red")
        graphics::axis(side = 4, at = pretty(range(D47)))
    }
    dat<-cbind(Tnew, D, d18Oc, D47) # Combine new data for export
    return(dat) # Return the new depth, d18Oc and D47 series
}