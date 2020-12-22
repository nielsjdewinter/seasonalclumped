#' Function used to linearly subsample data at new depth values
#'
#' @param data A vector of data to be interpolated
#' @param old_depth A vector containing the depth values belonging to \code{data}
#' @param new_depth A vector containing depth values at which the \code{data}
#' should be interpolated.
#' @param AV Should the subsampling take into account the mean value within the
#' sample interval? \code{TRUE/FALSE} If \code{FALSE}, the interpolated value
#' corresponding to the exact position is used instead of the mean of the interval
#' @param plot Should the result be plotted? \code{TRUE/FALSE}
#' @return A vector listing the values interpolated from \code{data} at the 
#' positions of \code{new_depth}
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
#' SR <- 0.75 # Set sampling resolution to 0.75 mm
#' # Create vector for all samples along entire shell length by applying constant
#' # sampling resolution
#' D <- seq(SR, sum(GR), SR)
#' D_cum <- cumsum(GR) # Create cumulative depth vector for all values
#' new_data <- subsample(SST, D_cum, D, AV = TRUE, plot = FALSE) # Interpolate
#' # SST values at the positions of D while calculating sample averages
#' @export
subsample <- function(data,
    old_depth,
    new_depth,
    AV = FALSE,
    plot = FALSE
    ){
    if(AV == TRUE){
        lwr <- new_depth - c(new_depth[1], diff(new_depth) / 2) # Find lower boundaries of sample intervals
        upr <- new_depth + c(diff(new_depth) / 2, new_depth[length(new_depth)]) # Find upper boundaries of sample intervals
        lwrpos <- vapply(lwr, function(x) which.min(abs(x - old_depth)), 1) # Find position of lower boundaries of sample intervals in old_data
        uprpos <- vapply(upr, function(x) which.min(abs(x - old_depth)), 1) # Find position of upper boundaries of sample intervals in old_data
        newdata <- vapply(1:length(lwrpos), function(x) mean(data[lwrpos[x]:uprpos[x]]), 1) # Calculate mean values between upper and lower boundaries for each sample interval
    }else{
        linterp <- stats::approx( # Calculate the new data value for each sample by linear interpolation
            x = old_depth,
            y = data,
            xout = new_depth,
            method = "linear"
        )
        newdata <- linterp$y
    }
    if(plot == TRUE){ # Create plot showing subsampling if requested
        grDevices::dev.new()
        plot(old_depth, data, type="l")
        graphics::points(new_depth, newdata, col="red")
    }
    return(newdata)
}