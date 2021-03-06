# Case 1: Control data (= ideal case)
# Set boundary conditions
Td <- seq(1, 12 * 365, 1) # Create timeline of 12 years in days
Ty <- Td / 365 # Convert to years
MAT <- 20 # Set mean annual temperature
Amp <- 10 # Set seasonal amplitude
Sext <- 2 * Amp # Calculate extent of seasonal variability
TSD <- 1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SST <- rnorm(length(Ty), MAT + Amp * sin(2 * pi * Ty), TSD) # Create virtual daily SST data
GR <- rep(10 / 365, length(Ty)) # Set growth rate to 10 mm/yr, create daily GR vector
DSD <- 0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw <- rnorm(length(Ty), rep(0, length(Ty)), DSD) # Set d18Osw to 0 permille VSMOW, create daily d18Osw vector

SR<-as.vector(c(0.1,0.2,0.45,0.75,1.55,3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case1 <- data.frame(column = rep(NA, sum(GR) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case1)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case1) - nrow(newdata)))
    }
    newdata <- cbind(Case1$column, newdata)
    # Add the new data to the storage dataframe
    Case1 <- cbind(Case1, newdata)
}
Case1$column <- NULL
colnames(Case1)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case1, file = "Case1.rda")