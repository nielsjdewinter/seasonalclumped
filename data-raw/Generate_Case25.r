# Case 25: Natural Case: Iceland
# Load data from provided dataset
GT <- read.csv(system.file("extdata", "Iceland_data.csv", package = "seasonalclumped"))
GT <- read.csv("inst/extdata/Iceland_data.csv") # From within package root directory
Ty <- GT[, 2]
GR <- GT[, 3]
SST <- GT[, 4]
d18Osw <- GT[, 5]

SR <- as.vector(c(0.1, 0.2, 0.45, 0.75, 1.55, 3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case25 <- data.frame(column = rep(NA, sum(GR, na.rm = TRUE) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR, na.rm = TRUE), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case25)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case25) - nrow(newdata)))
    }
    newdata <- cbind(Case25$column, newdata)
    # Add the new data to the storage dataframe
    Case25 <- cbind(Case25, newdata)
}
Case25$column <- NULL
colnames(Case25)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case25, file = "data/Case25.rda")