# Case 11: Increase in d18Osw (evaporation) in summer, d18Osw constant rest of year
# Set boundary conditions
Td <- seq(1,12*365,1) # Create timeline of 12 years in days
Ty <- Td/365 # Convert to years
MAT<-20 # Set mean annual temperature
Amp<-10 # Set seasonal amplitude
Sext<-2*Amp # Calculate extent of seasonal variability
TSD<-1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SST<-rnorm(length(Ty),MAT+Amp*sin(2*pi*Ty),TSD) # Create virtual daily SST data
GR<-rep(10/365,length(Ty)) # Set growth rate to 10 mm/yr, create daily GR vector
d18Oswmean<-0 # Set annual average d18O of seawater (0 permille VSMOW)
d18Oswpulse<-1 # Set increase in seawater d18O due to summer evaporation
# Create annual profile due to summer evaporation based on a Gaussian curve.
# Maximum value is equal to height of pulse.
# 6 sigma values span the entire summer season (sigma = 15 days, 3 sigma on each side of the center of the peak)
# The peak centered in the middle of the summer season (on day 137)
# For the rest of the year, d18Osw remains constant at 0 permille.
d18Oswyear<-c(rep(d18Oswmean,round(365/4,0)),d18Oswpulse*exp(-0.5*((seq(1,round(365/4,0))-round(365/(4*2),0))/round(365/(4*2*3),0))^2),rep(d18Oswmean,365-2*round(365/4,0)))
DSD<-0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw<-rnorm(length(Ty),rep(d18Oswyear,round(Ty[length(Ty)],0)),DSD) # Stitch years together to create d18Osw vector

SR<-as.vector(c(0.1,0.2,0.45,0.75,1.55,3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case11 <- data.frame(column = rep(NA, sum(GR) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case11)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case11) - nrow(newdata)))
    }
    newdata <- cbind(Case11$column, newdata)
    # Add the new data to the storage dataframe
    Case11 <- cbind(Case11, newdata)
}
Case11$column <- NULL
colnames(Case11)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case11, file = "Case11.rda")