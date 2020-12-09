# Case 15: Coastal case - Seasonal change in growth rate with fast growth in spring and linear growth decrease.
# Pulse of light d18Osw in spring and multi-annual cyclicity in SST
# Set boundary conditions
Td <- seq(1,12*365,1) # Create timeline of 12 years in days
Ty <- Td/365 # Convert to years
MAT<-20 # Set mean annual temperature
Amp<-10 # Set seasonal amplitude
Sext<-2*Amp # Calculate extent of seasonal variability

TSD<-1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SSTseas<-rnorm(length(Ty),MAT+Amp*sin(2*pi*Ty),TSD) # Create seasonal component of SST data
Pmulti<-10 # Set period of multi-annual variability (10 years)
Amulti<-1.5 # Set amplitude of multi-annual variability (5 degrees)
SSTmulti<-Amulti*sin(2*pi*1/Pmulti*Ty) # Create multi-annual component of SST data
SST<-SSTseas+SSTmulti # Combine components to create SST data

GRmeanstart<-10 # Set initial annual average growth rate (10 mm/yr)
GRmeanslope<--0.5 # Set annual reduction in growth rate
GRmean<-(GRmeanstart+GRmeanslope*Ty)
GRamp<-5 # Set seasonal amplitude of growth rate (5 mm/yr)
GR<-(GRmean+GRamp*sin(2*pi*Ty+0.5*pi))/365 # Calculate daily growth rates

d18Oswmean<-0 # Set annual average d18O of seawater (0 permille VSMOW)
d18Oswpulse<--1 # Set decrease in seawater d18O due to freshwater pulse
# Create annual profile due to freshwater pulse based on a Gaussian curve.
# Maximum value is equal to height of pulse.
# 6 sigma values span the entire spring season (sigma = 15 days, 3 sigma on each side of the center of the peak)
# The peak centered in the middle of the spring season (on day 46)
# For the rest of the year, d18Osw remains constant at 0 permille.
d18Oswyear<-c(d18Oswpulse*exp(-0.5*((seq(1,round(365/4,0))-round(365/(4*2),0))/round(365/(4*2*3),0))^2),rep(d18Oswmean,365-round(365/4,0)))
DSD<-0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw<-rnorm(length(Ty),rep(d18Oswyear,round(Ty[length(Ty)],0)),DSD) # Stitch years together to create d18Osw vector

SR<-as.vector(c(0.1,0.2,0.45,0.75,1.55,3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case15 <- data.frame(column = rep(NA, sum(GR) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case15)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case15) - nrow(newdata)))
    }
    newdata <- cbind(Case15$column, newdata)
    # Add the new data to the storage dataframe
    Case15 <- cbind(Case15, newdata)
}
Case15$column <- NULL
colnames(Case15)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case15, file = "Case15.rda")