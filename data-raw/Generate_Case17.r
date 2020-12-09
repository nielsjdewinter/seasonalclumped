# Case 17: Tropical case - Slight seasonal change in growth rate with slower growth in summer and linear growth decrease.
# Confined temperature seasonality, relatively strong multi-annual trend
# Strong d18Osw seasonality, light in summer, multi-annual trend in antiphase with multi-annual SST (ENSO-style; e.g. Iijima et al., 2005 - Geophysical Research Letters)
# Set boundary conditions
Td <- seq(1,12*365,1) # Create timeline of 12 years in days
Ty <- Td/365 # Convert to years
MAT<-28 # Set mean annual temperature (high/tropical)
Amp<-2 # Set seasonal amplitude (confined)
Sext<-2*Amp # Calculate extent of seasonal variability
TSD<-1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SSTseas<-rnorm(length(Ty),MAT+Amp*sin(2*pi*Ty),TSD) # Create virtual daily SST data
Pmulti<-10 # Set period of multi-annual variability (10 years)
Amulti<-1.5 # Set amplitude of multi-annual variability (5 degrees)
SSTmulti<-Amulti*sin(2*pi*1/Pmulti*Ty) # Create multi-annual component of SST data
SST<-SSTseas+SSTmulti # Combine components to create SST data

GRmeanstart<-10 # Set initial annual average growth rate (10 mm/yr)
GRmeanslope<--0.5 # Set annual reduction in growth rate
GRmean<-(GRmeanstart+GRmeanslope*Ty)
GRamp<-2 # Set seasonal amplitude of growth rate (confined, 2 mm/yr)
GR<-(GRmean-GRamp*sin(2*pi*Ty))/365 # Calculate daily growth rates, lowest growth rates in summer

d18Oswmean<-0 # Set annual average d18O of seawater (0 permille VSMOW)
d18Oswamp<-1 # Set seasonal amplitude of seawater d18O (1 permille)
d18Oswseas<-d18Oswmean-d18Oswamp*sin(2*pi*Ty)
d18OswPmulti<-10 # Set period of multi-annual d18Osw variability (10 years)
d18OswAmulti<-0.5 # Set amplitude of multi-annual d18Osw variability (1 permille)
d18Oswmulti<--1*d18OswAmulti*sin(2*pi*1/d18OswPmulti*Ty) # Create multi-annual d18Osw oscillation
DSD<-0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw<-rnorm(length(Ty),d18Oswseas+d18Oswmulti,DSD)

SR<-as.vector(c(0.1,0.2,0.45,0.75,1.55,3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case17 <- data.frame(column = rep(NA, sum(GR) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case17)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case17) - nrow(newdata)))
    }
    newdata <- cbind(Case17$column, newdata)
    # Add the new data to the storage dataframe
    Case17 <- cbind(Case17, newdata)
}
Case17$column <- NULL
colnames(Case17)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case17, file = "Case17.rda")