# Case 14: Marine case - Seasonal change in growth rate in phase with temperature with linear growth decrease and dynamic growth threshold.
# Multi-annual cyclicity in d18Osw (NAO-style; Sarafanov, 2009 - ICES Journal of Marine Science)
# Set boundary conditions
Td <- seq(1,12*365,1) # Create timeline of 12 years in days
Ty <- Td/365 # Convert to years
MAT<-20 # Set mean annual temperature
Amp<-10 # Set seasonal amplitude
Sext<-2*Amp # Calculate extent of seasonal variability
TSD<-1.5 # Set the degree of random non-seasonal noise on the SST curve ("weather")
SST<-rnorm(length(Ty),MAT+Amp*sin(2*pi*Ty),TSD) # Create virtual daily SST data

GRmeanstart<-10 # Set initial annual average growth rate (10 mm/yr)
GRmeanslope<--0.5 # Set annual reduction in growth rate
GRmean<-(GRmeanstart+GRmeanslope*Ty)
GRamp<-5 # Set seasonal amplitude of growth rate (5 mm/yr)
GR<-(GRmean+GRamp*sin(2*pi*Ty))/365 # Calculate daily growth rates
Gmax<-28-0.5*Ty # Create dynamic maximum growth temperature threshold
Gmin<-12+0.5*Ty # Create dynamic minimum growth temperature threshold
GR[(SST>Gmax)|(SST<Gmin)]<-0 # Apply dynamic thresholds to growth rate vector

d18OswPmulti<-10 # Set period of multi-annual d18Osw variability (10 years)
d18OswAmulti<-0.5 # Set amplitude of multi-annual d18Osw variability (0.5 permille)
DSD<-0.6 # Set the degree of random non-seasonal noise on the d18Osw curve ("salinity fluctuations")
d18Osw<-rnorm(length(Ty),d18OswAmulti*sin(2*pi*1/d18OswPmulti*Ty),DSD) # Create multi-annual d18Osw oscillation

SR<-as.vector(c(0.1,0.2,0.45,0.75,1.55,3.25)) # Set sampling resolutions at 3.3 mm (~3 yr-1), 1.55 mm (~6 yr-1; bimonthly), 0.75 mm (~12 yr-1; monthly), 0.45 mm (~25 yr-1), 0.2 mm (~50 yr-1) and 0.1 mm (~100 yr-1, maximum isotope sampling)
# Sampling resolutions for courser sampling are deliberately chosen as non-multiples of the growth rate (irregular numbers) to prevent bias against some months

# Loop through vector and calculate D, d18Oc and D47 data for all sampling densities
Case14 <- data.frame(column = rep(NA, sum(GR) / SR[1]))
for(i in 1:length(SR)){
    # Create vector for all samples along entire shell length by applying constant sampling resolution
    D <- seq(SR[i], sum(GR), SR[i])
    # Calculate virtual data
    newdata <- carbmodel(Ty, SST, GR, d18Osw, D, AV=TRUE)
    # Increase length of new data to match the storage dataframe
    if(nrow(newdata) < nrow(Case14)){
        newdata <- rbind(newdata, matrix(NA, ncol = ncol(newdata), nrow = nrow(Case14) - nrow(newdata)))
    }
    newdata <- cbind(Case14$column, newdata)
    # Add the new data to the storage dataframe
    Case14 <- cbind(Case14, newdata)
}
Case14$column <- NULL
colnames(Case14)[seq(1, 26, 5)] <- paste("SR_", SR)
save(Case14, file = "Case14.rda")