############################################################
# Main Simulation Model 
# 05/05/2020
# R.E.Baker
############################################################

require("plotrix")
require("raster")
require("rgdal")

# Load functions
source('~/Functions.R', encoding = 'UTF-8')

# load continental shapefile - to remove ocean based pixels
setwd("~/TM_WORLD_BORDERS")
overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))

# load humidity data, trim to within land-mass lines
setwd("~/ClimateDataPath")
q  <- brick("merra_weekly_climo.nc") # weekly average MERRA data
qmask <- mask(q,overlay)
qpt <- rasterToPoints(qmask)
qpt <- as.data.frame(qpt)
names(qpt)[1:2] <- c("lon","lat")

# NB/ coastal lat/lons shifted to nearest land-based pixel
locations <- data.frame(lat = c(40.7128, 51.5074, 28.7041, -37, -34.6037, -26.2041, -6.5, 10.6427, -4.4419), 
                        lon = c(-74.6, 0.1278, 77.1025, 144.9631, -58.7,28.0473, 106.8456, -71.6125,15.2663 ), 
                        names = c("NY","London","Delhi","Melbourne","Buenos Aires","Johannesburg","Jakarta","Maracaibo","Kinshasa")) # locations and lat lons

#Climate and immunity parameters
parms <- data.frame(model=c("Influenza","OC43","HKU1"), ClimVar = c(-180, -32.5,  -227.5), ImmLength = c(40, 62.5, 66.25) )

timeStart = 0 # intervention start time
timeEnd = 0 # intervention end time
R0fixed = 0 #reduction in R0 for intervention
R0max = 2.5 # coronavirus R0 max
R0min = 1.5 # coronavirus R0 min
Lead = 0 # start time of epidemic in year (weeks)
SSet = "Orig" # proportion of initial susceptibles - "Orig" = pop - 1
pop = 8e+06 # divide by pop to get I/N 
timeLengthSim = 728 # days to run Sim for


#### Run model for influenza parameters
flupreds <- NULL
for(i in 1:nrow(locations)){
  pred <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                                 R0max = R0max, Var = parms$ClimVar[1], Immunity =parms$ImmLength[1]*7, pop = pop, timestart = timeStart, timeend = timeEnd, 
                                 R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
flupreds[[i]] <- pred$I

}
predNOCLIM <- runModel(LatCity = locations$lat[1], LonCity =locations$lat[1], SSet = "Orig", Lead = Lead, R0min = R0min, R0max = R0max,
                                  Var = 0, Immunity = parms$ImmLength[1]*7, pop = pop, timestart = timeStart, timeend = timeEnd,
                       R0fixed = R0fixed, timeLengthSim = timeLengthSim, SHDAT = qpt)
flupreds[[10]] <- predNOCLIM$I


#### Run model for OC43 parameters
ocpreds <- NULL
for(i in 1:nrow(locations)){
  pred <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                   R0max = R0max, Var = parms$ClimVar[2], Immunity =parms$ImmLength[2]*7, pop = pop, timestart = timeStart, timeend = timeEnd, 
                   R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
  ocpreds[[i]] <- pred$I
  
}
predNOCLIM <- runModel(LatCity = locations$lat[1], LonCity =locations$lat[1], SSet = "Orig", Lead = Lead, R0min = R0min, R0max = R0max,
                       Var = 0, Immunity = parms$ImmLength[2]*7, pop = pop, timestart = timeStart, timeend = timeEnd,
                       R0fixed = R0fixed, timeLengthSim = timeLengthSim, SHDAT = qpt)
ocpreds[[10]] <- predNOCLIM$I


#### Run the model for HKU1 parameters
hkupreds <- NULL
for(i in 1:nrow(locations)){
  pred <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                   R0max = R0max, Var = parms$ClimVar[3], Immunity =parms$ImmLength[3]*7, pop = pop, timestart = timeStart, timeend = timeEnd, 
                   R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
  hkupreds[[i]] <- pred$I
  
}
predNOCLIM <- runModel(LatCity = locations$lat[1], LonCity =locations$lat[1], SSet = "Orig", Lead = Lead, R0min = R0min, R0max = R0max,
                       Var = 0, Immunity = parms$ImmLength[3]*7, pop = pop, timestart = timeStart, timeend = timeEnd,
                       R0fixed = R0fixed, timeLengthSim = timeLengthSim, SHDAT = qpt)
hkupreds[[10]] <- predNOCLIM$I


########## Plotting Figure 2
setwd("~/Plots")
pdf("Fig2.pdf",width=10,height=6)
par(mfrow=c(3,3))
par(mar=c(3,3,1,1))
for(i in c(1,4,7)){
    
    name = ""
    if(i==1){name = "Influenza"}
    year <- seq(1,5,1/364)[1:timeLengthSim]
    plot(year, flupreds[[i]]/pop, type="l", lwd = 2,ylab = "",xlab = "", main = name,
     ,xlim=c(1,3),ylim=c(0,1500000/pop),bty="n", col="black")
    lines(year, flupreds[[i+1]]/pop, type="l", lwd = 2, col="deepskyblue3")
    lines(year, flupreds[[i+2]]/pop, type="l", lwd = 2, col="firebrick")
    lines(year, flupreds[[10]]/pop, type="l", lwd = 1, col="grey64",lty =3)
    legend(2,1200000/pop, lwd=c(2,2,2), col=c("black","deepskyblue3","firebrick"), legend = locations$names[i:(i+2)])
    if(i ==7){title(xlab="Year", line = 2)}
    title(ylab="I/N", line = 2)
    year <- seq(1,5,1/364)[1:timeLengthSim]
    
    name = ""
    if(i==1){name = "OC43"}
    plot(year, ocpreds[[i]]/pop, type="l", lwd = 2,ylab = "",xlab = "", main = name,
         ,xlim=c(1,3),ylim=c(0,1500000/pop),bty="n", col="black")
    lines(year, ocpreds[[i+1]]/pop, type="l", lwd = 2, col="deepskyblue3")
    lines(year, ocpreds[[i+2]]/pop, type="l", lwd = 2, col="firebrick")
    lines(year, ocpreds[[10]]/pop, type="l", lwd = 1, col="grey64",lty =3)
    legend(2,1200000/pop, lwd=c(2,2,2), col=c("black","deepskyblue3","firebrick"), legend = locations$names[i:(i+2)])
    if(i ==7){title(xlab="Year", line = 2)}
    title(ylab="I/N", line = 2)
    
    name = ""
    if(i==1){name = "HKU1"}
    plot(year, hkupreds[[i]]/pop, type="l", lwd = 2,ylab = "",xlab = "", main = name,
         ,xlim=c(1,3),ylim=c(0,1500000/pop),bty="n", col="black")
    lines(year, hkupreds[[i+1]]/pop, type="l", lwd = 2, col="deepskyblue3")
    lines(year, hkupreds[[i+2]]/pop, type="l", lwd = 2, col="firebrick")
    lines(year, hkupreds[[10]]/pop, type="l", lwd = 1, col="grey64",lty =3)
    legend(2,1200000/pop, lwd=c(2,2,2), col=c("black","deepskyblue3","firebrick"), legend = locations$names[i:(i+2)])
    if(i ==7){title(xlab="Year", line = 2)}
    title(ylab="I/N", line = 2)
}
dev.off()

