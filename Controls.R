############################################################
# Running control scenarios
# 05/05/2020
# R.E.Baker
############################################################

require("plotrix")
require("raster")
require("rgdal")
require("plotfunctions")
require("pals")

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
locations <- data.frame(lat = c( 51.5074, 40.7128, 28.7041, 10.6427, -4.4419,-6.5, -26.2041,  -34.6037, -37), 
                        lon = c(0.1278, -74.6, 77.1025, -71.6125,15.2663,106.8456, 28.0473,-58.7,144.9631 ), 
                        names = c("London","NY","Delhi","Maracaibo","Kinshasa","Jakarta", "Johannesburg","Buenos Aires","Melbourne")) # locations and lat lons


#parms
parms <- data.frame(model=c("Influenza","OC43","HKU1"), ClimVar = c(-180, -32.5,  -227.5), ImmLength = c(40, 62.5, 66.25) )

timeStart = 0 # intervention start time
timeEnd = 0 # intervention end time
R0fixed = 1.1 #reduction in R0 for intervention
R0max = 2.5 # coronavirus R0 max
R0min = 1.5 # coronavirus R0 min
Lead = 0 # start tine of epidemic in year (weeks)
SSet = "Orig" # proportion of initial susceptibles - "Orig" = pop - 1
pop = 8e+06 # divide by pop to get I/N 
timeLengthSim = 728# days to run Sim for

## set control times
timeEndList <- c(0,seq(56,483,28))[1:17] 
timeStartList <- c(0, rep(28, times = (length(timeEndList) -1)))

matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[2], Immunity =parms$ImmLength[2]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop

#set min max range for plotting
minplot = 0
maxplot = 1.2
maxplotcal = 0



setwd("~/Plots")
pdf("controlocAllM1.pdf",width=8,height=8)
par(mfrow=c(2,2))
par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                       R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)
#dev.off()




matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[3], Immunity =parms$ImmLength[3]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop


par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)
#dev.off()




R0fixed = 1.3 #reduction in R0 for intervention
matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[2], Immunity =parms$ImmLength[2]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop



par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)




R0fixed = 1.3 #reduction in R0 for intervention
matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[3], Immunity =parms$ImmLength[3]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop



par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)
dev.off()







########################### second set with new control times

R0fixed = 1.1 #reduction in R0 for intervention
timeEndList <- c(0,seq(70,511,28))[1:17]
timeStartList <- c(0, rep(42, times = (length(timeEndList) -1))) # six weeks

matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[2], Immunity =parms$ImmLength[2]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop

setwd("~/Plots")
pdf("controlocAllM2.pdf",width=8,height=8)
par(mfrow=c(2,2))

#setwd("~/Dropbox/Coronavirus/Reviews/Plots")
#pdf("controlocR11M2.pdf",width=5,height=5)
par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)
#dev.off()




rbPal <- coolwarm(15)[3:13]

matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[3], Immunity =parms$ImmLength[3]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop



par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)





R0fixed = 1.3 #reduction in R0 for intervention

matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[2], Immunity =parms$ImmLength[2]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop




par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)





R0fixed = 1.3 #reduction in R0 for intervention

matoc <- matrix(NA, nrow = nrow(locations), ncol = length(timeEndList))
outdf <- NULL
for(j in 1:length(timeEndList)){
  for(i in 1:9){
    predoc43 <- runModel(LatCity =locations$lat[i], LonCity =locations$lon[i], SSet = SSet, Lead = Lead, R0min = R0min, 
                         R0max = R0max, Var = parms$ClimVar[3], Immunity =parms$ImmLength[3]*7, pop = pop, 
                         timestart = timeStartList[j], 
                         timeend = timeEndList[j], 
                         R0fixed = R0fixed, timeLengthSim =timeLengthSim, SHDAT = qpt)
    
    matoc[i,j] <- max(predoc43$I, na.rm=T)
    minidf <- data.frame(name = locations$names[i], control = timeEndList[j], max = max(predoc43$I, na.rm=T)/pop)
    outdf <- rbind(outdf, minidf)
  }
}
matoc <- matoc/pop



par(mar=c(7,6,1,4))
rbPal <- coolwarm(15)[3:13]
nuse = 17
plot(seq(1,9,1), seq(1,nuse,l=9), type="n", yaxt="n", xaxt="n", bty= "n",ylab="Length of control (months)",xlab="")
for(i in 1:nrow(locations)){
  sub <- outdf[outdf$name==locations$names[i],]
  sub <- sub[order(sub$control),]
  sub$maxrel <- sub$max/(max(sub$max[1]))
  sub$maxrel[sub$maxrel > maxplot] <- maxplot
  if(max(sub$maxrel) > maxplotcal){maxplotcal = max(sub$maxrel)}
  tes <- rbPal[as.numeric(cut(c(minplot,maxplot,sub$maxrel),breaks = 11))]
  sub$col2 <- tes[3:length(tes)]
  
  points(rep(i,length=nuse),seq(1,nuse,l=nuse), cex = 2*sub$maxrel + 0.5 , pch = 16, col = sub$col2)
  r0 <- grabR0Quick(LatCity =locations$lat[i], LonCity =locations$lon[i] , R0min = R0min, 
                    R0max = R0max, Var =  -227.5,  SHDAT = qpt)
  r02 <- c(r0,r0)
  r0log <- r02[timeEndList/7]
  minr0log <- 1 + which.max(r0log)
  if(i==3) {points(c(i,i,i), c(0,minr0log, minr0log + 12),pch = 4, col="white")}
  if(i != 3) {points(c(i,i,i), c(minr0log - 12,minr0log, minr0log + 12),pch = 4, col="white")}
}
axis(1, at=1:9, labels=locations$names,las =2, title="")
axis(2, at=1:nuse, c("No Control",seq(1,nuse-1,1)),las = 1)
gradientLegend(c(minplot,maxplot), color = rbPal, nCol = 30, pos = 0.5, side = 4, inside = F, fit.margin = T, 
               length = 0.5, n.seg=5, dec = 2)
dev.off()







