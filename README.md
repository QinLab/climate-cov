# Code for Baker et al. 2020

The main simulation model can be run using **SimMaster.R** which creates Figure 2 for nine different cities. 
The main function **runModel** can produce model runs under several different scenarios and acts as a wrapper for the climate-driven SIR model.
Inputs into runModel include:
* LatCity *Latitude of city* 
* LonCity *Longitude of city*
* SSet *Initial proportion susceptible* - either a number between 0 and 1 or "Orig" for N - 1
* Lead *Start time of the model in the year (days)* set to 0 to start at the beginning of the year
* R0min *Minimum R0*
* R0max *Maximum R0*
* Var *Climate dependency variable* 
* Immunity  
* pop = pop, 
* timestart = timeStart
* timeend = timeEnd, 
* R0fixed = R0fixed
* timeLengthSim =timeLengthSim
* SHDAT = qpt