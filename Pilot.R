#This script runs Maxent on various kelp forest species using
#a number of environmental map layers.
rm(list=ls())
require(raster)
require(corrplot)
require(dismo)
require(ENMeval)
require(dplyr)

setwd("/Users/levisimons/Desktop/Pilot/")

#Read in environmental map layers.  These are all of the .tif files.
environment.files <- list.files(pattern=".tif$",full.names=TRUE)

#Create a stack of environmental map layers.
env.data <- stack(c(environment.files))

#List out all species presence point files.
#Make sure the species file names contain the substring m_ and .csv.
species.files <- list.files(pattern="(.*?)m_(.*?).csv$",full.names=TRUE)

#Initiate the MaxEnt model package
maxent()

#Lopp through each species observation file.
for(species.file in species.files){
  #Read in each species presence points file.
  obs.points <- read.csv(file=species.file)
  
  #Get species name from file.  Store it as a character string.
  species.name <- as.character(unique(obs.points$species))
  
  #Remove species name column from obs.points.
  obs.points <- obs.points[,c("longitude","latitude")]
  
  #Initialize data containing environmental layer values at presence locations.
  presence.points <- obs.points
  
  #Loop through each environmental file layer.
  for(environment.file in environment.files){
    #Create a one column data frame with the values of an environmental layer at the
    #locations of a given species.
    temp <- as.data.frame(extract(raster(environment.file),obs.points))
    
    #Rename the column header to reflect the environmental layer's filename.
    colnames(temp) <- environment.file
    
    #Keep adding the environmental value columns to the original dataframe
    #with the location values.
    presence.points <- cbind(presence.points,temp)
  }
  #Convert NaN to NA.  Make sure the dataframe is forced into being numeric
  presence.points[is.na(presence.points)] <- NA
  
  #Add the species name column back in.
  presence.points$species <- species.name
  
  #Write out each dataframe containing environment values at species locations as a separate file.
  write.table(presence.points,paste(species.name,"EnvironmentalValues.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
  
  #Find species locations with complete environmental data.
  obs.points <- presence.points[complete.cases(presence.points),c("longitude","latitude")]
  
  #Generate a large set of random points where species were not observed (pseudo-absences).
  abs.points <- randomPoints(mask=env.data,n=20*nrow(obs.points),p=obs.points,tryf=50,lonlatCorrection=FALSE)
  
  #Run a MaxEnt model using your environmental data, presence points, and pseudo-absence points.
  xm <- maxent(x=env.data,p=obs.points,a=abs.points)
  
  #Determine the relative importance of environmental variables in predicting the presence of a species.
  XMImportance <- var.importance(xm)
  XMImportance <- XMImportance[,c(1,3)]
  
  #Write out the relative importance of environmental variables in predicting the presence of a species.
  write.table(XMImportance,paste(species.name,"MaxentRI.csv",sep=""),quote=FALSE,sep=",",row.names = FALSE)
}


#Plot a correlogram of environmental variables for each species data set.

#List out all files containing species locations and environmental data.
#Make sure the species file names contain the substring EnvironmentalValues.csv.
species.environment.files <- list.files(pattern="EnvironmentalValues.csv$",full.names=TRUE)

#In a loop read in all of files with species location data, and environmental values.
for(species.environment.file in species.environment.files){
  #Read in your file and store it as a dataframe.
  data.input <- read.csv(file=species.environment.file)
  
  #Get species name.
  species.name <- unique(data.input$species)
  
  #Calculate p-values for each correlation matrix for the environmental values
  #associated with the location of each species.
  temp <- cor.mtest(cor(data.input[,c("..Bathym100m.tif","..chl_2019.tif","..sss000_2019.tif","..sst000_2019.tif")],method="pearson",use="complete.obs"))
  p.mat <- temp$p
  
  #Open the png file to write the plot to.
  png(height=382, width=695, file=paste(species.name,".png",sep=""))
  
  #Create a correlation plot, assuming complete observations and a Pearson correlation.
  #Plot correlations with a p value less than 0.01.
  corrplot(cor(data.input[,c("..Bathym100m.tif","..chl_2019.tif","..sss000_2019.tif","..sst000_2019.tif")],method="pearson",use="complete.obs"),type="lower", p.mat=p.mat, sig=0.01, diag=FALSE, tl.srt=15,main=paste(species.name,"environmental correlogram"),mar=c(0,0,1,0))
  
  #Save out image file.
  dev.off()
}
