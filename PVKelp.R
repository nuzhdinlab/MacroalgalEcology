rm(list=ls()) #Clear memory
require(dplyr)

#Set the working directory.
wd <- "~/Desktop/PVKelp/"
setwd(wd)

#Read in kelp data.
KelpInput1 <- read.table("VRG_CRANE_Data_Density_2019-11-26.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
KelpInput2 <- read.table("VRG_CRANE_Data_NoStipes_2019-11-25.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
KelpInput3 <- read.table("VRG_CRANE_Data_Stipes_2019-11-25.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Rename variable name.
colnames(KelpInput1)[colnames(KelpInput1)=="Mean_Density_m2"] <-"MeanDensity"

#Calculate the mean abundance of kelp by sample date, site, depth zone, location, and species.
KelpInput2Summarized <- KelpInput2 %>%
  group_by(SampleDate, Site, DepthZone, BenthicReefSpecies, Latitude, Longitude) %>% 
  summarise_each(funs(mean))

#Calculate the mean density of kelp.
KelpInput2Summarized$MeanDensity <- KelpInput2Summarized$SumAbundance/KelpInput2Summarized$Area_m2

#Remove unused columns.
KelpInput2Summarized <- as.data.frame(KelpInput2Summarized[,c("SampleDate","Site","DepthZone","BenthicReefSpecies","Latitude","Longitude","MeanDensity")])

#Calculate the mean abundance of kelp by sample date, site, depth zone, location, and species.
KelpInput3Summarized <- KelpInput3 %>%
  group_by(SampleDate, Site, DepthZone, BenthicReefSpecies, Latitude, Longitude) %>% 
  summarise_each(funs(mean))

#Calculate the mean density of kelp.
KelpInput3Summarized$MeanDensity <- KelpInput3Summarized$Abundance/KelpInput3Summarized$Area_m2

#Remove unused columns.
KelpInput3Summarized <- as.data.frame(KelpInput3Summarized[,c("SampleDate","Site","DepthZone","BenthicReefSpecies","Latitude","Longitude","MeanDensity")])

#Create merged dataframe.
KelpData <- rbind(KelpInput1,KelpInput2Summarized,KelpInput3Summarized)

#Standardize date format.
KelpData$SampleDate <- as.Date(KelpData$SampleDate,format="%d-%b-%y")

#Remove duplicate lines.
KelpData <- KelpData[!duplicated(KelpData),]

#Export merged data frame to begin extracting map layer values at its coordinates.
write.table(KelpData,"KelpData.txt",quote=FALSE,sep="\t",row.names = FALSE)

#Split kelp data into separate files by species.
for(species in unique(KelpData$BenthicReefSpecies)){
  write.table(KelpData[KelpData$BenthicReefSpecies==species,],paste("KelpData",gsub(" ", "", species, fixed = TRUE),".txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
  print(paste(species,nrow(KelpData[KelpData$BenthicReefSpecies==species,])))
}
