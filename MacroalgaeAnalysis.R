require(dplyr)
require(zetadiv)
require(netassoc)
require(data.table)

setwd("~/Desktop/Macroalgae")
#Read in taxanomic count data from Santa Barbara LTER.
RawInput <- read.table("Annual_All_Species_Biomass_at_transect_20181127.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)
FilteredInput <- RawInput
#Set all non-data cells to NA.
FilteredInput[FilteredInput==-99999] <- NA

#Read in location data for all transects.
TransectLocations <- read.table("SBLTERSiteLocations.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Merge in location data with taxonomic data.
GISBioData <- left_join(FilteredInput,TransectLocations,by=c("SITE","TRANSECT"))

#Remove entries where taxa aren't resolved to genus.
GISBioData <- GISBioData[!is.na(GISBioData$TAXON_GENUS),]

#Set NA abundance values to 0 for downstream presence/absence operations.
GISBioData$PERCENT_COVER[is.na(GISBioData$PERCENT_COVER)] <- 0
GISBioData$DENSITY[is.na(GISBioData$DENSITY)] <- 0

#Create presence/absence column for taxa.  1 for present, 0 for absent.
#First sum the percent cover and density columns to do this.
GISBioData$PA <- GISBioData$PERCENT_COVER+GISBioData$DENSITY
#Now convert all values above 0 to 1.
GISBioData$PA[GISBioData$PA > 0] <- 1

#Create a unique ID column to represent each unique site-transect pairing.
GISBioData$UniqueID <- paste(GISBioData$SITE,GISBioData$TRANSECT,sep="")

#Get all unique taxa in the data set.
uniqueTaxa <- as.data.frame(unique(GISBioData$SCIENTIFIC_NAME))
colnames(uniqueTaxa) <- c("SCIENTIFIC_NAME")
#Coerce factor to character for species names.
uniqueTaxa$SCIENTIFIC_NAME <- as.character(uniqueTaxa$SCIENTIFIC_NAME)
uniqueTaxa <- arrange(uniqueTaxa,SCIENTIFIC_NAME)

#Read in substrate data for all transects.
#Substrate codes are here: https://portal.edirepository.org/nis/metadataviewer?packageid=knb-lter-sbc.15.27
#B = bedrock, BL = boulders larger than 1m diameter, BM = bolders with a diameter between 50 and 100cm.
#BS = boulders with a diameter between 25 and 50cm, C = rocks with a diameter less than 25cm.
#S = sand greater than 2.5cm deep, SH = shell debris, SS = sand less than 2.5cm deep.
SubstrateInput <- read.table("Annual_Substrate_All_Years_20181127.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,check.names=FALSE)

#Calculate average percent coverage by substrate type per unique year-site-transect
Substrate <- aggregate(SubstrateInput$PERCENT_COVER,by=list(SubstrateInput$YEAR,SubstrateInput$SITE,SubstrateInput$TRANSECT,SubstrateInput$SUBSTRATE_TYPE), FUN=mean, na.rm=TRUE)
#Rename columns
colnames(Substrate) <- c("YEAR","SITE","TRANSECT","SUBSTRATE_TYPE","PERCENT_SUBSTRATE_TYPE")

#Set the order out to which to calculate zeta diversity.
samplingNum <- 8

for(yr in unique(GISBioData$YEAR)){
  for(group in unique(GISBioData$COARSE_GROUPING)){
    selected <- subset(GISBioData,YEAR==yr & COARSE_GROUPING==group)
    #Create presence/absence matrix of taxa in samples.
    PresenceAbsence <- uniqueTaxa
    for(ID in unique(selected$UniqueID)){
      #Presence/Absence matrix for taxa.
      sampleDF <- subset(selected,UniqueID == ID)
      sampleDF <- sampleDF[,c("SCIENTIFIC_NAME","PA")]
      PresenceAbsence <- merge(PresenceAbsence,sampleDF,by="SCIENTIFIC_NAME",all=TRUE)
      colnames(PresenceAbsence)[which(names(PresenceAbsence)=="PA")] <- ID
    }
    #Reformat presence/absence matrix for use in zeta diversity package.
    #Rows for sample ID and columns 
    data.LTER <- as.data.frame(t(PresenceAbsence[,-c(1)]))
    colnames(data.LTER) <- uniqueTaxa$SCIENTIFIC_NAME
    #Set NA abundance values to 0 for downstream presence/absence operations.
    data.LTER[is.na(data.LTER)] <- 0
    #Calculate zeta diversity decay if enough species are present.
    if(nrow(sampleDF)>=samplingNum){
      #Calculate zeta diversity decay for taxa within each watershed set of samples.
      zetaDecay <- Zeta.decline.ex(data.LTER,orders=1:samplingNum,plot=FALSE,rescale=TRUE)
      ExpIntercept <- zetaDecay$zeta.exp$coefficients[1] #Zeta diversity exponential decay intercept.
      ExpExp <- zetaDecay$zeta.exp$coefficients[2] #Zeta diversity exponential decay exponent.
      ExpAIC <- zetaDecay$aic$AIC[1] #AIC coefficient Zeta diversity exponential decay.
      zeta_N <- Zeta.order.ex(data.LTER,order=samplingNum,rescale=TRUE)$zeta.val #Higher order zeta diversity measure.
      PLExp <- zetaDecay$zeta.pl$coefficients[2] #Zeta diversity power law decay exponent.
      PLAIC <- zetaDecay$aic$AIC[2] #AIC coefficient Zeta diversity power law decay.
      print(paste(yr,group,ID,nrow(sampleDF),ExpAIC,PLAIC))
    }
  }
}
