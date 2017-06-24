rm(list = ls())
# examine initial fitness values graphically
setwd("~/Documents/Research/SpeciationGenomics/fffits/")
fitness <- read.csv("InitialFitnessValues.csv")
boxplot(fitness$fitness ~ fitness$location)

#exame the deme index numbering
# all "any" tests here should return "FALSE"
demeIndexes <- read.csv("TestMakeDemesIndexes.txt", header = FALSE)
length(unique(demeIndexes$V1))
length(unique(demeIndexes$V2))
errs <- rep(FALSE, length(demeIndexes$V2))
for ( i in 2:length(demeIndexes$V2) ) {
  if ( demeIndexes$V2[i] < demeIndexes$V2[i-1] ) {
    print("Error!")
    errs[i] <- 1
  }
}
any(errs)
any(demeIndexes$V1 >= length(demeIndexes$V2))
any(demeIndexes$V1 < 0)

#check sums that should work out
segSitesTS <- read.csv("SegregatingSitesTS.csv")
nsamps <- dim(segSitesTS)[1]
for ( i in 1:nsamps ) {
  if ( segSitesTS$nSegregatingSites[i] != (segSitesTS$nNeutralSites[i] + segSitesTS$nBackgroundSelSites[i] + segSitesTS$nPositiveSelSites[i] + segSitesTS$nDivergentSelSites[i]) ) {
    print("Error in Seg Sites data!")
  } else {
    print(i)
  }
}

SFSdata <- read.csv("SFStimeSeries.csv")
nsamps <- dim(SFSdata)[1]
for ( i in 1:nsamps ) {
  if ( SFSdata$NumberOfSites[i] != (SFSdata$NumNeutralSites[i] + SFSdata$NumBGSsites[i] + SFSdata$NumPOSsites[i] + SFSdata$NumDIVsites[i]) ) {
    print("Error in SFSdata!")
  }
}
