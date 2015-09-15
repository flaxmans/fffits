# script to look at global allele frequencies over time

# allele frequency time series, global
setwd("~/fffitsdir/")
data <- read.csv("AlleleFreqTS.csv")
uniqueSites <- unique(data$SiteIndex)
firstCall <- TRUE
for ( i in 1:length(uniqueSites) ) {
    site <- uniqueSites[i]
    indexes <- intersect(which(data$SiteIndex %in% site), which(data$SiteClassCode == 0))
    if ( length(indexes) > 1 ) {
      time <- data$Time[indexes]
      counts <- data$DerivedAlleleFreq[indexes]
      if ( firstCall ) {
        plot(time, counts, type = "l", xlim = c(0,max(data$Time)), ylim = c(0,max(data$DerivedAlleleFreq)))
        firstCall <- FALSE
      } else {
        points(time, counts, type = "l")
      }
    }
}

# site frequency spectra, global
SFSdata <- read.csv("SFStimeSeries.csv")
uniqueTimes <- unique(SFSdata$Time)
firstCall <- TRUE
for ( i in 1:length(uniqueTimes) ) {
  time <- uniqueTimes[i]
  indexes <- which(SFSdata$Time %in% time)
  Count <- log10(SFSdata$DerivedAlleleCopyNumber[indexes])
  Frequency <- log10(SFSdata$NumberOfSites[indexes])
  if ( firstCall ) {
    plot(Count,Frequency, type = "l", 
         xlim = c(0,log10(max(SFSdata$DerivedAlleleCopyNumber))), 
         ylim = c(0,log10(max(SFSdata$NumberOfSites))), 
         xlab = "log10(count)", ylab = "log10(frequency)")
    firstCall <- FALSE
  } else {
    points(Count, Frequency, type = "l")
  }
}