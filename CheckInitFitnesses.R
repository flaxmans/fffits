rm(list = ls())
setwd("~/Documents/Research/SpeciationGenomics/fffits/")
fitness <- read.csv("InitialFitnessValues.csv")
boxplot(fitness$fitness ~ fitness$location)