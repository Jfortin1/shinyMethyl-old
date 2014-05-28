# setwd("./HudsonAnalysis/")
# load("RGSetHudson.Rda")


# library(minfi)
# library(IlluminaHumanMethylation450kmanifest)
# library(matrixStats)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(gmodels)


# RGSet <- updateObject(RGSetHudson)
	
# hudsonSummary <- shinySummarize(RGSet)


# load("/Users/Jean-Philippe/Desktop/hudsonSummary.Rda")
# source("~/shinyMethyl2/R/shinyMethylSet.R")


# source("~/shinyMethyl2/R/plotControlProbes.R")
# source("~/shinyMethyl2/R/plotDensities.R")
# source("~/shinyMethyl2/R/plotDesign.R")
# source("~/shinyMethyl2/R/plotPCA.R")
# source("~/shinyMethyl2/R/plotQC.R")
# source("~/shinyMethyl2/R/plotSex.R")
# source("~/shinyMethyl2/R/runShinyMethyl.R")


# # shinyMethylSet <<- orderByName(hudsonSummary)
# library(shiny)

# shinyMethylSet <<- orderByName()


 


# # sampleNames <- shinyMethylSet@sampleNames
# pd <- data.frame(plate=substr(sampleNames,1,6), position = substr(sampleNames,12,17))
# rownames(pd) <- sampleNames

# Gender <- returnPredictedGender(cnQuantiles=shinyMethylSet@cnQuantiles)
# Gender[c(1,4,8,20,2,18)] <- "M"
# pd$Gender <- Gender

# shinyMethylSet@phenotype <- pd

# runApp("~/Desktop/shinyMethyl") 
	
	
source("~/shinyMethyl2/R/shinyMethylSet.R")
source("~/shinyMethyl2/R/plotControlProbes.R")
source("~/shinyMethyl2/R/plotDensities.R")
source("~/shinyMethyl2/R/plotDesign.R")
source("~/shinyMethyl2/R/plotPCA.R")
source("~/shinyMethyl2/R/plotQC.R")
source("~/shinyMethyl2/R/plotSex.R")
source("~/shinyMethyl2/R/runShinyMethyl.R")

	
load("/Users/Jean-Philippe/Desktop/shinyMethyl/tcga.summary.Rda")
load("/Users/Jean-Philippe/Desktop/shinyMethyl/covariates.tcga.Rda")
load("/Users/Jean-Philippe/Desktop/norm.chip.summary.Rda")
tcga.summary@phenotype <- covariates.tcga
norm.chip.summary@phenotype <- covariates.tcga
shinyMethylSet <<- orderByName(tcga.summary)
shinyMethylSet.norm <<- orderByName(norm.chip.summary)
library(shiny)
runApp("~/Desktop/shinyMethyl") 
	
# TO DO LIST: 

	
# Change the controls so that they correspond to FunNorm
# Fix NA problem for design plot. 
# Make automated report using knitR 
# Make test functions. 
# Add cutoffs for quality control plot 
# Add a resume for the dataset:
# Add Type I/ Type II bias ..
# Add CpG island annotation.. 
# Clinical_M: Normalization destroys one curve. Intersting. Oh: only one sample for a group... to examine!
# For drawing densities, make sure all samples from group 1 draw first, all samples from group 2 draw second... etc
# For bad quality samples, add a download button.
# On-the-fly colors  
# Need to add legends to densities plots. 
# Position color gradient plot
 


# To put in documentation:
# CTRL to refresh the browser
# Mapping to genome.
# Using Illumina 450k manifest
# Raw processing
# Quality control plot: line +- 0.5 only works for blood samples
# rm(list=ls()) To make sure to clear the session when playing with new shinyMethylObject

# To do in knitr:
# Make a function so that real densities will be plotted.
# Number of pos: 20000


beta <- shinyMethylSet@betaQuantiles$II

returnFullDensity <- function(quantiles){
	n = 500
	for (i in 1:ncol(quantiles)){
		target.case <- sapply(1:499, function(j) {
            start <- quantiles[j,i]
            end <- quantiles[j+1,i]
            sequence <- seq(start, end,( end-start)/n)[-n]
            return(sequence)
        })
	print(i)
	}
}

library(knitr)
library(knitrBootstrap)
knit_bootstrap("makereport.Rmd",chooser=c('boot', 'code'))
  
        
        
        