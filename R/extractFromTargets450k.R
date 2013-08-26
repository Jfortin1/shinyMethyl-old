extractFromTargets450k <-
function(targets, file="extractedData.Rda", bySubsets = FALSE) {
	
	extractedData = NULL
	if (bySubsets==FALSE){
		extractedData <- returnSummary450k(targets)
	} else {
		start <- seq(1,nrow(targets),10)
		end   <- seq(10,nrow(targets),10)
		if (nrow(targets) %% 10 != 0){
			end <- c(end, nrow(targets))
		}
		for (k in 1:length(start)){
			print(k)
			currentTargets <- targets[start[k]:end[k],]
			extractedData <- mergeExtractedData(extractedData, returnSummary450k(currentTargets))
		}
	}
	
	save(extractedData, file=file)
	return(extractedData)	
}



mergeExtractedData <- function(extractedDataPart1, extractedDataPart2){
	if (is.null(extractedDataPart1)){
		return(extractedDataPart2)
	} else {
		newExtractedData <- vector("list", length(extractedDataPart2))
		names(newExtractedData) <- names(extractedDataPart2)
		for (i in 1:13){
		newExtractedData[[i]] <- mergeListOfMatrices(extractedDataPart1[[i]], extractedDataPart2[[i]], by= "column")
	}
		for (i in 14:14){
		 newExtractedData[[i]] <- mergeListOfVectors(extractedDataPart1[[i]], extractedDataPart2[[i]])
	 }
	for (i in 15:15){
		 newExtractedData[[i]] <- rbind(extractedDataPart1[[i]], extractedDataPart2[[i]])
	 }
	return(newExtractedData)
	}
}




mergeListOfVectors <- function(listPart1, listPart2){
	newList <- vector("list", length(listPart2))
	names(newList) <- names(listPart2)
		for (k in 1:length(listPart2)){
        	newList[[k]] <- c(listPart1[[k]], listPart2[[k]])
        }
	return(newList)
}

mergeListOfMatrices <- function(listPart1, listPart2, by ="column"){
	newList <- vector("list", length(listPart2))
	names(newList) <- names(listPart2)
	if (by== "column"){
        for (k in 1:length(listPart2)){
        	newList[[k]] <- cbind(listPart1[[k]], listPart2[[k]])
        }
	} else {
		for (k in 1:length(listPart2)){
        	newList[[k]] <- rbind(listPart1[[k]], listPart2[[k]])
        }
	}
	return(newList)
}


returnSummary450k <-
function(targets){
	
    library(minfi)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
	
	controlType <- c("BISULFITE CONVERSION I",
	"BISULFITE CONVERSION II",
	"EXTENSION",
	"HYBRIDIZATION",
	"NEGATIVE",
	"NON-POLYMORPHIC",
	"NORM_A",
	"NORM_C",
	"NORM_G",
	"NORM_T",
	"SPECIFICITY I",
	"SPECIFICITY II",
	"TARGET REMOVAL",
	"STAINING")
	
	RGSet <- read.450k.exp(targets = targets)
	r <- getRed(RGSet)
   g <- getGreen(RGSet)
	MSet.raw <- preprocessRaw(RGSet)
	meth <- getMeth(MSet.raw) 
	unmeth <-getUnmeth(MSet.raw)
	beta <- getBeta(MSet.raw)
	m <- getM(MSet.raw)
	cn <- meth + unmeth
	pd <- pData(RGSet)
	
	## Extraction of the controls
	greenControls=vector("list",length(controlType))
	redControls=vector("list",length(controlType))
	names(greenControls)=controlType
	names(redControls)=controlType
	
	for (i in 1:length(controlType)){
		if (controlType[i]!="STAINING"){
			ctrlAddress <- getControlAddress(
		        RGSet, controlType = controlType[i])
		} else {
			ctrlAddress <- getControlAddress(
		        RGSet, controlType = controlType[i])[c(2,3,4,6)]
		}
		redControls[[i]]=r[ctrlAddress,]
		greenControls[[i]]=g[ctrlAddress,]
		  
	}

	# Extraction of the undefined negative control probes
	locusNames <- getManifestInfo(RGSet, "locusNames")
	TypeI.Red <- getProbeInfo(RGSet, type = "I-Red")
   TypeI.Green <- getProbeInfo(RGSet, type = "I-Green")
    
   numberQuantiles <- 100
	probs <- 1:numberQuantiles/100
	
	
    greenNegativePooledA <- apply(getGreen(RGSet)[TypeI.Red$AddressA,], 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
    greenNegativePooledB <- apply(getGreen(RGSet)[TypeI.Red$AddressB,], 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
                            
    redNegativePooledA   <- apply(getRed(RGSet)[TypeI.Green$AddressA,], 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
    redNegativePooledB   <- apply(getRed(RGSet)[TypeI.Green$AddressB,], 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
    
    oobControls <- list(greenNegativePooledA,greenNegativePooledB,redNegativePooledA,redNegativePooledB) 
    names(oobControls) <- c("oobGreenA","oobGreenB","oobRedA","oobRedB")
                              

	# Defining the Type I, II Green and II Red probes:
	probesI <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "I")

	probesII <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "II")
	        
	### Chr probes:
	data(shinymethylAnnotation, package="shinyMethyl")
	# locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
	# autosomal <- names(locations[seqnames(locations) %in% paste0("chr", 1:22)])
	# chrY <- names(locations[seqnames(locations)=="chrY"])
	# chrX <- names(locations[seqnames(locations)=="chrX"])
   
	        
	probesIGrn <- intersect( probesI$Name[probesI$Color=="Grn"], autosomal )
	probesIRed   <- intersect( probesI$Name[probesI$Color=="Red"], autosomal )
	probesII     <- intersect( probesII$Name, autosomal )
	uProbeNames <- rownames(beta)
	uProbesIGrn <- intersect(uProbeNames, probesIGrn)
	uProbesIRed <- intersect(uProbeNames, probesIRed)
	uProbesII   <- intersect(uProbeNames, probesII)
	indicesIGrn <- match(uProbesIGrn, uProbeNames)
	indicesIRed <- match(uProbesIRed, uProbeNames)
	indicesII   <- match(uProbesII, uProbeNames)
	indicesX <- match(chrX, uProbeNames)
	indicesY <- match(chrY, uProbeNames)
	
	indList <- list(indicesIGrn, indicesIRed, indicesII, indicesX, indicesY)
	names(indList) <- c("IGrn", "IRed", "II","X","Y")
	
	
	# Extraction of the quantiles
	mQuantiles               <- vector("list",5)
	betaQuantiles            <- vector("list", 5)
	methQuantiles          <- vector("list", 5)
	unmethQuantiles        <- vector("list", 5)
	cnQuantiles     <- vector("list", 5)
	names(mQuantiles)        <- c("IGrn", "IRed", "II","X","Y")
	names(betaQuantiles)     <- c("IGrn", "IRed", "II","X","Y")
	names(methQuantiles)   <- c("IGrn", "IRed", "II","X","Y")
	names(unmethQuantiles) <- c("IGrn", "IRed", "II","X","Y")
	names(cnQuantiles) <- c("IGrn", "IRed", "II","X","Y")
	
	nq <- 500
	probs <- seq(0,1,1/(nq-1))
	
	for (i in 1:5){
		mQuantiles[[i]] <- 
		   apply(m[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		    )
		    
		betaQuantiles[[i]] <- 
		   apply(beta[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		    )
		    
		methQuantiles[[i]] <- 
		   apply(meth[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		    )
		    
		unmethQuantiles[[i]] <- 
		   apply(unmeth[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		    )
		    
		 cnQuantiles[[i]] <- 
		   apply(cn[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		    )
	}
	
    # Extraction of the densities
	mDensities               <- vector("list",5)
	betaDensities            <- vector("list", 5)
	methDensities          <- vector("list", 5)
	unmethDensities        <- vector("list", 5)
	cnDensities             <- vector("list", 5)
	names(mDensities)        <- c("IGrn", "IRed", "II","X","Y")
	names(betaDensities)     <- c("IGrn", "IRed", "II","X","Y")
	names(methDensities)   <- c("IGrn", "IRed", "II","X","Y")
	names(unmethDensities) <- c("IGrn", "IRed", "II","X","Y")
	names(cnDensities) <- c("IGrn", "IRed", "II","X","Y")

	for (i in 1:5){
		mDensities[[i]] <- 
		   apply(m[indList[[i]], ], 2, 
		       function(x) density(x, na.rm=T, bw=0.1, n = 500)$y
		    )
		    
		betaDensities[[i]] <- 
		   apply(beta[indList[[i]], ], 2, 
		       function(x) density(x, na.rm=T, bw=0.01, n = 500, from=0, to=1)$y
		 )
		    
		methDensities[[i]] <- 
		   apply(meth[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		)
		    
		unmethDensities[[i]] <- 
		   apply(unmeth[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		)
		
			cnDensities[[i]] <- 
		   apply(cn[indList[[i]], ], 2, 
		       function(x) quantile(x, probs=probs, na.rm=T)
		)
	}
	

   XYMedians <- extractXYMedians(RGSet)


	return(list(
		mQuantiles = mQuantiles,
		betaQuantiles = betaQuantiles,
		methQuantiles = methQuantiles,
		unmethQuantiles = unmethQuantiles,
		cnQuantiles = cnQuantiles, 
		mDensities = mDensities,
		betaDensities = betaDensities,
		methDensities = methDensities,
		unmethDensities = unmethDensities,
		cnDensities = cnDensities,
		greenControls = greenControls,
		redControls = redControls,
		oobControls = oobControls,
		XYMedians = XYMedians,
		pd = pd))
}






extractXYMedians <-
function(RGSet){
		library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
		library(minfi)
		locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
		chrY <- names(locations[seqnames(locations)=="chrY"])
		chrX <- names(locations[seqnames(locations)=="chrX"])
		rawSet <- preprocessRaw(RGSet)
		M <- getMeth(rawSet)
		U <- getUnmeth(rawSet)
		MX <- M[match(chrX,rownames(M)),]
		MY <- M[match(chrY,rownames(M)),]
		UX  <- U[match(chrX,rownames(U)),]
		UY  <- U[match(chrY,rownames(U)),]
		medianXU <-  apply(UX,2,function(x) median(x,na.rm=T))
		medianXM <-  apply(MX,2,function(x) median(x,na.rm=T))
		medianYU <- apply(UY,2,function(x) median(x,na.rm=T))
		medianYM <- apply(MY,2,function(x) median(x,na.rm=T))
		return(list(medianXU = medianXU,
		    medianXM = medianXM,
		    medianYU = medianYU,
		    medianYM = medianYM)
		)
	}


