extractFromNormData <-
function(normData, file="extractedNormData.Rda") {
	extractedNormData <- returnSummaryNormData(normData)
	save(extractedNormData, file=file)
	return(extractedNormData)
}



returnSummaryNormData <-    function(normData){
	
    require(minfi)
    require(IlluminaHumanMethylation450kmanifest)
    require(matrixStats)
    require(gmodels)


	meth <- normData$methMatrix
	unmeth <-normData$unmethMatrix
	beta <- meth / (meth+unmeth+100)
	m <- log2( (meth+1)/(unmeth+1))
                              

	# Defining the Type I, II Green and II Red probes:
	probesI <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "I")

	probesII <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "II")
	        
	### Chr probes:
	data(shinymethylAnnotation, package="shinyMethyl")

	        
	probesIGrn <- intersect( probesI$Name[probesI$Color=="Grn"], autosomal )
	probesIRed   <- intersect( probesI$Name[probesI$Color=="Red"], autosomal )
	probesII     <- intersect( probesII$Name, autosomal )

	uProbeNames <- rownames(beta)
	uProbesIGrn <- intersect(uProbeNames, probesIGrn)
	uProbesIRed <- intersect(uProbeNames, probesIRed)
	uProbesII   <- intersect(uProbeNames, probesII)
	uProbesX <- intersect(uProbeNames, chrX)
	uProbesY <-  intersect(uProbeNames, chrY)
	indicesIGrn <- match(uProbesIGrn, uProbeNames)
	indicesIRed <- match(uProbesIRed, uProbeNames)
	indicesII   <- match(uProbesII, uProbeNames)
    indicesX <- match(uProbesX, uProbeNames)
	indicesY <- match(uProbesY, uProbeNames)
	 
	 indList <- list(indicesIGrn, indicesIRed, indicesII, indicesX, indicesY)
	 names(indList) <- c("IGrn", "IRed", "II","X","Y")
	
	
	# Extraction of the quantiles
	mQuantiles               <- vector("list",5)
	betaQuantiles            <- vector("list", 5)
	methQuantiles          <- vector("list", 5)
	unmethQuantiles        <- vector("list", 5)
	names(mQuantiles)        <- c("IGrn", "IRed", "II","X","Y")
	names(betaQuantiles)     <- c("IGrn", "IRed", "II","X","Y")
	names(methQuantiles)   <- c("IGrn", "IRed", "II","X","Y")
	names(unmethQuantiles) <- c("IGrn", "IRed", "II","X","Y")
	
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
		    

	}
	




	pcaInfo <- returnPCScores(beta, 5000 )

	return(list(
		mQuantiles = mQuantiles,
		betaQuantiles = betaQuantiles,
		methQuantiles = methQuantiles,
		unmethQuantiles = unmethQuantiles,
		pcaInfo = pcaInfo,
        pd  = NULL))
}
