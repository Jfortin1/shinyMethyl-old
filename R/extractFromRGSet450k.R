extractFromRGSet450k <-
function(RGSet, file="extractedData.Rda") {
	extractedData <- returnSummaryRGSet450k(RGSet)
	save(extractedData, file=file)
	return(extractedData)
}



returnSummaryRGSet450k <-
function(RGSet){
	
    require(minfi)
    require(IlluminaHumanMethylation450kmanifest)
    require(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
    require(matrixStats)
    require(gmodels)

	
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
	
	MSet.raw <- preprocessRaw(RGSet)
	r <- getRed(RGSet)
    g <- getGreen(RGSet)
	meth <- getMeth(MSet.raw) 
	unmeth <-getUnmeth(MSet.raw)
	beta <- getBeta(MSet.raw)
	m <- getM(MSet.raw)
	cn <- meth + unmeth
	pd <- pData(RGSet)
	ann <- IlluminaHumanMethylation450kannotation.ilmn.v1.2

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
	
	greenOOB <- rbind(getGreen(RGSet)[TypeI.Red$AddressA,], getGreen(RGSet)[TypeI.Red$AddressB,])
    redOOB <- rbind(getRed(RGSet)[TypeI.Green$AddressA,], getRed(RGSet)[TypeI.Green$AddressB,])
    
    greenOOB <- apply(greenOOB, 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
   
    redOOB <- apply(redOOB, 2,
                               function(x)  quantile(x, probs=probs, na.rm=T)
                            )
                            
                            
     oob <- list(greenOOB = greenOOB,  redOOB = redOOB)            
   
   
   
                                 

	# Defining the Type I, II Green and II Red probes:
	probesI <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "I")

	probesII <- getProbeInfo(
	    IlluminaHumanMethylation450kmanifest,
	        type = "II")
	        
	### Chr probes:
	data(shinymethylAnnotation, package="shinyMethyl")
	# locations <- getLocations(ann)
	# autosomal <- names(locations[as.vector(seqnames(locations) %in% paste0("chr", 1:22))])
	# chrY <- names(locations[as.vector(seqnames(locations)=="chrY")])
	# chrX <- names(locations[as.vector(seqnames(locations)=="chrX")])
	        
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
	

	
	medianXU <-  unmethQuantiles$X[250,]
		medianXM <-  methQuantiles$X[250,]
		medianYU <- unmethQuantiles$Y[250,]
		medianYM <- methQuantiles$Y[250,]


	
	XYMedians <- list(medianXU = medianXU,
		    medianXM = medianXM,
		    medianYU = medianYU,
		    medianYM = medianYM)


	pcaInfo <- returnPCScores(beta, 5000 )

	return(list(
		mQuantiles = mQuantiles,
		betaQuantiles = betaQuantiles,
		methQuantiles = methQuantiles,
		unmethQuantiles = unmethQuantiles,
		cnQuantiles = cnQuantiles,
		greenControls = greenControls,
		redControls = redControls,
		XYMedians = XYMedians,
		oob = oob,
		pcaInfo = pcaInfo,
		pd = pd))
}
