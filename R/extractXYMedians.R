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
