returnPCScores <-
function(beta, numPositions, XYExcluded = TRUE){
		library(matrixStats)
		library(gmodels)
		library(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
		locations <- getLocations(IlluminaHumanMethylation450kannotation.ilmn.v1.2)
		chrY <- names(locations[seqnames(locations)=="chrY"])
		chrX <- names(locations[seqnames(locations)=="chrX"])
		XY <- union(chrY,chrX)
		
		if (XYExcluded & sum(is.na(match(XY,rownames(beta))))!=length(XY)  ){
			b <- beta[-match(XY,rownames(beta)),]
		} else {
			b <- beta
		}
		o <- order(-rowVars(b))[1:numPositions]
		pca <- fast.prcomp(t(b[o,]))
		scores <- pca$x
		PCpercs <- (pca$sdev^2)/(sum(pca$sdev^2))*100
		return(list(scores=scores, PCpercs = PCpercs ))
	}
