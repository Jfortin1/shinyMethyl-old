returnPCScores <-
function(beta, numPositions, XYExcluded = TRUE){
		library(matrixStats)
		library(gmodels)
		data(shinymethylAnnotation, package="shinyMethyl")
				
		uXY <- intersect(chrX, rownames(beta))
		
		if (XYExcluded & length(uXY)!=0 ){
			b <- beta[-match(uXY,rownames(beta)),]
		} else {
			b <- beta
		}
		o <- order(-rowVars(b))[1:numPositions]
		pca <- fast.prcomp(t(b[o,]))
		scores <- pca$x
		PCpercs <- (pca$sdev^2)/(sum(pca$sdev^2))*100
		return(list(scores=scores, PCpercs = PCpercs ))
	}
