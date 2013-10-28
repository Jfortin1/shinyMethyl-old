

plotDesign450k <- function(sampleNames, covariates, legendTitle){
	    palette(brewer.pal(8,"Set1"))
	     
	    boxColors <- buildColors(sampleNames, covariates)
	    nLevels <- length(unique(boxColors))
	    naPos <- boxColors=="NA"
	    boxColors <- factor(boxColors)
	    oldLevels <- levels(boxColors)
	    levels(boxColors) <- 1:(nLevels)
	    boxColors <- as.character(boxColors)
	    if (sum(naPos)>0){
	    	boxColors[naPos] <- "white"
	    }
	    
		numberOfChips <- length(unique(substr(sampleNames,1,10)))
		nrows <- ceiling(numberOfChips/8)
		par(pin=c(width=5,heigth=5),mai=c(0,0,1,0))
		plot(0:(nrows+1),seq(0,nrows,nrows/(nrows+1)),type="n",axes=FALSE,xlab="",ylab="", main = "Illumina HumanMethylation 450K Array Design")
		
		for (i in 1:nrows){
			plotRow(leftMargin= 0, plateytop=(i+0.75)-1, plateybottom=i-1, color = boxColors[(96*(i-1)+1):(96*i)])
		}
		legendColors <- 1: (nLevels)
		legendPch <- rep(15, nLevels)
		if (sum(naPos)>0){
	    	legendColors[oldLevels=="NA"] <- "black"
	    	legendPch[oldLevels=="NA"] <- 0
	    }
}
	



plotLegendDesign450k <- function(sampleNames, covariates, legendTitle){
	palette(brewer.pal(8,"Set1"))
	    
	    boxColors <- buildColors(sampleNames, covariates)
	    nLevels <- length(unique(boxColors))
	    naPos <- boxColors=="NA"
	    boxColors <- factor(boxColors)
	    oldLevels <- levels(boxColors)
	    levels(boxColors) <- 1:(nLevels)
	    boxColors <- as.character(boxColors)
	    if (sum(naPos)>0){
	    	boxColors[naPos] <- "white"
	    }

        legendColors <- 1: (nLevels)
		legendPch <- rep(15, nLevels)
        
        if ("NA" %in% oldLevels){
        	legendColors[oldLevels=="NA"] <- "white"
        }

        plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
        legend("topright",
		     pch= legendPch,
		     col = legendColors,
		     oldLevels,
		     cex=1.5,
		     title= legendTitle)  
}






## To build the color palette for the rectangles
buildColors <- function(sampleNames, colCovariate){
	        colCovariate[is.na(colCovariate)] <- "NA"
			colCovariate <- as.character(colCovariate)
		    chipInfo <- substr(sampleNames,1,10)
		    rowInfo  <- substr(sampleNames,14,14)
		    columnInfo <- substr(sampleNames,17,17)
		    numberOfChips <- length(unique(chipInfo))
		    boxColors <- rep("NA",numberOfChips*12)
		 
		for (i in 1:numberOfChips){
				for (columnIndex in 1:2){
					for (rowIndex in 1:6){
						indices <- intersect(which(chipInfo == unique(chipInfo)[i]), which(columnInfo == columnIndex))
						indices <- intersect(indices, which(rowInfo == rowIndex))
						if (!length(indices)==0){
							boxColors[12*(i-1)+1+6*(columnIndex-1)+(rowIndex-1)] <- colCovariate[indices]
						}
					}
				}	
			}
			
			return(boxColors)
}



    
    
    




### To plot a row of plates
plotRow <- function(leftMargin, plateytop, plateybottom, color){
		unit <- (plateytop - plateybottom)/6
		rightMargin = 1+23*unit
		platePosition <- seq(leftMargin,rightMargin, (rightMargin-leftMargin)/(7))
		for (j in 1:length(platePosition)){
		    plotPlate(plateybottom=plateybottom,plateytop=plateytop, platePosition[j], color[(12*(j-1)+1):(12*j)])
		}
}









### To plot a single array
plotPlate <- function(plateybottom, plateytop, platexleft, color){
		platexright <- platexleft + (plateytop - plateybottom)/3
		rect(ybottom = plateybottom,
		        ytop = plateytop, 
		            xleft= platexleft,
		                xright = platexright, 
		                    lwd=2,
		                        col=0)
		                        
		startsColumn1 <- sort(seq(plateybottom, plateytop, (plateytop-plateybottom)/6), decreasing=T)[-1]
		endsColumn1 <-  sort(seq(plateybottom, plateytop, (plateytop-plateybottom)/6), decreasing=T)[-7]
		
				### For plotting the left column:
			for (k in 1:length(startsColumn1)){
				rect(ybottom = startsColumn1[k], 
				       ytop=endsColumn1[k], 
				           xleft = platexleft, 
				               xright = (platexright+platexleft)/2, 
				                   col=color[k],
				                       lwd=3)
			  }
			
			### For plotting the right column:
			for (v in 1:length(startsColumn1)){
				rect(ybottom = startsColumn1[v], 
				       ytop=endsColumn1[v],
				          xleft = (platexleft+platexright)/2, 
				             xright = platexright, 
				                col = color[6+v],
				                   lwd=3)
			}
}












### To plot the internal controls
plotInternalControls <- function(y, col, main, sampleNames,
                                                 selectedSlides = NULL, 
                                                     slideNames = NULL){
	n <- length(y)
	
	plot(1:n, y, pch=20, cex=0.7,
	    main = main,
	       xlab = "Sample",
	           ylab = "Control intensity",
	               ylim = c(0, 2*max(y)),
	                   col = col)
	
	if (!is.null(selectedSlides)){
		m <- length(selectedSlides)
		for (i in 1:m){
			currentIndices = which( slideNames == selectedSlides[i] )
			
			points(currentIndices, 
			     y[currentIndices],
			         col = col[ currentIndices ],
			            cex = 1.5,
			                pch = 8)
		}
	}
	
	grid()
	abline(h=0, lty=3)
}

   
   
   
   
   
   
   
### To add the selected array to the internal controls plot
addHoverPoints <- function(y, sampleNames, xSelected, ySelected){
    	
    	if (!is.null(xSelected) && !is.null(ySelected)){
    		n = length(y)
		    xDiff <- (  (1:n) - rep(xSelected , n)   ) / n 
	        yDiff <- (y - ySelected) / (2*max(y))
	        clickedIndex <- which.min(xDiff^2 + yDiff^2)
	    
	    	points(clickedIndex, y[clickedIndex], pch = 1, cex = 4, lwd = 3)
		    points(clickedIndex, y[clickedIndex], pch = 20, cex = 1)
		    selectedName <- 
		        paste0( "Selected array:", as.character(sampleNames[clickedIndex]))
		    legend("topleft", cex = 1, pch = 1, legend = selectedName)
		    legend("topleft", cex = 1, pch= 20, legend = selectedName)
	    }
	    
	}












### Plot of the densities
densitiesPlot <- function(quantiles, main, xlab , xlim, ylim, bw, 
                            slideNames, solidLine = TRUE, mean = TRUE, 
                            slides = NULL, col, from, to ){
                            	
    batches= slides
    meanSample = apply(quantiles,1,mean)
    
    plot( density(meanSample, bw=bw) , 
                   main = main, 
                   ylab = "Density", 
                   xlab = xlab,
                   col = "white",
                   xlim = xlim, 
                   ylim = ylim
    )
                 
                 
      if (solidLine){
			lty = 1
		} else  {lty = 3}           
                                   
      if (mean){
          if (length(slides) != 0){
              names  = colnames(quantiles)
              batch  = substr(names, 1, 10)
              numberBatches = length(unique(batch))
              for (i in 1:length(batches)){
	  			    currentIndices=which(slideNames==batches[i])
	  			    currentColors <- col[which(slideNames==batches[i])]
	  				if (length(currentIndices==1)){
	  			    	currentMean=quantiles[,currentIndices]
	  			    } else {
	  			    	currentMean=apply(quantiles[,currentIndices],1,mean)
	  			    }
	  			    lines(density(currentMean, bw = bw, from= from, to = to), lwd=1.5, lty=lty,col=currentColors)
				}
				
				
          }
      } else  {
      	if (!is.null(slides)){
  				names = colnames(quantiles)
				batch = substr(names,1,10)
				numberBatches=length(unique(batch))	
	  			for (i in 1:length(batches)){
		  			currentIndices <- which(slideNames==batches[i])
		  			currentColors <- col[currentIndices]
		  			for (j in 1:length(currentIndices)){
		  				lines(density(quantiles[,currentIndices[j]], bw=bw, from=from, to=to),
		  				     col= currentColors[j], lwd = 1, lty=lty)
		  			}
	  			}
			  }
      }            
      grid()                                 	
}
	










	
	
plotQC <- function(unmethQuantiles, methQuantiles,
                       slides, slideNames, col){
	    palette(brewer.pal(8,"Set1"))
	    
	    mediansU <- unlist(unmethQuantiles[[2]][250,])        
	    mediansM <- unlist(methQuantiles[[2]][250,])
	    batches <- slides   
	    
	     xlim1 <- min(log2(mediansU))-0.2*(max(log2(mediansU))-min(log2(mediansU)))
         xlim2 <- max(log2(mediansU))+0.2*(max(log2(mediansU))-min(log2(mediansU)))
         ylim1 <- min(log2(mediansM))-0.2*(max(log2(mediansM))-min(log2(mediansM)))
         ylim2 <- max(log2(mediansM))+0.2*(max(log2(mediansM))-min(log2(mediansM)))
         
         plot(log2(mediansU), 
             log2(mediansM),
                xlim= c(xlim1, xlim2),
                    ylim=c(ylim1, ylim2),
                       cex=1,
                           pch=20,
                              col = col,
                                 main = "QC Plot")
         
         if (!is.null(slides)){
	    	numberSubsets=length(as.numeric(slides))
			 for (i in 1:numberSubsets){
	  		     currentIndices=which(slideNames == batches[i])
	  		     
	  			 points(log2(mediansU[currentIndices]), 
	  			                   log2(mediansM[currentIndices]),
	  			                       col=col[currentIndices], cex=1.5, pch=8
	  			 )
	  		}
	    }                     
  	
  		grid()
  	    
  		
}
	
	
	
	
	
	
	
	


addHoverQC <- function(y , 
                                      xSelected = NULL,
                                          ySelected = NULL,
                                               sampleNames, 
                                                   unmethQuantiles,
                                                       methQuantiles,
                                                            slides){
	     mediansU <- unlist(unmethQuantiles[[2]][250,]) 
	     mediansM <- unlist(methQuantiles[[2]][250,])
	     batches <- slides  
	     
	     n <- length(y)
	     xDiff <- ((1:n)-rep(xSelected, n)) / n
	     yDiff <- (y - ySelected)/(2*max(y))
		 clickedIndex <- which.min(xDiff^2 + yDiff^2)
		 
		 if (!is.null(xSelected) && !is.null(ySelected)){
			points(log2(mediansU[clickedIndex]), 
	  			            log2(mediansM[clickedIndex]),
	  			                col="black", 
	  			                     cex=4,
	  			                         pch=1,
	  			                             lwd=2
	  	  )
	  	  
	  	  points(log2(mediansU[clickedIndex]), 
	  			            log2(mediansM[clickedIndex]),
	  			                col="black",
	  			                    cex=2,
	  			                        pch=20
	  	  )
	  	  }
	  	 abline(v = log2(mediansU[clickedIndex]), lty = 3, lwd = 2, col = "black")
  	     abline(h = log2(mediansM[clickedIndex]),lty = 3, lwd = 2, col = "black")
  	  	
}
     
     
     
     
     
     














returnControlStat <- function(controlType, greenControls, redControls, controlNames){
	   index=match(controlType, controlNames)
	   greenControls.current <- greenControls[[index]]
	   redControls.current     <- redControls[[index]]
	   
		    if (controlType=="BISULFITE CONVERSION I"){
		  		controlStat <- (apply(greenControls.current[1:3,],2,mean) +
		  		       apply(redControls.current[7:9,],2,mean))/2 
		  	} else { 
		  	
		  	if (controlType=="BISULFITE CONVERSION II"){
		  		controlStat <- apply(redControls.current,2,mean)
		  	} else {
		  		
		  	if (controlType=="EXTENSION"){
		  		controlStat <- (apply(greenControls.current[3:4,],2,mean) +
		  		       apply(redControls.current[1:2,],2,mean))/2 
		    } else {
		    	
		  	if (controlType=="NEGATIVE"){
		  		controlStat <- (apply(greenControls.current,2,mean) +
		  		       apply(redControls.current,2,mean))/2 
		    }else {
		    	
		    if (controlType=="HYBRIDIZATION"){
		    	controlStat <- apply(greenControls.current,2,mean)	
		    } else {controlStat <- apply(greenControls.current,2,mean)
		    }}}}}
		    return(controlStat)
} 	











addHoverDensity <- function(y , 
                                      xSelected = NULL,
                                          ySelected = NULL,
                                               sampleNames, 
                                                            slides,
                                                                quantiles, bw, col){
		     batches <- slides  
		     
		     n <- length(y)
		     xDiff <- ((1:n)-rep(xSelected, n)) / n
		     yDiff <- (y - ySelected)/(2*max(y))
			 clickedIndex <- which.min(xDiff^2 + yDiff^2)
			 
			 if (!is.null(xSelected) && !is.null(ySelected)){
				lines(density(quantiles[,clickedIndex],bw=bw),
			  				      col=col[clickedIndex],lwd = 3)
		  	  }
}      











returnPredictedGender <- function(cutoff = (-3), XYMedians){
	
	
	x = XYMedians$medianXU + XYMedians$medianXM
	y = XYMedians$medianYM + XYMedians$medianYU
	diff <- log2(y)-log2(x)
	n <- length(diff)
	
	predictedGender = rep("M", n)
	predictedGender[which(diff < cutoff)] <- "F"     
	names(predictedGender) <- names(XYMedians$medianXU)   
	
	return(predictedGender)
}







plotPredictedGender <- function(cutoff = (-3),
                                                              XYMedians){
	palette("default")
	x <- XYMedians$medianXU + XYMedians$medianXM
	y <- XYMedians$medianYM + XYMedians$medianYU
	diff <- log2(y)-log2(x)
	
	n <- length(x)
	color = rep("lightskyblue", n)
	color[which(diff < cutoff)] <- "darkorange1"
	
	plot(diff, jitter(rep(0,n),factor=1.5), 
	          ylim = c(-1,1), pch=18, cex=2,col= color, 
	              yaxt="n", xlab="median CN(Y)  - median CN(X)", 
	                 ylab="")            
	abline(v=as.numeric(cutoff),lty=3,lwd=2)
	
	if (length(unique(color))!=1){
		legend("topright",c("Predicted male","Predicted female"),cex=2, pch=18, col=c("lightskyblue","darkorange1"))
	}

}








plotDiscrepancyGenders <- function(cutoff, predictedGender, covariates, XYMedians){
	x <- XYMedians$medianXU + XYMedians$medianXM
	y <- XYMedians$medianYM + XYMedians$medianYU
	diff <- log2(y)-log2(x)
	
	possibilities <- c("gender","Gender","sex","Sex","GENDER","SEX")
	sum <- sum(possibilities %in% colnames(covariates))
		if (sum>0){
		goodColumn <- possibilities[possibilities %in% colnames(covariates)][1]
		goodIndex <- match(goodColumn, colnames(covariates))
		givenGender <- as.character(covariates[,goodIndex])

        color <- predictedGender
        color[color == "M"] <- "lightskyblue"
        color[color == "F"] <- "darkorange1"

    	predictedGender[predictedGender=="lightskyblue"] <- "M"
	    predictedGender[predictedGender=="darkorange1"] <- "F"
	    diffGender <- rep(FALSE,length(predictedGender))
	    for (i in 1:length(predictedGender)){
	    	if (predictedGender[i]!=givenGender[i] && !is.na(givenGender[i])){
	    		diffGender[i] <- TRUE
	    	}
	    }
	    if (sum(diffGender)>0){
	    	color[diffGender] <- "black"
	    	 y <- jitter(rep(0,length(color)),factor=1.5)
		    cex <- rep(2,length(color))
		    y[diffGender] <- -0.2
		    cex[diffGender] <- 3
		    
		    plot(diff, y, 
			          ylim = c(-1,1), pch=18, cex=cex, col= color, 
			              yaxt="n", xlab="median CN(Y)  - median CN(X)",
			                  ylab="", main = "Gender prediction using the X and Y chromosomes intensities ",
			                      cex.main = 1.5, cex.lab = 1.5
			)
			                      
			                      
			abline(v=as.numeric(cutoff),lty=3,lwd=2)
			for (i in 1:length(diffGender)){
				if (diffGender[i]){
					text(diff[i],y[i]-0.2,names(XYMedians$medianXU[i]))
				}
			}
			if (length(unique(color))==3){
				legend("topright",c("Predicted male","Predicted female","Unmatching samples"),cex=1.5, pch=18, col=c("lightskyblue","darkorange1","black"))
			}
	    }    
	    
	}
}




