#### Created by Jean-Philippe Fortin
#### Aug 19  2013



sourceDir <- system.file("shinyMethyl",package = "shinyMethyl")
sourceDir <- paste0(sourceDir,"/","plotFunctions.R")
source(sourceDir)
library(RColorBrewer)




############################################################








shinyServer(function(input, output) {
	

betaQuantiles <- shinyMethylData$betaQuantiles
mQuantiles <- shinyMethylData$mQuantiles
methQuantiles <- shinyMethylData$methQuantiles
unmethQuantiles <- shinyMethylData$unmethQuantiles
greenControls <- shinyMethylData$greenControls
redControls <- shinyMethylData$redControls
# oobControls <- shinyMethylData$oobControls
XYMedians <- shinyMethylData$XYMedians
sampleDistance <- shinyMethylData$sampleDistance
covariates <- shinyMethylData$pd
covariates <- covariates[match(colnames(shinyMethylData$betaQuantiles$II),rownames(covariates)),]
pca<- shinyMethylData$pcaInfo$scores


	sampleNames <- colnames(betaQuantiles[[1]])
	slideNames     <- substr(sampleNames,1,10)
	arrayNames    <- substr(sampleNames,12,17)
	plateNames <- substr(sampleNames,1,7)
	 
	
	designInfo <- data.frame(sampleNames = sampleNames,
												slideNames  = slideNames,
												arrayNames = arrayNames,
												plateNames = plateNames)
	
	#### Ordering:
	designInfo <- designInfo[order(plateNames,slideNames,arrayNames), ]
	order <- match(designInfo$sampleNames, 
				                colnames(betaQuantiles[[1]]))


    sampleNames <- designInfo$sampleNames
	slideNames     <- designInfo$slideNames
	arrayNames    <- designInfo$arrayNames
    plateNames    <- designInfo$plateNames
	
	plate <- as.numeric(as.factor(plateNames))
	slides <- unique(slideNames)
	names(slideNames) <- slideNames
	names(slides) <- slides

	for (i in 1:3){
		betaQuantiles[[i]]         <- betaQuantiles[[i]][,order]
		mQuantiles[[i]]             <- mQuantiles[[i]][,order]
		methQuantiles[[i]]        <- methQuantiles[[i]][,order]
		unmethQuantiles[[i]]    <- unmethQuantiles[[i]][,order]
	}
	
	for (i in 1:12){
		greenControls[[i]]    <- greenControls[[i]][,order]
		redControls[[i]]        <- redControls[[i]][,order]
	}
	
	XYMedians$medianXU = XYMedians$medianXU[order]
	XYMedians$medianXM = XYMedians$medianXM[order]
	XYMedians$medianYU = XYMedians$medianYU[order]
	XYMedians$medianYM = XYMedians$medianYM[order]

	covariates <- covariates[order,]
	pca <- pca[order,]
   
    controlNames <- names(greenControls)


    sampleColors <<- as.numeric(as.factor(plate))

sampleColors <- reactive(
 color <- as.numeric(as.factor(covariates[,match(input$phenotype,colnames(covariates))]))
)


	   	


 #####---- Plot of the quality control using internal controls
 output$internalControls <- renderPlot({
   
   palette(brewer.pal(8,"Set1"))

   controlStat <- returnControlStat(input$controlType, 
                                      greenControls = greenControls,
                                          redControls = redControls, 
                                              controlNames = controlNames ) 	
                                              
    plotInternalControls(y = controlStat,
                                             col = sampleColors(),
                                                 main = input$controlType,
                                                     sampleNames = sampleNames,
                                                         selectedSlides = input$slides,
                                                             slideNames = slideNames)
                                                             
    addHoverPoints(y = controlStat, 
                                      xSelected = input$controlsHover$x, 
                                          ySelected = input$controlsHover$y, 
                                              sampleNames = sampleNames)

  })













  #####---- Plot of the densities
 output$rawDensities <- renderPlot({
      palette(brewer.pal(8,"Set1"))

	  index = match(input$probeType,c("I Green","I Red","II"))
	  
	  controlStat <- returnControlStat(input$controlType, 
                                      greenControls = greenControls,
                                          redControls = redControls, 
                                              controlNames = controlNames ) 		
                          
	   if (input$mOrBeta=="Beta-value"){
	   	from = 0
	   	to = 1
	   	densitiesPlot( quantiles = betaQuantiles[[index]],
	   	                              main = "BETA-VALUE DENSITIES", 
	   	                              xlab = "Beta-values",
	   	                              xlim  = c(-0.2,1.2),
	   	                              ylim = c(0,7),
	   	                              bw    <- input$bandwidth,
	   	                              slideNames = slideNames,
	   	                              solidLine = input$solidLine,
	   	                              mean = input$mean,
	   	                              slides = input$slides,
	   	                              col = sampleColors(), from=-4, to =4)
	   	
        abline(v=0,lty=3,lwd=3)
	    abline(v=1,lty=3,lwd=3)
	  
	    addHoverDensity(y = controlStat, 
                                xSelected = input$controlsHover$x, 
                                ySelected = input$controlsHover$y, 
                                sampleNames = sampleNames,
                                slides = input$slides, 
                                quantiles = betaQuantiles[[index]],
                                bw = input$bandwidth,
                                col = sampleColors()
      )    
      
	  } else {
	  	
	  from = -10
	  to = 10
	  	
	  densitiesPlot( quantiles = mQuantiles[[index]],
	   	                              main = "M-VALUE DENSITIES", 
	   	                              xlab = "M-values",
	   	                              xlim  = c(-8,8),
	   	                              ylim = c(0,0.5),
	   	                              bw    <- input$bandwidth2,
	   	                              slideNames = slideNames,
	   	                              solidLine = input$solidLine,
	   	                              mean = input$mean,
	   	                              slides = input$slides,
	   	                              col = sampleColors(), from=from, to=to)

	  abline(v=0,lty=3,lwd=3)
	  	
	  addHoverDensity(y = controlStat, 
                                      xSelected = input$controlsHover$x, 
                                      ySelected = input$controlsHover$y, 
                                      sampleNames = sampleNames,
                                      slides = input$slides, 
                                      quantiles = mQuantiles[[index]],
                                      bw  = input$bandwidth2,
                                      col = sampleColors()
      ) 
         
	 }

 })




      






#####---- Plot Quality control 
output$medianChannels <- renderPlot({
	
	plotQC(unmethQuantiles = unmethQuantiles, methQuantiles = methQuantiles,
	               slides = input$slides, 
	                    slideNames = slideNames, 
	                        col = sampleColors())
	 
   controlStat <- returnControlStat(input$controlType, 
                                      greenControls = greenControls,
                                          redControls = redControls, 
                                              controlNames = controlNames ) 	
    	
   addHoverQC(y = controlStat, 
                                      xSelected = input$controlsHover$x, 
                                          ySelected = input$controlsHover$y, 
                                              sampleNames = sampleNames,
                                                 unmethQuantiles = unmethQuantiles,
                                                     methQuantiles = methQuantiles,
                                                     slides = input$slides
     )
})











### Plot gender clustering
output$genderClustering <- renderPlot({

	cutoff <- -3
	if (!is.null(input$genderCutoff)){
		cutoff <- input$genderCutoff$x
	}
	
	plotPredictedGender(cutoff = cutoff, XYMedians = XYMedians)
	
	predictedGender <- returnPredictedGender(cutoff = cutoff, XYMedians = XYMedians)
	
	plotDiscrepancyGenders(cutoff = cutoff, 
	     predictedGender = predictedGender,
	          covariates = covariates, 
	              XYMedians = XYMedians
	  )
	
})














data <- reactive({
	
   cutoff <- -3
   if (!is.null(input$genderCutoff)){
		cutoff <- input$genderCutoff$x
   }
   
   return(returnPredictedGender(cutoff = cutoff, XYMedians))

 })









output$downloadClusters <- downloadHandler(

		 filename <- "predictedGender.csv",
		 content <- function(con){
		 	write.csv(data(),con)
		 }
		 
 )










 diffGenders <- reactive({
 	
 	predictedGender <- data()
 	diffs <- c()
 	possibilities <- c("gender","Gender","sex","Sex","GENDER","SEX")
	sum <- sum(possibilities %in% colnames(covariates))
	 if (sum>0){
	 	   goodColumn <- possibilities[possibilities %in% colnames(covariates)][1]
		   goodIndex <- match(goodColumn, colnames(covariates))
		   givenGender <- as.character(covariates[,goodIndex])
		   
		    diffGender <- rep(FALSE,length(predictedGender))
		  for (i in 1:length(predictedGender)){
	    	if (predictedGender[i]!=givenGender[i] && !is.na(givenGender[i])){
	    		diffGender[i] <- TRUE
	    	}
	    }
		    if (sum(diffGender)>0){
		    diffs <- names(XYMedians$medianXU)[diffGender]
		    }
		    
	}
    return(diffs)
    
 })
	
	






output$diffPrint <- renderPrint({

	possibilities <- c("gender","Gender","sex","Sex","GENDER","SEX")
	sum <- sum(possibilities %in% colnames(covariates))		
		
	if (length(diffGenders())!=0){
		print(diffGenders())
	} else if (sum>0) {
		print("The provided gender agrees with the predicted gender for all samples")
	} else {
		print("No gender was provided in the phenotype data")
	}
})
	 

	 







	   	
	   	
	   	

output$pcaPlot <- renderPlot({
	    palette(brewer.pal(8,"Set1"))

	    xMin <- min(pca[,as.numeric(input$pc1)])
	    xMax <- max(pca[,as.numeric(input$pc1)])
	    xRange <- xMax - xMin
	    xlim <- c(xMin-0.05*xRange, xMax+0.20*xRange)
	    
	    xlab <- paste("PC",as.numeric(input$pc1), " scores", sep="");
		ylab <- paste("PC",as.numeric(input$pc2), " scores", sep="");
	    plot(pca[,as.numeric(input$pc1)],
	         pca[,as.numeric(input$pc2)],
	             col = sampleColors(),   pch = 18,  cex = 2,
	                 xlab=  xlab, ylab = ylab, xlim = xlim,
	                     main= "Principal component analysis (PCA)",
	                         cex.main = 1.5, cex.lab = 1.5
	     )
	     
	     uColor  <- unique(sampleColors())
         uCov    <- unique( covariates[,match(input$phenotype,colnames(covariates))])
	     legend("bottomright",
	         legend = uCov, pch= 18, col = uColor, cex = 1.5,
		         title = input$phenotype
		  )
	      grid()

})







 
 
 
 
### For the moment assuming that covariates are factors...
output$modelPrint <- renderPrint({
	x <- factor(plate)
    y <- shinyMethylData$pcaInfo$scores[,as.numeric(input$pcToExplore)]
    
     if (exists("covariates")){
	   	cov <- covariates[match(rownames(shinyMethylData$pcaInfo$scores),rownames(covariates)),]
	   	x <- (as.factor(cov[,match(input$covToRegress,colnames(cov))]))
	}
    
    
   model <- lm(y~ x)
    summary(model)
	 })








output$arrayDesign <- renderPlot({
	

	color <- covariates[,match(input$phenotype,colnames(covariates))]
    plotDesign450k(as.character(sampleNames), covariates = color , legendTitle = input$phenotype)
 

})





output$arrayDesignLegend <- renderPlot({
	
	color <- covariates[,match(input$phenotype,colnames(covariates))]

    plotLegendDesign450k(as.character(sampleNames), covariates =color , legendTitle = input$phenotype)
 })

})




	
