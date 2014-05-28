setGeneric("getGreenControls", function(object, ...) standardGeneric("getGreenControls"))
setGeneric("getRedControls", function(object, ...) standardGeneric("getRedControls"))
setGeneric("getPCA", function(object, ...) standardGeneric("getPCA"))
setGeneric("orderByName", function(object, ...) standardGeneric("orderByName"))

setClass("shinyMethylSet", 
         representation(sampleNames = "character",
                        phenotype = "data.frame",
                        mQuantiles ="list",
                        betaQuantiles = "list",
                        methQuantiles = "list",
                        unmethQuantiles = "list",
                        cnQuantiles = "list",
                        greenControls = "list",
                        redControls = "list",
                        pca = "list",
                        originObject = "character"
                        )
         )
	
shinyMethylSet <- function(sampleNames = new("character"), 
                           phenotype = new("data.frame"),  
                           mQuantiles = new(vector("list",5)),
                           betaQuantiles = new(vector("list",5)),
                           methQuantiles = new(vector("list",5)),
                           unmethQuantiles = new(vector("list",5)),
                           cnQuantiles = new(vector("list",5)),
                           greenControls = new(vector("list",12)),
                           redControls = new(vector("list",12)),
                           pca = new("list"),
                           originObject = new("character")
                           ) {
    set <- new("shinyMethylSet", 
               sampleNames = sampleNames,
               phenotype = phenotype,
               mQuantiles = mQuantiles,
               betaQuantiles = betaQuantiles,
               methQuantiles = methQuantiles,
               unmethQuantiles = unmethQuantiles,
               cnQuantiles = cnQuantiles,
               greenControls = greenControls,
               redControls = redControls,
               pca = pca,
               originObject = originObject
               )
    set
}

setMethod("orderByName", signature(object = "shinyMethylSet"), 
          function(object){
              sampleNames    <- object@sampleNames
              slideNames     <- substr(sampleNames,1,10)
              arrayNames     <- substr(sampleNames,12,17)
              plateNames     <- substr(sampleNames,1,6)
              
              designInfo <- data.frame(sampleNames = sampleNames,
                                       slideNames  = slideNames,
                                       arrayNames  = arrayNames,
                                       plateNames  = plateNames
                                       )
              designInfo <- designInfo[order(plateNames,slideNames,arrayNames), ]
              o <- match(designInfo$sampleNames, sampleNames)
              for (i in 1:length(object@betaQuantiles)){
                  object@betaQuantiles[[i]] <- object@betaQuantiles[[i]][,o]
              }
              for (i in 1:length(object@mQuantiles)){
                  object@mQuantiles[[i]] <- object@mQuantiles[[i]][,o]
              }
              for (i in 1:length(object@methQuantiles)){
                  object@methQuantiles[[i]] <- object@methQuantiles[[i]][,o]
              }
              for (i in 1:length(object@unmethQuantiles)){
                  object@unmethQuantiles[[i]] <- object@unmethQuantiles[[i]][,o]
              }
              for (i in 1:length(object@cnQuantiles)){
                  object@cnQuantiles[[i]] <- object@cnQuantiles[[i]][,o]
              }

              object@sampleNames <- object@sampleNames[o]
		object@pca$scores  <- object@pca$scores[o,]
		object@pca$percs  <- object@pca$percs[o]
              
              object@phenotype <- object@phenotype[o,]
              
              for (i in 1:length(object@greenControls)){
                  object@greenControls[[i]] <- object@greenControls[[i]][,o]
                  object@redControls[[i]]   <- object@redControls[[i]][,o]
              }
              
              object
          })

setMethod("show",signature(object="shinyMethylSet"),
          function(object){
              cat(" Object: ",class(object),"\n")
              n <- length(object@sampleNames)
              cat(" Number of samples: ",n, "\n")
              nCovs <- ncol(object@phenotype)
              cat(" Phenotype: ",nCovs, "covariates \n")
              cat(" Method: ", object@originObject, "\n")
          })

setMethod("getBeta",signature(object="shinyMethylSet"),
          function(object){
              object@betaQuantiles
          })

setMethod("pData",signature(object="shinyMethylSet"),
          function(object){
              object@phenotype
          })

setMethod("getMeth",signature(object="shinyMethylSet"),
          function(object){
              object@methQuantiles
          })

setMethod("getUnmeth",signature(object="shinyMethylSet"),
          function(object){
              object@unmethQuantiles
          })

setMethod("getM",signature(object="shinyMethylSet"),
          function(object){
              object@mQuantiles
          })

setMethod("getCN",signature(object="shinyMethylSet"),
          function(object){
              object@cnQuantiles
          })

setMethod("getGreenControls",signature(object="shinyMethylSet"),
          function(object){
              object@greenControls
          })

setMethod("getRedControls",signature(object="shinyMethylSet"),
          function(object){
              object@redControls
          })

setMethod("getPCA",signature(object="shinyMethylSet"),
          function(object){
              object@pca
          })

setMethod("sampleNames",signature(object="shinyMethylSet"),
          function(object){
              object@sampleNames
          })


## setMethod("getSex", signature(object="shinyMethylSet"),
## function(object, cutoff = -2){
## nQuantiles <- nrow(object@cnQuantiles)
## mid <- as.integer(nQuantiles/2)
## xMed <- object@cnQuantiles$X[mid,]
## yMed <- object@cnQuantiles$Y[mid,]
## dd <- yMed - xMed
## k <- kmeans(dd, centers = c(min(dd), max(dd)))

## sex0 <- ifelse(dd < cutoff, "F", "M")
## sex0 <- .checkSex(sex0)
## sex1 <- ifelse(k$cluster == which.min(k$centers), "F", "M")
## sex1 <- .checkSex(sex1)

## if(!identical(sex0,sex1))
## warning("An inconsistency was encountered while determining sex. One possibility is that only one sex is present. 
## We recommend further checks, for example with the plotSex function.")
## df <- DataFrame(xMed = xMed, yMed = yMed, predictedSex = sex0)
## rownames(df) <- object@sampleNames
## df
## })

.checkSex <- function(sex) {
    if(! (is.character(sex) && !any(is.na(sex)) && all(sex %in% c("M", "F"))))
        stop("'sex' seems wrong (needs to be a character, without missing values, of 'M' and 'F'")
    sex
}

## Function to combine two shinyMethylSet. 	
shinyCombine <- 	function(shinyMethylSet1, shinyMethylSet2){
    if (!is(shinyMethylSet1, "shinyMethylSet") | !is(shinyMethylSet2, "shinyMethylSet") ){
        stop("Both objects must be shinyMethylSet")
    }
    c.shinyMethylSet <- new("shinyMethylSet")
    c.shinyMethylSet@sampleNames <- c(shinyMethylSet1@sampleNames, shinyMethylSet2@sampleNames)
    c.shinyMethylSet@phenotype   <- rbind(shinyMethylSet1@phenotype, shinyMethylSet2@phenotype)
    
    ## Merging the quantile matrices
    for (i in 1:5){
        c.shinyMethylSet@mQuantiles[[i]] <- cbind(shinyMethylSet1@mQuantiles[[i]],shinyMethylSet2@mQuantiles[[i]])
        c.shinyMethylSet@betaQuantiles[[i]] <- cbind(shinyMethylSet1@betaQuantiles[[i]],shinyMethylSet2@betaQuantiles[[i]])
        c.shinyMethylSet@methQuantiles[[i]] <- cbind(shinyMethylSet1@methQuantiles[[i]],shinyMethylSet2@methQuantiles[[i]])
        c.shinyMethylSet@unmethQuantiles[[i]] <- cbind(shinyMethylSet1@unmethQuantiles[[i]],shinyMethylSet2@unmethQuantiles[[i]])
        c.shinyMethylSet@cnQuantiles[[i]] <- cbind(shinyMethylSet1@cnQuantiles[[i]],shinyMethylSet2@cnQuantiles[[i]])
    }
    names(c.shinyMethylSet@mQuantiles) <- names(c.shinyMethylSet@betaQuantiles) <- 
        names(c.shinyMethylSet@methQuantiles) <- names(c.shinyMethylSet@unmethQuantiles) <-
            names(c.shinyMethylSet@cnQuantiles) <- names(shinyMethylSet1@mQuantiles)
    
    ## Merging the control matrices
    for (i in 1:14){
        c.shinyMethylSet@greenControls[[i]] <- cbind(shinyMethylSet1@greenControls[[i]],shinyMethylSet2@greenControls[[i]])
        c.shinyMethylSet@redControls[[i]] <- cbind(shinyMethylSet1@redControls[[i]],shinyMethylSet2@redControls[[i]])
    }
    names(c.shinyMethylSet@greenControls) <- names(c.shinyMethylSet@redControls) <- names(shinyMethylSet1@redControls)
    c.shinyMethylSet@pca <- list() # PCA information is lost. 
    c.shinyMethylSet@originObject <- "Merging"
    c.shinyMethylSet
}

