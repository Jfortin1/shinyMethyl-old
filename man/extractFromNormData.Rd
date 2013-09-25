\name{extractFromNormData}
\alias{extractFromNormData}
\title{
Extract necessary information from a normalized 450k dataset
}
\description{
Function that extracts the necessary data from a normalized 450k dataset containing only the methylated and unmethylated intensities. 
}
\usage{
extractFromNormData(normData, file = "extractedNormData.Rda")
}
\arguments{
  \item{normData}{
  A list containing two elements: \code{$methMatrix} is the matrix of intensities in the methylated channel, and  \code{$unmethMatrix} is the matrix of intensities in the unmethylated channel
}
  \item{file}{
  A \code{string} specifying the name of the file to be saved
}
}
\value{
Returns a data frame containing the necessary data for the \code{shinyMethyl} application from a normalized dataset. In addition, the extracted data is saved in the current directory under the provided name \code{file}.
}
\author{
Jean-Philippe Fortin
}
\seealso{
\code{\link{extractFromRGSet450k}}
\code{\link{minfi}}
}
\examples{
if (require(minfiData)){
	baseDir <- system.file("extdata", package = "minfiData")
    targets <- read.450k.sheet(baseDir)
    RGSet <- read.450k.exp(base = baseDir, targets = targets)
    rawData <- preprocessRaw(RGSet)
    methMatrix <- getMeth(rawData)
    unmethMatrix <- getUnmeth(rawData)
    normData <- list(methMatrix = methMatrix, unmethMatrix = unmethMatrix)
    
    extractedNormData <- extractFromNormData(normData)
}   
}
