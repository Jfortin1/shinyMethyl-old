\name{extractFromRGSet450k}
\alias{extractFromRGSet450k}
\title{
Extract necessary information for shinyMethyl from an RGMethylSet
}
\description{
Function that extracts the necessary data for the shinyMethyl application from an \code{RGMethylSet}
}
\usage{
extractFromRGSet450k(RGSet, file = "extractedData.Rda")
}
\arguments{
  \item{RGSet}{
  An \code{RGChannelSet} object (see the \code{minfi} package)
}
  \item{file}{
  A \code{string} specifying the name of the file to be saved
}
}
\value{
Returns a data frame containing the necessary data for the \code{shinyMethyl} application. In addition, the extracted data is saved in the current directory under the provided name \code{file}.
}
\author{
Jean-Philippe Fortin
}
\seealso{
\code{\link{extractFromTargets450k}}
\code{\link{minfi}}
}
\examples{
if (require(minfiData)){
	baseDir <- system.file("extdata", package = "minfiData")
    targets <- read.450k.sheet(baseDir)
    RGSet <- read.450k.exp(base = baseDir, targets = targets)
    
    extractedData <- extractFromRGSet450k(RGSet)
}   
}
