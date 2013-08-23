\name{runShinyMethyl}
\alias{runShinyMethyl}
\title{
Run the interactive shinyMethyl session
}
\description{
	Function to run the interactive shinyMethyl session after the necessary data has been extracted. 
}
\usage{
runShinyMethyl(extractedData)
}
\arguments{
  \item{extractedData}{
data extracted from an \code{RGChannelSet} or a targets \code{data.frame} using \code{extractFromRGSet450k} or \code{extractFromTargets450k} 
}
}
\value{
No value returned. Instead the shinyMethyl interactive session is launched.
}
\author{
Jean-Philippe Fortin
}
\seealso{
\code{\link{extractFromTargets450k}}
\code{\link{extractFromRGSet450k}}
}
\examples{
\dontrun{
if (require(minfiData)){
	baseDir <- system.file("extdata", package = "minfiData")
    targets <- read.450k.sheet(baseDir)
    targets$Sample_Plate <- substr(targets$Slide,1,7)
    RGSet <- read.450k.exp(base = baseDir, targets = targets)
    
    extractedData <- extractFromRGSet450k(RGSet)
    
    runShinyMethyl(extractedData)
}
}

}

