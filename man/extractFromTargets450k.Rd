\name{extractFromTargets450k}
\alias{extractFromTargets450k}
\title{
Extract necessary information for shinyMethyl from a data frame
}
\description{
Extract necessary information for shinyMethyl from a targets \code{data.frame}.
}
\usage{
extractFromTargets450k(targets, file = "extractedData.Rda", bySubsets = FALSE)
}

\arguments{
  \item{targets}{
  \code{data.frame} used to read the IDAT files, as described in the \code{minfi} vignette. See details 
}
  \item{file}{
A \code{string} specifying the name of the file to be saved
}
  \item{bySubsets}{
  If \code{TRUE}, the function will extract the data by loading 10 samples at a time. 
}
}
\value{
Returns a data frame containing the necessary data for the \code{shinyMethyl} application. In addition, the extracted data is saved in the current directory under the provided name \code{file}.
}
\details{
\code{targets} is assumed to have a column named \code{Basename} pointing to the base name of a two color IDAT file, ie. a name that can be made into a real IDAT file by appending either \code{_Red.idat} or \code{_Grn.idat}. Other columns will be used as default covariates for the \code{shinyMethyl} application
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
    
    extractedData <- extractFromTargets450k(targets = targets)
}
}

