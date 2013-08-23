\name{runExampleTCGA}
\alias{runExampleTCGA}
\title{
Run the interactive shinyMethyl session from the shinyMethylData example dataset
}
\description{
Function to run the interactive shinyMethyl session with previously extracted data from 369 Head and Neck Cancer public samples from The Cancer Genome Atlas (TCGA). The exampled dataset is located in the shinyMethylData package. The dataset contains 310 tumor samples, 50 matched normal and 9 replicates of a control cell line. There are 8 different plates (batches). 
}
\usage{
runExampleTCGA()
}
\value{
No value returned. Instead the shinyMethyl interactive session is launched.
}
\references{
The Cancer Genome Atltas (TCGA) Head and Neck Cancer dataset:
http://cancergenome.nih.gov/cancersselected/headandneck
}
\author{
Jean-Philippe Fortin
}
\seealso{
\code{\link{runShinyMethyl}}
\code{\link[shinyMethylData]{TCGAExtractedData}}
}
\examples{
\dontrun{
if (require(shinyMethylData)){
	 runExampleTCGA()
}
}
}

