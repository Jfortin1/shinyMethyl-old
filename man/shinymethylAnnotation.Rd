\name{shinymethylAnnotation}
\alias{shinymethylAnnotation}
\alias{autosomal}
\alias{chrY}
\alias{chrX}
\docType{data}
\title{
Chromosome annotation for the 450k array probes
}
\description{
Chromosome annotation of the 450k array probes. Annotation was obtained with the \code{minfi} package with the function \code{minfi:::getAnnotations}. The annotation is used to extract the median intensity of the methylated and unmethylated channels for the X and Y chromosomes probes separately in order to predict the sample genders. 
}
\usage{data(shinymethylAnnotation)}
\format{
\code{chrX}, \code{chrY} and \code{autosomal}: vectors containing the probe names of the probes respectively mapped to the X-chromosome, the Y-chromosome and to autosomal chromosomes
}
\keyword{datasets}
