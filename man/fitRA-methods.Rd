\name{fitRA-methods}
\docType{methods}
\alias{fitRA-methods}
\alias{fitRA,AffyBatch-method}
\title{fit RA object}
\description{
This function fits a Residual-Array (RA) object from an object of class AffyBatch, using a robust linear model at 
the probe level.
}
\usage{
\S4method{fitRA}{AffyBatch}(object,verbose=FALSE)
}
\arguments{
  \item{object}{on object of class AffyBatch}
  \item{verbose}{print some information. Default to FALSE}
}
\details{
}
\value{
An onbject of class RA containing:
  \item{data}{a dataframe containing array effect estimates and residual variance}
  \item{ProbeID}{character. Probe ids}
}
\references{Stefano Calza, Wolfgang Raffelsberger, Alexander Ploner, Jose Sahel, Thierry Leveillard
and Yudi Pawitan. Filtering genes to improve sensitivity in oligonucleotide
microarray data analysis - 2007 NAR}
\author{Stefano Calza}
\note{
}
\seealso{\code{\link{RAplot}}, \code{\link{Flush}},
  \code{\link{FlushSet}},
  \code{\link[affyPLM]{fitPLM}}}
\examples{

\dontrun{
data(choe)
choe.RA <- fitRA(choe)
}
}
\keyword{methods}
\keyword{models}
\keyword{manip}


