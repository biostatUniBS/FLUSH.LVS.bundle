\name{FlushSet-methods}
\docType{methods}
\alias{FlushSet-methods}
\alias{FlushSet,RA-method}
\alias{FlushSet}
\title{Compute Flushing criteria}
\description{
  This function fits a quantile regression on an object of class RA
}
\usage{
\S4method{FlushSet}{RA}(object, proportion = 0.6, delta = NULL, 
        lambda = delta, df=10,...)
}
\arguments{
  \item{object}{An object of class "RA"}
  \item{proportion}{Proportion of genes to filter out, i.e. the quantile
    to estimate in the quantile regression}
  \item{delta}{weighting parameter}
  \item{lambda}{weighting parameter}
  \item{df}{Degrees of freedom of the B-spline}
  \item{\dots}{Not used}
}
\details{
}
\value{
And object of class 'RA' with a new slot
}
\references{Stefano Calza, Wolfgang Raffelsberger, Alexander Ploner, Jose Sahel, Thierry Leveillard
and Yudi Pawitan. Filtering genes to improve sensitivity in oligonucleotide
microarray data analysis - 2007 NAR}
\author{Stefano Calza <calza@med.unibs.it> }
\note{
}
\seealso{}
\examples{

\dontrun{demo(FLUSH.tour)
FlushSet(choe.RA)
FlushSet(choe.RA, delta = .5)
}
}
\keyword{methods}
\keyword{manip}
